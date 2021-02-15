#!/usr/bin/env python

from Bio import SeqIO
import click
from click_option_group import optgroup
import copy
import gzip
from ignite.contrib.engines import common
from ignite.engine import (Events, create_supervised_trainer,
    create_supervised_evaluator)
from ignite.handlers import (ModelCheckpoint, EarlyStopping,
    global_step_from_engine)
from ignite.metrics import Accuracy, Loss
import json
# import math
# import matplotlib.pyplot as plt
import numpy as np
import os
# import pandas as pd
# import seaborn as sns
# from sklearn.metrics import (
#     average_precision_score, precision_recall_curve,
#     roc_auc_score, roc_curve,
#     matthews_corrcoef
# )
# from sklearn.model_selection import train_test_split
# from time import time
import torch
# import torch.nn as nn
# # from torch.optim.lr_scheduler import ReduceLROnPlateau
from torch.utils.data import DataLoader, TensorDataset

# Models
from models import one_hot_encode, reverse_complement_one_hot_encoding
from models.danq import DanQ, get_criterion, get_optimizer

# Globals
scripts_dir = os.path.dirname(os.path.realpath(__file__))

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "training_file", type=click.Path(exists=True, resolve_path=True)
)
@click.argument(
    "validation_file", type=click.Path(exists=True, resolve_path=True)
)
@click.option(
    "-d", "--dummy-dir",
    help="Dummy directory.",
    type=click.Path(resolve_path=True),
    default="/tmp/",
    show_default=True
)
@click.option(
    "-o", "--output-dir",
    help="Output directory.",
    type=click.Path(resolve_path=True),
    default="./",
    show_default=True
)
@click.option(
    "-l", "--length",
    help="Sequence length (in bp).",
    type=int,
    default=201,
    show_default=True
)
@click.option(
    "-r", "--reverse-complement",
    help="Reverse complement training sequences.",
    is_flag=1
)
@click.option(
    "-t", "--threads",
    help="Threads to use.",
    type=int,
    default=1,
    show_default=True
)
@optgroup.group("Trainer")
@optgroup.option(
    "--epochs",
    help="Number of epochs to train.",
    type=int,
    default=100,
    show_default=True
)
@optgroup.option(
    "--log-interval",
    help="Number of batches to wait before logging training status.",
    type=int,
    default=10,
    show_default=True
)
@optgroup.option(
    "--patience",
    help="Number of epochs to wait before stopping training if no improvement.",
    type=int,
    default=10,
    show_default=True
)
@optgroup.option(
    "--train-batch-size",
    help="Batch size for training.",
    type=int,
    default=100,
    show_default=True
)
@optgroup.option(
    "--val-batch-size",
    help="Batch size for validation.",
    type=int,
    default=1000,
    show_default=True
)
@optgroup.group("Optimizer")
@optgroup.option(
    "--lr",
    help="Learning rate.",
    type=float,
    default=1e-4,
    show_default=True
)
# @optgroup.option(
#     "--momentum",
#     help="Momentum factor.",
#     type=float,
#     default=0.9,
#     show_default=True
# )
# @optgroup.option(
#     "--weight-decay",
#     help="Weight decay.",
#     type=float,
#     default=1e-6,
#     show_default=True
# )

def main(**params):

    # Initialize
    statistics = {}
    device = "cpu"
    if torch.cuda.is_available():
        device = "cuda"
    model = DanQ(params["length"]).to(device)
    criterion = get_criterion()
    # optimizer = get_optimizer(model.parameters(), params["lr"],
    #     params["momentum"], params["weight_decay"])
    optimizer = get_optimizer(model.parameters(), params["lr"])

    # Create output dir
    if not os.path.exists(params["output_dir"]):
        os.makedirs(params["output_dir"])

    # Get DataLoaders
    train_loader, val_loader = __get_data_loaders(params["training_file"],
        params["validation_file"], params["train_batch_size"],
        params["val_batch_size"], params["reverse_complement"],
        params["threads"])

    # Trainer
    trainer = create_supervised_trainer(model, optimizer, criterion,
        device=device)
    @trainer.on(Events.ITERATION_COMPLETED(every=params["log_interval"]))
    def log_training_loss(trainer):
        print(f"Epoch[{trainer.state.epoch}] Loss: {trainer.state.output:.4f}")
    @trainer.on(Events.EPOCH_COMPLETED)
    def log_training_results(trainer):
        train_evaluator.run(train_loader)
        metrics = train_evaluator.state.metrics
        print((
            f"Training Results - Epoch: {trainer.state.epoch}  "
            f"Avg accuracy: {metrics['accuracy']:.4f} "
            f"Avg loss: {metrics['loss']:.4f}"
        ))
        statistics.setdefault(trainer.state.epoch, {})
        statistics[trainer.state.epoch].setdefault("train",
            copy.deepcopy(metrics))
    @trainer.on(Events.EPOCH_COMPLETED)
    def log_validation_results(trainer):
        val_evaluator.run(val_loader)
        metrics = val_evaluator.state.metrics
        print((
            f"Validation Results - Epoch: {trainer.state.epoch}  "
            f"Avg accuracy: {metrics['accuracy']:.4f} "
            f"Avg loss: {metrics['loss']:.4f}"
        ))
        statistics.setdefault(trainer.state.epoch, {})
        statistics[trainer.state.epoch].setdefault("train",
            copy.deepcopy(metrics))

    # Evaluators
    metrics = {"accuracy": Accuracy(output_transform=output_transform),
        "loss": Loss(criterion)}
    train_evaluator = create_supervised_evaluator(model, metrics=metrics,
        device=device)
    val_evaluator = create_supervised_evaluator(model, metrics=metrics,
        device=device)

    # Checkpoint
    model_checkpoint = ModelCheckpoint(params["output_dir"], "best", n_saved=1,
        score_function=score_function, score_name="val_acc",
        global_step_transform=global_step_from_engine(trainer))
    val_evaluator.add_event_handler(Events.COMPLETED, model_checkpoint,
        {"model": model})

    # Eearly stop
    early_stop = EarlyStopping(patience=params["patience"],
        score_function=score_function, trainer=trainer)
    val_evaluator.add_event_handler(Events.COMPLETED, early_stop)

    # Train
    trainer.run(train_loader, max_epochs=params["epochs"])

    # Write
    # json_file = os.path.join(params["output_dir"])
    print(statistics)

def __get_data_loaders(training_file, val_file, train_batch_size=100, 
    val_batch_size=1000, reverse_complement=False, threads=1):

    # Initialize
    train_Xs = []
    train_ys = []
    val_Xs = []
    val_ys = []

    # Xs / ys
    with gzip.open(training_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            _, y = record.description.split()
            train_Xs.append(one_hot_encode(str(record.seq).upper()))
            train_ys.append([float(y)])
    with gzip.open(val_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            _, y = record.description.split()
            val_Xs.append(one_hot_encode(str(record.seq).upper()))
            val_ys.append([float(y)])

    # Reverse complement
    if reverse_complement:
        n = len(train_Xs)
        for i in range(n):
            encoded_seq = reverse_complement_one_hot_encoding(train_Xs[i])
            train_Xs.append(encoded_seq)
            train_ys.append(train_ys[i])

    # TensorDatasets
    train_set = TensorDataset(torch.Tensor(train_Xs),
        torch.Tensor(train_ys))
    val_set = TensorDataset(torch.Tensor(val_Xs),
        torch.Tensor(val_ys))

    # DataLoaders
    kwargs = dict(batch_size=train_batch_size, shuffle=True,
        num_workers=threads)
    train_loader = DataLoader(train_set, **kwargs)
    kwargs = dict(batch_size=val_batch_size, shuffle=True,
        num_workers=threads)
    val_loader = DataLoader(val_set, **kwargs)

    return(train_loader, val_loader)

def output_transform(output):
    y_pred, y = output
    y_pred = y_pred.greater_equal(0.5)
    return(y_pred, y)

def score_function(engine):
    return(engine.state.metrics["accuracy"])

if __name__ == "__main__":
    main()