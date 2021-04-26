#!/usr/bin/env python

from Bio import SeqIO
import click
from click_option_group import optgroup
import copy
import gzip
from ignite.contrib.handlers.tqdm_logger import ProgressBar
from ignite.contrib.metrics import GpuInfo, ROC_AUC
from ignite.engine import (Events, create_supervised_trainer,
    create_supervised_evaluator)
from ignite.handlers import (ModelCheckpoint, EarlyStopping,
    global_step_from_engine)
from ignite.metrics import Loss
import json
# import math
# import matplotlib.pyplot as plt
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

# Local imports
from utils.architectures.cam import Model, get_criterion, get_optimizer
from utils.metrics import PearsonR
from utils.sequence import one_hot_encode, reverse_complement_one_hot_encoding

# Globals
scripts_dir = os.path.dirname(os.path.realpath(__file__))
bar_format="{desc}[{n_fmt}/{total_fmt}] {percentage:3.0f}%|{bar:20}{postfix} [{elapsed}<{remaining}]"

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
    "-r", "--reverse-complement",
    help="Reverse complement train sequences.",
    is_flag=1
)
@click.option(
    "-t", "--threads",
    help="Number of threads to use.",
    type=int,
    default=1,
    show_default=True
)
@optgroup.group("CAM")
@optgroup.option(
    "--seq-length",
    help="Sequence length (in bp).",
    type=int,
    required=True
)
@optgroup.option(
    "--n-cnn-units",
    help="Number of individual CNN units.",
    type=int,
    default=10,
    show_default=True
)
@optgroup.option(
    "--output",
    help="Output type.",
    type=click.Choice(["binary", "linear"]),
    default="binary",
    show_default=True
)
@optgroup.group("Ignite")
@optgroup.option(
    "--epochs",
    help="Number of epochs to train.",
    type=int,
    default=100,
    show_default=True
)
@optgroup.option(
    "--patience",
    help="If no improvement, number of epochs to wait before stopping training.",
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
    default=1e-03,
    show_default=True
)

def main(**params):

    # Initialize
    statistics = {}
    device = "cpu"
    if torch.cuda.is_available():
        device = "cuda"
    model = Model(params["seq_length"], params["n_cnn_units"]).to(device)
    criterion = get_criterion(params["output"])
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
    @trainer.on(Events.EPOCH_COMPLETED)
    def log_results(trainer):
        train_evaluator.run(train_loader)
        metrics = train_evaluator.state.metrics
        log_message = f"Training - Epoch: {trainer.state.epoch}"
        for metric in metrics:
            log_message += f", {metric}: {metrics[metric]:.4f}"
        pbar.log_message(log_message)
        statistics.setdefault(trainer.state.epoch, {})
        statistics[trainer.state.epoch].setdefault("training",
            copy.deepcopy(metrics))
        val_evaluator.run(val_loader)
        metrics = val_evaluator.state.metrics
        log_message = f"Validation - Epoch: {trainer.state.epoch}"
        for metric in metrics:
            log_message += f", {metric}: {metrics[metric]:.4f}"
        pbar.log_message(log_message)
        statistics[trainer.state.epoch].setdefault("validation",
            copy.deepcopy(metrics))

    # Progress bar
    GpuInfo().attach(trainer, name="gpu")
    pbar = ProgressBar(bar_format=bar_format, persist=True)
    pbar.attach(trainer, ["gpu:0 mem(%)", "gpu:0 util(%)"])

    # Evaluators
    metrics = {"loss": Loss(criterion)}
    if params["output"] == "binary":
        metrics.setdefault("aucROC", ROC_AUC(output_transform))
        # score_function = score_function_binary
    else:
        metrics.setdefault("R", PearsonR())
        # score_function = score_function_linear
    train_evaluator = create_supervised_evaluator(model, metrics=metrics,
        device=device)
    val_evaluator = create_supervised_evaluator(model, metrics=metrics,
        device=device)

    # Checkpoint
    model_checkpoint = ModelCheckpoint(params["output_dir"], "best", n_saved=1,
        global_step_transform=global_step_from_engine(trainer),
        score_function=score_function)
    val_evaluator.add_event_handler(Events.COMPLETED, model_checkpoint,
        {"model": model})

    # Eearly stop
    early_stop = EarlyStopping(patience=params["patience"],
        score_function=score_function, trainer=trainer)
    val_evaluator.add_event_handler(Events.COMPLETED, early_stop)

    # Train
    trainer.run(train_loader, max_epochs=params["epochs"])

    # Write
    json_file = os.path.join(params["output_dir"], "statistics.json")
    with open(json_file, "wt") as handle:
        handle.write(json.dumps(statistics, indent=4))

def __get_data_loaders(training_file, val_file, train_batch_size=100, 
    val_batch_size=1000, reverse_complement=False, threads=1):

    # Initialize
    train_Xs = []
    train_ys = []
    val_Xs = []
    val_ys = []

    # Xs / ys
    handle = __get_handle(training_file)
    for record in SeqIO.parse(handle, "fasta"):
        _, y = record.description.split()
        train_Xs.append(one_hot_encode(str(record.seq).upper()))
        train_ys.append([float(y)])
    handle.close()
    handle = __get_handle(val_file)
    for record in SeqIO.parse(handle, "fasta"):
        _, y = record.description.split()
        val_Xs.append(one_hot_encode(str(record.seq).upper()))
        val_ys.append([float(y)])
    handle.close()

    # Reverse complement
    if reverse_complement:
        n = len(train_Xs)
        for i in range(n):
            encoded_seq = reverse_complement_one_hot_encoding(train_Xs[i])
            train_Xs.append(encoded_seq)
            train_ys.append(train_ys[i])
        # n = len(val_Xs)
        # for i in range(n):
        #     encoded_seq = reverse_complement_one_hot_encoding(val_Xs[i])
        #     val_Xs.append(encoded_seq)
        #     val_ys.append(val_ys[i])

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

def __get_handle(file_name):
    if file_name.endswith("gz"):
        handle = gzip.open(file_name, "rt")
    else:
        handle = open(file_name, "rt")
    return(handle)

def output_transform(output):
    y_pred, y = output
    y_pred = y_pred.greater_equal(0.5)
    return(y_pred, y)

def score_function(engine):
    return(engine.state.metrics["loss"]*-1)

def score_function_binary(engine):
    return(engine.state.metrics["aucROC"])

def score_function_linear(engine):
    return(engine.state.metrics["R"])

if __name__ == "__main__":
    main()