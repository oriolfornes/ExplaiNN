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
from ignite.metrics import Loss
from io import StringIO
from itertools import cycle, islice
import json
import logging
import numpy as np
import os
from scipy.stats import pearsonr, spearmanr
from sklearn.model_selection import PredefinedSplit
import sys
import torch
from torch.utils.data import DataLoader, TensorDataset

# Local imports
from utils.architectures import CAM as Model, get_criterion, get_optimizer
from utils.metrics import PearsonR, SpearmanR
from utils.sequence import one_hot_encode, reverse_complement_one_hot_encoding

# CUDA
device = "cpu"
if torch.cuda.is_available():
    device = "cuda"

# Globals
scripts_dir = os.path.dirname(os.path.realpath(__file__))
bar_format = "{desc}[{n_fmt}/{total_fmt}] {percentage:3.0f}%"
bar_format += "|{bar:20}{postfix} [{elapsed}<{remaining}]"

# Logger
streamer = StringIO()
logging.basicConfig(format="%(asctime)-15s %(message)s", level=logging.INFO,
    stream=streamer)
logger = logging.getLogger("CAM")

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "training_file", type=click.Path(exists=True, resolve_path=True)
)
@click.argument(
    "test_file", type=click.Path(exists=True, resolve_path=True)
)
@click.option(
    "-d", "--debug-mode",
    help="Debug mode.",
    is_flag=True
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
    "-s", "--splits",
    help="Number of calibration splits.",
    default=3,
    show_default=True
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
    "--cnn-units",
    help="Number of individual CNN units.",
    type=int,
    default=16,
    show_default=True
)
@optgroup.group("Ignite")
@optgroup.option(
    "--batch-size",
    help="Batch size for training.",
    type=int,
    default=64,
    show_default=True
)
@optgroup.option(
    "--epochs",
    help="Number of epochs to train.",
    type=int,
    default=128,
    show_default=True
)
@optgroup.option(
    "--patience",
    help="If no improvement, number of epochs to wait before stopping training.",
    type=int,
    default=20,
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

    # Create output dir
    if not os.path.exists(params["output_dir"]):
        os.makedirs(params["output_dir"])

    # Get Xs/ys
    Xs, ys = __get_Xs_ys(params["training_file"])

    # Get sequence length
    seq_length = max(Xs[0].shape)

    # Model calibration
    json_file = os.path.join(params["output_dir"], "calibration.json")
    if not os.path.exists(json_file):

        # Initialize
        calibration = {}

        # Get calibration splits
        splits = __get_calibration_splits(Xs, ys, params["splits"])

        for split in range(params["splits"]):

            # Logger
            streamer.truncate(0)
            logger.info(f"Calibration - split: {split+1}, start")
            sys.stdout.write(streamer.getvalue())

            # Get DataLoaders
            Xs_train, ys_train, Xs_val, ys_val = splits[split]
            train_loader, val_loader = __get_data_loaders(list(Xs_train),
                list(ys_train), list(Xs_val), list(ys_val),
                params["batch_size"], params["reverse_complement"],
                params["threads"])

            # Calibrate
            results = calibrate(train_loader, val_loader, seq_length,
                params["cnn_units"], params["lr"], params["epochs"],
                params["patience"])
            calibration.setdefault(split, results)

            # Logger
            streamer.truncate(0)
            logger.info(f"Calibration - split: {split+1}, finished")
            sys.stdout.write(streamer.getvalue())

        # Write
        with open(json_file, "wt") as handle:
            handle.write(json.dumps(calibration, indent=4))

    else:

        # Load
        with open(json_file) as handle:
            calibration = json.load(handle)

    # Model evaluation
    json_file = os.path.join(params["output_dir"], "evaluation.json")
    if not os.path.exists(json_file):

        # Get max. epochs
        max_epochs = __get_max_epochs(calibration, params["splits"])

        # Get Xs/ys
        Xs_test, ys_test = __get_Xs_ys(params["test_file"])

        # Get DataLoaders
        train_loader, test_loader = __get_data_loaders(list(Xs), list(ys),
            list(Xs_test), list(ys_test), params["batch_size"],
            params["reverse_complement"], params["threads"])

        # Train
        model_file = os.path.join(params["output_dir"], "model.pt")
        if not os.path.exists(model_file):

            # Logger
            streamer.truncate(0)
            logger.info(f"Training: start")
            sys.stdout.write(streamer.getvalue())

            train(train_loader, model_file, seq_length, params["cnn_units"],
                params["lr"], round(max_epochs), params["output_dir"])

            # Logger
            streamer.truncate(0)
            logger.info(f"Training: finished")
            sys.stdout.write(streamer.getvalue())

        # Evaluate
        evaluation = evaluate(test_loader, model_file, seq_length,
            params["cnn_units"])

        # Write
        with open(json_file, "wt") as handle:
            handle.write(json.dumps(evaluation, indent=4))

def calibrate(train_loader, val_loader, seq_length, cnn_units=10, lr=1e-03,
    epochs=128, patience=20):

    # Initialize
    calibration = dict({"best_val_loss": None})
    model = Model(seq_length, cnn_units).to(device)
    criterion = get_criterion(output="linear")
    optimizer = get_optimizer(model.parameters(), lr)

    # Trainer
    trainer = create_supervised_trainer(model, optimizer, criterion,
        device=device)
    @trainer.on(Events.EPOCH_COMPLETED)
    def __log(trainer):
        train_evaluator.run(train_loader)
        train_metrics = train_evaluator.state.metrics
        val_evaluator.run(val_loader)
        val_metrics = val_evaluator.state.metrics
        streamer.truncate(0)
        log_message = f"Epoch: {trainer.state.epoch}"
        train_loss = train_metrics["loss"]
        log_message += f", Train loss: {train_loss:.4f}"
        val_loss = val_metrics["loss"]
        log_message += f", Val. loss: {val_loss:.4f}"
        logger.info(log_message)
        pbar.log_message(streamer.getvalue().strip())
        # Calibration
        calibration.setdefault(trainer.state.epoch, {})
        calibration[trainer.state.epoch].setdefault("training",
            copy.deepcopy(train_metrics))
        calibration[trainer.state.epoch].setdefault("validation",
            copy.deepcopy(val_metrics))
        if calibration["best_val_loss"] is None:
            calibration["best_train_loss"] = train_loss
            calibration["best_val_loss"] = val_loss
            calibration["best_epoch"] = trainer.state.epoch
        elif val_loss < calibration["best_val_loss"]:
            calibration["best_train_loss"] = train_loss
            calibration["best_val_loss"] = val_loss
            calibration["best_epoch"] = trainer.state.epoch
        # Early stopping criteria
        if trainer.state.epoch == calibration["best_epoch"] + patience:
            streamer.truncate(0)
            logger.info("EarlyStopping: Stop training")
            pbar.log_message(streamer.getvalue().strip())
            trainer.terminate()

    # Progress bar
    GpuInfo().attach(trainer, name="gpu")
    pbar = ProgressBar(bar_format=bar_format, persist=True)
    pbar.attach(trainer, ["gpu:0 mem(%)", "gpu:0 util(%)"])

    # Evaluators
    metrics = {"loss": Loss(criterion)}
    train_evaluator = create_supervised_evaluator(model, metrics=metrics,
        device=device)
    val_evaluator = create_supervised_evaluator(model, metrics=metrics,
        device=device)

    # Train
    trainer.run(train_loader, max_epochs=epochs)

    return(calibration)

def train(train_loader, model_file, seq_length, cnn_units=10, lr=1e-03,
    epochs=128, output_dir="./"):

    # Initialize
    training = dict({"best_train_loss": None})
    model = Model(seq_length, cnn_units).to(device)
    criterion = get_criterion(output="linear")
    optimizer = get_optimizer(model.parameters(), lr)

    # Models directory
    models_dir = os.path.join(output_dir, "pts")
    if not os.path.exists(models_dir):
        os.makedirs(models_dir)

    # Trainer
    trainer = create_supervised_trainer(model, optimizer, criterion,
        device=device)
    @trainer.on(Events.EPOCH_COMPLETED)
    def __log(trainer):
        train_evaluator.run(train_loader)
        train_metrics = train_evaluator.state.metrics
        streamer.truncate(0)
        log_message = f"Epoch: {trainer.state.epoch}"
        train_loss = train_metrics["loss"]
        log_message += f", Train loss: {train_loss:.4f}"
        logger.info(log_message)
        pbar.log_message(streamer.getvalue().strip())
        # Save model
        if training["best_train_loss"] is None:
            training["best_train_loss"] = train_loss
            best_model = os.path.join(models_dir,
                "best_model_epoch=%s_loss=%s.pt" % (str(trainer.state.epoch),
                str(train_loss)))
            training.setdefault("best_model", best_model)
            torch.save(model.state_dict(), training["best_model"])
        elif train_loss < training["best_train_loss"]:
            training["best_train_loss"] = train_loss
            best_model = os.path.join(models_dir,
                "best_model_epoch=%s_loss=%s.pt" % (str(trainer.state.epoch),
                str(train_loss)))
            training["best_model"] = best_model
            torch.save(model.state_dict(), training["best_model"])

    # Progress bar
    GpuInfo().attach(trainer, name="gpu")
    pbar = ProgressBar(bar_format=bar_format, persist=True)
    pbar.attach(trainer, ["gpu:0 mem(%)", "gpu:0 util(%)"])

    # Evaluator
    metrics = {"loss": Loss(criterion)}
    train_evaluator = create_supervised_evaluator(model, metrics=metrics,
        device=device)

    # Train
    trainer.run(train_loader, max_epochs=epochs)

    # Symlink to best model
    os.symlink(training["best_model"], model_file)

def evaluate(test_loader, model_file, seq_length, cnn_units=10):

    # Initialize
    model = Model(seq_length, cnn_units)
    model.load_state_dict(torch.load(model_file))
    model.to(device)
    model.eval() # set the model in evaluation mode

    for Xs, ys in test_loader:
        Xs = Xs.to(device)
        ys = ys.numpy().flatten()
        with torch.no_grad():
            outputs = model(Xs)
        scores = outputs.data.cpu().numpy().flatten()
        print(pearsonr(scores, ys))
        print(spearmanr(scores, ys))

    exit(0)

    # # Evaluator
    # metrics = dict({"Pearson": PearsonR(), "Spearman": SpearmanR()})
    # test_evaluator = create_supervised_evaluator(model, metrics=metrics,
    #     device=device)
    # test_evaluator.run(test_loader)
    # test_metrics = test_evaluator.state.metrics
    # pearson = test_metrics["Pearson"]
    # spearman = test_metrics["Spearman"]
    # streamer.truncate(0)
    # log_message = f"Evaluation:"
    # log_message += f" Pearson: {pearson:.4f}, Spearman: {spearman:.4f}"
    # logger.info(log_message)
    # sys.stdout.write(streamer.getvalue())

    # return(test_metrics)

def __get_Xs_ys(training_file):

    # Initialize
    Xs = []
    ys = []

    # Xs / ys
    handle = __get_handle(training_file)
    for record in SeqIO.parse(handle, "fasta"):
        _, y_list = record.description.split()
        Xs.append(one_hot_encode(str(record.seq).upper()))
        ys.append([float(y) for y in y_list.split(";")])

    return(np.array(Xs), np.array(ys))

def __get_calibration_splits(Xs, ys, n_splits=3):

    # Initialize
    splits = []
    ixs = list(range(len(ys)))

    # Sort indices by value
    ixs = np.array(sorted(ixs, key=lambda x: ys[x]))

    # Predefined splits
    predefined_splits = list(islice(cycle(list(range(n_splits))), ixs.size))
    ps = PredefinedSplit(predefined_splits)
    for train_index, val_index in ps.split():
        ixs_train, ixs_val = ixs[train_index], ixs[val_index]
        Xs_train, Xs_val = Xs[ixs_train], Xs[ixs_val]
        ys_train, ys_val = ys[ixs_train], ys[ixs_val]
        splits.append((Xs_train, ys_train, Xs_val, ys_val))

    return(splits)

def __get_data_loaders(Xs_train, ys_train, Xs_val=None, ys_val=None,
    batch_size=64, reverse_complement=False, threads=1):

    # Reverse complement
    if reverse_complement:
        n = len(Xs_train)
        for i in range(n):
            encoded_seq = reverse_complement_one_hot_encoding(Xs_train[i])
            Xs_train.append(encoded_seq)
            ys_train.append(ys_train[i])

    # TensorDatasets
    train_set = TensorDataset(torch.Tensor(Xs_train),
        torch.Tensor(ys_train))
    if Xs_val is not None and ys_val is not None:
        val_set = TensorDataset(torch.Tensor(Xs_val),
            torch.Tensor(ys_val))
    else:
        val_set = None

    # DataLoaders
    kwargs = dict(batch_size=batch_size, shuffle=True,
        num_workers=threads)
    train_loader = DataLoader(train_set, **kwargs)
    if val_set is not None:
        kwargs = dict(batch_size=len(Xs_val), shuffle=True,
            num_workers=threads)
        val_loader = DataLoader(val_set, **kwargs)
    else:
        val_loader = None

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

def __get_max_epochs(calibration, n_splits=3):

    epochs = []

    for split in sorted(calibration):

        epochs.append(calibration[split]["best_epoch"])

    return(np.mean(epochs)*(float(n_splits-1)/n_splits))

if __name__ == "__main__":
    main()