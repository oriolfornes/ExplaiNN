#!/usr/bin/env python

from Bio import SeqIO
import click
from click_option_group import optgroup
import gzip
import math
import numpy as np
import os
import random
import torch
from torch.utils.data import DataLoader, TensorDataset

# Local imports
from utils.architectures import DanQ, NonStrandSpecific, get_loss_criterion, \
                                get_metrics, get_optimizer
from utils.selene import Trainer
from utils.sequence import one_hot_encode, rc_one_hot_encoding

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "training_file",
    type=click.Path(exists=True, resolve_path=True),
)
@click.argument(
    "validation_file",
    type=click.Path(exists=True, resolve_path=True),
    nargs=-1,
)
@click.option(
    "-b", "--batch-size",
    help="Batch size.",
    type=int,
    default=2**6,
    show_default=True,
)
@click.option(
    "-d", "--debugging",
    help="Debugging mode.",
    is_flag=True,
)
@click.option(
    "-o", "--output-dir",
    help="Output directory.",
    type=click.Path(resolve_path=True),
    default="./",
    show_default=True,
)
@click.option(
    "-t", "--threads",
    help="Number of threads to use.",
    type=int,
    default=1,
    show_default=True,
)
@click.option(
    "-r", "--rev-complement",
    help="Reverse complement training sequences.",
    is_flag=1
)
@click.option(
    "-v", "--val-samples",
    help="Number of validation samples to use.",
    type=int,
    default=640000,
    show_default=True,
)
@optgroup.group("DanQ")
@optgroup.option(
    "--input-data",
    help="Type of input data.",
    type=click.Choice(["binary", "linear"]),
    required=True,
)
@optgroup.option(
    "--strand-specific",
    help="By default, score each input sequence and its " + \
         "reverse-complement and return the mean.",
    is_flag=True,
)
@optgroup.group("Optimizer")
@optgroup.option(
    "--lr",
    help="Learning rate for the `torch.optim.Adam` optimizer.",
    type=float,
    default=1e-03,
    show_default=True,
)
@optgroup.group("Selene")
@optgroup.option(
    "--checkpoint-resume",
    help="Resume training from model file.",
    type=click.Path(exists=True, resolve_path=True),
)
# @optgroup.option(
#     "--report-steps",
#     help="Report stats `n` steps.",
#     type=int,
#     default=1000,
#     show_default=True,
# )
@optgroup.option(
    "--max-epochs",
    help="Number of epochs to train.",
    type=int,
    default=128,
    show_default=True,
)
@optgroup.option(
    "--patience",
    help="If no improvement, wait `n` epochs before stopping training.",
    type=int,
    default=16,
    show_default=True,
)

def main(**params):

    # Create output dir
    if not os.path.exists(params["output_dir"]):
        os.makedirs(params["output_dir"])

    # Get DataLoaders
    Xs_train, ys_train = _get_Xs_ys(params["training_file"],
        params["debugging"])
    Xs_val, ys_val = [], []
    for validation_file in params["validation_file"]:
        Xs, ys = _get_Xs_ys(validation_file, params["debugging"])
        Xs_val.extend(Xs); ys_val.extend(ys)
    train_loader, val_loader = _get_data_loaders(list(Xs_train),
        list(ys_train), list(Xs_val), list(ys_val), params["batch_size"],
        params["rev_complement"], params["threads"], params["val_samples"])
    data_loaders = dict({"train": train_loader, "validation": val_loader})

    # Get steps per epoch
    steps_per_epoch = math.ceil(
        len(train_loader.dataset) / float(params["batch_size"])
    )

    # Train
    model = DanQ(max(Xs_train[0].shape))
    if not params["strand_specific"]:
        model = NonStrandSpecific(model)
    loss_criterion = get_loss_criterion(input_data=params["input_data"])
    metrics = get_metrics(input_data=params["input_data"])
    optimizer = get_optimizer(model.parameters(), params["lr"])
    trainer = Trainer(
        model,
        data_loaders,
        loss_criterion,
        metrics,
        optimizer,
        max_steps=steps_per_epoch*params["max_epochs"],
        patience=steps_per_epoch*params["patience"],
        report_stats_every_n_steps=steps_per_epoch,
        output_dir=params["output_dir"],
        cpu_n_threads=params["threads"],
        use_cuda=torch.cuda.is_available(),
        checkpoint_resume=params["checkpoint_resume"],
    )
    trainer.train_and_validate()

def _get_Xs_ys(fasta_file, debugging=False):

    # Initialize
    Xs = []
    ys = []

    # Xs / ys
    handle = __get_handle(fasta_file)
    for record in SeqIO.parse(handle, "fasta"):
        _, y_list = record.description.split()
        Xs.append(one_hot_encode(str(record.seq).upper()))
        ys.append([float(y) for y in y_list.split(";")])

    # Return 1,000 sequences
    if debugging:
        return(np.array(Xs)[:1000], np.array(ys)[:1000])

    return(np.array(Xs), np.array(ys))

def _get_data_loaders(Xs_train, ys_train, Xs_val=None, ys_val=None,
    batch_size=2**6, reverse_complement=False, threads=1, val_samples=640000):

    # Initialize
    ix = 0
    Xs_val_nested = []
    ys_val_nested = []
    random.seed(123)

    # Reverse complement
    if reverse_complement:
        n = len(Xs_train)
        for i in range(n):
            encoded_seq = rc_one_hot_encoding(Xs_train[i])
            Xs_train.append(encoded_seq)
            ys_train.append(ys_train[i])

    # Shuffle validation data (i.e. if multiple were provided)
    while ix < len(Xs_val):
        Xs_val_nested.append(Xs_val[ix:ix+2])
        ys_val_nested.append(ys_val[ix:ix+2])
        ix += 2
    z = list(zip(Xs_val_nested, ys_val_nested))
    random.shuffle(z)
    Xs_val = [j for i in z for j in i[0]]
    ys_val = [j for i in z for j in i[1]]

    # TensorDatasets
    train_set = TensorDataset(torch.Tensor(Xs_train), torch.Tensor(ys_train))
    val_set = TensorDataset(torch.Tensor(Xs_val[:val_samples]),
        torch.Tensor(ys_val[:val_samples]))

    # DataLoaders
    kwargs = dict(batch_size=batch_size, num_workers=threads)
    train_loader = DataLoader(train_set, **kwargs)
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

if __name__ == "__main__":
    main()