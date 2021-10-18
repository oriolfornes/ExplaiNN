#!/usr/bin/env python

from Bio import SeqIO
import click
from click_option_group import optgroup
import gzip
import math
import numpy as np
import pandas as pd
import os
import random
import shutil
import torch
from torch.utils.data import DataLoader, TensorDataset

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
    # nargs=-1,
)
@click.option(
    "-b", "--batch-size",
    help="Batch size.",
    type=int,
    default=2**6,
    show_default=True,
)
@click.option(
    "-i", "--initialize",
    help="Number of models to initialize.",
    type=int,
    default=1,
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
    "-r", "--rev-complement",
    help="Reverse complement training sequences.",
    is_flag=True,
)
@click.option(
    "-t", "--threads",
    help="Number of CPU threads to use.",
    type=int,
    default=1,
    show_default=True,
)
@click.option(
    "-v", "--val-samples",
    help="Number of validation samples to use. " + \
         "By default, it uses all samples.",
    type=int,
)
@optgroup.group("CAM")
@optgroup.option(
    "--clamp-weights",
    help="Constrain weights to being non-negative.",
    is_flag=True,
)
@optgroup.option(
    "--cnn-units",
    help="Number of CNN units to use.",
    type=int,
    default=2**4,
    show_default=True,
)
@optgroup.option(
    "--kernel-size",
    help="Kernel size (i.e. motif width).",
    type=int,
    default=19,
    show_default=True,
)
@optgroup.option(
    "--no-padding",
    help="Do not apply padding to input sequences. " + \
         "By default, it applies a padding equal to the kernel size.",
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

    ##############
    # Load Data  #
    ##############

    # Get data
    Xs_train, ys_train = _get_Xs_ys(params["training_file"],
        params["debugging"], params["rev_complement"])
    Xs_val, ys_val = _get_Xs_ys(params["training_file"],
        params["debugging"], params["rev_complement"])

    # Get DataLoaders
    train_loader, val_loader = _get_data_loaders(list(Xs_train),
        list(ys_train), list(Xs_val), list(ys_val), params["batch_size"],
        # params["threads"], params["val_samples"])
        params["val_samples"])
    data_loaders = dict({"train": train_loader, "validation": val_loader})

    ##############
    # Train      #
    ############## 

    # Initialize
    sequence_length = max(Xs_train[0].shape)
    n_features = max(ys_train[0].shape)
    if np.unique(ys_train[:, 0]).size == 2:
        input_data = "binary"
    else:
        input_data = "linear"
    steps_per_epoch = math.ceil(
        len(train_loader.dataset) / float(params["batch_size"])
    )

    # For each initialization...
    best_init_loss = None
    best_init_model = None
    for i in range(params["initialize"]):

        # Skip if already done
        init_dir = os.path.join(params["output_dir"], f"init.{i}")
        if not os.path.exists(init_dir):

            # Create output dir
            os.makedirs(init_dir)

            # Train
            _train(sequence_length, n_features, data_loaders, input_data,
                steps_per_epoch, params["cnn_units"], params["kernel_size"],
                params["clamp_weights"], params["no_padding"], params["lr"],
                int(params["patience"] * .5), params["patience"],
                params["threads"], init_dir)

        # Get best initialized
        df = pd.read_csv(os.path.join(init_dir, "validation.txt"), sep="\t")
        init_loss = df.loss.min()
        init_model = os.path.join(init_dir, "best_model.pth.tar")
        if best_init_model is None:
            best_init_loss = init_loss
            best_init_model = init_model
        elif init_loss < best_init_loss:
            best_init_loss = init_loss
            best_init_model = init_model

    # Train
    shutil.copy(best_init_model, params["output_dir"])
    _train(sequence_length, n_features, data_loaders, input_data,
        steps_per_epoch, params["cnn_units"], params["kernel_size"],
        params["clamp_weights"], params["no_padding"], params["lr"],
        params["max_epochs"], params["patience"], params["threads"],
        params["output_dir"], best_init_model)

def _get_Xs_ys(fasta_file, debugging=False, reverse_complement=False):

    from sequence import one_hot_encode, rc_one_hot_encoding

    # Initialize
    Xs = []
    ys = []

    # Xs / ys
    handle = __get_handle(fasta_file)
    for record in SeqIO.parse(handle, "fasta"):
        _, y_list = record.description.split()
        Xs.append(one_hot_encode(str(record.seq).upper()))
        ys.append([float(y) for y in y_list.split(";")])
    handle.close()

    # Reverse complement
    if reverse_complement:
        n = len(Xs)
        for i in range(n):
            Xs.append(rc_one_hot_encoding(Xs[i]))
            ys.append(ys[i])

    # Return 1,000 sequences
    if debugging:
        return(np.array(Xs, dtype=object)[:1000],
            np.array(ys, dtype=object)[:1000])

    return(np.array(Xs, dtype=object), np.array(ys, dtype=object))

def _get_data_loaders(Xs_train, ys_train, Xs_val=None, ys_val=None,
    # batch_size=2**6, threads=1, val_samples=None):
    batch_size=2**6, val_samples=None):

    # Initialize
    ix = 0
    Xs_val_nested = []
    ys_val_nested = []
    random.seed(123)

    # Shuffle validation data (i.e. if multiple validation sets were provided)
    if Xs_val is not None and ys_val is not None:
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
    if Xs_val is not None and ys_val is not None:
        if val_samples is None:
            val_set = TensorDataset(torch.Tensor(Xs_val[:val_samples]),
                torch.Tensor(ys_val[:val_samples]))
        else:
            val_set = TensorDataset(torch.Tensor(Xs_val), torch.Tensor(ys_val))

    # DataLoaders
    # kwargs = dict(batch_size=batch_size, num_workers=threads)
    kwargs = dict(batch_size=batch_size)
    train_loader = DataLoader(train_set, **kwargs)
    if Xs_val is not None and ys_val is not None:
        val_loader = DataLoader(val_set, **kwargs)

    if Xs_val is not None and ys_val is not None:
        return(train_loader, val_loader)
    else:
        return(train_loader)

def __get_handle(file_name, mode="rt"):
    if file_name.endswith("gz"):
        handle = gzip.open(file_name, mode)
    else:
        handle = open(file_name, mode)

    return(handle)

def _train(sequence_length, n_features, data_loaders, input_data,
           steps_per_epoch, cnn_units=8, kernel_size=19, clamp_weights=False,
           no_padding=False, lr=1e-03, max_epochs=128, patience=16, threads=1,
           output_dir="./", checkpoint_resume=None):

    from architectures import CAM, get_loss, get_metrics, get_optimizer
    from selene import Trainer

    # Train
    model = CAM(cnn_units, kernel_size, sequence_length, n_features,
        clamp_weights=clamp_weights, no_padding=no_padding)
    loss_criterion = get_loss(input_data=input_data)
    metrics = get_metrics(input_data=input_data)
    optimizer = get_optimizer(model.parameters(), lr)
    trainer = Trainer(
        model,
        data_loaders,
        loss_criterion,
        metrics,
        optimizer,
        max_steps=steps_per_epoch*max_epochs,
        patience=steps_per_epoch*patience,
        report_stats_every_n_steps=steps_per_epoch,
        output_dir=output_dir,
        cpu_n_threads=threads,
        use_cuda=torch.cuda.is_available(),
        checkpoint_resume=checkpoint_resume,
    )
    trainer.train_and_validate()

if __name__ == "__main__":
    main()