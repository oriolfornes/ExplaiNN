#!/usr/bin/env python

import click
from click_option_group import optgroup
import math
import numpy as np
import pandas as pd
import os
import shutil
import torch
from torch.utils.data import DataLoader, TensorDataset

# Local imports
from architectures import DanQ, get_loss, get_metrics, get_optimizer
from selene import Trainer
from sequence import one_hot_encode_many, rc_one_hot_encoding_many
from train import _get_Xs_ys_ids, _get_data_loader

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
    "-i", "--initialize",
    help="Number of models to initialize.",
    type=int,
    default=1,
    show_default=True,
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
@optgroup.group("Optimizer")
@optgroup.option(
    "--lr",
    help="Learning rate for the `torch.optim.Adam` optimizer.",
    type=float,
    default=1e-04,
    show_default=True,
)
@optgroup.group("Selene")
@optgroup.option(
    "--checkpoint-resume",
    help="Resume training from model file.",
    type=click.Path(exists=True, resolve_path=True),
)
@optgroup.option(
    "--cpu-threads",
    help="Number of CPU threads to use.",
    type=int,
    default=1,
    show_default=True,
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

def main(**args):

    ##############
    # Load Data  #
    ##############

    # Get data
    Xs_train, ys_train, _ = _get_Xs_ys_ids(args["training_file"],
        args["debugging"], args["rev_complement"])
    Xs_val, ys_val, _ = _get_Xs_ys_ids(args["validation_file"],
        args["debugging"], args["rev_complement"])

    # Get DataLoaders
    dl_train = _get_data_loader(Xs_train, ys_train, args["batch_size"], True)
    dl_val = _get_data_loader(Xs_val, ys_val, args["batch_size"], True)
    data_loaders = dict({"train": dl_train, "validation": dl_val})

    ##############
    # Train      #
    ##############

    # Initialize
    sequence_length = Xs_train[0].shape[1]
    n_features = ys_train[0].shape[0]
    if np.unique(ys_train[:, 0]).size == 2:
        input_data = "binary"
    else:
        input_data = "linear"
    steps_per_epoch = math.ceil(len(dl_train.dataset) / dl_train.batch_size)

    # For each initialization...
    best_init_loss = None
    best_init_model = None
    for i in range(args["initialize"]):

        # Skip if already done
        d = os.path.join(args["output_dir"], f"init.{i}")
        if not os.path.isdir(d):

            # Create output dir
            os.makedirs(d)

            # Train for a few epoch
            _train(sequence_length, n_features, data_loaders, input_data,
                steps_per_epoch, args["lr"], 5, args["patience"],
                args["cpu_threads"], d)

        # Get best initialized
        df = pd.read_csv(os.path.join(d, "validation.txt"), sep="\t")
        init_loss = df.loss.min()
        init_model = os.path.join(d, "best_model.pth.tar")
        if best_init_model is None:
            best_init_loss = init_loss
            best_init_model = init_model
        elif init_loss < best_init_loss:
            best_init_loss = init_loss
            best_init_model = init_model

    # Train
    shutil.copy(best_init_model, args["output_dir"])
    _train(sequence_length, n_features, data_loaders, input_data,
        steps_per_epoch, args["lr"], args["max_epochs"], args["patience"],
        args["cpu_threads"], args["output_dir"], best_init_model)

def _train(sequence_length, n_features, data_loaders, input_data,
           steps_per_epoch, lr=1e-03, max_epochs=128, patience=16,
           cpu_threads=1, output_dir="./", checkpoint_resume=None):

    # Train
    model = DanQ(sequence_length, n_features)
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
        cpu_n_threads=cpu_threads,
        use_cuda=torch.cuda.is_available(),
        checkpoint_resume=checkpoint_resume,
    )
    trainer.train_and_validate()

if __name__ == "__main__":
    main()