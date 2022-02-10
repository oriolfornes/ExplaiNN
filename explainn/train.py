#!/usr/bin/env python

import click
from click_option_group import optgroup
import math
import numpy as np
import os
import pandas as pd
import pickle
import shutil
import time
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

# Local imports
from architectures import ExplaiNN, get_loss, get_metrics, get_optimizer
from selene import Trainer
from sequence import one_hot_encode_many, rc_one_hot_encoding_many
from utils import get_file_handle

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
    default=100,
    show_default=True,
)
@click.option(
    "-c", "--cpu-threads",
    help="Number of CPU threads to use.",
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
    help="Reverse complement sequences.",
    is_flag=True,
)
@click.option(
    "-t", "--time-me",
    help="Return the program's running execution time in seconds.",
    is_flag=True,
)
# @click.option(
#     "-v", "--val-samples",
#     help="Number of validation samples to use. " + \
#          "By default, use all samples.",
#     type=int,
# )
@optgroup.group("ExplaiNN")
# @optgroup.option(
#     "--clamp-weights",
#     help="Constrain weights to being non-negative.",
#     is_flag=True,
# )
@optgroup.option(
    "--cnn-units",
    help="Number of CNN units to use.",
    type=int,
    default=100,
    show_default=True,
)
@optgroup.option(
    "--kernel-size",
    help="Kernel size (i.e. max. length for learning motifs).",
    type=int,
    default=19,
    show_default=True,
)
# @optgroup.option(
#     "--no-padding",
#     help="Do not apply padding to input sequences. " + \
#          "By default, it applies a padding equal to the kernel size.",
#     is_flag=True,
# )
@optgroup.group("Optimizer")
@optgroup.option(
    "--lr",
    help="Learning rate for the `torch.optim.Adam` optimizer.",
    type=float,
    default=0.003,
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
    default=100,
    show_default=True,
)
@optgroup.option(
    "--patience",
    help="If no improvement, wait `n` epochs before stopping training.",
    type=int,
    default=10,
    show_default=True,
)
@optgroup.group("Transfer Learning")
@optgroup.option(
    "--filter-weights",
    help="Filter weights from the pre-training step.",
    type=click.Path(exists=True, resolve_path=True),
    multiple=True,
)
@optgroup.option(
    "--freeze",
    help="Do not update the pre-training filter weights during the fine-tuning step.",
    is_flag=True,
)

def main(**args):

    # Start execution
    start_time = time.time()

    # Initialize
    if not os.path.exists(args["output_dir"]):
        os.makedirs(args["output_dir"])

    ##############
    # Load Data  #
    ##############

    # Get data
    seqs_train, labels_train, _ = _get_seqs_labels_ids(args["training_file"],
        args["debugging"], args["rev_complement"])
    seqs_val, labels_val, _ = _get_seqs_labels_ids(args["validation_file"],
        args["debugging"], args["rev_complement"])

    # Get DataLoaders
    data_loader_train = _get_data_loader(seqs_train, labels_train,
        args["batch_size"], shuffle=True)
    data_loader_val = _get_data_loader(seqs_val, labels_val,
        args["batch_size"], shuffle=True)
    data_loaders = dict({"train": data_loader_train,
        "validation": data_loader_val})

    ##########################
    # (if) Transfer Learning #
    ##########################

    # Get pre-training filter IDs and weights
    filter_ids = []
    filter_weights = []
    for pickle_file in args["filter_weights"]:
        fh = get_file_handle(pickle_file, "rb")
        w = pickle.load(fh)
        fh.close()
        for k, v in w.items():
            if isinstance(k, str):
                filter_ids.append(k)
            elif isinstance(k, tuple):
                filter_ids.append(";".join(k))
            # Save filter weights as tensors
            filter_weights.append(torch.tensor(v.T))

    # Save IDs of filters from the pre-training step
    # as they will be used in the fine-tunning step
    if len(filter_ids) > 0 and len(filter_ids) == len(filter_weights):
        txt_file = os.path.join(args["output_dir"],
            "filter-ids-from-pre-training-step.txt")
        fh = get_file_handle(txt_file, "wt")
        fh.write("\n".join(filter_ids))
        fh.close()

    ##############
    # Train      #
    ##############

    # Initialize
    sequence_length = seqs_train[0].shape[1]
    n_features = labels_train[0].shape[0]
    if np.unique(labels_train[:, 0]).size == 2:
        input_data = "binary"
    else:
        input_data = "linear"
    steps_per_epoch = math.ceil(
        len(data_loader_train.dataset) / data_loader_train.batch_size
    )

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
                steps_per_epoch, args["cnn_units"], args["kernel_size"],
                args["lr"], 5, args["patience"], args["cpu_threads"], d,
                filter_weights, args["freeze"])

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
        steps_per_epoch, args["cnn_units"], args["kernel_size"],
        args["lr"], args["max_epochs"], args["patience"], args["cpu_threads"],
        args["output_dir"], filter_weights, args["freeze"], best_init_model)

    # Finish execution
    if args["time_me"]:
        seconds = format(time.time() - start_time, ".2f")
        txt_file = os.path.join(args["output_dir"],
            f"time-{os.path.basename(__file__)}.txt")
        fh = get_file_handle(txt_file, "wt")
        fh.write(f"{seconds} seconds")
        fh.close()

def _get_seqs_labels_ids(tsv_file, debugging=False, reverse_complement=False):

    # Seqs / labels / ids
    df = pd.read_csv(tsv_file, sep="\t", header=None)
    ids = df.pop(0).values
    seqs = one_hot_encode_many(df.pop(1).values)
    labels = df.values

    # Reverse complement
    if reverse_complement:
        seqs = np.append(seqs, rc_one_hot_encoding_many(seqs), axis=0)
        labels = np.append(labels, labels, axis=0)
        ids = np.append(ids, ids, axis=0)

    # Return 1,000 sequences
    if debugging:
        return(seqs[:1000], labels[:1000], ids[:1000])

    return(seqs, labels, ids)

def _get_data_loader(seqs, labels, batch_size=100, shuffle=False):

    # TensorDatasets
    dataset = TensorDataset(torch.Tensor(seqs), torch.Tensor(labels))

    # Avoid Error: Expected more than 1 value per channel when training
    batch_size = __avoid_expect_more_than_1_value_per_channel(len(dataset),
        batch_size)
        
    return(DataLoader(dataset, batch_size, shuffle=shuffle))

def __avoid_expect_more_than_1_value_per_channel(n, batch_size):

    if n % batch_size == 1:
        return(__avoid_expect_more_than_1_value_per_channel(n, batch_size - 1))

    return(batch_size)

def _train(sequence_length, n_features, data_loaders, input_data,
           steps_per_epoch, cnn_units=100, kernel_size=19, lr=0.003,
           max_epochs=100, patience=10, cpu_threads=1, output_dir="./",
           filter_weights=[], freeze=False, checkpoint_resume=None):

    # Initialize
    freeze_top_n_filters = 0
    exp_model = ExplaiNN(cnn_units, kernel_size, sequence_length, n_features)
    loss_criterion = get_loss(input_data=input_data)
    metrics = get_metrics(input_data=input_data)
    optimizer = get_optimizer(exp_model.parameters(), lr)

    # Transfer learning
    if len(filter_weights) > 0:
        finetuning_filter_weights = exp_model.linears[0].weight.data
        for i in range(finetuning_filter_weights.shape[0]):
            finetuning_filter_weights[i] = filter_weights[i]
            if freeze:
                freeze_top_n_filters += 1
        finetuning_filter_weights.requires_grad = True
        exp_model.linears[0].weight = nn.Parameter(finetuning_filter_weights)

    # Train
    trainer = Trainer(
        exp_model,
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
        freeze_top_n_filters=freeze_top_n_filters
    )
    trainer.train_and_validate()

if __name__ == "__main__":
    main()