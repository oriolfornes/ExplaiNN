#!/usr/bin/env python

from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import click
import gzip
import numpy as np
import os
import pandas as pd
import torch
from tqdm import tqdm
bar_format = "{percentage:3.0f}%|{bar:20}{r_bar}"

# Local imports
from architectures import DanQ, get_metrics
from jaspar import get_figure, reformat_motif
from sequence import rc_many
from train import _get_Xs_ys_ids, _get_data_loader

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "model_file",
    type=click.Path(exists=True, resolve_path=True),
)
@click.argument(
    "test_file",
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

def main(**args):

    ##############
    # Load Data  #
    ##############

    # Get data
    Xs, ys, _ = _get_Xs_ys_ids(args["test_file"], args["debugging"],
        args["rev_complement"])

    # Get DataLoader
    data_loader = _get_data_loader(Xs, ys, args["batch_size"])

    # Load model
    model = _load_model(args["model_file"])

    ##############
    # Test       #
    ############## 

    # Initialize
    if np.unique(ys[:, 0]).size == 2:
        input_data = "binary"
    else:
        input_data = "linear"

    # Create output dirs
    if not os.path.isdir(args["output_dir"]):
        os.makedirs(args["output_dir"])

    # Get labels and predictions
    labels, predictions = _get_labs_preds(model, data_loader)

    # Get performance metrics
    metrics = get_metrics(input_data=input_data)
    tsv_file = os.path.join(args["output_dir"], "performance-metrics.tsv")
    if not os.path.exists(tsv_file):
        data = []
        for m in metrics:
            p = _get_performances(labels, predictions, input_data, metrics[m],
                args["rev_complement"])
            data.append([m] + p)
        column_names = ["metric", "global"] + list(range(ys.shape[1]))
        df = pd.DataFrame(data, columns=column_names)
        df.to_csv(tsv_file, sep="\t", index=False)

def _load_model(model_file):

    # Initialize
    device = "cuda" if torch.cuda.is_available() else "cpu"

    # Get model
    selene_dict = torch.load(model_file)
    model = DanQ(
        selene_dict["options"]["sequence_length"],
        selene_dict["options"]["n_features"],
        selene_dict["options"]["weights_file"],
    )
    model.load_state_dict(selene_dict["state_dict"])
    model.to(device)

    return(model)

def _get_labs_preds(model, data_loader):

    # Counter
    idx = 0
    x = len(data_loader.dataset)
    n_features = model._options["n_features"]
    n_channels = model._DanQ__channels_after_bilstm
    labels = torch.zeros((x, n_features), dtype=torch.float32)
    predictions = torch.zeros((x, n_features), dtype=torch.float32)
    device = "cuda" if torch.cuda.is_available() else "cpu"

    with torch.no_grad():
        for seqs, labs in tqdm(iter(data_loader), total=len(data_loader),
                bar_format=bar_format):

            # Prepare input
            seqs = seqs.to(device)
            labs = labs.to(device)

            # Get labels
            labels[idx:idx+seqs.shape[0]] = labs.cpu()

            # Get predictions
            outs = model.Conv1(seqs)
            outs = torch.nn.functional.relu(outs)
            outs = model.Maxpool(outs)
            outs = model.Drop1(outs)
            outs , _ = model.BiLSTM(torch.transpose(outs, 1, 2))
            outs = model.Linear1(outs.contiguous().view(-1, n_channels*640))
            preds = model.Linear2(outs)
            predictions[idx:idx+seqs.shape[0]] = preds.cpu()

            # Index increase
            idx += seqs.shape[0]

    return(labels.numpy(), predictions.numpy())

def _get_performances(labels, predictions, input_data, metric,
                      rev_complement=False):

    # Initialize
    performances = []

    if rev_complement:
        fwd = __get_fwd_rev(predictions, "fwd")
        rev = __get_fwd_rev(predictions, "rev")
        p = np.empty(fwd.shape)
        ys = __get_fwd_rev(labels, "fwd")
        # Average predictions from forward and reverse strands
        for i in range(p.shape[1]):
            p[:, i] = np.mean([fwd[:, i], rev[:, i]], axis=0)
            if input_data == "binary":
                p[:, i] = torch.sigmoid(torch.from_numpy(p[:, i])).numpy()
    else:
        if input_data == "binary":
            p = torch.sigmoid(torch.from_numpy(predictions)).numpy()
        else:
            p = predictions
        ys = labels

    # For each class...
    performances.append(metric(ys, p))
    for i in range(ys.shape[1]):
        performances.append(metric(ys[:, i], p[:, i]))

    return(performances)

def __get_fwd_rev(arr, strand):

    if strand == "fwd" or strand == "+":
        return(arr[:len(arr)//2])
    elif strand == "rev" or strand == "-":
        return(arr[len(arr)//2:])

if __name__ == "__main__":
    main()