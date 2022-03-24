#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import click
import copy
from io import StringIO
import numpy as np
import os
import pandas as pd
import random
import re
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder
import subprocess as sp
import sys

# Globals
scripts_dir = os.path.dirname(os.path.realpath(__file__))

# Local imports
sys.path.insert(0, os.path.abspath(os.path.join(scripts_dir, os.pardir))) 
from interpret import _get_sequences
from train import _get_seqs_labels_ids
from utils import get_file_handle

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "tsv_file", type=click.Path(exists=True, resolve_path=True),
)
@click.option(
    "-e", "--experiment-type",
    help="Experiment type.",
    type=click.Choice(["AFFI-SEQ", "CHIP-SEQ", "HT-SELEX", "PBM", "SMILE-SEQ"],
        case_sensitive=False),
    required=True,
)
@click.option(
    "-f", "--format",
    help="Output format.",
    type=click.Choice(["EXPLAINN", "FASTA"], case_sensitive=False),
    default="EXPLAINN",
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
    "-p", "--prefix",
    help="Output prefix.",
)

def main(**args):

    # Create output dir
    if not os.path.exists(args["output_dir"]):
        os.makedirs(args["output_dir"])

    # Get positive/negative sequences from ExplaiNN TSV file
    _to_sequences(args["tsv_file"], args["experiment_type"].upper(),
        args["format"].upper(), args["output_dir"], args["prefix"])

def _to_sequences(tsv_file, experiment_type, oformat="EXPLAINN",
                  output_dir="./", prefix=None):

    # Get sequences
    seqs = _get_sequences(tsv_file)

    # Get labels and sequence ids
    _, labels, ids = _get_seqs_labels_ids(tsv_file)

    # Get positives and negatives
    positives, negatives = __get_pos_neg_seqs(seqs, labels, ids,
        experiment_type)

    # Save sequences
    if oformat == "EXPLAINN":
        # Get DataFrame
        data = []
        if prefix is None:
            tsv_file = os.path.join(output_dir, "sequences.tsv.gz")
        else:
            tsv_file = os.path.join(output_dir, f"{prefix}.tsv.gz")
        for record in positives + negatives:
            data.append([str(record.id), str(record.seq),
                float(record.description)])
        df = pd.DataFrame(data)
        df.to_csv(tsv_file, sep="\t", header=False, index=False,
            compression="gzip")
    else:
        # Positive sequences
        if prefix is None:
            fasta_file = os.path.join(output_dir, "pos.fa")
        else:
            fasta_file = os.path.join(output_dir, f"{prefix}.pos.fa")
        handle = get_file_handle(fasta_file, "wt")
        SeqIO.write(positives, handle, "fasta")
        # Negative sequences
        if prefix is None:
            fasta_file = os.path.join(output_dir, "neg.fa")
        else:
            fasta_file = os.path.join(output_dir, f"{prefix}.neg.fa")
        handle = get_file_handle(fasta_file, "wt")
        SeqIO.write(negatives, handle, "fasta")

def __get_pos_neg_seqs(seqs, labels, ids, experiment_type):

    # Initialize
    n = .05
    positives = []
    negatives = []

    # Affi-seq/HT-SELEX
    if experiment_type == "AFFI-SEQ" or experiment_type == "HT-SELEX":
        # Positives = 1s in nth class
        pos_idxs = np.where(labels.T[labels.shape[1] - 1] == 1.)[0]
        # Negatives = 1s in class 0
        neg_idxs = np.where(labels.T[0] == 1.)[0]
    # ChIP-seq/SMiLE-seq
    elif experiment_type == "CHIP-SEQ" or experiment_type == "SMILE-SEQ":
        # Positives = 1s
        pos_idxs = np.where(labels.flatten() == 1.)[0]
        # Negatives = 0s
        neg_idxs = np.where(labels.flatten() == 0.)[0]
    # PBM
    elif experiment_type == "PBM":       
        # Positives = top `n` percentile probes
        pos_idxs = np.argsort(-labels.flatten())[:int(max(labels.shape)*n)]
        # Negatives = bottom `n` percentile probes
        neg_idxs = np.argsort(labels.flatten())[:int(max(labels.shape)*n)]

    # Positives
    for s, i in zip(seqs[pos_idxs], ids[pos_idxs]):
        positives.append(SeqRecord(Seq(s), id=str(i), description="1."))

    # Negatives
    for s, i in zip(seqs[neg_idxs], ids[neg_idxs]):
        negatives.append(SeqRecord(Seq(s), id=str(i), description="0."))

    return(positives, negatives)

if __name__ == "__main__":
    main()