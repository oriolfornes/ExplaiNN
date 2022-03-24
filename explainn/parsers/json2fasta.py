#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import click
import json
import gzip
import os
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "json_file", type=click.Path(exists=True, resolve_path=True)
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
@click.option(
    "-t", "--test",
    help="Create a test set.",
    is_flag=True,
)

def main(**args):

    # Create output dir
    if not os.path.exists(args["output_dir"]):
        os.makedirs(args["output_dir"])

    # Initialize
    train = []
    validation = []
    test = []

    # Load JSON
    with open(args["json_file"]) as handle:
        sequences = json.load(handle)
    sequences.pop(0)

    # Ys
    enc = OneHotEncoder()
    arr = [[i] for i in range(len(sequences[0]) - 1)]
    enc.fit(arr)
    ys = [";".join(map(str, i)) for i in enc.transform(arr).toarray()]

    # Get records
    records = [[] for _ in range(len(ys))]
    for i in range(len(sequences)):
        for j in range(1, len(sequences[i])):
            record = SeqRecord(Seq(sequences[i][j][1]), id=sequences[i][j][0],
                description=ys[j-1])
            records[j-1].append(record)

    # Get data splits
    for i in range(len(records)):

        # Split records 80/20
        t, rest = train_test_split(records[i], test_size=.2, random_state=1)
        train.append(t)

        if args["test"]:
            # Split the 20 split 50/50
            v, t = train_test_split(rest, test_size=.5, random_state=1)
            validation.append(v)
            test.append(t)
        else:
            validation.append(rest)

    # Save sequences
    if args["prefix"] is None:
        f = os.path.join(args["output_dir"], "train.fa.gz")
    else:
        f = os.path.join(args["output_dir"], "%s.train.fa.gz" % args["prefix"])
    with gzip.open(f, "wt") as handle:
        sequences = [j for i in map(list, zip(*train)) for j in i]
        SeqIO.write(sequences, handle, "fasta")
    if args["prefix"] is None:
        f = os.path.join(args["output_dir"], "validation.fa.gz")
    else:
        f = os.path.join(args["output_dir"], "%s.validation.fa.gz" % \
            args["prefix"])
    with gzip.open(f, "wt") as handle:
        sequences = [j for i in  map(list, zip(*validation)) for j in i]
        SeqIO.write(sequences, handle, "fasta")
    if args["test"]:
        if args["prefix"] is None:
            f = os.path.join(args["output_dir"], "test.fa.gz")
        else:
            f = os.path.join(args["output_dir"], "%s.test.fa.gz" % \
                args["prefix"])
        with gzip.open(f, "wt") as handle:
            sequences = [j for i in  map(list, zip(*test)) for j in i]
            SeqIO.write(sequences, handle, "fasta")

if __name__ == "__main__":
    main()