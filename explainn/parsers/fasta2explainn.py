#!/usr/bin/env python

from Bio import SeqIO
import click
import copy
from io import StringIO
import os
import pandas as pd
import random
import re
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder
import subprocess as sp

from utils import get_file_handle

# Globals
scripts_dir = os.path.dirname(os.path.realpath(__file__))

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "fasta_file", type=click.Path(exists=True, resolve_path=True),
)
@click.option(
    "-k", "--kmer",
    help="K-mer size to shuffle.",
    type=int,
    default=2,
    show_default=True
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
    "-s", "--shuffle-lowercase",
    help="Shuffle lowercase letters.",
    is_flag=True
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

    # Get TSV files for CAMs
    _to_CAM(args["fasta_file"], args["kmer"], args["output_dir"],
        args["prefix"], args["shuffle_lowercase"], args["test"])

def _to_CAM(fasta_file, kmer=2, output_dir="./", prefix=None, shuffle_lc=False,
    create_test=False):

    # Initialize
    data = []
    regexp = re.compile(r"[a-z]+")

    # Positive sequences
    fh = get_file_handle(fasta_file, "rt")
    for record in SeqIO.parse(fh, "fasta"):
        s = str(record.seq)
        if shuffle_lc:
            # 1) extract blocks of lowercase letters;
            # 2) shuffle the letters; and
            # 3) put the shuffled letters back
            l = list(s)
            for m in re.finditer(regexp, s):
                sublist = l[m.start():m.end()]
                random.shuffle(sublist)
                l[m.start():m.end()] = copy.copy(sublist)
            s = "".join(l)
        data.append([record.id, s.upper(), 1.])
    fh.close()

    # Negative sequences
    cmd = "biasaway k -f %s -k %s -e 1" % (fasta_file, kmer)
    p = sp.run([cmd], shell=True, stdout=sp.PIPE, stderr=sp.DEVNULL)
    fh = StringIO(p.stdout.decode())
    for record in SeqIO.parse(fh, "fasta"):
        data.append([record.id, str(record.seq).upper(), 0.])
    fh.close()

    # Get DataFrame
    df = pd.DataFrame(data)
    df = df.groupby(1).max().reset_index()
    df = df.reindex(sorted(df.columns), axis=1)

    # Get data splits
    train, validation = train_test_split(df, test_size=0.2)
    if create_test:
        validation, test = train_test_split(validation, test_size=0.5)

    # Save sequences
    if prefix is None:
        tsv_file = os.path.join(output_dir, "train.tsv.gz")
    else:
        tsv_file = os.path.join(output_dir, f"{prefix}.train.tsv.gz")
    train.to_csv(tsv_file, sep="\t", header=False, index=False,
        compression="gzip")
    if prefix is None:
        tsv_file = os.path.join(output_dir, "validation.tsv.gz")
    else:
        tsv_file = os.path.join(output_dir, f"{prefix}.validation.tsv.gz")
    validation.to_csv(tsv_file, sep="\t", header=False, index=False,
        compression="gzip")
    if create_test:
        if prefix is None:
            tsv_file = os.path.join(output_dir, "test.tsv.gz")
        else:
            tsv_file = os.path.join(output_dir, f"{prefix}.test.tsv.gz")
        test.to_csv(tsv_file, sep="\t", header=False, index=False,
            compression="gzip")

if __name__ == "__main__":
    main()