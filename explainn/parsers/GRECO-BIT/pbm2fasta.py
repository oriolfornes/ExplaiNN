#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import click
from functools import partial
import gzip
from multiprocessing import Pool
import os
import pandas as pd
import random
import re
from sklearn.model_selection import train_test_split
from tqdm import tqdm
bar_format = "{percentage:3.0f}%|{bar:20}{r_bar}"

# Globals
scripts_dir = os.path.dirname(os.path.realpath(__file__))
CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "intensities_dir", type=click.Path(exists=True, resolve_path=True)
)
@click.option(
    "-n", "--no-linker",
    help="Exclude the linker sequence.",
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
    "-t", "--threads",
    help="Threads to use.",
    type=int,
    default=1,
    show_default=True
)

def main(**args):

    # Initialize
    intensity_files = []

    # Create output dir
    if not os.path.exists(args["output_dir"]):
        os.makedirs(args["output_dir"])

    # Get intensity files
    d = args["intensities_dir"]
    for intensity_file in os.listdir(d):
        intensity_files.append(os.path.join(d, intensity_file))

    # Get FASTA sequences
    kwargs = {"total": len(intensity_files), "bar_format": bar_format}
    pool = Pool(args["threads"])
    p = partial(__get_FASTA_sequences, no_linker=args["no_linker"],
        output_dir=args["output_dir"])
    for _ in tqdm(pool.imap(p, intensity_files), **kwargs):
        pass

def __get_FASTA_sequences(intensity_file, no_linker=False, output_dir="./"):

    # Initialize
    sequences = []
    cols = [4, 5, 6, 7]
    col_names = ["name", "sequence", "linker_sequence", "signal"]
    m = re.search("^(\S+)\.\w+@\S+@\S+@\w+\.(\S+)\.\S+\.tsv$",
        os.path.split(intensity_file)[1])
    prefix = "%s@%s" % (m.group(1), m.group(2))

    # Intensities as pandas DataFrame
    df = pd.read_csv(intensity_file, sep="\t", skiprows=1, usecols=cols,
        names=col_names)

    # Sequences
    for _, row in df.iterrows():
        if no_linker:
            seq = Seq(row["sequence"])
        else:
            seq = Seq(row["sequence"] + row["linker_sequence"])
        sequences.append(
            SeqRecord(seq, row["name"], description=str(row["signal"]))
        )

    # Split sequences 80/20 for training/validation
    train, validation = train_test_split(sequences, test_size=.2)

    # Save sequences
    random.shuffle(train)
    sequences_file = os.path.join(output_dir, "%s.train.fa.gz" % prefix)
    with gzip.open(sequences_file, "wt") as handle:
        SeqIO.write(train, handle, "fasta")
    random.shuffle(validation)
    sequences_file = os.path.join(output_dir, "%s.validation.fa.gz" % prefix)
    with gzip.open(sequences_file, "wt") as handle:
        SeqIO.write(validation, handle, "fasta")

if __name__ == "__main__":
    main()