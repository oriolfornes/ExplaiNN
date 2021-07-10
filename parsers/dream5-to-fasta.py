#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import click
from functools import partial
import gzip
from multiprocessing import Pool
import numpy as np
import os
import pandas as pd
import random
import re
from tqdm import tqdm
bar_format = "{percentage:3.0f}%|{bar:20}{r_bar}"

# Globals
scripts_dir = os.path.dirname(os.path.realpath(__file__))

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "challenge_file", type=click.Path(exists=True, resolve_path=True)
)
@click.argument(
    "samples_file", type=click.Path(exists=True, resolve_path=True)
)
@click.argument(
    "raw_dir", type=click.Path(exists=True, resolve_path=True)
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

def main(**params):

    # Create output dir
    if not os.path.exists(params["output_dir"]):
        os.makedirs(params["output_dir"])
    for data_split in ["Train", "Test"]:
        if not os.path.exists(os.path.join(params["output_dir"], data_split)):
            os.makedirs(os.path.join(params["output_dir"], data_split))

    # Load files as pandas DataFrames
    challenge = pd.read_csv(params["challenge_file"], sep="\t", header=0)
    samples = pd.read_csv(params["samples_file"], sep="\t", header=0)

    # Get FASTA sequences
    kwargs = {"total": len(challenge), "bar_format": bar_format}
    pool = Pool(params["threads"])
    p = partial(__get_FASTA_sequences, challenge=challenge, samples=samples,
        raw_dir=params["raw_dir"], no_linker=params["no_linker"],
        output_dir=params["output_dir"])
    for _ in tqdm(pool.imap(p, range(len(challenge))), **kwargs):
        pass

def __get_FASTA_sequences(ix, challenge, samples, raw_dir, no_linker=False,
    output_dir="./"):

    # Skip
    if pd.isnull(challenge.iloc[ix,]["Training Array"]):
        return

    # Initialize
    plasmid = challenge.iloc[ix,]["Plasmid"]
    train_array = challenge.iloc[ix,]["Training Array"]
    test_array = "ME" if train_array == "HK" else "HK"
    source_name = samples[samples["title"] == "%s_%s_8mer" % \
        (plasmid, train_array)]["source name"].values[0]
    m = re.search("(\S+), %s design" % train_array, source_name)
    tf = m.group(1)
    names=["name", "sequence", "linker_sequence", "signal"]

    # Save train sequences
    sequences_file = os.path.join(output_dir, "Train", "%s@%s.fa.gz" % \
        (tf, plasmid))
    if not os.path.exists(sequences_file):

        # Initialize
        sequences = []

        # Train intensities as pandas DataFrame
        intensity_file = os.path.join(raw_dir, "%s_%s.raw.tar.gz" % \
            (plasmid, train_array))
        df = pd.read_csv(intensity_file, sep="\t", skiprows=1,
            usecols=[4, 5, 6, 7], names=names)

        # Sequences
        for _, row in df.iterrows():
            if not pd.isnull(row["sequence"]):
                if no_linker:
                    s = Seq(row["sequence"])
                else:
                    s = Seq(row["sequence"] + row["linker_sequence"])
                signal = np.log10(row["signal"])
                record = SeqRecord(s, row["name"], description=str(signal))
                sequences.append(record)
        random.shuffle(sequences)

        # Save sequences
        with gzip.open(sequences_file, "wt") as handle:
            SeqIO.write(sequences, handle, "fasta")

    # Save test sequences
    sequences_file = os.path.join(output_dir, "Test", "%s@%s.fa.gz" % \
        (tf, plasmid))
    if not os.path.exists(sequences_file):

        # Initialize
        sequences = []

        # Test intensities as pandas DataFrame
        intensity_file = os.path.join(raw_dir, "%s_%s.raw.tar.gz" % \
            (plasmid, test_array))
        df = pd.read_csv(intensity_file, sep="\t", skiprows=1,
            usecols=[4, 5, 6, 7], names=names)

        # Sequences
        for _, row in df.iterrows():
            if not pd.isnull(row["sequence"]):
                if no_linker:
                    s = Seq(row["sequence"])
                else:
                    s = Seq(row["sequence"] + row["linker_sequence"])
                signal = np.log10(row["signal"])
                record = SeqRecord(s, row["name"], description=str(signal))
                sequences.append(record)
        random.shuffle(sequences)

        # Save sequences
        with gzip.open(sequences_file, "wt") as handle:
            SeqIO.write(sequences, handle, "fasta")

if __name__ == "__main__":
    main()