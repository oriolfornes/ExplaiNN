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

def main(**params):

    # Initialize
    global intensity_files
    intensity_files = {}

    # Create output dir
    if not os.path.exists(params["output_dir"]):
        os.makedirs(params["output_dir"])

    # Get intensity files
    d = params["intensities_dir"]
    for intensity_file in os.listdir(d):
        m = re.search("^(\S+)@\S+@\S+@\w+\.\S+\.\S+\.tsv$", intensity_file)
        intensity_files.setdefault(m.group(1), [])
        intensity_files[m.group(1)].append(os.path.join(d, intensity_file))

    # Get FASTA sequences
    kwargs = {"total": len(intensity_files), "bar_format": bar_format}
    pool = Pool(params["threads"])
    p = partial(__get_FASTA_sequences, no_linker=params["no_linker"],
        output_dir=params["output_dir"])
    for _ in tqdm(pool.imap(p, intensity_files), **kwargs):
        pass

def __get_FASTA_sequences(tf, no_linker=False, output_dir="./"):

    # Initialize
    records = []
    words = []
    cols = [4, 5, 6, 7]
    names = ["name", "sequence", "linker_sequence", "signal"]

    # For each intensity file...
    for intensity_file in intensity_files[tf]:

        # Initialize
        m = re.search("^\S+@\S+@\S+@\w+\.(\S+)\.\S+\.tsv$",
            os.path.split(intensity_file)[1])
        #prefix = "%s@%s" % (m.group(1), m.group(2))
        words.append(m.group(1))

        # Intensities as pandas DataFrame
        df = pd.read_csv(intensity_file, sep="\t", skiprows=1, usecols=cols,
            names=names)
        df["name"] = [f"{words[-1]}::{n}" for n in df["name"].tolist()]

        # Save records
        for _, row in df.iterrows():
            if no_linker:
                s = Seq(row["sequence"])
            else:
                s = Seq(row["sequence"] + row["linker_sequence"])
            r = SeqRecord(s, row["name"], description=str(row["signal"]))
            records.append(r)

    # Save sequences
    random.shuffle(records)
    prefix = "%s@%s" % (tf, "+".join(sorted(words)))
    sequences_file = os.path.join(output_dir, "%s.fa.gz" % prefix)
    if os.path.exists(sequences_file):
        return
    with gzip.open(sequences_file, "wt") as handle:
        SeqIO.write(records, handle, "fasta")

if __name__ == "__main__":
    main()