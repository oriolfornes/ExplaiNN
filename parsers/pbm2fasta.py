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

    # Create output dir
    if not os.path.exists(params["output_dir"]):
        os.makedirs(params["output_dir"])

    # Get intensity files
    d = params["intensities_dir"]
    intensity_files = [os.path.join(d, f) for f in os.listdir(d)]

    # Get FASTA sequences
    kwargs = {"total": len(intensity_files), "bar_format": bar_format}
    pool = Pool(params["threads"])
    p = partial(__get_FASTA_sequences, no_linker=params["no_linker"],
        output_dir=params["output_dir"])
    for _ in tqdm(pool.imap(p, intensity_files), **kwargs):
        pass

def __get_FASTA_sequences(intensity_file, no_linker=False, output_dir="./"):

    # Initialize
    records = []
    # YY1.NA@PBM.ME@PBM14342.5GTGAAATTGTTATCCGCTCT@QNZS.wimpy-thistle-catfish.Train.tsv
    # YY1.NA@PBM.HK@PBM14358.5GTGAAATTGTTATCCGCTCT@QNZS.breezy-magnolia-dane.Val.tsv
    # YY1.NA@PBM.ME@PBM14342.5GTGAAATTGTTATCCGCTCT@SDQN.squirrely-blue-quokka.Train.tsv
    # YY1.NA@PBM.HK@PBM14358.5GTGAAATTGTTATCCGCTCT@SDQN.scanty-xanthic-bulldog.Val.tsv
    m = re.search("^(\S+)@\S+@\S+@\w+\.(\S+)\.\S+\.tsv$",
        os.path.split(intensity_file)[1])
    prefix = "%s@%s" % (m.group(1), m.group(2))
    sequences_file = os.path.join(output_dir, "%s.fa.gz" % prefix)

    # Intensities as pandas DataFrame
    df = pd.read_csv(intensity_file, sep="\t", skiprows=1, usecols=[4, 5, 6, 7],
        names=["name", "sequence", "linker_sequence", "signal"])
    # df.sort_values(by="signal", ascending=False, inplace=True)
    # positives = df.iloc[:int(df.shape[0]/2.), :].copy()
    # negatives = df.iloc[int(df.shape[0]/2.):, :].copy()
    # negatives.sort_values(by="signal", inplace=True)

    # Save sequences
    for _, row in df.iterrows():
        if no_linker:
            s = Seq(row["sequence"])
        else:
            s = Seq(row["sequence"] + row["linker_sequence"])
        r = SeqRecord(s, row["name"], description=str(row["signal"]))
        records.append(r)
    random.shuffle(records)
    # records.append([])
    # for _, row in positives.iterrows():
    #     if no_linker:
    #         s = Seq(row["sequence"])
    #     else:
    #         s = Seq(row["sequence"] + row["linker_sequence"])
    #     record = SeqRecord(s, row["name"], description=str(row["signal"]))
    #     records[-1].append(record)
    # records.append([])
    # for _, row in negatives.iterrows():
    #     if no_linker:
    #         s = Seq(row["sequence"])
    #     else:
    #         s = Seq(row["sequence"] + row["linker_sequence"])
    #     record = SeqRecord(s, row["name"], description=str(row["signal"]))
    #     records[-1].append(record)
    with gzip.open(sequences_file, "wt") as handle:
        SeqIO.write(records, handle, "fasta")
        # sequences = []
        # while True:
        #     exit_loop = True
        #     for i in range(len(records)):
        #         if len(records[i]) > 0:
        #             exit_loop = False
        #             sequences.append(records[i].pop())
        #         else:
        #             pass
        #     if exit_loop:
        #         break
        # SeqIO.write(sequences, handle, "fasta")

if __name__ == "__main__":
    main()