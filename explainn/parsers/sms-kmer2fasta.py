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
import subprocess as sp
from tqdm import tqdm
bar_format = "{percentage:3.0f}%|{bar:20}{r_bar}"

# Globals
scripts_dir = os.path.dirname(os.path.realpath(__file__))

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "reads_dir", type=click.Path(exists=True, resolve_path=True)
)
@click.option(
    "-d", "--dummy-dir",
    help="Dummy directory.",
    type=click.Path(resolve_path=True),
    default="/tmp/",
    show_default=True
)
@click.option(
    "-k", "--kmer",
    help="K-mer size for KMC.",
    type=int,
    default=8,
    show_default=True
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

    # Get intervals files
    d = params["reads_dir"]
    reads_files = [os.path.join(d, f) for f in os.listdir(d)]

    # Get FASTA sequences
    kwargs = {"total": len(reads_files), "bar_format": bar_format}
    pool = Pool(params["threads"])
    p = partial(__get_FASTA_sequences, dummy_dir=params["dummy_dir"],
        kmer=params["kmer"], output_dir=params["output_dir"])
    for _ in tqdm(pool.imap(p, reads_files), **kwargs):
        pass

def __get_FASTA_sequences(reads_file, dummy_dir="/tmp/", kmer=10,
    output_dir="./"):

    # Initialize
    # positives = []
    # negatives = []
    # YY1.NA@SMS@SRR3405098.5GCTCTTCCGATCTTTGGATC.3GATCCAAGATCGGAAGAGCT@Reads.tasty-beige-bonobo.Train.fastq.gz
    # YY1.NA@SMS@SRR3405097.5GCTCTTCCGATCTGATGCCA.3TGGCATCGATCGGAAGAGCT@Reads.skinny-olivine-dunker.Val.fastq.gz
    sequences = []
    m = re.search("^(\S+)\.\S+@\S+@\S+@\w+\.(\S+)\.\w+\.fastq.gz$",
        os.path.split(reads_file)[1])
    prefix = "%s@%s" % (m.group(1), m.group(2))
    sequences_file = os.path.join(output_dir, "%s.fa.gz" % prefix)
    if os.path.exists(sequences_file):
        return

    # Count k-mers with KMC3
    database = os.path.join(dummy_dir, prefix)
    cmd = f"kmc -k{kmer} -ci{0} -cs{10**9} {reads_file} {database} {dummy_dir}"
    _ = sp.run([cmd], shell=True, cwd=scripts_dir, stderr=sp.DEVNULL,
        stdout=sp.DEVNULL)
    dummy_file = os.path.join(dummy_dir, "%s.txt" % prefix)
    cmd = f"kmc_dump {database} {dummy_file}"
    _ = sp.run([cmd], shell=True, cwd=scripts_dir, stderr=sp.DEVNULL,
        stdout=sp.DEVNULL)
    df = pd.read_csv(dummy_file, sep="\t", names=["Kmer", "Freq."])
    df["Prob."] = df["Freq."].div(df["Freq."].sum())
    kmers = dict(zip(df["Kmer"].to_list(), df["Prob."].to_list()))
    for k in frozenset(kmers):
        s = Seq(k)
        kmers.setdefault(str(s.reverse_complement()), kmers[k])

    # Sequences
    for s in SeqIO.parse(gzip.open(reads_file, "rt"), "fastq"):
        # Skip read
        if s.seq.count("N") > 0:
            continue
        score = 0.
        for i in range(len(s) - kmer):
            score += kmers[str(s.seq[i:i+kmer])]
        score = str(np.log10(score))
        sequences.append(SeqRecord(s.seq, id=s.id, description=score))
    random.shuffle(sequences)

    # Save sequences
    with gzip.open(sequences_file, "wt") as handle:
        SeqIO.write(sequences, handle, "fasta")

    # Delete tmp files
    os.remove("%s.kmc_pre" % database)
    os.remove("%s.kmc_suf" % database)
    os.remove(dummy_file)

if __name__ == "__main__":
    main()