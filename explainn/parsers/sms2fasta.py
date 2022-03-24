#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import click
from functools import partial
import gzip
from multiprocessing import Pool
import os
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
    global reads_files
    reads_files = {}

    # Create output dir
    if not os.path.exists(params["output_dir"]):
        os.makedirs(params["output_dir"])

    # Get reads files
    d = params["reads_dir"]
    for reads_file in os.listdir(d):
        m = re.search("^(\S+)@\S+@\S+@\w+\.\S+\.\S+\.fastq.gz$", reads_file)
        reads_files.setdefault(m.group(1), [])
        reads_files[m.group(1)].append(os.path.join(d, reads_file))

    # Get FASTA sequences
    kwargs = {"total": len(reads_files), "bar_format": bar_format}
    pool = Pool(params["threads"])
    p = partial(__get_FASTA_sequences, dummy_dir=params["dummy_dir"],
        kmer=params["kmer"], output_dir=params["output_dir"])
    for _ in tqdm(pool.imap(p, reads_files), **kwargs):
        pass

def __get_FASTA_sequences(tf, dummy_dir="/tmp/", kmer=2,
    output_dir="./"):

    # Initialize
    sequences = []
    positives = []
    negatives = []
    words = []

    # For each intervals file...
    for reads_file in reads_files[tf]:

        # Initialize
        m = re.search("^\S+@\S+@\S+@\w+\.(\S+)\.\S+\.fastq.gz$",
            os.path.split(reads_file)[1])
        words.append(m.group(1))

        # For each sequence...
        for s in SeqIO.parse(gzip.open(reads_file, "rt"), "fastq"):
            if s.seq.count("N") > 0:
                continue
            s.id = f"{words[-1]}::{s.id}"
            sequences.append(s)

    # Positive sequences
    prefix = "%s@%s" % (tf, "+".join(sorted(words)))
    dummy_file = os.path.join(dummy_dir, "%s.fa" % prefix)
    sequences_file = os.path.join(output_dir, "%s.fa.gz" % prefix)
    if os.path.exists(sequences_file):
        return
    with open(dummy_file, "wt") as handle:
        for s in sequences:
            handle.write(s.format("fasta"))
            positives.append(SeqRecord(s.seq, id=s.id, description="1."))
    random.shuffle(positives)

    # Negative sequences by dinucleotide shuffling
    cmd = "biasaway k -f %s -k %s -e 1 > %s.biasaway" % (dummy_file, kmer,
        dummy_file)
    _ = sp.run([cmd], shell=True, cwd=scripts_dir, stderr=sp.DEVNULL)
    for s in SeqIO.parse("%s.biasaway" % dummy_file, "fasta"):
        negatives.append(SeqRecord(s.seq, description="0.",
            id="%s::shuf" % s.description.split(" ")[-1]))
    random.shuffle(negatives)

    # Save sequences
    with gzip.open(sequences_file, "wt") as handle:
        sequences = [j for i in zip(positives, negatives) for j in i]
        SeqIO.write(sequences, handle, "fasta")

    # Delete tmp files
    os.remove(dummy_file)

if __name__ == "__main__":
    main()