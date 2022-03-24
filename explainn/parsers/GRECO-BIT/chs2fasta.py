#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import click
import copy
from functools import partial
import gzip
from multiprocessing import Pool
import os
import pandas as pd
from pybedtools import BedTool
from pybedtools.helpers import cleanup
import random
import re
from sklearn.model_selection import train_test_split
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
    "intervals_dir", type=click.Path(exists=True, resolve_path=True)
)
@click.argument(
    "genome_file", type=click.Path(exists=True, resolve_path=True)
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
    default=3,
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

def main(**args):

    # Initialize
    intervals_files = []

    # Create output dir
    if not os.path.exists(args["output_dir"]):
        os.makedirs(args["output_dir"])

    # Get chrom sizes
    chrom_sizes = __get_chrom_sizes(args["genome_file"])

    # Get intervals files
    d = args["intervals_dir"]
    for intervals_file in os.listdir(d):
        intervals_files.append(os.path.join(d, intervals_file))

    # Get FASTA sequences
    kwargs = {"total": len(intervals_files), "bar_format": bar_format}
    pool = Pool(args["threads"])
    p = partial(__get_FASTA_sequences, genome_file=args["genome_file"],
        chrom_sizes=chrom_sizes, dummy_dir=args["dummy_dir"],
        kmer=args["kmer"], output_dir=args["output_dir"])
    for _ in tqdm(pool.imap(p, intervals_files), **kwargs):
        pass

def __get_chrom_sizes(genome_file):

    # Initialize
    chrom_sizes = {}

    with open("%s.sizes" % genome_file) as handle:
        for line in handle:
            chrom, size = line.strip("\n").split("\t")
            chrom_sizes.setdefault(chrom, tuple([0, size]))

    return(chrom_sizes)

def __get_FASTA_sequences(intervals_file, genome_file, chrom_sizes,
    dummy_dir="/tmp/", kmer=3, output_dir="./"):

    # Initialize
    positives = []
    cols = [0, 3, 8]
    col_names = ["chrom", "summit", "name"]
    m = re.search("^(\S+)\.\w+@\S+@\S+@\w+\.(\S+)\.\S+\.peaks$",
        os.path.split(intervals_file)[1])
    prefix = "%s@%s" % (m.group(1), m.group(2))
    dummy_file = os.path.join(dummy_dir, "%s.fa" % prefix)

    # Intervals as pandas DataFrame
    df = pd.read_csv(intervals_file, sep="\t", skiprows=1, usecols=cols,
        names=col_names)

    # Get start & end
    df["start"] = df["summit"] - 1
    df["end"] = df["summit"]
    df = df[["chrom", "start", "end", "name"]]

    # Get BedTool object
    s = "\n".join(["\t".join(map(str, r.tolist())) for _, r in df.iterrows()])
    a = BedTool(s, from_string=True)
    a.set_chromsizes(chrom_sizes)
    a = a.slop(b=100)
    a.sequence(fi=genome_file, name=True)

    # Positive sequences
    for s in SeqIO.parse(a.seqfn, "fasta"):
        positives.append(SeqRecord(s.seq, id=s.id, description="1.0"))

    # Split sequences 80/20 for training/validation
    train_split, validation_split = train_test_split(positives, test_size=.2)
    train = [train_split, []]
    validation = [validation_split, []]

    # Negative sequences by dinucleotide shuffling
    with open(dummy_file, "wt") as handle:
        SeqIO.write(train[0], handle, "fasta")
    cmd = "biasaway k -f %s -k %s -e 1 > %s.biasaway" % (dummy_file, kmer,
        dummy_file)
    _ = sp.run([cmd], shell=True, cwd=scripts_dir, stderr=sp.DEVNULL)
    for s in SeqIO.parse("%s.biasaway" % dummy_file, "fasta"):
        id = "%s::shuf" % s.description.split(" ")[-1]
        train[1].append(SeqRecord(s.seq, description="0.0", id=id))
    with open(dummy_file, "wt") as handle:
        SeqIO.write(validation[0], handle, "fasta")
    cmd = "biasaway k -f %s -k %s -e 1 > %s.biasaway" % (dummy_file, kmer,
        dummy_file)
    _ = sp.run([cmd], shell=True, cwd=scripts_dir, stderr=sp.DEVNULL)
    for s in SeqIO.parse("%s.biasaway" % dummy_file, "fasta"):
        id = "%s::shuf" % s.description.split(" ")[-1]
        validation[1].append(SeqRecord(s.seq, description="0.0", id=id))
    random.shuffle(train[0])
    random.shuffle(train[1])
    random.shuffle(validation[0])
    random.shuffle(validation[1])

    # Save sequences
    sequences_file = os.path.join(output_dir, "%s.train.fa.gz" % prefix)
    with gzip.open(sequences_file, "wt") as handle:
        sequences = [j for i in map(list, zip(*train)) for j in i]
        SeqIO.write(sequences, handle, "fasta")
    sequences_file = os.path.join(output_dir, "%s.validation.fa.gz" % prefix)
    with gzip.open(sequences_file, "wt") as handle:
        sequences = [j for i in map(list, zip(*validation)) for j in i]
        SeqIO.write(sequences, handle, "fasta")

    # Delete tmp files
    cleanup()
    os.remove(dummy_file)
    os.remove("%s.biasaway" % dummy_file)

if __name__ == "__main__":
    main()