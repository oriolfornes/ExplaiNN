#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import click
from functools import partial
import gzip
from multiprocessing import Pool
import os
import pandas as pd
from pybedtools import BedTool
from pybedtools.helpers import cleanup
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

    # Create output dir
    if not os.path.exists(params["output_dir"]):
        os.makedirs(params["output_dir"])

    # Get chrom sizes
    chrom_sizes = __get_chrom_sizes(params["genome_file"])

    # Get intervals files
    d = params["intervals_dir"]
    intervals_files = [os.path.join(d, f) for f in os.listdir(d)]

    # Get FASTA sequences
    kwargs = {"total": len(intervals_files), "bar_format": bar_format}
    pool = Pool(params["threads"])
    p = partial(__get_FASTA_sequences, genome_file=params["genome_file"],
        chrom_sizes=chrom_sizes, dummy_dir=params["dummy_dir"],
        kmer=params["kmer"], output_dir=params["output_dir"])
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
    dummy_dir="/tmp/", kmer=2, output_dir="./"):

    # Initialize
    positives = []
    negatives = []
    # ZNF189.FL@CHS@SI0201@Peaks.wiggy-green-bombay.Train.peaks
    # ZNF189.FL@CHS@SI0201@Peaks.flabby-thistle-numbat.Val-B.peaks
    m = re.search("^(\S+)\.\S+@\S+@\S+@\w+\.(\S+)\.\S+\.peaks$",
        os.path.split(intervals_file)[1])
    prefix = "%s@%s" % (m.group(1), m.group(2))
    dummy_file = os.path.join(dummy_dir, "%s.fa" % prefix)
    sequences_file = os.path.join(output_dir, "%s.fa.gz" % prefix)

    # Intervals as pandas DataFrame
    df = pd.read_csv(intervals_file, sep="\t", skiprows=1, usecols=[0, 3, 8],
        names=["chrom", "summit", "name"])
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
    with open(dummy_file, "wt") as handle:
        for s in SeqIO.parse(a.seqfn, "fasta"):
            handle.write(s.format("fasta"))
            positives.append(SeqRecord(s.seq, id=s.id, description="1."))
    random.shuffle(positives)

    # Negative sequences by dinucleotide shuffling
    cmd = "biasaway k -f %s -e 1" % dummy_file
    output = sp.run([cmd], shell=True, cwd=scripts_dir, stdout=sp.PIPE,
        stderr=sp.DEVNULL, encoding="utf8")
    with open(dummy_file, "wt") as handle:
        for line in output.stdout.splitlines():
            handle.write("%s\n" % line)
    for s in SeqIO.parse(dummy_file, "fasta"):
        negatives.append(SeqRecord(s.seq, description="0.",
            id="%s::shuf" % s.description.split(" ")[-1]))
    random.shuffle(negatives)

    # Save sequences
    with gzip.open(sequences_file, "wt") as handle:
        sequences = [j for i in zip(positives, negatives) for j in i]
        SeqIO.write(sequences, handle, "fasta")

    # Delete tmp files
    cleanup()
    os.remove(dummy_file)

if __name__ == "__main__":
    main()