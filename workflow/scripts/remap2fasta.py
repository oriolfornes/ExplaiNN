#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
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
from sklearn.model_selection import train_test_split
import subprocess as sp
from tqdm import tqdm

# Globals
scripts_dir = os.path.dirname(os.path.realpath(__file__))

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "remap_dir", type=click.Path(exists=True, resolve_path=True)
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
    if not os.path.exists(os.path.join(params["output_dir"], "train")):
        os.makedirs(os.path.join(params["output_dir"], "train"))
    if not os.path.exists(os.path.join(params["output_dir"], "validation")):
        os.makedirs(os.path.join(params["output_dir"], "validation"))

    # Get chrom sizes
    chrom_sizes = __get_chrom_sizes(params["genome_file"])

    # Get intervals files
    d = params["remap_dir"]
    remap_files = [os.path.join(d, f) for f in sorted(os.listdir(d)) \
        if f.startswith("remap2020")]

    __get_FASTA_sequences(remap_files[1], genome_file=params["genome_file"],
        chrom_sizes=chrom_sizes, dummy_dir=params["dummy_dir"],
        kmer=params["kmer"], output_dir=params["output_dir"])
    exit(0)

    # Get FASTA sequences
    bar_format = "{percentage:3.0f}%|{bar:20}{r_bar}"
    kwargs = {"total": len(remap_files), "bar_format": bar_format}
    pool = Pool(params["threads"])
    p = partial(__get_FASTA_sequences, genome_file=params["genome_file"],
        chrom_sizes=chrom_sizes, dummy_dir=params["dummy_dir"],
        kmer=params["kmer"], output_dir=params["output_dir"])
    for _ in tqdm(pool.imap(p, remap_files), **kwargs):
        pass

def __get_chrom_sizes(genome_file):

    # Initialize
    chrom_sizes = {}

    with open("%s.sizes" % genome_file) as handle:
        for line in handle:
            chrom, size = line.strip("\n").split("\t")
            chrom_sizes.setdefault(chrom, tuple([0, size]))

    return(chrom_sizes)

def __get_FASTA_sequences(remap_file, genome_file, chrom_sizes,
    dummy_dir="/tmp/", kmer=2, output_dir="./"):

    # Initialize
    prefix = re.search("^remap2020\_(\S+)\_nr\_macs2\_hg38\_v1\_0\.bed\.gz$",
        os.path.split(remap_file)[1]).group(1)
    dummy_file = os.path.join(dummy_dir, "%s.fa" % prefix)

    # Intervals as pandas DataFrame
    df = pd.read_csv(remap_file, sep="\t", usecols=[0, 4, 6, 7],
        names=["chrom", "count", "start", "end"])
    df = df[df["chrom"].isin(chrom_sizes.keys())].sort_values("count",
        ascending=False)
    df = df[["chrom", "start", "end", "name"]]
    print(df)
    exit(0)

    # Get BedTool object
    b = BedTool("\n".join(["\t".join(map(str, row.tolist())) \
        for _, row in df.iterrows()]), from_string=True)
    b.set_chromsizes(chrom_sizes)
    b = b.slop(b=100)
    m = b.sort().merge()
    print(len(b), len(m))
    exit(0)
    b.sequence(fi=genome_file, name=True)

    # Positive sequences
    sequences = []
    # with open(dummy_file, "wt") as handle:
    for count, s in enumerate(SeqIO.parse(b.seqfn, "fasta")):
        # Filter non-ACGTs
        if re.search("[^acgtACGT]", str(s.seq)):
            continue
        # handle.write(seq.format("fasta"))
        sequences.append(SeqRecord(s.seq, id="%s%s" % (count + 1, s.id),
            description="1."))

    # Split sequences into train and validation
    sequences = __split_seqs_into_train_validation(sequences)

    # Negative sequences by dinucleotide shuffling
    for i in range(len(sequences)):
        with open(dummy_file, "wt") as handle:
            for s in sequences[i]:
                handle.write(s.format("fasta"))
        cmd = "./fasta-shuffle-letters -kmer 2 -dna -tag '::shuf' -seed 1 %s" \
            % dummy_file
        output = sp.run([cmd], shell=True, cwd=scripts_dir, stdout=sp.PIPE,
            stderr=sp.DEVNULL, encoding="utf8")
        with open(dummy_file, "wt") as handle:
            for line in output.stdout.splitlines():
                handle.write("%s\n" % line)
        for s in SeqIO.parse(dummy_file, "fasta"):
            sequences[i].append(SeqRecord(s.seq, id=s.id, description="0."))

    # Save sequences
    seqs_file = os.path.join(output_dir, "train", "%s.fa.gz" % prefix)
    with gzip.open(seqs_file, "wt") as handle:
        SeqIO.write(sequences[0], handle, "fasta")
    seqs_file = os.path.join(output_dir, "validation", "%s.fa.gz" % prefix)
    with gzip.open(seqs_file, "wt") as handle:
        SeqIO.write(sequences[1], handle, "fasta")

    # Delete tmp files
    cleanup()
    os.remove(dummy_file)

def __split_seqs_into_train_validation(sequences):

    # Initialize
    gc_contents = {}
    train_seqs = []
    val_seqs = []

    # For each record...
    for record in sequences:
        gc_content = round(GC(str(record.seq)))
        gc_contents.setdefault(gc_content, [])
        gc_contents[gc_content].append(record)

    # Get train (80%) and test (20%) sets:
    # Raschka S. Model Evaluation, Model Selection, and Algorithm Selection
    # in Machine Learning. arXiv; 2018.
    # Typically, we assign 2/3 to the training set and 1/3 of the data to the
    # test set. Other common training/test splits are 60/40, 70/30, or 80/20,
    # or even 90/10 if the dataset is relatively large.
    for gc_content in sorted(gc_contents):
        if len(gc_contents[gc_content]) < 5:
            random.seed(1)
            choices = [c for c in random.choices(["t", "v"], weights=[4, 1],
                k=len(gc_contents[gc_content]))]
            for c in range(len(choices)):
                if choices[c] == "t":
                    train_seqs.append(gc_contents[gc_content][c])
                else:
                    val_seqs.append(gc_contents[gc_content][c])
        else:
            data_splits = train_test_split(gc_contents[gc_content],
                test_size=.2, random_state=1)
            train_seqs += data_splits[0]
            val_seqs += data_splits[1]

    return(train_seqs, val_seqs)

if __name__ == "__main__":
    main()