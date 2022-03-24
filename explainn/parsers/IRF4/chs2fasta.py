#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import click
import gzip
import os
from pybedtools import BedTool
from pybedtools.helpers import cleanup
import random
import re
from sklearn.model_selection import train_test_split
import subprocess as sp

from utils import get_chrom_sizes, get_file_handle

# Globals
scripts_dir = os.path.dirname(os.path.realpath(__file__))

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "chip_seq_dir", type=click.Path(exists=True, resolve_path=True)
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

def main(**args):

    # Create output dir
    if not os.path.exists(args["output_dir"]):
        os.makedirs(args["output_dir"])

    # Get chrom sizes
    chrom_sizes = get_chrom_sizes("%s.sizes" % args["genome_file"])

    # Get FASTA sequences
    d = args["chip_seq_dir"]
    for bed_file in ["WT.bed", "T95R.bed"]:
        __get_FASTA_sequences(os.path.join(args["chip_seq_dir"], bed_file),
            args["genome_file"], chrom_sizes, args["dummy_dir"], args["kmer"],
            args["output_dir"])

def __get_FASTA_sequences(bed_file, genome_file, chrom_sizes,
    dummy_dir="/tmp/", kmer=3, output_dir="./"):

    # Initialize
    positives = []
    m = re.search("^(\S+)\.bed$", os.path.split(bed_file)[1])
    prefix = m.group(1)
    dummy_file = os.path.join(dummy_dir, "%s.fa" % prefix)

    # Get BedTool object
    a = BedTool(bed_file)
    a.set_chromsizes(chrom_sizes)
    a.sequence(fi=genome_file, name=True)

    # Positive sequences
    for s in SeqIO.parse(a.seqfn, "fasta"):
        positives.append(SeqRecord(s.seq, id=s.id, description="1.0"))

    # Split sequences 80/20 for training/validation
    train_split, validation_split = train_test_split(positives, test_size=.2)
    train = [train_split, []]
    validation = [validation_split, []]

    # Negative sequences by dinucleotide shuffling
    fh = get_file_handle(dummy_file, "wt")
    SeqIO.write(train[0], fh, "fasta")
    fh.close()
    cmd = "biasaway k -f %s -k %s -e 1 > %s.biasaway" % (dummy_file, kmer,
        dummy_file)
    _ = sp.run([cmd], shell=True, cwd=scripts_dir, stderr=sp.DEVNULL)
    for s in SeqIO.parse("%s.biasaway" % dummy_file, "fasta"):
        id = "%s::shuf" % s.description.split(" ")[-1]
        train[1].append(SeqRecord(s.seq, description="0.0", id=id))
    fh = get_file_handle(dummy_file, "wt")
    SeqIO.write(validation[0], fh, "fasta")
    fh.close()
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
    fh = get_file_handle(sequences_file, "wt")
    sequences = [j for i in map(list, zip(*train)) for j in i]
    SeqIO.write(sequences, fh, "fasta")
    fh.close()
    sequences_file = os.path.join(output_dir, "%s.validation.fa.gz" % prefix)
    fh = get_file_handle(sequences_file, "wt")
    sequences = [j for i in map(list, zip(*validation)) for j in i]
    SeqIO.write(sequences, handle, "fasta")
    fh.close()

    # Delete tmp files
    cleanup()
    os.remove(dummy_file)
    os.remove("%s.biasaway" % dummy_file)

if __name__ == "__main__":
    main()