#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import click
from copy import copy
import gzip
import os
import random
import re
from sklearn.model_selection import train_test_split

from utils import get_file_handle

# Globals
CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}
prefixes = {
    "YWIB1_A07_AffiSeq_Read1.fastq.gz": "WT",
    "YWIB1_A07_AffiSeq_Read2.fastq.gz": "WT",
    "YWIB2_A07_AffiSeq_Read1.fastq.gz": "WT",
    "YWIB2_A07_AffiSeq_Read2.fastq.gz": "WT",
    "YWIB3_A07_AffiSeq_Read1.fastq.gz": "WT",
    "YWIB3_A07_AffiSeq_Read2.fastq.gz": "WT",
    "YWIB1_B07_AffiSeq_Read1.fastq.gz": "T95R",
    "YWIB1_B07_AffiSeq_Read2.fastq.gz": "T95R",
    "YWIB2_B07_AffiSeq_Read1.fastq.gz": "T95R",
    "YWIB2_B07_AffiSeq_Read2.fastq.gz": "T95R",
    "YWIB3_B07_AffiSeq_Read1.fastq.gz": "T95R",
    "YWIB3_B07_AffiSeq_Read2.fastq.gz": "T95R"
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "reads_dir", type=click.Path(exists=True, resolve_path=True)
)
@click.option(
    "-o", "--output-dir",
    help="Output directory.",
    type=click.Path(resolve_path=True),
    default="./",
    show_default=True
)

def main(**args):

    # Initialize
    grouped_reads_files = {}

    # Create output dir
    if not os.path.exists(args["output_dir"]):
        os.makedirs(args["output_dir"])

    # Group reads files from different cycles
    d = args["reads_dir"]
    for reads_file in sorted(prefixes):
        grouped_reads_files.setdefault(prefixes[reads_file], list())
        grouped_reads_files[prefixes[reads_file]].append(
            os.path.join(d, reads_file)
        )

    # Get FASTA sequences
    for prefix, reads_files in grouped_reads_files.items():
        __get_FASTA_sequences(prefix, sorted(reads_files), args["output_dir"])

def __get_FASTA_sequences(prefix, reads_files, output_dir="./"):

    # Initialize
    train = []
    validation = []
    test = []
    sequences = {}
    cycles2sequences = {}
    cycles = ["0." for _ in range(0, len(reads_files), 2)]

    # For each read file...
    for i, j in enumerate(range(0, len(reads_files), 2)):

        # Sequences
        fh = get_file_handle(reads_files[j], "rt")
        for s in SeqIO.parse(fh, "fastq"):
            if s.seq.count("N") > 0: continue
            seq = str(s.seq)
            sequences.setdefault(seq, copy(cycles))
            sequences[seq][i] = "1."
        fh.close()
        fh = get_file_handle(reads_files[j+1], "rt")
        for s in SeqIO.parse(fh, "fastq"):
            if s.seq.count("N") > 0: continue
            seq = str(s.seq)
            sequences.setdefault(seq, copy(cycles))
            sequences[seq][i] = "1."
        fh.close()

    # For each sequence...
    for i, s in enumerate(sorted(sequences)):
        cycles2sequences.setdefault(tuple(sequences[s]), [])
        cycles2sequences[tuple(sequences[s])].append((i, s))  

    # For each cycle...
    for c in cycles2sequences:

        # Initialize
        records = []
        d = ";".join(c)

        # Sequences as records
        for i, s in cycles2sequences[c]:
            records.append(SeqRecord(Seq(s), id=str(i), description=d))

        # Split sequences 80/20 for training/validation
        try:
            train_split, validation_split = train_test_split(
                records, test_size=.2
            )
        except:
            train_split = records
            validation_split = []
        train.extend(train_split)
        validation.extend(validation_split)

    # Random shuffle
    random.shuffle(train)
    random.shuffle(validation)

    # Save sequences
    sequences_file = os.path.join(output_dir, "%s.train.fa.gz" % prefix)
    fh = get_file_handle(sequences_file, "wt")
    SeqIO.write(train, fh, "fasta")
    fh.close()
    sequences_file = os.path.join(output_dir, "%s.validation.fa.gz" % prefix)
    fh = get_file_handle(sequences_file, "wt")
    SeqIO.write(validation, fh, "fasta")
    fh.close()

if __name__ == "__main__":
    main()