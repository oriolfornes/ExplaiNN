#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import click
from copy import copy
from functools import partial
import gzip
from multiprocessing import Pool
import os
import random
import re
from sklearn.model_selection import train_test_split
from tqdm import tqdm
bar_format = "{percentage:3.0f}%|{bar:20}{r_bar}"

# Globals
CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
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
@click.option(
    "-t", "--threads",
    help="Threads to use.",
    type=int,
    default=1,
    show_default=True
)

def main(**args):

    # Initialize
    reads_files = []
    grouped_reads_files = {}

    # Create output dir
    if not os.path.exists(args["output_dir"]):
        os.makedirs(args["output_dir"])

    # Get reads files
    d = args["reads_dir"]
    for reads_file in os.listdir(d):
        reads_files.append(reads_file)

    # Group reads files from different cycles
    for reads_file in sorted(reads_files):
        m = re.search("^(\S+@\S+@\S+)\.C\d\.5\w+\.3\w+@Reads", reads_file)
        grouped_reads_files.setdefault(m.group(1), list())
        grouped_reads_files[m.group(1)].append(os.path.join(d, reads_file))
    grouped_reads_files = list(grouped_reads_files.values())

    # Get FASTA sequences
    kwargs = {"total": len(grouped_reads_files), "bar_format": bar_format}
    pool = Pool(args["threads"])
    p = partial(__get_FASTA_sequences, output_dir=args["output_dir"])
    for _ in tqdm(pool.imap(p, grouped_reads_files), **kwargs):
        pass

def __get_FASTA_sequences(reads_files, output_dir="./"):

    # Initialize
    tf = None
    prefixes = []
    train = []
    validation = []
    sequences = {}
    cycles2sequences = {}
    cycles = ["0." for _ in reads_files]

    # For each read file...
    for i, reads_file in enumerate(sorted(reads_files)):

        # Get prefixes
        m = re.search("^(\S+)\.\w+@\S+@\S+@\w+\.(\S+)\.\S+\.fastq.gz$",
            os.path.basename(reads_file))
        tf = m.group(1)
        prefixes.append(m.group(2))

        # Sequences
        for s in SeqIO.parse(gzip.open(reads_file, "rt"), "fastq"):
            if s.seq.count("N") > 0: continue
            seq = str(s.seq)
            sequences.setdefault(seq, copy(cycles))
            sequences[seq][i] = "1."

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
    prefix = "%s@%s" % (tf, "+".join(prefixes))
    sequences_file = os.path.join(output_dir, "%s.train.fa.gz" % prefix)
    with gzip.open(sequences_file, "wt") as handle:
        SeqIO.write(train, handle, "fasta")
    sequences_file = os.path.join(output_dir, "%s.validation.fa.gz" % prefix)
    with gzip.open(sequences_file, "wt") as handle:
        SeqIO.write(validation, handle, "fasta")

if __name__ == "__main__":
    main()