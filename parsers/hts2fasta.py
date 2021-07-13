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
# import pandas as pd
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

def main(**params):

    # Create output dir
    if not os.path.exists(params["output_dir"]):
        os.makedirs(params["output_dir"])

    # Get reads files
    d = params["reads_dir"]
    reads_files = [os.path.join(d, f) for f in os.listdir(d)]

    # Group reads files from different cycles
    grouped_reads_files = {}
    for f in sorted(reads_files):
        # ZNF8.FL@HTS.Lys@AAT_A_TC40NACGGGA.C1.5ACGACGCTCTTCCGATCTTC.3ACGGGAAGATCGGAAGAGCA@Reads.flabby-yellow-insect.Train.fastq.gz
        # ZNF8.FL@HTS.Lys@AAT_A_TC40NACGGGA.C2.5ACGACGCTCTTCCGATCTTC.3ACGGGAAGATCGGAAGAGCA@Reads.squirrely-pumpkin-wrasse.Train.fastq.gz
        # ZNF8.FL@HTS.Lys@AAT_A_TC40NACGGGA.C3.5ACGACGCTCTTCCGATCTTC.3ACGGGAAGATCGGAAGAGCA@Reads.gamy-chartreuse-sparrow.Train.fastq.gz
        m = re.search("^(\S+@\S+@\S+)\.C\d\.5\w+\.3\w+", os.path.basename(f))
        grouped_reads_files.setdefault(m.group(1), list())
        grouped_reads_files[m.group(1)].append(f)
    grouped_reads_files = [v for v in sorted(grouped_reads_files.values())]

    # Get FASTA sequences
    kwargs = {"total": len(grouped_reads_files), "bar_format": bar_format}
    pool = Pool(params["threads"])
    p = partial(__get_FASTA_sequences, output_dir=params["output_dir"])
    for _ in tqdm(pool.imap(p, grouped_reads_files), **kwargs):
        pass

def __get_FASTA_sequences(reads_files, output_dir="./"):

    # Initialize
    tf = None
    prefixes = []
    sequences = {}
    cycles = ["0." for _ in reads_files]
    # ZNF8.FL@HTS.Lys@AAT_A_TC40NACGGGA.C1.5ACGACGCTCTTCCGATCTTC.3ACGGGAAGATCGGAAGAGCA@Reads.flabby-yellow-insect.Train.fastq.gz
    # ZNF8.FL@HTS.Lys@AAT_A_TC40NACGGGA.C2.5ACGACGCTCTTCCGATCTTC.3ACGGGAAGATCGGAAGAGCA@Reads.squirrely-pumpkin-wrasse.Train.fastq.gz
    # ZNF8.FL@HTS.Lys@AAT_A_TC40NACGGGA.C3.5ACGACGCTCTTCCGATCTTC.3ACGGGAAGATCGGAAGAGCA@Reads.gamy-chartreuse-sparrow.Train.fastq.gz

    # For each read file...
    for i, reads_file in enumerate(reads_files):

        # Initialize
        m = re.search("^(\S+)@\S+@\S+@\w+\.(\S+)\.\S+\.fastq.gz$",
            os.path.basename(reads_file))
        tf = m.group(1)
        prefixes.append(m.group(2))

        # Sequences
        for s in SeqIO.parse(gzip.open(reads_file, "rt"), "fastq"):
            if s.seq.count("N") > 0: continue
            seq = str(s.seq)
            sequences.setdefault(seq, copy(cycles))
            sequences[seq][i] = "1."

    # Sequences as records
    records = []
    for id, seq in enumerate(sorted(sequences)):
        desc = ";".join(sequences[seq])
        records.append(SeqRecord(Seq(seq), id=str(id), description=desc))
    random.shuffle(records)

    # Save sequences
    prefix = "%s@%s" % (tf, "+".join(prefixes))
    sequences_file = os.path.join(output_dir, "%s.fa.gz" % prefix)
    with gzip.open(sequences_file, "wt") as handle:  
        SeqIO.write(records, handle, "fasta")

if __name__ == "__main__":
    main()