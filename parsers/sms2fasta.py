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

    # Create output dir
    if not os.path.exists(params["output_dir"]):
        os.makedirs(params["output_dir"])

    # Get reads files
    d = params["reads_dir"]
    reads_files = [os.path.join(d, f) for f in os.listdir(d)]

    # Get FASTA sequences
    kwargs = {"total": len(reads_files), "bar_format": bar_format}
    pool = Pool(params["threads"])
    p = partial(__get_FASTA_sequences, dummy_dir=params["dummy_dir"],
        kmer=params["kmer"], output_dir=params["output_dir"])
    for _ in tqdm(pool.imap(p, reads_files), **kwargs):
        pass

def __get_FASTA_sequences(reads_file, dummy_dir="/tmp/", kmer=2,
    output_dir="./"):

    # Initialize
    positives = []
    negatives = []
    # YY1.NA@SMS@SRR3405098.5GCTCTTCCGATCTTTGGATC.3GATCCAAGATCGGAAGAGCT@Reads.tasty-beige-bonobo.Train.fastq.gz
    # YY1.NA@SMS@SRR3405097.5GCTCTTCCGATCTGATGCCA.3TGGCATCGATCGGAAGAGCT@Reads.skinny-olivine-dunker.Val.fastq.gz
    m = re.search("^(\S+)@\S+@\S+@\w+\.(\S+)\.\S+\.fastq.gz$",
        os.path.basename(reads_file))
    prefix = "%s@%s" % (m.group(1), m.group(2))
    dummy_file = os.path.join(dummy_dir, "%s.fa" % prefix)
    sequences_file = os.path.join(output_dir, "%s.fa.gz" % prefix)
    if os.path.exists(sequences_file):
        return

    # Positive sequences
    with open(dummy_file, "wt") as handle:
        for s in SeqIO.parse(gzip.open(reads_file, "rt"), "fastq"):
            if s.seq.count("N") > 0: continue
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