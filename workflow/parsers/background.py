#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import click
from functools import partial
import gzip
from multiprocessing import Pool
import os
import pandas as pd
from pybedtools import BedTool
from pybedtools.helpers import cleanup
import re
from tqdm import tqdm
from urllib.request import urlretrieve
bar_format = "{percentage:3.0f}%|{bar:20}{r_bar}"

# Globals
scripts_dir = os.path.dirname(os.path.realpath(__file__))

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "input_file", type=click.Path(exists=True, resolve_path=True)
)
@click.option(
    "-g", "--genome",
    help="Genome.",
    type=click.Choice(["hg38", "mm10"]),
    required=True
)
@click.option(
    "-i", "--input",
    help="Input type.",
    type=click.Choice(["bed", "fasta"]),
    required=True
)
@click.option(
    "-l", "--length",
    help="Sequence length (in bp).",
    type=int,
    required=True
)
@click.option(
    "-o", "--output-file",
    help="Output file.  [default: stdout]",
    type=click.File(mode="w"),
    default="-"
)
@click.option(
    "-t", "--threads",
    help="Threads to use.",
    type=int,
    default=1,
    show_default=True
)

def main(**params):

    seqs = get_background_seqs(params["input_file"], params["genome"],
        params["input"], params["length"], params["threads"])

    # # Get FASTA sequences
    # kwargs = {"total": len(intensity_files), "bar_format": bar_format}
    # pool = Pool(params["threads"])
    # p = partial(__get_FASTA_sequences, no_linker=params["no_linker"],
    #     output_dir=params["output_dir"])
    # for _ in tqdm(pool.imap(p, intensity_files), **kwargs):
    #     pass

def get_background_seqs(input_file, genome, input_type, length, threads=1):

    # Get cis-regulatory regions (CREs)
    cres = __get_CREs(genome, length)

def __get_CREs(genome, length):

    # Initialize
    bed_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), ".%s-ccREs.bed" % genome
    )
    url = "https://api.wenglab.org/screen_v13/fdownloads/"
    if genome == "hg38":
        url += "GRCh38-ccREs.bed"
    else:
        url += "mm10-ccREs.bed"

    # Download CREs
    if not os.path.exists(bed_file):
        urlretrieve(url, bed_file)

    # BedTool
    b = BedTool(bed_file)

    # Get centers
    centers = []
    for i in b:
        center = int((i.start+i.end)/2)
        centers.append([i.chrom, str(center - 1), str(center)])
    c = BedTool("\n".join(["\t".join(c) for c in centers]), from_string=True)
    print(c)
    exit(0)



def __get_FASTA_sequences(intensity_file, no_linker=False, output_dir="./"):

    # Initialize
    sequences = []
    prefix = re.search("^(\S+).spatialDetrend_quantNorm.pbm.\S+.txt$",
        os.path.split(intensity_file)[1]).group(1)
    sequences_file = os.path.join(output_dir, "%s.fa.gz" % prefix)

    # Intensities as pandas DataFrame
    df = pd.read_csv(intensity_file, sep="\t", skiprows=1, usecols=[4, 5, 6, 7],
        names=["name", "sequence", "linker_sequence", "signal"])

    # Save sequences
    for _, row in df.iterrows():
        if no_linker:
            s = Seq(row["sequence"])
        else:
            s = Seq(row["sequence"] + row["linker_sequence"])
        record = SeqRecord(s, row["name"], description=str(row["signal"]))
        sequences.append(record)
    with gzip.open(sequences_file, "wt") as handle:
        SeqIO.write(sequences, handle, "fasta")

if __name__ == "__main__":
    main()