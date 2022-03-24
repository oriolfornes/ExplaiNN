#!/usr/bin/env python

import click

import os
import pandas as pd
from pybedtools import BedTool
from pybedtools.helpers import cleanup
import re
import subprocess as sp
from tqdm import tqdm
bar_format = "{percentage:3.0f}%|{bar:20}{r_bar}"

from utils import get_chrom_sizes, get_file_handle

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
@click.argument(
    "dhs_file", type=click.Path(exists=True, resolve_path=True)
)
# @click.argument(
#     "tfs_file", type=click.Path(exists=True, resolve_path=True)
# )
@click.option(
    "-d", "--dummy-dir",
    help="Dummy directory.",
    type=click.Path(resolve_path=True),
    default="/tmp/",
    show_default=True
)
@click.option(
    "-l", "--length",
    help="Sequence length.",
    type=click.Choice(["201", "501", "1001"]),
    default="201",
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

    # Create output dir
    if not os.path.exists(args["output_dir"]):
        os.makedirs(args["output_dir"])

    # Get chrom sizes
    chrom_sizes = get_chrom_sizes("%s.sizes" % args["genome_file"])

    # Get already processed TFs
    tfs = set()
    for tsv_file in os.listdir(args["output_dir"]):
        m = re.search("^(\S+).(train|validation|test).tsv.gz$", tsv_file)
        tfs.add(m.group(1))

    # Get BED files
    remap_files = []
    for bed_file in os.listdir(args["remap_dir"]):
        m = re.search("^remap2022_(\S+)_nr_macs2_\w+_v1_0.bed.gz$", bed_file)
        if m.group(1) not in tfs:
            remap_files.append(os.path.join(args["remap_dir"], bed_file))
    # d = args["remap_dir"]
    # fh = get_file_handle(args["tfs_file"], "rt")
    # for line in fh:
    #     tf = line.strip()
    #     f = os.path.join(d, f"remap2022_{tf}_nr_macs2_mm10_v1_0.bed.gz")
    #     if os.path.exists(f):
    #         remap_files.append(f)
    # fh.close()

    # Get FASTA sequences
    kwargs = {"total": len(remap_files), "bar_format": bar_format}
    pool = Pool(args["threads"])
    p = partial(__get_FASTA_sequences, genome_file=args["genome_file"],
        dhs_file=args["dhs_file"], chrom_sizes=chrom_sizes,
        dummy_dir=args["dummy_dir"], length=int(args["length"]),
        output_dir=args["output_dir"])
    for _ in tqdm(pool.imap(p, remap_files), **kwargs):
        pass

def __get_FASTA_sequences(remap_file, genome_file, dhs_file, chrom_sizes,
    dummy_dir="/tmp/", length=201, output_dir="./"):

    # Initialize
    prefix = re.search("^remap2022_(\S+)_nr_macs2_\w+_v1_0.bed.gz$",
        os.path.split(remap_file)[1]).group(1)

    # Intervals as pandas DataFrame
    df = pd.read_csv(remap_file, sep="\t", usecols=[0, 6, 7],
        names=["chrom", "start", "end"])
    df = df[df["chrom"].isin(chrom_sizes.keys())]

    # Get non-redundant BedTool object (i.e. positive sequences)
    b = BedTool("\n".join(["\t".join(map(str, row.tolist())) \
            for _, row in df.iterrows()]), from_string=True).sort()
    b.set_chromsizes(chrom_sizes)
    b = b.slop(b=int((length-1)/2)).filter(lambda x: len(x) == length).sort()
    m = b.merge()
    nr = BedTool("\n".join(["\t".join(map(str, i.fields[-3:])) \
            for i in m.closest(b, t="first")]), from_string=True)
    nr.sequence(fi=genome_file)
    positive_file = os.path.join(dummy_dir, "%s_pos.fa" % prefix)
    nr.save_seqs(positive_file)

    # Get non-overlapping BedTool object (i.e. negative sequences)
    b = BedTool(dhs_file)
    nov = b.intersect(nr, v=True)
    nov.sequence(fi=genome_file)
    negative_file = os.path.join(dummy_dir, "%s_neg.fa" % prefix)
    nov.save_seqs(negative_file)

    # Subsample negative sequences by %GC
    json_file = os.path.join(dummy_dir, "%s.json" % prefix)
    cmd = "./match-seqs-by-gc.py -f -o %s %s %s" % \
        (json_file, negative_file, positive_file)
    _ = sp.run([cmd], shell=True, cwd=scripts_dir, stdout=sp.DEVNULL,
        stderr=sp.DEVNULL)

    # Subsample negative sequences by %GC
    cmd = "./json2cam.py -o %s -p %s --test %s" % \
        (output_dir, prefix, json_file)
    _ = sp.run([cmd], shell=True, cwd=scripts_dir, stdout=sp.DEVNULL,
        stderr=sp.DEVNULL)

    # Delete tmp files
    cleanup()
    os.remove(positive_file)
    os.remove(negative_file)
    os.remove(json_file)

if __name__ == "__main__":
    main()