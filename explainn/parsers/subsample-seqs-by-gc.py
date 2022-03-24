#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqUtils import GC
import click
import math
import random
import sys

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "fasta_file",
    type=click.Path(exists=True, resolve_path=True),
)
@click.option(
    "-f", "--filter-masked",
    help="Filter masked DNA sequences.",
    is_flag=True,
)
@click.option(
    "-s", "--subsample",
    help="Number of sequences to subsample.",
    type=int,
    default=1000,
    show_default=True,
)
@click.option(
    "-o", "--output-file",
    help="Output file.  [default: STDOUT]",
    type=click.Path(writable=True, readable=False, resolve_path=True,
        allow_dash=True),
)

def main(**args):

    # Group sequences by %GC content
    gc_groups = {}
    for record in SeqIO.parse(args["fasta_file"], "fasta"):
        if args["filter_masked"]:
            if record.seq.count("N") or record.seq.count("n"):
                continue
        gc = round(GC(record.seq))
        gc_groups.setdefault(gc, [])
        gc_groups[gc].append(record)

    # Subsampling
    sampled = []
    random_seed = 123
    norm_factor = args["subsample"] / \
        sum([len(v) for v in gc_groups.values()])
    for i in sorted(gc_groups):
        random.Random(random_seed).shuffle(gc_groups[i])
        sampled.extend(gc_groups[i][:math.ceil(len(gc_groups[i])*norm_factor)])
    random.Random(random_seed).shuffle(sampled)

    # Write
    if args["output_file"] is not None:
        handle = open(args["output_file"], "wt")
    else:
        handle = sys.stdout
    SeqIO.write(sampled[:args["subsample"]], handle, "fasta")
    handle.close()

if __name__ == "__main__":
    main()
