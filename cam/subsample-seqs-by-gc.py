#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqUtils import GC
import click
import math
import random
import sys

# Local imports
from train import __get_handle

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

def cli(**params):

    # Group sequences by %GC content
    gc_groups = {}
    for record in SeqIO.parse(params["fasta_file"], "fasta"):
        if params["filter_masked"]:
            if record.seq.count("N") or record.seq.count("n"):
                continue
        gc = round(GC(record.seq))
        gc_groups.setdefault(gc, [])
        gc_groups[gc].append(record)

    # Subsampling
    sampled = []
    random_seed = 123
    norm_factor = params["subsample"] / \
        sum([len(v) for v in gc_groups.values()])
    for i in sorted(gc_groups):
        random.Random(random_seed).shuffle(gc_groups[i])
        sampled.extend(gc_groups[i][:math.ceil(len(gc_groups[i])*norm_factor)])
    random.Random(random_seed).shuffle(sampled)

    # Write
    if params["output_file"] is not None:
        handle = __get_handle(params["output_file"], "wt")
    else:
        handle = sys.stdout
    SeqIO.write(sampled[:params["subsample"]], handle, "fasta")
    handle.close()

if __name__ == "__main__":
    cli()
