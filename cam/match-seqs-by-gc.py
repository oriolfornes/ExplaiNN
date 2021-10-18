#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqUtils import GC
import click
import json
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
    nargs=-1,
)
@click.option(
    "-f", "--filter-masked",
    help="Filter masked DNA sequences.",
    is_flag=True,
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
    for i in range(len(params["fasta_file"])):
        fasta_file = params["fasta_file"][i]
        for record in SeqIO.parse(fasta_file, "fasta"):
            if params["filter_masked"]:
                if record.seq.count("N") or record.seq.count("n"):
                    continue
            gc = round(GC(record.seq))
            gc_groups.setdefault(gc, [[] for i in range(len(params["fasta_file"]))])
            gc_groups[gc][i].append(record)

    # Matching
    matched = []
    random_seed = 123
    for i in sorted(gc_groups):
        for j in range(len(gc_groups[i])):
            random.Random(random_seed).shuffle(gc_groups[i][j])
        min_len = min(len(gc_groups[i][j]) for j in range(len(gc_groups[i])))
        for j in range(min_len):
            matched.append([i])
            for k in range(len(gc_groups[i])):
                record = gc_groups[i][k][j]
                matched[-1].extend([record.id, str(record.seq)])

    # Write
    if params["output_file"] is not None:
        handle = __get_handle(params["output_file"], "wt")
    else:
        handle = sys.stdout
    json.dump(matched, handle, indent=4, sort_keys=True)
    handle.close()

if __name__ == "__main__":
    cli()
