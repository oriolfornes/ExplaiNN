#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqUtils import GC
import click
import pandas as pd
import random

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.option(
    "-a",
    type=click.Path(exists=True, resolve_path=True),
    multiple=True,
    help="FASTA file A.",
    required=True
)
@click.option(
    "-b",
    type=click.Path(exists=True, resolve_path=True),
    multiple=True,
    help="FASTA file B.",
    required=True
)
@click.option(
    "-f", "--filter-masked",
    help="Filter masked DNA sequences.",
    is_flag=True
)
@click.option(
    "-o", "--output-file",
    help="Output file.  [default: STDOUT]",
    type=click.Path(writable=True, readable=False, resolve_path=True,
        allow_dash=True),
    default="-",
)

def cli(**params):

    # Group sequences by %GC content
    fasta_files = ["a", "b"]
    gc_groups = {}
    for i in range(len(fasta_files)):
        for fasta_file in params[fasta_files[i]]:
            for record in SeqIO.parse(fasta_file, "fasta"):
                if params["filter_masked"]:
                    if "N" in record.seq:
                        continue
                gc = round(GC(record.seq))
                gc_groups.setdefault(gc, [[], []])
                gc_groups[gc][i].append(record)

    # Downsampling
    sampled = []
    random_seed = 123
    for i in sorted(gc_groups):
        random.Random(random_seed).shuffle(gc_groups[i][0])
        random.Random(random_seed).shuffle(gc_groups[i][1])
        for j in range(min([len(gc_groups[i][0]), len(gc_groups[i][1])])):
            a = gc_groups[i][0][j]
            b = gc_groups[i][1][j]
            sampled.append([i, a.id, str(a.seq), b.id, str(b.seq)])

    # Write
    column_names = ["%GC", "ASeqId", "ASeq", "BSeqId", "BSeq"]
    df = pd.DataFrame(sampled, columns=column_names)
    with click.open_file(params["output_file"], "wt") as f:
        df.to_csv(f, sep="\t")

if __name__ == "__main__":
    cli()
