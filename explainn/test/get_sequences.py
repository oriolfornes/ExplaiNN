#!/usr/bin/env python

from Bio import SeqIO
import click
import gzip
import os

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "fasta_file", type=click.Path(exists=True, resolve_path=True), nargs=-1
)
@click.option(
    "-o", "--output-dir",
    help="Output directory.",
    type=click.Path(resolve_path=True),
    default="./",
    show_default=True
)
@click.option(
    "-p", "--positives",
    help="Description of positive sequences.",
    required=True
)
@click.option(
    "-n", "--negatives",
    help="Description of negative sequences.",
    required=True
)

def main(**args):

    # Initialize
    positives = []
    negatives = []

    # Create output dir
    if not os.path.exists(args["output_dir"]):
        os.makedirs(args["output_dir"])

    # Get FASTA sequences
    for f in args["fasta_file"]:
        with gzip.open(f, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.description.split(" ")[1] == args["positives"]:
                    positives.append(record)
                elif record.description.split(" ")[1] == args["negatives"]:
                    negatives.append(record)

    # Save sequences
    sequences_file = os.path.join(args["output_dir"], "positives.fa")
    with open(sequences_file, "wt") as handle:
        SeqIO.write(positives, handle, "fasta")
    sequences_file = os.path.join(args["output_dir"], "negatives.fa")
    with open(sequences_file, "wt") as handle:
        SeqIO.write(negatives, handle, "fasta")

if __name__ == "__main__":
    main()