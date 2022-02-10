#!/usr/bin/env python

from Bio import motifs
from Bio.motifs import jaspar
import click
import os

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "jaspar_file", type=click.Path(exists=True, resolve_path=True)
)
@click.argument(
    "damo_file", type=click.Path(resolve_path=True)
)

def main(**args):

    # Get motif
    with open(args["jaspar_file"]) as handle:
        m = motifs.read(handle, "jaspar")
        m.pseudocounts = jaspar.calculate_pseudocounts(m)

    # Save motif
    with open(args["damo_file"], "w") as handle:
        handle.write(f">{m.matrix_id} {m.name}" + "\n")
        for letter in "ACGT":
            s = " ".join([format(round(i, 3), ".3f") for i in m.pwm[letter]])
            handle.write(f"{letter} | {s}" + "\n")

if __name__ == "__main__":
    main()