#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import click
import copy
import json
import random
import re
import sys

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
    "-s", "--shuffle-lowercase",
    help="Shuffle lowercase letters.",
    is_flag=True
)
@click.option(
    "-o", "--output-file",
    help="Output file.  [default: STDOUT]",
    type=click.Path(writable=True, readable=False, resolve_path=True,
        allow_dash=True),
)

def main(**args):

    # Initialize
    regexp = re.compile(r"[a-z]+")

    # Group sequences by %GC content
    gc_groups = {}
    for i in range(len(args["fasta_file"])):
        fasta_file = args["fasta_file"][i]
        for record in SeqIO.parse(fasta_file, "fasta"):
            if args["filter_masked"]:
                if record.seq.upper().count("N"):
                    continue
            gc = round(GC(record.seq))
            gc_groups.setdefault(gc, [[] for i in range(len(args["fasta_file"]))])
            gc_groups[gc][i].append(record)

    # Matching
    matched = [["labels"] + list(args["fasta_file"])]
    random_seed = 123
    for i in sorted(gc_groups):
        for j in range(len(gc_groups[i])):
            random.Random(random_seed).shuffle(gc_groups[i][j])
        min_len = min(len(gc_groups[i][j]) for j in range(len(gc_groups[i])))
        for j in range(min_len):
            matched.append([i])
            for k in range(len(gc_groups[i])):
                record = gc_groups[i][k][j]
                s = str(record.seq)
                if args["shuffle_lowercase"]:
                    # 1) extract blocks of lowercase letters;
                    # 2) shuffle the letters; and
                    # 3) put the shuffled letters back
                    l = list(s)
                    for m in re.finditer(regexp, s):
                        sublist = l[m.start():m.end()]
                        random.shuffle(sublist)
                        l[m.start():m.end()] = copy.copy(sublist)
                    s = "".join(l)
                matched[-1].extend([[record.id, s]])

    # Write
    if args["output_file"] is not None:
        handle = open(args["output_file"], "wt")
    else:
        handle = sys.stdout
    json.dump(matched, handle, indent=4, sort_keys=True)
    handle.close()

if __name__ == "__main__":
    main()
