#!/usr/bin/env python

import argparse
from Bio import motifs
import numpy as np
import sys

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    Parses arguments provided through the command line.
    """

    # Initialize
    parser = argparse.ArgumentParser()

    # Mandatory args
    parser.add_argument("motif_file", metavar="motif.jaspar")
    parser.add_argument(
        "--format",
        choices=["meme", "pcm", "pssm"],
        default="pssm",
        help="output format (default: pssm)"
    )
    parser.add_argument(
        "--out-file",
        help="output file (default: stdout)"
    )

    return(parser.parse_args())

def reformat_jaspar_motif(motif_file, oformat, out_file):

    # Initialize
    l = []

    # From https://biopython.readthedocs.io/en/latest/chapter_motifs.html
    for m in motifs.parse(open(motif_file), "jaspar"):
        l.append(m)

    reformat_motifs(l, oformat, out_file)

def reformat_motifs(l, oformat, out_file):

    # Initialize
    s = ""

    # Print
    if oformat == "meme":
        s += "MEME version 4\n\n"
        s += "ALPHABET= ACGT\n\n"
        s += "strands: + -\n\n"
        s += "Background letter frequencies (from uniform background):\n"
        s += "A 0.25000 C 0.25000 G 0.25000 T 0.25000\n"
    for m in l:
        if oformat == "meme":
            name = m.name
            consensus = m.pssm.consensus
            w = len(consensus)
            nsites = int(sum([m.counts[n][0] for n in "ACGT"]))
            s += "\nMOTIF %s %s\n" % (name, consensus)
            s += "letter-probability matrix: alength= 4 w= %s nsites= %s E= 0\n" % \
                (w, nsites)
            for row in np.transpose(np.array(list(m.pwm.values()))):
                s+= " ".join([str(round(r, 8)).rjust(11) for r in row]) + "\n"
        elif oformat == "pcm":
            m.pseudocounts = motifs.jaspar.calculate_pseudocounts(m)
            for row in np.transpose(np.array(list(m.counts.values()))):
                s+= "\t".join(list(map(str, list(map(int, row))))) + "\n"
        else:
            for row in np.transpose(np.array(list(m.pssm.values()))):
                s+= " ".join([str(round(r, 8)).rjust(11) for r in row]) + "\n"

    # Write
    if out_file is None:
        sys.stdout.write(s)
    else:
        with open(out_file, "w") as handle:
            handle.write(s)

def main():

    # Parse arguments
    args = parse_args()

    reformat_jaspar_motif(args.motif_file, args.format, args.out_file)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()