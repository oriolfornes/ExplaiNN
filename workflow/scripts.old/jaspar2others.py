#!/usr/bin/env python

import argparse
from Bio import motifs
import numpy as np

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

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # From https://biopython.readthedocs.io/en/latest/chapter_motifs.html
    m = motifs.read(open(args.motif_file), "jaspar")
    if args.format != "pcm":
        m.pseudocounts = motifs.jaspar.calculate_pseudocounts(m)

    # Print
    if args.format == "meme":
        name = m.name
        consensus = m.pssm.consensus
        w = len(consensus)
        nsites = int(sum([m.counts[n][0] for n in "ACGT"]))
        print("MEME version 4")
        print()
        print("ALPHABET= ACGT")
        print()
        print("strands: + -")
        print()
        print("Background letter frequencies (from uniform background):")
        print("A 0.25000 C 0.25000 G 0.25000 T 0.25000")
        print()
        print("MOTIF %s %s" % (name, consensus))
        print()
        print(
            "letter-probability matrix: alength= 4 w= %s nsites= %s E= 0" \
            % (w, nsites)
        )
        for row in np.transpose(np.array(list(m.pwm.values()))):
            print(" ".join([str(round(r, 8)).rjust(11) for r in row]))
    elif args.format == "pcm":
        for row in np.transpose(np.array(list(m.counts.values()))):
            print("\t".join(list(map(str, list(map(int, row))))))
    else:
        for row in np.transpose(np.array(list(m.pssm.values()))):
            print(" ".join([str(round(r, 8)).rjust(11) for r in row]))

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()