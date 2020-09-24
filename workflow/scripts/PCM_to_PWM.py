#!/usr/bin/env python

import argparse
from Bio import motifs

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

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # From https://biopython.readthedocs.io/en/latest/chapter_motifs.html
    m = motifs.read(open(args.motif_file), "jaspar")
    m.pseudocounts = motifs.jaspar.calculate_pseudocounts(m)

    # Print PWM
    for value in m.pssm.values():
        print(" ".join([str(round(v, 8)).rjust(11) for v in value]))

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()