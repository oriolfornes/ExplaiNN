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

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # From https://biopython.readthedocs.io/en/latest/chapter_motifs.html
    m = motifs.read(open(args.motif_file), "jaspar")
    m.pseudocounts = motifs.jaspar.calculate_pseudocounts(m)

    # Print PWM
    for row in np.transpose(np.array(m.pssm.values())):
        print(" ".join([str(round(e, 8)).rjust(11) for e in row]))

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()