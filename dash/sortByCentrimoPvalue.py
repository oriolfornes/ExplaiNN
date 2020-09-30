#!/usr/bin/env python

import argparse
from decimal import Decimal
import os

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
    parser.add_argument("centrimo_dir", metavar="centrimo-dir")

    return(parser.parse_args())

def get_sorted_motifs(centrimo_dir):

    # Initialize
    motifs = []
    pvalues = []

    for d in os.listdir(centrimo_dir):
        if os.path.isdir(os.path.join(centrimo_dir, d)):
            centrimo_file = os.path.join(centrimo_dir, d, "centrimo.tsv")
            with open(centrimo_file, "r") as f:
                for line in f:
                    if line.startswith("   1"):
                        line = line.split("\t")
                        motifs.append(line[1])
                        pvalues.append(Decimal(line[5]))

    return(list(sorted(zip(motifs, pvalues), key=lambda x: x[1])))

def main():

    # Parse arguments
    args = parse_args()

    # Sort motifs by centrimo p-value
    for motif, pvalue in get_sorted_motifs(args.centrimo_dir):
        print(motif, str(pvalue))

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()