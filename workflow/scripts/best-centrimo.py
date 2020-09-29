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

def main():

    # Parse arguments
    args = parse_args()

    # Initialize
    motif = ""
    pvalue = 1.

    for d in os.listdir(args.centrimo_dir):
        if os.path.isdir(os.path.join(args.centrimo_dir, d)):
            centrimo_file = os.path.join(args.centrimo_dir, d, "centrimo.tsv")
            with open(centrimo_file, "r") as f:
                for line in f:
                    if line.startswith("   1"):
                        line = line.split("\t")
                        d = Decimal(line[5])
                        if d < pvalue:
                            motif = line[1]
                            pvalue = d

    print("%s\t%s" % (motif, str(pvalue)))

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()