#!/usr/bin/env python

import argparse
from Bio import motifs
from math import log2, sqrt
import numpy as np

# Example of PFM:
# A  [ 7  6  0 20  0  0  0  0  5  5 ]
# C  [ 6  5 20  0 20 20 20 20  4  7 ]
# G  [ 2  3  0  0  0  0  0  0  1  3 ]
# T  [ 5  6  0  0  0  0  0  0 10  5 ]

# Example of PSSM:
#  0.40806226   0.2184107  -2.4521041   1.7873355  -2.4521041  -2.4521041  ...
#   0.2184107           0   1.7873355  -2.4521041   1.7873355   1.7873355  ...
# -0.97243147 -0.57111238  -2.4521041  -2.4521041  -2.4521041  -2.4521041  ...
#           0   0.2184107  -2.4521041  -2.4521041  -2.4521041  -2.4521041  ...

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns an
    {argparse} object.
    """

    # Initialize
    parser = argparse.ArgumentParser()

    # Mandatory args
    parser.add_argument("motif_file", metavar="motif.jaspar")

    # Optional args
    bgp = .25
    parser.add_argument(
        "-A", default=bgp, metavar="FLOAT",
        help="background prob for A (default: %s)" % bgp
    )
    parser.add_argument(
        "-C", default=bgp, metavar="FLOAT",
        help="background prob for C (default: %s)" % bgp
    )
    parser.add_argument(
        "-G", default=bgp, metavar="FLOAT",
        help="background prob for G (default: %s)" % bgp
    )
    parser.add_argument(
        "-T", default=bgp, metavar="FLOAT",
        help="background prob for T (default: %s)" % bgp
    )

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # From TFBS module for TFBS::Matrix::PFM
    # sub to_PWM  {
    #     my ($self, %args) = @_;
    #     my $bg = ($args{'-bg_probabilities' } || $self->{'bg_probabilities'});
    #     my $bg_pdl = 
    #     transpose pdl ($bg->{'A'}, $bg->{'C'}, $bg->{'G'}, $bg->{'T'});
    #     my $nseqs = $self->pdl_matrix->sum / $self->length;
    #     my $q_pdl = ($self->pdl_matrix +$bg_pdl*sqrt($nseqs))
    #         / 
    #         ($nseqs + sqrt($nseqs));
    #     my $pwm_pdl = log2(4*$q_pdl);
    #     my $PWM = TFBS::Matrix::PWM->new
    #     ( (map {("-$_", $self->{$_}) } keys %$self),
    #         # do not want tags to point to the same arrayref as in $self:
    #     -tags => \%{ $self->{'tags'}}, 
    #     -bg_probabilities => \%{ $self->{'bg_probabilities'}}, 
    #     -matrix    => $pwm_pdl
    #     );
    #     return $PWM;
    # }
    motif = motifs.read(open(args.motif_file), "jaspar")
    pfm = np.array(list(motif.counts.values()))
    bg = np.array([[args.A], [args.C], [args.G], [args.T]])
    n = np.sum(pfm, axis=0)[0]
    pwm = 4 * (pfm + bg * sqrt(n)) / (n + sqrt(n))
    for row in pwm:
        print(" ".join([str(round(x, 8)).rjust(11) for x in map(log2, row)]))

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()