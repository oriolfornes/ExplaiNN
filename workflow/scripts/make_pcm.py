#!/usr/bin/python
#*-* coding: utf-8 *-*

""" Produces a PCM from hits file obtained with HOMER findMotifs.pl. """

import sys
import getopt


###############################################################################
#                               MAIN
###############################################################################
if __name__ == "__main__":
    usage = '''
    %s -i <input file>
    ''' % (sys.argv[0])

    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:h")
    except getopt.GetoptError:
        sys.exit(str(getopt.GetoptError) + usage)

    input_file = None
    for o, a in opts:
        if o == "-i":
            input_file = a
        else:
            sys.exit(usage)
    if not input_file:
        sys.exit(usage)

    pcm = {}
    with open(input_file) as stream:
        indx = 0
        for line in stream:
            if not line.startswith('FASTA ID'):
                indx += 1
                hit = line.split('\t')[2]
                strand = line.split('\t')[4]
                from Bio.Seq import Seq
                sequence = Seq(hit)
                if strand == '-':
                    sequence = sequence.reverse_complement()
                if indx == 1:
                    length = len(sequence)
                    pcm = {'A': [0] * length, 'C': [0] * length, 'G': [0] *
                           length, 'T': [0] * length}
                for index, char in enumerate(sequence):
                    pcm[char][index] += 1

    for letter in "ACGT":
        line = " ".join([value.__str__() for value in pcm[letter]])
        print line
