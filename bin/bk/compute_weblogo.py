#!/usr/bin/python
#*-* coding: utf-8 *-*

""" """

import sys
import getopt


###############################################################################
#                               MAIN
###############################################################################
if __name__ == "__main__":
    usage = '''
    %s -j <jaspar file> -o <output file name>
    ''' % (sys.argv[0])

    try:
        opts, args = getopt.getopt(sys.argv[1:], "j:o:h")
    except getopt.GetoptError:
        sys.exit(str(getopt.GetoptError) + usage)

    jaspar_file = None
    output = None
    for o, a in opts:
        if o == '-j':
            jaspar_file = a
        elif o == '-o':
            output = a
        else:
            sys.exit(usage)
    if not(jaspar_file or output):
        sys.exit(usage)

    from Bio import motifs
    with open(jaspar_file) as stream:
        motif = motifs.read(stream, 'jaspar')
        motif.weblogo(output)
