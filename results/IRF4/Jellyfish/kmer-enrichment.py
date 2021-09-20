#!/usr/bin/env python

from Bio import SeqIO

wt = "./HT-SELEX/WT.mer_counts.fa"
t95r = "./HT-SELEX/T95R.mer_counts.fa"

kmer_counts = {}

for seq_record in SeqIO.parse(wt, "fasta"):
    kmer = sorted([str(seq_record.seq), str(seq_record.seq.complement())])[0]
    kmer_counts.setdefault(kmer, [0, 0])
    kmer_counts[kmer][0] += int(seq_record.id)

for seq_record in SeqIO.parse(t95r, "fasta"):
    kmer = sorted([str(seq_record.seq), str(seq_record.seq.complement())])[0]
    kmer_counts.setdefault(kmer, [0, 0])
    kmer_counts[kmer][1] += int(seq_record.id)

for kmer in sorted(kmer_counts):
    print("%s\t%s\t%s" % (kmer, kmer_counts[kmer][0], kmer_counts[kmer][1]))
