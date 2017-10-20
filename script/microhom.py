#!/usr/bin/python

import pysam
import sys
import re
from difflib import SequenceMatcher

genome = pysam.Fastafile("/Users/Nick_curie/Documents/Curie/Data/Genomes/Dmel_v6.12/Dmel_6.12.fasta")

def reversed_seq(x):
    return x[::-1]

def get_parts(lookup):
    (chrom1, bp1, bp2) = re.split(':|-',lookup)
    bp2 = (int(bp2) -1)
    return(chrom1, int(bp1), int(bp2))

def getMechanism(homlen, inslen):
    if inslen >= 10:
        mechanism="false_positives"
    elif homlen <= 2 or inslen >=1 and inslen <= 10:
        mechanism="NHEJ"
    elif homlen >= 3 and homlen <= 100:
        mechanism="Alt-EJ"
    elif homlen >= 100:
        mechanism="NAHR"
    else:
        mechanism="NA"

    return(mechanism)

def microhomology(seq1, seq2):
    mh = 0
    pos_count = 0
    seq = ""
    position =""
    longest_hom=""
    for i in range(len(upstream_seq)):
        pos_count += 1
        if upstream_seq[i] == downstream_seq[i]:
            mh += 1
            seq += upstream_seq[i]
            position = pos_count
            longest_hom = mh

        else:
            mh = 0
            longest_hom=0
            break
    return(position , longest_hom, seq)

def longestMatch(seq1, seq2):
    s = SequenceMatcher(None, seq1, seq2)
    match = s.find_longest_match(0, len(seq1), 0, len(seq2))
    upstream_start = match[0]
    upstream_end = match[0]+match[2]
    seq = upstream_seq[match[0]:(match[0]+match[2])]
    downstream_start = match[1]
    downstream_end = match[1]+match[2]
    return(upstream_start, upstream_end, downstream_start, downstream_end, seq)


pos = sys.argv[1]
split_read = sys.argv[2]

(chrom1, bp1, bp2) = get_parts(pos)

upstream = bp1 - 200
downstream = bp2 + 200

upstream_n = genome.fetch(chrom1, upstream, bp1)
downstream_seq = genome.fetch(chrom1, bp2, downstream)
upstream_seq = reversed_seq(upstream_n)


# upstream_seq = 'ATACATTGGCCTTGGCTTAGACTTAGATCTAGACCTGAAAATAACCTGCCGAAAAGACCCGCCCGACTGTTAATACTTTACGCGAGGCTCACCTTTTTGTTGTGCTCCC'
# downstream_seq = 'ATACACGAAAAGCGTTCTTTTTTTGCCACTTTTTTTTTATGTTTCAAAACGGAAAATGTCGCCGTCGTCGGGAGAGTGCCTCCTCTTAGTTTATCAAATAAAGCTTTCG'

##################
## Microomology ##
##################

(position, longest_hom, mhseq) =  microhomology(upstream_seq, downstream_seq)

if(longest_hom>=1):
    print("")
    print("* Microhomology at breakpoint: %s (%s bp)") % (mhseq, longest_hom)

    upstream_base_end = bp1 - len(mhseq)
    downstream_base_end = bp2 + len(mhseq)
    marker = "^"*len(mhseq)
    print(" Upstream:      %s") % (upstream_seq)
    print(" Downstream:    %s") % (downstream_seq)
    print(" Microhomology: %s") % (marker)
    n=200

else:
    print("\n* No microhomology found\n")
    n=10


##############
## Homology ##
##############

if n >= len(upstream_seq):
    n -= 1

(upstream_start, upstream_end, downstream_start, downstream_end, seq) = longestMatch(upstream_seq[0:n], downstream_seq[0:n])

print("* Longest homologous sequence +/- %s bp from breakpoint: %s (%s bp)") % (n, seq, len(seq))
umarker = " "*upstream_start + "^"*len(seq)
dmarker = " "*downstream_start + "^"*len(seq)
print(" Upstream:      %s") % (upstream_seq)
print(" Homology:      %s") % (umarker)
print(" Downstream:    %s") % (downstream_seq)
print(" Homology:      %s\n") % (dmarker)
# print("Upstream: %s:%s-%s") % (chrom1, upstream_start, upstream_end)
# print("Downstream: %s:%s-%s") % (chrom1, downstream_start, downstream_end)

###############
## Inserions ##
###############

# Align split read to upstream seq (normal orientation)
(upstream_start, upstream_end, split_start, split_end, upseq) = longestMatch(upstream_n, split_read)
print(" Upstream:      %s") % (upstream_n[upstream_start:upstream_end])
print(" Split read:    %s--/--%s\n") % (split_read[split_start:split_end], split_read[split_end:len(split_read)])

bp_start = split_end # for extracting the inserted seq

# Align split read to downstream seq (normal orientation)
(downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, split_read)

# Calculate length of aligned sequences
aligned_up = len(upseq)
aligned_down = len(downseq)
seqbuffer = " "*len(split_read[0:split_start])

print(" Split read:    %s--/--%s") % (split_read[0:split_start], split_read[split_start:split_end])
print(" Downstream:    %s     %s\n") % (seqbuffer, downstream_seq[downstream_start:downstream_end])

# Split read length - aligned portion = insetion size
insertion_size = len(split_read) - aligned_up - aligned_down
print("* %s bp insertion '%s' at breakpoint\n") % (insertion_size, split_read[bp_start:bp_start+insertion_size])


###############
## Mechanism ##
###############

# Calculate mechanism
mechanism=getMechanism(longest_hom, insertion_size)
print("* Mechanism: %s\n") % (mechanism)









#
