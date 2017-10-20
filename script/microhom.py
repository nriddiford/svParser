#!/usr/bin/python

import pysam
import sys
import re
from difflib import SequenceMatcher

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-p', action="store", help="Breakpoint position (chr:bp1-bp2) [Required]", dest='location', required=True)
parser.add_argument('-s', action="store", help="Split read sequence for detection of insertions", dest="split")
parser.add_argument('-n', action="store", help="Number of bases to look for extended homology. [Default: 200 if microhomology found, 10 if not]", dest='homspace',type=int)

args = parser.parse_args()
pos = args.location
split_read = args.split
n=args.homspace

genome = pysam.Fastafile("/Users/Nick_curie/Documents/Curie/Data/Genomes/Dmel_v6.12/Dmel_6.12.fasta")

def reversed_seq(x):
    return x[::-1]

def get_parts(lookup):
    (chrom1, bp1, bp2) = re.split(':|-',lookup)
    return(chrom1, int(bp1), int(bp2))

def getMechanism(homlen, inslen):
    if inslen >= 10:
        mechanism="FoSTeS"
    elif (homlen <= 2) or (inslen >=1 and inslen <= 10):
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
    longest_hom=0
    for i in range(len(upstream_seq)):
        pos_count += 1
        if upstream_seq[i] == downstream_seq[i]:
            mh += 1
            seq += upstream_seq[i]
            position = pos_count
            longest_hom = mh

        else:
            mh = 0
            # longest_hom=0
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

(chrom1, bp1, bp2) = get_parts(pos)

upstream = bp1 - 200
downstream = bp2 + 200

upstream_n = genome.fetch(chrom1, upstream, bp1)
downstream_seq = genome.fetch(chrom1, bp2, downstream)
upstream_seq = reversed_seq(upstream_n)

# upstream_seq = 'ATACATTGGCCTTGGCTTAGACTTAGATCTAGACCTGAAAATAACCTGCCGAAAAGACCCGCCCGACTGTTAATACTTTACGCGAGGCTCACCTTTTTGTTGTGCTCCC'
# downstream_seq = 'ATACACGAAAAGCGTTCTTTTTTTGCCACTTTTTTTTTATGTTTCAAAACGGAAAATGTCGCCGTCGTCGGGAGAGTGCCTCCTCTTAGTTTATCAAATAAAGCTTTCG'


# split_read = reversed_seq(split_read)

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
    if args.homspace is None:
        if longest_hom<=3:
            n=10
        else:
            n=200

else:
    print("\n* No microhomology found\n")
    if args.homspace is None:
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

if n <= 10 and len(seq) >= 3:
    longest_hom=len(seq)

if args.split is not None:
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
    print(len(split_read) - aligned_up - aligned_down)

    if insertion_size >= 1 and insertion_size < 10:
        inserted_seq = split_read[bp_start:bp_start+insertion_size]
        print("* %s bp insertion '%s' at breakpoint\n") % (insertion_size, inserted_seq)
    elif insertion_size >= 10:
        inserted_seq = split_read[bp_start:bp_start+insertion_size]
        print("* %s bp insertion '%s' at breakpoint\n") % (insertion_size, inserted_seq)

        (upstream_start, upstream_end, inserted_start, inserted_end, aligned) = longestMatch(upstream_n, inserted_seq)
        if len(aligned) > insertion_size/3:
            print(" Upstream:      %s") % (upstream_n[upstream_start:upstream_end])
            print(" Insertion:     %s\n") % (inserted_seq)
        (downstream_start, downstream_end, inserted_start, inserted_end, aligned) = longestMatch(downstream_seq, inserted_seq)
        if len(aligned) > insertion_size/3:
            print(" Downstream:      %s--/--%s") % (downstream_seq[downstream_start:downstream_end], downstream_seq[downstream_end:len(downstream_seq)])
            print(" Insertion:       %s\n") % (inserted_seq)
            print("%s bp of inserted sequence found on downstream sequence") % (len(aligned))


else:
    print("No split read sequence provided. Unable to find insetion at bp")
    insertion_size = 0

###############
## Mechanism ##
###############

# Calculate mechanism
mechanism=getMechanism(longest_hom, insertion_size)
print("* Mechanism: %s") % (mechanism)
print("  * %s bp insertion") % (insertion_size)
print("  * %s bp homology at breakpoints") % (longest_hom)










#
