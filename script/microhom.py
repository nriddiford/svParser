#!/usr/bin/python

import pysam
import sys
import re
from difflib import SequenceMatcher

from Bio.Seq import Seq

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-p', action="store", help="Breakpoint position (chr:bp1-bp2) [Required]", dest='location', required=True)
parser.add_argument('-s', action="store", help="Split read sequence for detection of insertions", dest="split")
parser.add_argument('-n', action="store", help="Number of bases to look for extended homology. [Default: 200 if microhomology found, 10 if not]", dest='homspace',type=int)
#
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

def microhomology(seq1, seq2):
    mh = 0
    pos_count = 0
    seq = ""
    position =""
    longest_hom=0
    for i in range(len(seq1)):
        pos_count += 1
        # print("Up   %s:     %s") % (i, seq1[-i-1])
        # print("Down %s: %s") % (i, seq2[i])
        if seq1[-i:] == seq2[:i]:
            # print(seq1[-i:])
            mh += 1
            seq += seq2[i]
            position = pos_count
            longest_hom = mh
    return(position , longest_hom, seq)


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

def longestMatch(seq1, seq2):
    s = SequenceMatcher(None, seq1, seq2)
    match = s.find_longest_match(0, len(seq1), 0, len(seq2))

    block = s.get_matching_blocks()
    seq1_start = match[0]
    seq1_end = match[0]+match[2]
    seq = seq1[match[0]:(match[0]+match[2])]
    seq2_start = match[1]
    seq2_end = match[1]+match[2]

    # print("Seq1:%s\nSeq2:%s") % (seq1, seq2)
    # print(match)
    return(seq1_start, seq1_end, seq2_start, seq2_end, seq)

(chrom1, bp1, bp2) = get_parts(pos)

# bp1 -= 1
bp2 -= 1
upstream = bp1 - 100
downstream = bp2 + 100

upstream_seq = genome.fetch(chrom1, upstream, bp1)
downstream_seq = genome.fetch(chrom1, bp2, downstream)
# upstream_seq = reversed_seq(upstream_n)

upstream_seq = 'GGGGGTTTTTTTTTTCCCCC'
downstream_seq = 'CCCCCTTTTTTTTTTGGGGG'

# split_read = reversed_seq(split_read)

##################
## Microomology ##
##################

print("Upstream:   %s") % (upstream_seq)
print("Downstream: %s") % (downstream_seq)

(position, longest_hom, mhseq) =  microhomology(upstream_seq, downstream_seq)

if(longest_hom>=1):
    print("")
    print("* Microhomology at breakpoint: %s (%s bp)") % (mhseq, longest_hom)
    downstream_spacer = (len(upstream_seq) - len(mhseq))
    dmarker = " "*(downstream_spacer) + "^"*len(mhseq)
    downstream_buff = " "*downstream_spacer

    # upstream_base_end = bp1 - len(mhseq)
    # downstream_base_end = bp2 + len(mhseq)
    marker = "^"*len(mhseq)
    print(" Upstream:      %s") % (upstream_seq)
    print(" Downstream:    %s%s") % (downstream_buff, downstream_seq)
    print(" Microhomology: %s\n") % (dmarker)

else:
    print("\n* No microhomology found\n")

if args.homspace is None:
    if longest_hom<=3:
        n=10
    else:
        n=200

##############
## Homology ##
##############

if n >= len(upstream_seq):
     n = len(upstream_seq)

# upstream_seq[-n:] takes n characters from the end of the string
(downstream_start, downstream_end, upstream_start, upstream_end, seq) = longestMatch(downstream_seq[0:n],upstream_seq[-n:])

downstream_spacer = len(upstream_seq)  - n  +  upstream_start - downstream_start
print("* Longest homologous sequence +/- %s bp from breakpoint: %s (%s bp)") % (n, seq, len(seq))

seq_diff = len(upstream_seq) - len(upstream_seq[-n:])

umarker = " "*(upstream_start+seq_diff) + "^"*len(seq)
dmarker = " "*(downstream_start+downstream_spacer) + "^"*len(seq)

padded_downstream = " "*downstream_spacer + downstream_seq

print(" Upstream:      %s") % (upstream_seq)
print(" Homology:      %s") % (umarker)
print(" Downstream:    %s") % (padded_downstream)
print(" Homology:      %s\n") % (dmarker)


if args.split is not None:
    ###############
    ## Inserions ##
    ###############

    print("\n* Split read aligned to upstream sequence:")
    # Align split read to upstream seq (normal orientation)
    (upstream_start, upstream_end, split_start, split_end, upseq) = longestMatch(upstream_seq, split_read)
    print(" Upstream:      %s") % (upstream_seq[upstream_start:upstream_end])
    print(" Split read:    %s--/--%s\n") % (split_read[split_start:split_end], split_read[split_end:len(split_read)])

    bp_start = split_end # for extracting the inserted seq

    # Align split read to downstream seq (normal orientation)
    # seq = Seq(split_read)
    ###
    #split_read = seq.reverse_complement()

    # split_read = split_read[split_end:]
    (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, split_read[bp_start:])
    # Calculate length of aligned sequences
    split_len = len(split_read)
    aligned_up = len(upseq)
    aligned_down = len(downseq)
    # Split read length - aligned portion = insetion size
    insertion_size = int(split_len - aligned_up - aligned_down)
    inserted_seq = split_read[bp_start:bp_start+insertion_size]

    print("\n* Split read aligned to downstream sequence:")
    if insertion_size <= 0:
        deleted_bases = "."*(downstream_start)
        seqbuffer = " "*(5+len(split_read[bp_start:]))

        print(" Split read:    %s--/--%s") % (split_read[0:split_end], split_read[split_end:])
        print(" Downstream:    %s%s\n") % (seqbuffer, downstream_seq[0:downstream_end])

    else:
        seqbuffer = " "*(downstream_start+5+insertion_size)

        print(" Split read:    %s--/--%s") % (split_read[0:split_end-1], split_read[split_end-1:])
        print(" Downstream:    %s     %s\n") % (seqbuffer, downstream_seq[downstream_start:])

    print("* %s bp insertion '%s' at breakpoint\n") % (insertion_size, inserted_seq)

    if insertion_size >= 3:
        deletion = 0
        (inserted_start, inserted_end, upstream_start, upstream_end, aligned) = longestMatch(inserted_seq,upstream_seq)
        if len(aligned) > 3:
            splitbuffer = " "*upstream_start
            print(" Upstream:      %s--/--") % (upstream_seq)
            print(" Insertion:     %s%s\n") % (splitbuffer, inserted_seq[0:len(aligned)])
            insertion_pos = (len(upstream_seq) - len(aligned))
            print("* %s bp of inserted sequence -%s bps from breakpoint on upstream sequence\n") % (len(aligned),insertion_pos)

        (downstream_start, downstream_end, inserted_start, inserted_end, aligned) = longestMatch(downstream_seq, inserted_seq)
        if len(aligned) > 3:
            print(downstream_start)
            splitbuffer = " "*(downstream_start+5)
            print(" Downstream:     --/--%s") % (downstream_seq)
            print(" Insertion:      %s%s\n") % (splitbuffer, inserted_seq)
            print("* %s bp of inserted sequence +%s bps from breakpoint on downstream sequence\n") % (len(aligned),downstream_start)

    elif insertion_size < 0:
        insertion_size -= downstream_start

        deletion = 1
        insertion_size = abs(insertion_size)
        deleted_seq = split_read[bp_start-1:bp_start]
        print("* %s bp deletion '%s' at breakpoint\n") % (insertion_size, deleted_seq)

    else:
        deletion = 0

else:
    print("No split read sequence provided. Unable to find insetion at bp")
    insertion_size = 0

###############
## Mechanism ##
###############

# Calculate mechanism
mechanism=getMechanism(longest_hom, insertion_size)
print("* Mechanism: %s") % (mechanism)
if deletion:
    print("  * %s bp deletion") % (insertion_size)
else:
    print("  * %s bp insertion") % (insertion_size)
print("  * %s bp homology at breakpoints") % (longest_hom)



#
