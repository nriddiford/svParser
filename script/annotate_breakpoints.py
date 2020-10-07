#!/usr/bin/env python
import os, re, sys
import ntpath
import pandas as pd
import numpy as np
import ntpath
from optparse import OptionParser
# import xlwt
from difflib import SequenceMatcher
import pysam


def parse_svParser(options):
    breakpoints = pd.read_csv(options.variants, delimiter="\t", index_col=False, na_filter=False)
    breakpoints['type'] = breakpoints['type'] + ":" + breakpoints['event'].map(str) + ":" + breakpoints['configuration'].map(str)

    return breakpoints


def get_ref_seq(c, start, end, options):
    genome = pysam.Fastafile(options.genome)
    reference_seq = genome.fetch(c, start, end)

    return reference_seq



def parse_splitvision(options):
    annotated = pd.read_excel(options.annotation, sheet_name=0)

    annotated = annotated.rename(index=str, columns={"variant_type": "type",
                                                     "sampleID": "sample",
                                                     "breakpoint_microhomology(sequence)": "microhomology",
                                                     "insertions(sequence)": "inserted_seq"})

    annotated = annotated.fillna('-')

    for index, row in annotated.iterrows():
        c1 = row['ChrA']
        c2 = row['ChrB']
        b1 = row['PosA']
        b2 = row['PosB']

        mh, ins, upstream_seq, downstream_seq = row[['microhomology', 'inserted_seq', 'regionA_sequence', 'regionB_sequence']]
        upstream_ref = get_ref_seq(c1, b1-100, b1, options)
        downstream_ref = get_ref_seq(c1, b2, b2+100, options)

        print(row['type'])
        mh_len = len(mh)
        if mh == '-':
            mh_len = 0
        ins_len = len(ins)
        if ins == '-':
            ins_len = 0

        templated_up = 0
        templated_down = 0
        if ins_len >= 3:
            print("Insered seq %s (length: %s)" % (ins, ins_len))
            templated_up, t_up_seq = templated_search(ins, upstream_ref, 'up')
            templated_down, t_down_seq = templated_search(ins, downstream_ref, 'down')

        templated_length = max(templated_up, templated_down)
        annotated.loc[index, 'mechanism'] = getMechanism(mh_len, ins_len, templated_length)
        mh_string = "mh:%s" % mh_len
        ins_string = "ins:%s" % ins_len
        temp_up_string = "tup:%s" % templated_up
        temp_down_string = "tdown:%s" % templated_down

        annotated.loc[index, 'notes'] = '; '.join([mh_string, ins_string, temp_up_string, temp_down_string])

        # print("mh_seq: %s [%s]") % (mh, mh_len)
        # print("inserted seq: %s [%s]") % (ins, ins_len)
        print("* mechanism: %s" % annotated.loc[index, 'mechanism'])

    annotated = annotated.filter(items=['type', 'microhomology', 'mechanism', 'inserted_seq', 'contig_sequence', 'notes'])
    annotated = annotated.fillna('-')

    return annotated


def getMechanism(homlen, inslen, templen):
    if homlen >= 20:
        mechanism = "NAHR"
    elif inslen >= 5 or templen >= 5:
        mechanism = "FoSTeS"
    elif homlen >= 1 and inslen <= 5:
        mechanism = "Alt-EJ"
    else:
        mechanism = "NHEJ"

    return mechanism


def longestMatch(seq1, seq2):
    s = SequenceMatcher(None, seq1, seq2)
    match = s.find_longest_match(0, len(seq1), 0, len(seq2))

    block = s.get_matching_blocks()
    seq1_start = match[0]
    seq1_end = match[0]+match[2]
    seq = seq1[match[0]:(match[0]+match[2])]
    seq2_start = match[1]
    seq2_end = match[1]+match[2]

    return seq1_start, seq1_end, seq2_start, seq2_end, seq


def templated_search(inserted_seq, reference, direction):
    (inserted_start, inserted_end, reference_start, reference_end, aligned) = longestMatch(inserted_seq, reference)
    if len(aligned) >= 2:
        templated = aligned
        templated_insertion_size = len(templated)

        if direction == 'up':
            insertion_pos = (len(reference) - reference_end)
            print("* %s bp of inserted sequence -%s bps from breakpoint on %s sequence") % (len(templated), insertion_pos, direction)
            splitbuffer = " "*reference_start
            print(" Upstream:      %s--/--") % (reference)
            print(" Insertion:     %s%s\n") % (splitbuffer, templated)

        else:
            insertion_pos = (reference_start)
            print "* %s bp of inserted sequence +%s bps from breakpoint on %s sequence" % (len(templated), insertion_pos, direction)
            splitbuffer = " "*(reference_start+5)
            print(" Downstream:     --/--%s") % (reference)
            print(" Insertion:      %s%s\n") % (splitbuffer, templated)
    else:
        templated_insertion_size = 0
        templated = ''

    return templated_insertion_size, templated


def annotate_vars(options, vars, ann):
    if not options.out_file:
        sample = ntpath.basename(options.variants).split("_")[0]
        options.out_file = os.path.join(sample + '_microhomology.txt')

    if ann.empty:
        print("No annotations")
        annotated = vars
        annotated = annotated.assign(inserted_seq = np.nan,
                         contig_sequence = np.nan)

    else:
        cols = ['microhomology']
        vars.drop(cols, inplace=True, axis=1)

        annotated = pd.merge(vars, ann,  how='left', on='type')
        annotated = annotated.fillna('-')

        for index, row in annotated.iterrows():
            notes, new_notes = row[['notes_x', 'notes_y']]
            if not notes:
                notes = new_notes
            elif not new_notes == '-':
                notes = '; '.join([notes, new_notes])
            annotated.loc[index, 'notes'] = notes

        annotated = annotated.rename(index=str, columns={"microhomology_y": "microhomology",
                                                          "mechanism_y": "mechanism"})

    annotated['type'] = annotated['type'].str.partition(':')[[0,2]].rename({0: 'type', 2: 'event'}, axis=1)

    new_order = ['event', 'source', 'type', 'chromosome1', 'bp1', 'chromosome2', 'bp2', 'split_reads', 'disc_reads','genotype',
                 'id', 'length(Kb)', 'position', 'configuration', 'allele_frequency', 'log2(cnv)', 'microhomology', 'inserted_seq',
                 'consensus', 'contig_sequence', 'mechanism', 'bp1_locus', 'bp2_locus', 'affected_genes', 'status', 'notes']

    annotated = annotated[new_order]
    annotated.to_csv(options.out_file, sep="\t", index=False)

    return True


def get_args():
    parser = OptionParser()

    parser.add_option("-v",
                      "--variants",
                      dest="variants",
                      action="store",
                      help="svParser format file")

    parser.add_option("-a",
                      "--annotations",
                      dest="annotation",
                      action="store",
                      help="Breakpoints annotated by SplitVision")

    parser.add_option("-g",
                      "--genome",
                      dest="genome",
                      action="store",
                      help="Genome Fasta file")

    parser.add_option("-o",
                      "--out_file",
                      dest="out_file",
                      action="store",
                      help="File to write annotated variants to")

    parser.set_defaults(genome = "/Users/Nick_curie/Documents/Curie/Data/Genomes/Dmel_v6.12/Dmel_6.12.fasta")

    options, args = parser.parse_args()

    if not options.variants or not options.annotation:
        parser.print_help()
        print
        sys.exit("[!] Both a variants file and annotation file are required. Exiting.")
    else:
        return options, args


def main():
        options, args = get_args()

        try:
            vars = parse_svParser(options)
            ann = parse_splitvision(options)
            annotate_vars(options, vars, ann)

        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n")
            return


if __name__ == "__main__":
    sys.exit(main())