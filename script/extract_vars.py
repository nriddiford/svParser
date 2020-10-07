#!/usr/bin/env python
import os, re, sys
import pandas as pd
import ntpath
import json
from optparse import OptionParser


def assign_confidence(row):
    if row['source'] == "CNV-Seq" and row['chromosome1'] == 'X' and row['bp1'] >= 2500000 and row['bp2'] <= 3500000:
        return 'precise'
    if row['split_reads'] == '-' and row['disc_reads'] == '-':
        return 'imprecise'
    return 'precise'


def extract_vars(options):
    breakpoints = pd.read_csv(options.variants, delimiter="\t", index_col=False, na_filter=False)

    sample = ntpath.basename(options.variants).split("_")[0]

    breakpoints['sample'] = sample

    bp1 = []
    bp2 = []
    hit_genes = []

    breakpoints[['gene1', 'locus1']] = breakpoints['bp1_locus'].str.split(', ', n=2, expand=True)
    breakpoints[['gene2', 'locus2']] = breakpoints['bp2_locus'].str.split(', ', n=2, expand=True)

    breakpoints['locus1'].replace([None], 'intergenic', inplace=True)
    breakpoints['locus2'].replace([None], 'intergenic', inplace=True)
    breakpoints['microhomology'].replace(['-'], 0, inplace=True)

    breakpoints['confidence'] = breakpoints.apply(assign_confidence, axis=1)

    breakpoints.drop_duplicates(['event', 'chromosome1', 'bp1', 'chromosome2', 'bp2'], inplace=True)

    false_calls = ['F', 'aF']
    if options.normal:
        print("This sample is non-tumour. Only removing calls marked as 'F'")
        false_calls = ['F']

    for index, row in breakpoints.iterrows():
        if row['status'] not in false_calls:
            if options.write_genes:
                affected_genes = row['affected_genes'].split(", ")
                for g in affected_genes:
                    if 'DEL' in row['type'] and row['log2(cnv)'] < -0.2 or 'DUP' in row['type'] and row['log2(cnv)'] > 0.2:
                        hit_genes.append([row['event'], row['sample'], row['genotype'], row['type'], row['allele_frequency'], row['length(Kb)'], row['log2(cnv)'], row['chromosome1'], g])

            if options.write_breakpoints:
                bp1.append([row['event'], 'bp1', row['sample'], row['genotype'], row['chromosome1'], row['bp1'], row['gene1'], row['locus1'], row['chromosome2'], row['bp2'],  row['gene2'], row['locus2'], row['type'],  row['length(Kb)'], row['allele_frequency'], row['confidence'], row['microhomology'], row['mechanism']])
                bp2.append([row['event'], 'bp2', row['sample'], row['genotype'], row['chromosome2'], row['bp2'], row['gene2'], row['locus2'], row['chromosome1'], row['bp1'],  row['gene1'], row['locus1'], row['type'],  row['length(Kb)'], row['allele_frequency'], row['confidence'], row['microhomology'], row['mechanism']])

    if options.write_genes:
        genes_out = os.path.join(sample + "_hit_genes.txt")
        with open(genes_out, 'w') as genes:
            for l in hit_genes:
                genes.write('\t'.join(map(str, l)) + '\n')

    if options.write_breakpoints:
        # options.out_file = os.path.join(sample + "_" + options.out_file)
        with open(options.out_file, 'a') as extracted_vars:
            for l, m in zip(bp1, bp2):
                extracted_vars.write('\t'.join(map(str, l)) + '\n')
                extracted_vars.write('\t'.join(map(str, m)) + '\n')


def get_args():
    parser = OptionParser()

    parser.add_option("-v",
                      "--variants",
                      dest="variants",
                      action="store",
                      help="svParser format file")

    parser.add_option("--file_name",
                      dest="file_name",
                      action="store",
                      help="generic file name for output. [Default: 'microhomology.txt'")

    parser.add_option("-b",
                      "--write_breakpoints",
                      dest="write_breakpoints",
                      action="store_true",
                      help="Extract breakpoints.txt?")

    parser.add_option("-g",
                      "--write_genes",
                      dest="write_genes",
                      action="store_true",
                      help="Also extract hit_genes.txt?")

    parser.add_option("-o",
                      "--out_file",
                      dest="out_file",
                      action="store",
                      help="File to write all_breakpoints.txt to")

    parser.add_option("-n", "--normal", dest="normal", action="store_true", help="This is a normal tissue")

    parser.set_defaults(file_name='microhomology.txt')
    options, args = parser.parse_args()

    if not options.variants:
        parser.print_help()
        print
        sys.exit("[!] Must provide a variants file. Exiting.")
    else:
        return options, args


def main():
    """Extract files (all_bps.txt / all_genes.txt format) for a variants file"""
    options, args = get_args()
    try:
        extract_vars(options)
    except IOError as err:
        sys.stderr.write("IOError " + str(err) + "\n")
        return


if __name__ == "__main__":
    sys.exit(main())