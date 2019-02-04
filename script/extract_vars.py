#!/usr/bin/env python
import os, re, sys
import pandas as pd
import ntpath
from optparse import OptionParser


def assign_confidence (row):
   if row['split_reads'] == '-' and row['disc_reads'] == '-':
       return 'imprecise'
   else:
       return 'precise'


def extract_vars(options):
    breakpoints = pd.read_csv(options.variants, delimiter="\t", index_col=False, na_filter=False)

    sample = ntpath.basename(options.variants).split("_")[0]

    breakpoints['sample'] = sample

    bp1 = []
    bp2 = []

    breakpoints[['gene1', 'locus1']] = breakpoints['bp1_locus'].str.split(', ', n=2, expand=True)
    breakpoints[['gene2', 'locus2']] = breakpoints['bp2_locus'].str.split(', ', n=2, expand=True)

    breakpoints['locus1'].replace([None], 'intergenic', inplace=True)
    breakpoints['locus2'].replace([None], 'intergenic', inplace=True)
    breakpoints['microhomology'].replace(['-'], 0, inplace=True)

    breakpoints['confidence'] = breakpoints.apply(assign_confidence, axis=1)

    for index, row in breakpoints.iterrows():
        if row['status'] != 'F':
            bp1.append([row['event'], 'bp1', row['sample'], row['genotype'], row['chromosome1'], row['bp1'], row['gene1'], row['locus1'], row['chromosome2'], row['bp2'],  row['gene2'], row['locus2'], row['type'],  row['length(Kb)'], row['allele_frequency'], row['confidence'], row['microhomology'], row['mechanism']])
            bp2.append([row['event'], 'bp2', row['sample'], row['genotype'], row['chromosome2'], row['bp2'], row['gene2'], row['locus2'], row['chromosome1'], row['bp1'],  row['gene1'], row['locus1'], row['type'],  row['length(Kb)'], row['allele_frequency'], row['confidence'], row['microhomology'], row['mechanism']])


    if not options.out_file:
        sample = ntpath.basename(options.variants).split("_")[0]
        options.out_file = os.path.join(sample + '_microhomology.txt')

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

    parser.add_option("-o",
                      "--out_file",
                      dest="out_file",
                      action="store",
                      help="File to write all_breakpoints.txt to")

    options, args = parser.parse_args()

    if not options.variants:
        parser.print_help()
        print
        sys.exit("[!] Must provide a variants file. Exiting.")
    else:
        return options, args


def main():
        options, args = get_args()

        try:
            extract_vars(options)
        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n")
            return


if __name__ == "__main__":
    sys.exit(main())