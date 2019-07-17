#!/usr/bin/python

from optparse import OptionParser
import pandas as pd
import os, sys

##
## Find manually annoated CNV calls in a variant summary file produced by svParser (e.g. 'sample.annotated.txt')
## Outputs CNV position to specified directory
##

def getCNVs(calls):
    cnvs = []

    breakpoints = pd.read_csv(calls, delimiter="\t", index_col=False, na_filter=False)
    breakpoints = breakpoints[breakpoints['source'].str.contains('cnv', case=False)]

    for index, row in breakpoints.iterrows():
        length = round((row['bp2']-row['bp1'])/1e3, 1)
        line = list( map(str, [row['chromosome1'], row['bp1'], row['bp2'], row['type'], length, row['log2(cnv)']]) )
        line = '\t'.join(line)
        cnvs.append(line)


    # with open(calls, 'U') as sv_calls_file:
    #     for l in sv_calls_file:
    #         parts = l.rstrip().split('\t')
    #
    #         if "cnv" in parts[1].lower():
    #             try:
    #                 depth = parts[18]
    #             except IndexError:
    #                 depth = '-'
    #
    #             line = [parts[3], parts[4], parts[6], parts[2], parts[11], depth]
    #
    #             line = '\t'.join(line)
    #             cnvs.append(line)

    return(cnvs)


def printCNVs(cnvs, options):
    base_name = (os.path.basename(options.annotated_variants))
    outfile = base_name.split('_')[0] + "_" + "CNVs.txt"
    output = os.path.join(options.cnv_dir, outfile)

    print("Writing CNVs to: '%s'" % output)

    with open(output, 'w') as cnv_out:
        for l in cnvs:
            cnv_out.write("%s\n" % l)


def get_args():
    parser = OptionParser()

    parser.add_option("-f",
                      "--annotated_variants",
                      dest="annotated_variants",
                      help="annotated variants file")

    parser.add_option("-d",
                      "--cnv_dir",
                      dest="cnv_dir",
                      action="store",
                      help="Directory to write cnvs to")

    options, args = parser.parse_args()

    if not options.annotated_variants or not options.cnv_dir:
        parser.print_help()
        exit("err: Must provide both '-f' and '-d' options. Exiting")
    return options, args


def main():
        options, args = get_args()

        cnvs = getCNVs(options.annotated_variants)

        print(cnvs)

        if cnvs:
            printCNVs(cnvs, options)


if __name__ == "__main__":
    sys.exit(main())