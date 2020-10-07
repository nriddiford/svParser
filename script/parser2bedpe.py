#!/usr/bin/env python
import os, re, sys
import ntpath
import pandas as pd
import ntpath
from optparse import OptionParser


def get_fields(options):
    breakpoints = pd.read_csv(options.in_file, delimiter="\t", index_col=False, na_filter=False)

    breakpoints = breakpoints.loc[breakpoints['status'] == '']
    sources = ['lumpy', 'delly', 'novobreak']
    breakpoints = breakpoints.loc[breakpoints['source'].isin(sources)]

    breakpoints['type'] = breakpoints['type'] + ":" + breakpoints['event'].map(str) + ":" + breakpoints['configuration'].map(str)
    breakpoints = breakpoints.drop_duplicates(subset='type', keep="last")

    breakpoints = breakpoints.filter(items=['type', 'chromosome1', 'bp1', 'chromosome2', 'bp2'])

    return breakpoints


def write_bedpe(df, options):

    if not options.out_file:
        name = ntpath.basename(options.in_file).split("_")[0]
        options.out_file = '_'.join([name, 'breakpoints.bedpe'])

    print("Writing bedpe to '%s'" % options.out_file)

    new_order = ['chromosome1', 'bp1', 'chromosome2', 'bp2', 'type']

    df = df.reindex(columns=new_order)

    df.to_csv(options.out_file, sep="\t", index=False, header=False)


def get_args():
    parser = OptionParser()

    parser.add_option("-i",
                      "--in_file",
                      dest="in_file",
                      action="store",
                      help="Tab separated output file from svParser",)

    parser.add_option("-o",
                      "--out_file",
                      dest="out_file",
                      action="store",
                      help="bedectory to write bed files to")

    options, args = parser.parse_args()

    if not options.in_file:
        parser.print_help()
        print
    return options, args


def main():
        options, args = get_args()

        if options.in_file:
            try:
                fields = get_fields(options)
                write_bedpe(fields, options)
            except IOError as err:
                sys.stderr.write("IOError " + str(err) + "\n");
                return

if __name__ == "__main__":
    sys.exit(main())