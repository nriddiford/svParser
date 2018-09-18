#!/usr/bin/env python
import sys, os
import pandas as pd
from pprint import pprint

from optparse import OptionParser

def get_file_name(variants_file):
    base_name = (os.path.splitext(variants_file)[0])
    out_base = base_name.split('_')[0]
    outfile = out_base + "_" + "stitched_SVs.txt"
    return(outfile)

def parse_vars(options):
    out_file = get_file_name(options.inFile)

    print("Stitching together variants for: %s" % options.inFile)
    print("Writing stitched variants to: %s" % out_file)

    df = pd.read_csv(options.inFile, delimiter="\t")

    events = []

    df = df[['event', 'chromosome1', 'bp1', 'chromosome2', 'bp2']]

    df_dict = df.to_dict()

    for i in df_dict:
        print i

    return(df_dict)

    # for i in df.index:
    #     event = df.loc[i, 'event']
    #     c1 = df.loc[i, 'chromosome1']
    #     b1 = int(df.loc[i, 'bp1'])
    #     c2 = df.loc[i, 'chromosome2']
    #     b2 = int(df.loc[i, 'bp2'])
    #
    #     e = Events(event, c1, b1, c2, b2)
    #
    # events.append(e.gatherEvents())
    #
    # return(events)


def stitch(options):
    variants = parse_vars(options)
    pprint(variants)


class Events(object):
    def __init__(self, event, c1, b1, c2, b2):
        self.event = event
        self.c1 = c1
        self.b1 = b1
        self.c2 = c2
        self.b2 = b2


    def gatherEvents(self):
        events = {}
        events[self.event] = [self.c1, self.b1, self.c2, self.b2]

        return events


def main():
    parser = OptionParser()

    parser.add_option("-i",
        "--inFile",
        dest="inFile",
        help="An annotated variants file produced by sv2gene " + \
             "accepts both '_annotated_SVs.txt' and " + \
             "'_reannotated_SVs.txt' files",
        metavar="FILE")

    options, args = parser.parse_args()

    if options.inFile is None:
        parser.print_help()
        print
    else:
        try:
            stitch(options)

        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n");
            return

if __name__ == "__main__":
    sys.exit(main())
