#!/usr/bin/env python
import sys, os
import pandas as pd
from pprint import pprint

from optparse import OptionParser
from collections import defaultdict
import json


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

    right_end = defaultdict(lambda: defaultdict(dict))
    left_end = defaultdict(lambda: defaultdict(dict))
    vars = {}

    df = df[['event', 'chromosome1', 'bp1', 'chromosome2', 'bp2']]

    df = df.sort_values(['event', 'chromosome1', 'bp1', 'chromosome2', 'bp2'])

    seen_l = defaultdict(lambda: defaultdict(dict))
    seen_r = defaultdict(lambda: defaultdict(dict))

    # TODO : Change back to real data names
    for index in df.index:
        event = df.loc[index, 'event']
        c1 = df.loc[index, 'chromosome1']
        b1 = int(df.loc[index, 'bp1'])
        c2 = df.loc[index, 'chromosome2']
        b2 = int(df.loc[index, 'bp2'])

        for i in range(b1-10, b1+10):
            for j in range(b2-10, b2+10):
                if i in seen_l[c1] and j in seen_r[c2]:

                    index = sameIndex(seen_l[c1][i], seen_r[c2][j])

                    if index is not None:
                        if seen_l[c1][i][index] == seen_r[c2][j][index]:
                            if event != seen_l[c1][i][index]:
                                print"Seen: %s %s in event %s" % (event, [c1, b1, c2, b2], seen_l[c1][i][index])
                                event = seen_l[c1][i][index]

        seen_l[c1].setdefault(b1, []).append(event)
        seen_r[c2].setdefault(b2, []).append(event)

        vars.setdefault(event, []).append([c1, b1, c2, b2])

        left_end[c1].setdefault(b1, []).append([event, c1, b1, c2, b2])
        right_end[c2].setdefault(b2, []).append([event, c1, b1, c2, b2])

        # left_end[c1][b1] = [event, c1, b1, c2, b2]
        # right_end[c2][b2] = [event, c1, b1, c2, b2]


    return(vars, right_end, left_end)


def sameIndex(l1, l2):
    for i, v in enumerate(l1):
        if l1[i] == l2[i]: return i



def stitch(options):
    vars, right_end, left_end = parse_vars(options)

    # print(json.dumps(left_end, indent=4))

    complex_events = {}

    for c in sorted(right_end.keys()):
        for b1 in sorted(right_end[c].keys()):
            for i in range(b1-250, b1+250):
                if i in left_end[c]:
                    e1 = right_end[c][b1][0][0]
                    e2 = left_end[c][i][0][0]
                    if e1 == e2:
                        continue
                    if e1 not in complex_events:
                        complex_events.setdefault(e1, []).append([e1, e2])
                        print("Overlap between right event %s [%s] and left event %s [%s]") % (e1, b1, e2, i)



    print(json.dumps(complex_events, indent=4, sort_keys=True))
    return right_end, left_end



def rec_dd():
    return defaultdict(rec_dd)


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
