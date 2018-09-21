#!/usr/bin/env python
from __future__ import print_function
import sys, os
import pandas as pd

from optparse import OptionParser
from collections import defaultdict
import json


def get_file_name(variants_file):
    base_name = (os.path.splitext(variants_file)[0])
    out_base = base_name.split('_')[0]
    outfile = out_base + "_" + "stitched_SVs.txt"
    return outfile


def parse_vars(options):
    out_file = get_file_name(options.inFile)

    print("Stitching together variants for: %s" % options.inFile)
    print("Writing stitched variants to: %s" % out_file)

    df = pd.read_csv(options.inFile, delimiter="\t")
    df = df[['event', 'chromosome1', 'bp1', 'chromosome2', 'bp2']]
    df = df.sort_values(['event', 'chromosome1', 'bp1', 'chromosome2', 'bp2'])

    right_end = defaultdict(lambda: defaultdict(dict))
    left_end = defaultdict(lambda: defaultdict(dict))

    seen_l = defaultdict(lambda: defaultdict(dict))
    seen_r = defaultdict(lambda: defaultdict(dict))

    variants = {}

    for row in df.itertuples():
        idx, event, c1, b1, c2, b2 = row
        for i in range(b1-10, b1+10):
            for j in range(b2-10, b2+10):
                if i in seen_l[c1] and j in seen_r[c2]:
                    index = same_index(seen_l[c1][i], seen_r[c2][j])
                    if index is not None:
                        if seen_l[c1][i][index] == seen_r[c2][j][index]:
                            if event != seen_l[c1][i][index]:
                                print("Seen: %s %s in event %s" % (event, [c1, b1, c2, b2], seen_l[c1][i][index]))
                                event = seen_l[c1][i][index]

        seen_l[c1].setdefault(b1, []).append(event)
        seen_r[c2].setdefault(b2, []).append(event)

        variants.setdefault(event, []).extend([c1, b1, c2, b2])

        left_end[c1].setdefault(b1, []).extend([event, c1, b1, c2, b2])
        right_end[c2].setdefault(b2, []).extend([event, c1, b1, c2, b2])

    return variants, right_end, left_end


def same_index(l1, l2):
    for i, v in enumerate(l1):
        if l1[i] == l2[i]:
            return i


def stitch(options):
    """Tie together events where right edge and left edge are within 250 bases"""
    variants, right_end, left_end = parse_vars(options)
    complex_events = {}

    for c in sorted(right_end.keys()):
        for b1 in sorted(right_end[c].keys()):
            for i in range(b1-250, b1+250):
                if i in left_end[c]:
                    e1 = right_end[c][b1][0]
                    e2 = left_end[c][i][0]
                    if e1 == e2:
                        continue
                    if e1 not in complex_events:
                        complex_events.setdefault(e1, []).extend([e1, e2])
                        print("Overlap between right event %s [%s] and left event %s [%s]" % (e1, b1, e2, i))


    pack_complex_vars(complex_events)

    return complex_events


def rec_dd():
    return defaultdict(rec_dd)


def pack_complex_vars(complex):
    """Join complex events where the value is also a key
       Delete keys where values combined into prior event"""
    # print(json.dumps(complex, indent=4, sort_keys=True))
    extended = join_events(complex)

    for old in sorted(complex.keys()): # For old keys
        for new, joined in sorted(extended.iteritems()): # And new keys/values
            if new != old:
                for j in joined:
                    if j in complex[old] and old in extended:
                        print("Event %s and %s have value %s" % (old, new, j))
                        extended.setdefault(old, []).extend(complex[new])

    # extended = join_events(complex)


        # for new in extended[k]:
        #     if new in extended:
        #         extended.setdefault(k, []).extend(extended[new])
        #         # del extended[new]
        #         break
        # extended[k] = list(set(extended[k]))

    # extended = flatten(extended)

    print(json.dumps(extended, indent=4, sort_keys=True))


def join_events(d):
    extended = {}
    seen = []
    for k, l in sorted(d.iteritems()):
        for event in l:
            if event in d: # Is this event also a key?
                if event not in seen:
                    extended.setdefault(k, []).extend(d[event])
                    if k != event: # Delete key if we're adding its contents to another key
                        extended.pop(event, None)
                    else: # else the event and key are the same, so we can remove it from list
                        extended[k].remove(event)
                seen.append(event)

    return extended


def main():
    parser = OptionParser()

    parser.add_option("-i",
        "--inFile", dest="inFile",
        help="An annotated variants file produced by sv2gene "
             "accepts both '_annotated_SVs.txt' and "
             "'_reannotated_SVs.txt' files",
        metavar="FILE")

    options, args = parser.parse_args()

    if options.inFile is None:
        parser.print_help()
        print()
    else:
        try:
            stitch(options)

        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n")
            return


if __name__ == "__main__":
    sys.exit(main())