#!/usr/bin/env python
from __future__ import print_function
import sys, os
import pandas as pd

from optparse import OptionParser
from collections import defaultdict
import json
import ntpath


def parse_vars(options):
    print("Stitching together variants for: %s" % options.inFile)

    sample = ntpath.basename(options.inFile).split("_")[0]
    temp = os.path.join(options.out_dir, sample + '_temp.txt')

    df = pd.read_csv(options.inFile, delimiter="\t")

    if not 'event' in df:
        print("No event column")
        df['event'] = df.index
        df['notes'] = ""

    df = df.sort_values(['event', 'chromosome1', 'bp1', 'chromosome2', 'bp2'])

    right_end = defaultdict(lambda: defaultdict(dict))
    left_end = defaultdict(lambda: defaultdict(dict))

    seen_l = defaultdict(lambda: defaultdict(dict))
    seen_r = defaultdict(lambda: defaultdict(dict))
    seen = []

    for idx, row in df.iterrows():
        event, c1, b1, c2, b2, notes, genotype = row[['event', 'chromosome1', 'bp1', 'chromosome2', 'bp2', 'notes', 'genotype']]
        if genotype != 'somatic_tumour': continue

        for i in range(b1 - 10, b1 + 10):
            for j in range(b2 - 10, b2 + 10):
                if i in seen_l[c1] and j in seen_r[c2]:
                    index = same_index(seen_l[c1][i], seen_r[c2][j])
                    if index is not None:
                        if seen_l[c1][i][index] == seen_r[c2][j][index]:
                            if event != seen_l[c1][i][index]:
                                if event not in seen: print("Seen: %s %s in event %s" % (event, [c1, b1, c2, b2], seen_l[c1][i][index]))
                                seen.append(event)
                                event = seen_l[c1][i][index]
                                df.loc[idx, 'event'] = event

        seen_l[c1].setdefault(b1, []).append(event)
        seen_r[c2].setdefault(b2, []).append(event)

        left_end[c1].setdefault(b1, []).append([event, c1, b1, c2, b2])
        right_end[c2].setdefault(b2, []).append([event, c1, b1, c2, b2])

    df = df.sort_values(['event', 'chromosome1', 'bp1', 'chromosome2', 'bp2'])
    df.to_csv(temp, sep="\t", index=False)

    return right_end, left_end, temp


def same_index(l1, l2):
    for i, v in enumerate(l1):
        try:
            if l1[i] == l2[i]:
                return i
        except IndexError:
            pass


def stitch(right_end, left_end, options):
    """Tie together events where right edge and left edge are within 250 bases"""
    window = options.window
    complex_events = {}
    seen = []

    complex_events = join_same_end(left_end, window, complex_events)
    complex_events = join_same_end(right_end, window, complex_events)

    for c in sorted(right_end.keys()):
        for b1 in sorted(right_end[c].keys()):
            for i in range(b1 - window, b1 + window):
                if i in left_end[c]:
                    e1 = right_end[c][b1][0][0]
                    e2 = left_end[c][i][0][0]
                    if e1 == e2:
                        continue
                    seen_key = '_'.join(map(str, [e1, e2]))
                    if seen_key not in seen: print("Overlap between right event %s [%s] and left event %s [%s]" % (e1, b1, e2, i))
                    seen.append(seen_key)
                    complex_events.setdefault(e1, []).extend([e1, e2])

    # print(json.dumps(complex_events, indent=4, sort_keys=True))
    return pack_complex_vars(complex_events)


def join_same_end(d, window, complex_events):
    seen = []
    for c in sorted(d.keys()):
        for b1 in sorted(d[c].keys()):
            seen.append(b1)
            first_event = d[c][b1][0][0]
            for i, var in enumerate(sorted(d[c][b1])):
                # if first_event == var[0]: continue
                # complex_events.setdefault(first_event, []).extend([first_event, var[0]])
                for j in range(b1 - window, b1 + window):
                    if j in d[c] and j not in seen:
                        complex_events.setdefault(first_event, []).extend([first_event, d[c][j][0][0]])

    return complex_events


def pack_complex_vars(complex_events):
    """Join complex events where the value is also a key
       Delete keys where values combined into prior event"""

    extended = join_events(complex_events)

    for old in sorted(complex_events.keys()):
        for new, joined in sorted(extended.iteritems()):  # And new keys/values
            if new != old:
                for j in joined:
                    if j in complex_events[old] and old in extended:
                        # print("Event %s and %s have value %s" % (old, new, j))
                        extended.setdefault(old, []).extend(complex_events[new])
                        extended.pop(new, None)

    for k, v in sorted(extended.iteritems()):
        extended[k] = list(set(v))
        extended[k].sort()

    return extended


def join_events(d):
    """If a value in d is also a key, add values from older event into
      later event"""
    extended = {}
    seen = []
    for k, l in sorted(d.iteritems()):
        for event in l:
            if event in d:  # Is this event also a key?
                if event not in seen:
                    extended.setdefault(k, []).extend(d[event])
                    if k != event:  # Delete key if we're adding its contents to another key
                        extended.pop(event, None)
                seen.append(event)

    return extended


def print_complex(complex_events, options, temp):
    sample = ntpath.basename(options.inFile).split("_")[0]
    out_file = os.path.join(options.out_dir, sample + '_stitched.txt')

    print("Writing stitched variants to: %s" % out_file)
    df = pd.read_csv(temp, delimiter="\t")
    df = df.sort_values(['event', 'chromosome1', 'bp1', 'chromosome2', 'bp2'])
    events = defaultdict(lambda: defaultdict(int))

    for index, row in df.iterrows():
        event, sv_type, notes = row[['event', 'type', 'notes']]
        for e, j in sorted(complex_events.iteritems()):
            # if event == e:
            #     linked_events = "_".join(map(str, complex_events[event]))
            #     configuration = ":".join(map(str, [e, linked_events]))
            #     df.loc[index, 'configuration'] = configuration
            #     df.loc[index, 'type'] = '_'.join(["COMPLEX", sv_type])
            for joined in j:
                if event == joined:
                    # print("Event %s is in %s - key: %s" % (event, joined, e))
                    linked_events = "_".join(map(str, complex_events[e]))
                    configuration = ":".join(map(str, [event, linked_events]))
                    if len(complex_events[e]) > 1:
                        df.loc[index, 'configuration'] = configuration
                    df.loc[index, 'type'] = '_'.join(["COMPLEX", sv_type])
                    df.loc[index, 'event'] = e

        events[df.loc[index, 'event']][df.loc[index, 'type']] += 1

    # Finally, substitute "COMPELX_" for "" in events that are the same event chained together (DEL/DUPs)
    for event in events:
        if len(events[event]) == 1:
            # print("%s Not complex" % (event))
            val = events[event].keys()
            # print("Old value = %s" % (val))
            new = val[0].replace('COMPLEX_', '')
            # print("New value = %s" % (new) )
            df.loc[df['event'] == event, ['type']] = new

    # print(json.dumps(events, indent=4, sort_keys=True))

    df = df.sort_values(['event', 'chromosome1', 'bp1', 'chromosome2', 'bp2'])
    os.remove(temp)
    df.to_csv(out_file, sep="\t", index=False)


def main():
    parser = OptionParser()
    parser.add_option("-i", "--inFile", dest="inFile", help="An annotated variants file produced by sv2gene "
                                                            "accepts both '_annotated_SVs.txt' and "
                                                            "'_reannotated_SVs.txt' files", metavar="FILE")
    parser.add_option("-w", "--window", dest="window", action="store", type=int, help="The distance to search for connected breakpoints")
    parser.add_option("-o", "--out_dir", dest="out_dir", action="store", help="Directory to write output to " + "[Default: '.']")
    parser.set_defaults(window=1000, out_dir=os.getcwd())
    options, args = parser.parse_args()

    print("Tying variants +/- %s" % options.window)

    if options.inFile is None:
        parser.print_help()
        print()
    else:
        try:
            right_end, left_end, temp_out = parse_vars(options)
            complex_events = stitch(right_end, left_end, options)
            print_complex(complex_events, options, temp_out)
        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n")
            return


if __name__ == "__main__":
    sys.exit(main())
