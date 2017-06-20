#!/usr/bin/env python

import sys
import os

from optparse import OptionParser

def index_exists(ls, i):
    return (0 <= i < len(ls)) or (-len(ls) <= i < 0)

parser = OptionParser()

parser.add_option("-f",
    "--in_file",
    dest="in_file",
    help="Input")

options, args = parser.parse_args()

if not options.in_file:
    parser.error('No input file provided')

base_name = (os.path.splitext(options.in_file)[0])
out_base = base_name.split('.')[0]
outfile = out_base + "." + "cleaned_SVs.txt"
print("Writing true positives to '%s'" % outfile)

false_calls_file = 'all_samples_false_calls.txt'

def remove_false_positives(false_calls_file, input_file, clean_output):
    """Remove entries in input_file if the 20th column equals F."""
    i = 1
    filtered_calls = 0
    with open(false_calls_file, 'a+') as false_calls, open(input_file,'U') as infile, open(clean_output, 'w') as clean_files:
        false_calls.seek(0)
        seen_lines = {line.rstrip() for line in false_calls}
        for l in infile:
            parts = l.rstrip().split('\t')

            if i == 1 and parts[0] == 'event':
                clean_files.write(l)
                continue

            match = "%s_" % out_base + "_".join(parts[3:7])
            notes = '-'

            if index_exists(parts, 20):
                notes = parts[20]

            if notes == 'F':
                filtered_calls += 1
                if not match in seen_lines:
                    false_calls.write("%s\n" % (match))

            elif match in seen_lines:
                filtered_calls += 1
            else:
                clean_files.write(l)
            i += 1
    print "Removed %s false positives from %s" % (filtered_calls, input_file)

remove_false_positives(false_calls_file=false_calls_file, input_file=options.in_file, clean_output=outfile)
