#!/usr/bin/env python

import sys
import os

from optparse import OptionParser

def get_file_name(variants_file):
    base_name = (os.path.splitext(variants_file)[0])
    out_base = base_name.split('_')[0]
    outfile = out_base + "_" + "stitched_SVs.txt"
    return(outfile)

def parse_vars(variants_file):
    out_file = get_file_name(variants_file)

    print("Stitching together variants for: %s" % variants_file)
    print("Writing stitched variants to: %s" % out_file)

    with open(variants_file, 'U') as inFile:
        variant = []
        for l in inFile:
            parts = l.rstrip().strip('\"').split('\t')

            if parts[0] == 'event':
                header = l
                continue

            try:
                notes = parts[21]
            except IndexError:
                notes = 0

            if notes == 'F':
                continue

            if parts[9] != 'somatic_tumour':
                continue


        return(parts, header, out_file)

def stitch(variants_file):
    variants, head, outFile = parse_vars(variants_file)

    # with open(outFile, 'w') as stitched_SVs:
    #     stitched_SVs.write("%s\n" % head)

    




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
            stitch(options.inFile)

        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n");
            return

if __name__ == "__main__":
    sys.exit(main())
