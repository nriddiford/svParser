#!/usr/bin/python

import sys
from optparse import OptionParser
import os

parser = OptionParser()

parser.add_option("-f",
                  "--annotated_files",
                  dest="SV_calls",
                  help="SV calls file"
                  )

options, args = parser.parse_args()

if not options.SV_calls:
    parser.error('No input file provided')

base_name = (os.path.splitext(options.SV_calls)[0])
sample = base_name.split('_')[0]
outfile = sample + "_" + "CNVs.txt"

dir = os.path.dirname(__file__)
# path = os.path.join(dir, '../data/cnv')
path = os.path.join(dir, '..')

output = os.path.join(path, outfile)

def getCNVs(SV_calls):
    cnvs = []
    with open(SV_calls, 'U') as sv_calls_file:
        for l in sv_calls_file:
            parts = l.rstrip().split('\t')

            if "cnv" in parts[1].lower():
                try:
                    depth = parts[16]
                except IndexError:
                    depth = '-'

                line = [parts[3], parts[4], parts[6], parts[2], parts[10], depth]

                line = '\t'.join(line)
                cnvs.append(line)
    return(cnvs)

def printCNVs(cnvs, output):
    with open(output, 'w') as cnv_out:
        for l in cnvs:
            cnv_out.write("%s\n" % l)

cnvs = getCNVs(options.SV_calls)
if cnvs:
    printCNVs(cnvs, output)
