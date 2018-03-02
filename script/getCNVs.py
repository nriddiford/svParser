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

script_dir = os.path.dirname(__file__)
path = os.path.join(script_dir, '../data/cnv')
output = os.path.join(path, outfile)

def getCNVs(calls):
    cnvs = []
    with open(calls, 'U') as sv_calls_file:
        for l in sv_calls_file:
            parts = l.rstrip().split('\t')

            if "cnv" in parts[1].lower():
                try:
                    depth = parts[17]
                except IndexError:
                    depth = '-'

                line = [parts[3], parts[4], parts[6], parts[2], parts[11], depth]

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
