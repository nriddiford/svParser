#!/usr/bin/env python

import sys

genome_bps = 'bps_accross_genome.txt'
all_bps = 'all_bps.txt'

false_positves_file = 'all_samples_false_calls.txt'

false_calls = set()
with open(false_positves_file, 'U') as fp_file:
    for l in fp_file:
        parts = l.rstrip().split('_')
        key1 =  "%s_%s_%s" % (parts[0], parts[1], parts[2])
        key2 =  "%s_%s_%s" % (parts[0], parts[3], parts[4])
        false_calls.add(key1)
        false_calls.add(key2)

# Filter if bp coordinates found in 'all_samples_false_calls.txt'
def filter(coords, false_calls, outfile, line):
    if not coords in false_calls:
        outfile.write(line)

def remove_false_positives_from_bps(genome_bps_file, all_bps_file, all_bps_file_clean, genome_bps_file_clean):
    """Removes false positives from the files 'all_bps_new.txt' and 'bps_accross_genome_new.txt'"""

    with open(genome_bps,'U') as gen_in, open(all_bps,'U') as bps_in, open('all_bps_new.txt', 'w') as all_bps_out, open('bps_accross_genome_new.txt', 'w') as genome_bps_out:
        for l in bps_in:
            parts = l.rstrip().split('\t')
            # print(parts)

            foo = [parts[2], parts[3], parts[4] ]
            coords = "_".join(foo)

            filter(coords, false_calls, all_bps_out, l)

        for l in gen_in:
            parts = l.rstrip().split('\t')
            coords = "%s_%s_%s" % (parts[0], parts[1], parts[2])

            filter(coords, false_calls, genome_bps_out, l)

remove_false_positives_from_bps(genome_bps_file=genome_bps, all_bps_file=all_bps, all_bps_file_clean='all_bps_new.txt', genome_bps_file_clean='bps_accross_genome_new.txt')
