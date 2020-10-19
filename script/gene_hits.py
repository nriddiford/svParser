import fnmatch
import ntpath
import numpy as np
import os
import pandas as pd
import sys
from operator import itemgetter
from optparse import OptionParser


def read_bed(bed_dir):
    files = fnmatch.filter(os.listdir(bed_dir), '*.bed')
    all_bed = {}
    for f in sorted(files):
        name = ntpath.basename(f).split(".")[0]
        file = os.path.join(bed_dir, f)
        df = pd.read_csv(file, delimiter="\t", index_col=False, na_filter=False, header=None)
        df.columns = ['chrom', 'start', 'stop', 'gene', 'id']
        all_bed[name] = list(df['gene'])

    return all_bed

def read_file(bed_file):
    df = pd.read_csv(bed_file, delimiter="\t", index_col=False, na_filter=False, header=None)
    df.columns = ['chrom', 'start', 'stop', 'gene', 'id', 'groups']

    genes = list(df['gene'])
    all_bed = dict(zip(list(df['gene']), list(df['groups'])))

    return all_bed, genes


def get_vars(options, all_files, genes):
    all_genes = pd.read_csv(options.variants, delimiter="\t", index_col=False, na_filter=False)
    all_genes.columns = ['event', 'sample', 'genotype', 'type', 'allele_frequency', 'length', 'cnv', 'chromosome', 'gene']

    all_genes['cell_fraction'] = np.where(all_genes['chromosome'] == "X", all_genes['allele_frequency'], all_genes['allele_frequency']*2)

    sam = "all samples"
    if options.sample:
        sam = options.sample
        all_genes = all_genes[(all_genes["sample"] == options.sample)]

    print("Looking for variants with cell fraction >= %s <= %s in %s") % (options.min_cf, options.max_cf, sam)

    hit_genes = all_genes[(all_genes['cell_fraction'] >= options.min_cf) &
                          (all_genes['cell_fraction'] <= options.max_cf)]

    # hit_genes = hit_genes[(hit_genes["type"] == "DEL")]

    # all_files = read_bed(options.bed_dir)
    ranked_hits = []

    for idx, row in hit_genes.iterrows():
        hit_gene = row['gene'].lower()
        if hit_gene in all_files:
            ranked_hits.append([row['sample'], row['event'], row['gene'], row['type'], row['cell_fraction'], all_files[hit_gene]])

    for h in sorted(ranked_hits, key=itemgetter(4), reverse=True):
        print(h)
    return ranked_hits


def get_muts(options, all_files, genes):
    muts = pd.read_csv(options.muts, index_col=False, na_filter=False)

    muts['event'] = 'NA'
    sam = "all samples"

    if options.sample:
        sam = options.sample
        muts = muts[(muts["sample_old"] == options.sample)]

    print("Looking for protein coding mutations with cell fraction >= %s <= %s in %s") % (options.min_cf, options.max_cf, sam)

    hit_genes = muts[(muts['cell_fraction'] >= options.min_cf) &
                          (muts['cell_fraction'] <= options.max_cf)]


    ranked_hits = []
    for idx, row in hit_genes.iterrows():
        hit_gene = row['gene'].lower()
        if hit_gene in all_files:
            ranked_hits.append([row['sample_old'], ':'.join([row['chr'], str(row['pos'])]), row['gene'], row['impact'], row['cell_fraction'], all_files[hit_gene]])

    for h in sorted(ranked_hits, key=itemgetter(4), reverse=True):
        print(h)

    return ranked_hits


def get_snvs(options, all_files, genes):
    all_genes = pd.read_csv(options.snvs, delimiter="\t", index_col=False, na_filter=False)
    all_genes['cell_fraction'] = np.where((all_genes['chromosome'] == "X") | (all_genes['chromosome'] == "Y"), all_genes['allele_frequency'], all_genes['allele_frequency']*2)
    all_genes['event'] = 'NA'
    sam = "all samples"
    if options.sample:
        sam = options.sample
        all_genes = all_genes[(all_genes["sample"] == options.sample)]

    print("Looking for SNVs with cell fraction >= %s <= %s in %s") % (options.min_cf, options.max_cf, sam)

    hit_genes = all_genes[(all_genes['cell_fraction'] >= options.min_cf) &
                          (all_genes['cell_fraction'] <= options.max_cf)]
    if options.status:
        hit_genes = all_genes[all_genes['status'] == options.status]

    # all_files = read_bed(options.bed_dir)
    ranked_hits = []
    # for b in all_files.iterkeys():
    #     genes = all_files[b]
    for idx, row in hit_genes.iterrows():
        hit_gene = row['gene'].lower()
        snpEff_gene = row['snpEff_anno'].lower()
        if hit_gene in all_files:
            ranked_hits.append([row['sample'], ':'.join([row['chromosome'], str(row['pos'])]), row['gene'], row['variant_type'], row['status'], row['cell_fraction'], all_files[hit_gene]])
        # elif row['snpEff_anno'].lower() in all_files:
        #     print("found in snpEff annotated gene")
        #     ranked_hits.append([row['sample'], ':'.join([row['chromosome'], str(row['pos'])]), row['snpEff_anno'], row['variant_type'], row['status'], row['cell_fraction'], all_files[hit_gene]])

    for h in sorted(ranked_hits, key=itemgetter(5), reverse=True):
        print(h)

    return ranked_hits

def print_hits(options, hits):

    df = pd.DataFrame(hits)
    df.to_csv(options.out, sep="\t", index=False)

    # for h in sorted(hits, key=itemgetter(5), reverse=True):
    #     options.out.write("\t".join(map(str, [h])))
    # hits = hits.sort_values(['cell_fraction', 'sample'])
    # hits.to_csv(options.out, sep="\t", index=False)

def main():
    parser = OptionParser()
    parser.add_option("-v", "--variants", dest="variants", help="An annotated variants file produced by sv2gene ", metavar="FILE")
    parser.add_option("-s", "--snvs", dest="snvs", help="An annotated snvs file produced by snv2gene ", metavar="FILE")
    parser.add_option("-m", "--muts", dest="muts", help="Protein coding mutations, annotated by dnds", metavar="FILE")
    parser.add_option("--sample", dest="sample", action="store", help="The sample to inspect")
    parser.add_option("-b", "--bed_dir", dest="bed_dir", help="The directory containing bed files")
    parser.add_option("-f", "--bed_file", dest="bed_file", help="A bed file to scan")
    parser.add_option("--min", dest="min_cf", action="store", type=float, help="Min cell fraction")
    parser.add_option("--max", dest="max_cf", action="store", type=float, help="Max cell fraction")
    parser.add_option("--status", dest="status", action="store", type=str, help="Status for SNVS/indels")
    parser.add_option("-o", "--out_file", dest="out", action="store", type=str, help="File to write hits to")

    parser.set_defaults(min_cf=0.1, max_cf=1)
    options, args = parser.parse_args()

    if (options.variants is None and options.snvs is None and options.muts is None) or (options.bed_dir is None and options.bed_file is None):
        parser.print_help()
        print
    else:
        try:
            if options.bed_dir:
                genes = read_bed(options.bed_dir)
            elif options.bed_file:
                df, genes = read_file(options.bed_file)

            if options.variants:
                hits = get_vars(options, df, genes)
            elif options.snvs:
                hits = get_snvs(options, df, genes)
            else:
                hits = get_muts(options, df, genes)

            print_hits(options, hits)

        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n")
            return


if __name__ == "__main__":
    sys.exit(main())
