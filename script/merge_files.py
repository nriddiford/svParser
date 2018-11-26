import os
import fnmatch
import pandas as pd
import ntpath


pattern = '*_reannotated_SVs.txt'
path = '.'
files = fnmatch.filter(os.listdir('.'), '*_reannotated_SVs.txt')

header = list(pd.read_csv(files[0], delimiter="\t", nrows=1).columns)

header.insert(0, 'sample')

with open('all_samples.txt', 'w') as all_out:
    all_out.write('\t'.join(header) + '\n')
    for f in sorted(files):
        sample = ntpath.basename(f).split("_")[0]
        df = pd.read_csv(f, delimiter="\t", index_col=False, na_filter=False)
        for idx, row in df.iterrows():
            all_out.write(sample + '\t')
            all_out.write('\t'.join(map(str, row)) + '\n')
