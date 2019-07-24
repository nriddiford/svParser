[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/nriddiford/svParser/blob/master/LICENSE)

# svParser

Filter, genotype, annotate and combine VCF files from structural variant callers LUMPY, DELLY and novobreak.
Run without options or with `--help` or `-h` to print usage statement.

This tool is under constant development. Please feel free to [contact me](mailto:nick.riddiford@curie.fr), or [raise an issue](https://github.com/nriddiford/svParser/issues) if you encounter any problems.


# Table of Contents
* [Installation](#installation)
* [Summarise variants](#summarise-variants)
* [Genotyping variants](#genotyping-variants)
* [Filtering variants](#filtering-variants)
* [Test filters against true positives](#test-filters-on-true-positive-variant-calls)
* [Inspecting variants](#inspect-specific-variant-call)
* [Browse variants](#browse-variants)
* [Print variants](#print-variants)
* [Combine variants from multiple sources](#combining-calls-from-multiple-sources)

## Installation
```
git clone https://github.com/nriddiford/svParser.git
```

Install cpanm for easy installation of requirements:

```
wget -O- http://cpanmin.us | perl - -l ~/perl5 App::cpanminus local::lib
```

To avoid issues with root access, I recommended using a local Perl library, which can be setup as follows (following [this](https://stackoverflow.com/questions/2980297/how-can-i-use-cpan-as-a-non-root-user)):

```
eval `perl -I ~/perl5/lib/perl5 -Mlocal::lib`
```
Add the following lines to your `bash_profile` (or`.profile/.bashrc/`):

```
eval `perl -I ~/perl5/lib/perl5 -Mlocal::lib`
export MANPATH=$HOME/perl5/man:$MANPATH
```

Once you have cpanm installed, then installing the dependencies should be as simple as:

```
cd svParser
bash install_deps.sh
```


## Summarise variants
A good place to start is with a summary of variants called in VCF file:


#### Read vcf file from lumpy (or delly/novobreak) and see summary of variants called:
`perl script/svParse.pl -v data/Droso_R7.lumpy.vcf`
It's reccomended to explicity tag files on the caller that has produced them using the `--method` flag:   
`perl script/svParse.pl -v data/Droso_R7.lumpy.vcf -m l`


## Genotyping variants
The defualt behaviour of svParser is to consider all variants present in a VCF file, and classify these according to genotype. Genotype will be assigned as follows:
* germline_recurrent - variants with high quality read support in both a tumour and normal sample and at least one other sample in a Panel of Normals (PON)
* germline_private - variants with high quality read support in both a tumour and normal sample but NOT in a any samples in a PON
* somatic_normal - variants with high quality read support in normal sample, but not in tumour sample
* somatic_tumour - variants with high quality read support in tumour sample, but not in normal sample

Genotypes can also be filtered in using the genotype filters outlined in the [filter section](#filters-available)

**IMPORTANT:**
*  Variants called by LUMPY must also be genotyped by [SVTyper](https://github.com/hall-lab/svtyper)
*  Genotyping requires that the samples in the VCF are ordered as follows: `TUMOUR NORMAL PON.1 ... PON.N`


## Filtering variants
The real power of svParser comes from its ability to easily filter variant calls on anumber of different criteria, and quickly assess how this affects your callset.
It is highly recommended to play around with differnt combinations of filters that suit your needs. Filter flags can be used with any of the other options aid in fine tuning


#### Filters available:
```
su   -   number of tumour reads supporting variant. Expects integer e.g. `-f su=5`
dp   -   minimum depth for both tumour normal at variant site. Expects integer e.g. `-f dp=10`
rdr  -   supporting reads/tumour depth - a value of 1 would mean all reads support variant. Expects integer/float e.g. `-f rdr=0.2`
sq   -   phred-scaled variant likelihood. Expects integer e.g. `-f sq=10`
chr  -   filter out chromosomes not in `chroms.txt` (if not provided, defaults to chromosomes 2L 2R 3L 3R 4 X Y). Expects binary e.g. `-f chr=1`
st   -   only keep somatic TUMOUR variants. Expects binary e.g. `-f st=1`
sn   -   only keep somatic NORMAL variants. Expects binary e.g. `-f sn=1`
gp   -   only keep germline PRIVATE variants. Expects binary e.g. `-f gp=1`
gr   -   only keep germline RECURRENT variants. Expects binary e.g. `-f gr=1`
a    -   apply default combination of filters. Equivalent to:

perl script/svParse.pl \
-v data/Droso_R7.lumpy.vcf \
-f su=4 \     # min 4 reads supporting event in tumour
-f dp=10 \    # min read depth of 10 in both tumour/normal
-f rdr=0.1 \  # min 10% of reads at breakpoint supporting variant
-f sq=10 \    # min Log10 likelihood of 10
-f chr=1      # filter out calls on chromosomes not in 'chroms.txt'
```

* In addition, users can provide a bed file containing regions to exclude by using `-e [path/to/exclude.bed]`

### Filter-check-filter strategy to hone filters
Filters can be used in combination with other features of svParser to experiment with filters that remove obvious false positives while retaining true positives. E.g:


#### See summary of variants that have >= 4 reads supporting call in tumour (`-f su=4`)
`perl script/svParse.pl -v data/Droso_R7.lumpy.vcf -m l -f su=4`


#### See summary of variants that have >= 4 reads supporting call in tumour and have both breakpoints on chromosomes contained in `chroms.txt`
`perl script/svParse.pl -v data/Droso_R7.lumpy.vcf -m l -f su=4 -f chr=1`


## Test filters on true positive variant calls
If you have a set of true positive calls that you want to protect from filtering, you can provide a set of sample-specfic true positives:

```
sample  caller   type   chrom1   bp1      chrom2    bp2
R7      lumpy    DEL    X        456393   X         4588700
```
`perl script/svParse.pl -v data/Droso_R7.lumpy.vcf -m l -f su=4 -f chr=1 -t tests/`

This allows you to see how adjusting filters affects your true positives.


## Inspect specific variant call
Sometimes it's useful to see the information behind a specific variant call. svParser makes it easy to dig out a variant using it's `id` (taken from the VCF column `ID`)
If the variant is affected by supplied filters this will be reported in the summary


#### Investigate a specific variant (by ID) using the `--id` flag:
`perl script/svParse.pl -v data/Droso_R7.lumpy.vcf -m l -i 1`


## Browse variants
Browse variants using the `--dump` flag. This cycles through each line of the VCF file, printing an easy-to-read breakdown for each variant.
Press any key to move to the next variant, or `q` to quit

`perl script/svParse.pl -v data/Droso_R7.lumpy.vcf -d`


#### Browse all variants with su>=5 on X chromosome (`-c X -f su=5`):
`perl script/svParse.pl -v data/Droso_R7.lumpy.vcf -m l -f su=5 -d -c X`


#### Browse all variants within a specific window on X chromosome (`-c X:3000000-3500000`):
`perl script/svParse.pl -v data/Droso_R7.lumpy.vcf -m l -d -c X:3000000-3500000`


#### Browse all variants that passed read depth filter (`-f dp=20`) filter within a specific window on X chromosome:
`perl script/svParse.pl -v data/Droso_R7.lumpy.vcf -m l -f dp=20 -d -c X:3000000-3500000`


## Print variants
Write a new VCF file (`filtered/input.filtered.vcf`) and a summary txt file containing useful info (`filtered/summary/input.filtered.summary.txt`).


#### Print (`-p`) all variants that passed all defualt filters (`-f a`):
`perl script/svParse.pl -v data/Droso_R7.lumpy.vcf -f a -p`


# Combining calls from multiple sources
Another useful feature of svParser is its ability to combine calls made by multiple callers into a unified call set per sample


### Run on all samples
#### Generate filtered vcf and summary file for calls made by delly, lumpy and novoBreak:

```
for file in data/lumpy/*.vcf
do
  perl script/svParse.pl -v $file -f a -m l -p
done
```

```
for file in data/delly/*.vcf
do
  perl script/svParse.pl -v $file -f a -m d -p
done
```

```
for file in data/novobreak/*.vcf
do
  perl script/svParse.pl -v $file -f a -m n -p
done
```

### Merge calls together in VCF using merge_vcf
#### Run [merge_vcf](https://github.com/ljdursi/mergevcf/tree/master/mergevcf) on all filtered files. (Run from `filtered` directory):

```
cd filtered
perl ../script/merge_vcf.pl
```

### Merge summary output foreach sample
#### Run svMerger.pl for all files in sample group. E.g. `input1.delly.filtered.summary.txt` `input1.lumpy.filtered.summary.txt` `input1.novobreak.filtered.summary.txt`

```
cd filtered/summary
perl ../../script/svMerger.pl -f input1.*
```

### Cluster sv calls that are close together (+/- 50 bp):

```
for f in *_merged_SVs.txt
do
  perl ../../script/svClusters.pl $f -d 50
  rm $f
done
```

### Annotate calls from a .gtf file
#### This takes a user-provided .gtf file and annotates variants with what feature/gene they affect:

```
features=/path/to/annotations.gtf

for clustered_file in *clustered_SVs.txt
do
  perl ../../script/sv2gene.pl -f $features -i $clustered_file
  rm $clustered_file
done
```


### Cleaning results
The recommended next step from here is to inspect calls that remain. Any obvious false positives can be marked as false positives by entering `F` in the `T/F` column for each `_annotated.txt` file.

These can then be removed by using:

```
for annofile in *_annotated_SVs.txt
do
  python ../../script/clean.py -f $annofile
done
```

# Plotting results
See [svBreaks](https://github.com/nriddiford/svBreaks) for plotting functions. This takes as input the annotated SV calls per sample
