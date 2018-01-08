# svParser

Parse VCF output from structural variant callers LUMPY and DELLY (also supports input from novoBreak).
Run without options or with `--help` or `-h` to print usage statement

## Parsing structural variants from VCF files called by LUMPY DELLY and novoBreak

## Summarise variants

A good place to start is with a summary of variants called in VCF file:

#### Read vcf file from lumpy (or delly/novobreak) and see summary of variants called:
`perl script/svParse.pl -v data/NA12878.NA12891.NA12892.vcf`

#### It's reccomended to explicity tag files on the caller that has produced them using the `--type` flag:
`perl script/svParse.pl -v data/NA12878.NA12891.NA12892.vcf -t l`

## Filtering variants
The defualt behaviour of svParser is to consider all variants present in a VCF file, and classify these according to genotype. Genotype will be assigned as follows:
* germline_recurrent - variants with high quality read support in both a tumour and normal sample and at least one other sample in a Panel of Normals (PON)
* germline_private - variants with high quality read support in both a tumour and normal sample but NOT in a any samples in a PON
* somatic_normal - variants with high quality read support in normal sample, but not in tumour sample
* somatic_tumour - variants with high quality read support in tumour sample, but not in normal sample

Variants called by LUMPY must also be genotyped by [SVTyper](https://github.com/hall-lab/svtyper)

The real power of svParser comes from its ability to easily filter variant calls on anumber of different criteria, and quickly assess how this affects your callset.
It is highly recommended to play around with differnt combinations of filters that suit your needs. Filter flags can be used with any of the other options aid in fine tuning

#### Filters available:
* `su`     - number of tumour reads supporting var. Expects integer e.g. `-f su=5`
* `dp`     - minimum depth for both tumour normal at variant site. Expects integer e.g. `-f dp=10`
* `rdr`    - supporting reads/tumour depth - a value of 1 would mean all reads support variant. Expects integer/float e.g. `-f rdr=0.2`
* `sq`     - phred-scaled variant likelihood. Expects integer e.g. `-f sq=10`
* `chr`    - filter out chromosomes not in `chroms.txt` (if not provided, defaults to chromosomes 2L 2R 3L 3R 4 X Y). Expects binary e.g. `-f chr=1`
* `s=1`    - only keep somatic events. Expects binary e.g. `-f s=1`
* `a`      - apply default combination of filters. Equivalent to (e.g.):

 ```
 perl script/svParse.pl \
 -v data/Droso_R7.lumpy.vcf \
 -f su=4 \     # min 4 reads supporting event in tumour
 -f dp=10 \    # min read depth of 10 in both tumour/normal
 -f rdr=0.1 \  # min 10% of reads at breakpoint supporting variant
 -f sq=10 \    # min Log10 likelihood of 10
 -f chr=1      # filter out calls on chromosomes not in 'chroms.txt'
 ```

* In addition, users can provide a bed file containing regions to exclude by using `-e [path/to/exclude.bed]`

# Filter-check-filter strategy to hone filters
Filters can be used in combination with other features of svParser to experiment with filters that remove obvious false positives while retaining true positives. E.g:

#### See summary of variants that have >= 4 reads supporting call in tumour (`-f su=4`)
`perl script/svParse.pl -v data/NA12878.NA12891.NA12892.vcf -t l -f su=4`

#### See summary of variants that have >= 4 reads supporting call in tumour and have both breakpoints on chromosomes contained in `chroms.txt`
`perl script/svParse.pl -v data/NA12878.NA12891.NA12892.vcf -t l -f su=4 -f chr=1`


## Inspect specific variant call

#### Investigate a specific variant (by ID) using the `--id` flag:
`perl script/svParse.pl -v data/NA12878.NA12891.NA12892.vcf -t l -i 1`

## Browse variants (-d)

Browse variants using the `--dump` flag. This cycles through each line of the VCF file, printing an easy-to-read breakdown for each variant.
Press any key to move to the next variant, or `q` to quit

`perl script/svParse.pl -v data/Droso_R7.lumpy.vcf -d`

#### Browse all variants with su>=5 on X chromosome (`-c X -f su=5`):
`perl script/svParse.pl -v data/Droso_R7.lumpy.vcf -t l -f su=5 -d -c X`

#### Browse all variants within a specific window on X chromosome (`-c X:3000000-3500000`):
`perl script/svParse.pl -v data/Droso_R7.lumpy.vcf -t l -d -c X:3000000-3500000`

#### Browse all variants that passed read depth filter (`-f dp=20`) filter within a specific window on X chromosome:
`perl script/svParse.pl -v data/Droso_R7.lumpy.vcf -t l -f dp=20 -d -c X:3000000-3500000`


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
  perl script/svParse.pl -v $file -f a -t l -p
done
```

```
for file in data/delly/*.vcf
do
  perl script/svParse.pl -v $file -f a -t d -p
done
```

```
for file in data/novobreak/*.vcf
do
  perl script/svParse.pl -v $file -f a -t n -p
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
  perl ../../script/svClusters.pl $f
  rm $f
done
```

### Annotate calls from a .gtf file

#### This takes a user-provided .gtf file and annotates variants with what feature/gene they affect:

```
features=/path/to/annotations.gtf

for clustered_file in *clustered_SVs.txt
do
    perl perl ../../script/sv2gene.pl -f $features -i $clustered_file
  fi
  rm $clustered_file
done
```


### Cleaning results

The recommended next step from here is to inspect calls that remain. Any obvious false positives can be marked as false positives by entering `F` in the `T/F` field for each `_annotated.txt` file.

These can then be removed by using:

```
for annofile in *_annotated_SVs.txt
do
  python ../../script/clean.py -f $annofile
done
```

# Plotting results
See [svBreaks](https://github.com/nriddiford/svBreaks) for plotting functions. This takes as input the annotated SV calls per sample

# To do
- [x] Exclude file
- [x] User provided chromosomes on whcih to filter
- [ ] Germline tagging
