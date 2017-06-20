#!/bin/sh

## usage
usage() {
    echo "
usage:   run_parser.sh [options]
options:
  -f    filter
  -m    merge
  -a    annotate
  -c    clean-up false positives
  -s    stats
  -h    show this message
"
}

filter=0
merge=0
annotate=0
clean=0
stats=0

while getopts 'fmacsh' flag; do
  case "${flag}" in
    f)  filter=1 ;;
    m)  merge=1 ;;
    a)  annotate=1 ;;
    c)  clean=1 ;;
    s)  stats=1 ;;
    h)  usage
        exit 0 ;;
  esac
done


if [[ $filter -eq 1 ]]
then
  for delly_file in data/lumpy/*.vcf
  do
    echo "perl script/sv_parse.1.0.pl -v $delly_file -f a -t l -p"
    perl script/sv_parse.1.0.pl -v $delly_file -f a -t l -p
  done

  for lumpy_file in data/delly/*.vcf
  do
    echo "perl script/sv_parse.1.0.pl -v $lumpy_file -f a -t d -p"
    perl script/sv_parse.1.0.pl -v $lumpy_file -f a -t d -p
  done

  for novo_file in data/novobreak/*.vcf
  do
    echo "perl script/sv_parse.1.0.pl -v $novo_file -f a -t n -p"
    perl script/sv_parse.1.0.pl -v $novo_file -f a -t n -p
  done

fi

cd filtered

mergeVCF=`which mergevcf || true`

if [[ -z "$mergeVCF" ]]
then
  usage
  echo -e "Error: mergevcf was not found. Please set in path\n"
  exit 1
fi

if [[ $merge -eq 1 ]]
then
  echo "perl ../script/merge_vcf.pl"
  perl ../script/merge_vcf.pl
fi

cd summary

if [[ $merge -eq 1 ]]
then
  samples+=( $(ls -1 *.txt | cut -d '.' -f 1 | sort -u ) )

  for ((i=0;i<${#samples[@]};++i))
  do
    echo "perl ../../script/svMerger.pl -f ${samples[i]}*.txt"
    perl ../../script/svMerger.pl -f ${samples[i]}*.txt
  done

fi

cd merged

if [[ $merge -eq 1 ]]
then
  for f in *_merged_SVs.txt
  do
    echo "perl ../../../script/svClusters.pl $f"
    perl ../../../script/svClusters.pl $f
    rm $f
  done
fi

#features=/Users/Nick/Documents/Curie/Data/Genomes/Dmel_v6.12/Features/dmel-all-r6.12.gtf
features=/Users/Nick_curie/Documents/Curie/Data/Genomes/Dmel_v6.12/Features/dmel-all-r6.12.gtf

if [[ $annotate -eq 1 ]]
then

  if [ -f all_genes.txt ]
  then
    rm all_genes.txt
  fi

  if [ -f all_bps.txt ]
  then
    rm all_bps.txt
  fi

  for clustered_file in *clustered_SVs.txt
  do
    echo "Annotating $clustered_file"
    echo "perl ../../../script/sv2gene.pl $features $clustered_file"
    perl ../../../script/sv2gene.pl $features $clustered_file
    rm $clustered_file
  done

fi

if [[ $clean -eq 1 ]]
then
  for annofile in *_SVs.txt
  do
    echo "Removing calls marked as flase positives in 'all_samples_false_calls.txt'"
    echo "python ../../../script/clean.py -f $annofile"
    python ../../../script/clean.py -f $annofile
    # rm $annofile
  done
fi

if [[ $stats -eq 1 ]]
then

  if [ -z all_genes.txt ]
  then
    echo "'all_genes' not found! Exiting"
    exit 1
  fi

  echo "Calculating breakpoint stats..."
  echo "perl script/bpstats.pl all_bps.txt"
  perl ../../../script/bpstats.pl all_bps.txt

fi

exit 0
