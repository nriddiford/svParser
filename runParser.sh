#!/bin/sh

## usage
usage() {
    echo "
usage:   run_parser.sh [options]
options:
  -f    filter
  -m    merge
  -a    annotate
  -h    show this message
"
}

filter=0
merge=0
annotate=0

while getopts 'fmah' flag; do
  case "${flag}" in
    f)  filter=1 ;;
    m)  merge=1 ;;
    a)  annotate=1 ;;
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

features=/Users/Nick/Documents/Curie/Data/Genomes/Dmel_v6.12/Features/dmel-all-r6.12.gtf

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

  for merged_file in *merged_SVs.txt
  do
    echo "Annotating $merged_file"
    echo "perl ../../../script/sv2gene.v0.1.pl $features $merged_file"
    perl ../../../script/sv2gene.v0.1.pl $features $merged_file
  done

fi

exit 0
