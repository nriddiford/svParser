#!/bin/sh

## usage
usage() {
    echo "
usage:   run_parser.sh [options]
options:
  -f    filter
  -m    merge
  -a    annotate
  -c    clean-up false positives anad reannotate
  -s    stats
  -ns   stats for notch-excluded hits
  -h    show this message
"
}

filter=0
merge=0
annotate=0
clean=0
stats=0
nstats=0

while getopts 'fmacsnh' flag; do
  case "${flag}" in
    f)  filter=1 ;;
    m)  merge=1 ;;
    a)  annotate=1 ;;
    c)  clean=1 ;;
    s)  stats=1 ;;
    n)  nstats=1 ;;
    h)  usage
        exit 0 ;;
  esac
done

#script_bin=/Users/Nick/iCloud/Desktop/script_test/SV_Parser/script # home
script_bin=/Users/Nick_curie/Desktop/script_test/SV_Parser/script # work

if [[ $filter -eq 1 ]]
then
  for lumpy_file in data/lumpy/*.vcf
  do
    echo "perl script/svParse.pl -v $lumpy_file -f a -t l -p"
    perl script/svParse.pl -v $lumpy_file -f a -t l -p
  done

  for delly_file in data/delly/*.vcf
  do
    echo "perl script/svParse.pl -v $delly_file -f a -t d -p"
    perl script/svParse.pl -v $delly_file -f a -t d -p
  done

  for novo_file in data/novobreak/*.vcf
  do
    echo "perl script/svParse.pl -v $novo_file -f a -t n -p"
    perl script/svParse.pl -v $novo_file -f a -t n -p
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
    echo "perl ../../script/svMerger.pl -f ${samples[i]}.*.txt"
    perl ../../script/svMerger.pl -f ${samples[i]}.*.txt
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

#features=/Users/Nick/Documents/Curie/Data/Genomes/Dmel_v6.12/Features/dmel-all-r6.12.gtf # home
features=/Users/Nick_curie/Documents/Curie/Data/Genomes/Dmel_v6.12/Features/dmel-all-r6.12.gtf # work

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
    echo "perl $script_bin/script/sv2gene.pl -f $features -i $clustered_file"
    perl $script_bin/script/sv2gene.pl -f $features -i $clustered_file
    rm $clustered_file
  done

fi

if [[ $clean -eq 1 ]]
then
  echo "Removing calls marked as flase positives in 'all_samples_false_calls.txt'"
  for annofile in *_annotated_SVs.txt
  do
    python $script_bin/clean.py -f $annofile
  done

  for clean_file in *cleaned_SVs.txt
  do
    if [[ ! -s $clean_file ]]
    then
      rm $clean_file
    fi

    if [[ -f $clean_file ]]
    then
      perl $script_bin/sv2gene.pl -r -f $features -i $clean_file
    fi
  done

  echo "Writing bp info for cleaned, reannotated SV calls to 'all_bps_cleaned.txt')"

  if [[ -f 'all_bps_cleaned.txt' ]]
  then
    rm 'all_bps_cleaned.txt'
  fi

  for reanno_file in *reannotated_SVs.txt
  do
    python $script_bin/getbps.py -f $reanno_file
  done

  echo "Removing false positives from bp file 'all_bps_cleaned.txt', writing new bp file to 'all_bps_new.txt'"
  python $script_bin/update_bps.py

fi

if [[ $stats -eq 1 ]]
then

  if [ -z all_genes.txt ]
  then
    echo "'all_genes' not found! Exiting"
    exit 1
  fi

  echo "Calculating breakpoint stats..."
  perl $script_bin/bpstats.pl all_bps_new.txt

fi

if [[ $nstats -eq 1 ]]
then

  if [ -z all_genes.txt ]
  then
    echo "'all_genes' not found! Exiting"
    exit 1
  fi

  echo "Calculating breakpoint stats..."
  perl $script_bin/bpstats.pl -n all_bps_new.txt

fi

exit 0
