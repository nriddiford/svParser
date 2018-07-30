#!/bin/sh
set -euo pipefail
## usage
usage() {
    echo "
usage:   run_parser.sh [options]
options:
  -f    filter
  -m    merge
  -o    output directory
  -a    annotate
  -c    clean-up false positives anad reannotate
  -s    stats
  -n    stats for notch-excluded hits
  -h    show this message
"
}

filter=0
merge=0
annotate=0
replace=0
out_dir='filtered/'
data_dir='data/'
cnv_dir=
exclude_file='/Users/Nick_curie/Documents/Curie/Data/Genomes/Dmel_v6.12/Mappability/dmel6_unmappable_100.bed'

while getopts 'fmarho:d:c:e:' flag; do
  case "${flag}" in
    f)  filter=1 ;;
    m)  merge=1 ;;
    a)  annotate=1 ;;
    r)  replace=1 ;;
    d)  data_dir="$OPTARG" ;;
    o)  out_dir="$OPTARG" ;;
    c)  cnv_dir="$OPTARG" ;;
    e)  exclude_file="$OPTARG" ;;
    h)  usage
        exit 0 ;;
  esac
done

if [[ $# -eq 0 ]]
then
  usage
  exit 0
fi

# Change to the 'script' dir in svParser
dir=$(dirname "$0")
script_bin="$dir/script"
script_bin="$( cd "$script_bin" ; pwd -P )"

echo "Reading data from '$data_dir'"
echo "Writing data to '$out_dir'"
echo "Exclude file set to '$exclude_file'"

mkdir -p "$out_dir/summary"

if [[ $filter -eq 1 ]]
then

  for lumpy_file in "$data_dir/lumpy/*.vcf"
  do
    if [ ! -f $lumpy_file ]
    then
      echo "No files in $data_dir/lumpy/"
    else
      echo "perl $script_bin/svParse.pl -v $lumpy_file -m l -f chr=1 -f su=4 -f dp=10 -f rdr=0.1 -f sq=10 -e $exclude_file -o $out_dir -p"
      perl $script_bin/svParse.pl -v $lumpy_file -m l -f chr=1 -f su=4 -f dp=10 -f rdr=0.1 -f sq=10 -e $exclude_file -p -o $out_dir
    fi
  done

  for delly_file in "$data_dir/delly/*.vcf"
  do
    if [ ! -f $delly_file ]
    then
      echo "No files in $data_dir/delly/"
    else
      echo "perl $script_bin/svParse.pl -v $delly_file -m d -f chr=1 -f su=4 -f dp=10 -f rdr=0.1 -f sq=10 -e $exclude_file -o $out_dir -p"
      perl $script_bin/svParse.pl -v $delly_file -m d -f chr=1 -f su=4 -f dp=10 -f rdr=0.1 -f sq=10 -e $exclude_file -p -o $out_dir
    fi
  done

  for novo_file in "$data_dir/novobreak/*.vcf"
  do
    if [ ! -f $novo_file ]
    then
      echo "No files in $data_dir/novobreak/"
    else
      echo "perl $script_bin/svParse.pl -v $novo_file -m n -f chr=1 -f su=4 -f dp=10 -f rdr=0.1 -f sq=10 -e $exclude_file -o $out_dir -p"
      perl $script_bin/svParse.pl -v $novo_file -m n -f chr=1 -f su=4 -f dp=10 -f rdr=0.1 -f sq=10 -e $exclude_file -p -o $out_dir
    fi
  done

  for cnv_file in "$data_dir/cnv/*.txt"
  do
    if [ ! -f $cnv_file ]
    then
      echo "No files in $data_dir/cnv/"
    else
      echo "perl $script_bin/parseCNV.pl -c $cnv_file -o $out_dir/summary"
      perl $script_bin/parseCNV.pl -c $cnv_file -o $out_dir/summary
    fi
  done

fi

function getBase(){
  stem=$(basename "$1" )
  name=$(echo $stem | cut -d '.' -f 1)
  echo $name
}

if [ -n "$cnv_dir" ]
then
  cd $out_dir/summary
  samples+=( $(ls -1 *.filtered.summary.txt | cut -d '.' -f 1 | sort -u ) )

  for ((i=0;i<${#samples[@]};++i))
  do
    # echo $out_dir/summary/${samples[i]}*.txt
    # echo $cnv_dir/${samples[i]}.*.cnv
    perl $script_bin/findCNV.pl -c $cnv_dir/${samples[i]}.*.cnv -v $out_dir/summary/${samples[i]}*.txt
  done
fi

#   echo "Getting read depth information from $cnv_dir"
#   for cnv_file in $(ls -1 $cnv_dir/*.cnv)
#   do
#     cname=$(getBase $cnv_file)
#     for sum_file in $(ls -1 $out_dir/summary/*.txt)
#     do
#       sname=$(getBase $sum_file)
#       if [[ $cname == $sname ]]
#       then
#         echo" perl $script_bin/findCNV.pl -c $cnv_dir -v $sum_file"
#       fi
#     done
#   done
# fi


cd $out_dir


if [[ $merge -eq 1 ]]
then

  mergeVCF=`which mergevcf || true`

  if [[ -z "$mergeVCF" ]]
  then
    usage
    echo -e "Error: mergevcf was not found. Please set in path\n`pip install mergevcf`"
    exit 1
  fi

  echo "perl $script_bin/merge_vcf.pl"
  #perl $script_bin/merge_vcf.pl
fi

cd summary


if [[ $merge -eq 1 ]]
then
  mkdir -p "$out_dir/summary/merged/"

  samples+=( $(ls -1 *.txt | cut -d '.' -f 1 | sort -u ) )

  for ((i=0;i<${#samples[@]};++i))
  do
    echo "perl $script_bin/svMerger.pl -f ${samples[i]}.*.txt"
    perl $script_bin/svMerger.pl -f ${samples[i]}.*.txt
  done

fi

cd merged

if [[ $merge -eq 1 ]]
then
  for f in *_merged_SVs.txt
  do
    echo "perl $script_bin/svClusters.pl $f"
    perl $script_bin/svClusters.pl $f
    rm $f
  done
fi

#features=/Users/Nick/Documents/Curie/Data/Genomes/Dmel_v6.12/Features/dmel-all-r6.12.gtf # home
features=/Users/Nick_curie/Documents/Curie/Data/Genomes/Dmel_v6.12/Features/dmel-all-r6.12.gtf # work

if [[ $annotate -eq 1 ]]
then

  if [ -f "all_genes.txt" ] && [ -f "all_bps.txt" ]
  then
    rm "all_genes.txt"
    rm "all_bps.txt"
  fi

  if [ -f all_samples_false_calls.txt ]
  then
    # if [ -f *_annotated_SVs.txt ]
    # then
    # for annofile in *_annotated_SVs.txt
    # do
    #   if [ -e "$annofile" ]
    #   then
    #     echo "Updating 'all_samples_false_calls.txt' with false positive calls from $annofile"
    #     echo "Updating 'all_samples_whitelist.txt' with whitelisted calls from $annofile"
    #     python $script_bin/clean.py -f $annofile
    #   fi
    # done
      rm *cleaned_SVs.txt
    # fi
  fi

  for clustered_file in *clustered_SVs.txt
  do
    echo "Annotating $clustered_file"
    # Should check both files individually
    if [ -f all_samples_false_calls.txt ] && [ -f all_samples_whitelist.txt ]
    then
      echo "perl $script_bin/sv2gene.pl -f $features -i $clustered_file -b all_samples_false_calls.txt -w all_samples_whitelist.txt"
      perl $script_bin/sv2gene.pl -f $features -i $clustered_file -b all_samples_false_calls.txt -w all_samples_whitelist.txt
    elif [ -f all_samples_false_calls.txt ]
    then
      echo "perl $script_bin/sv2gene.pl -f $features -i $clustered_file -b all_samples_false_calls.txt"
      perl $script_bin/sv2gene.pl -f $features -i $clustered_file -b all_samples_false_calls.txt
    elif [ -f all_samples_whitelist.txt ]
    then
      echo "perl $script_bin/sv2gene.pl -f $features -i $clustered_file -w all_samples_whitelist.txt"
      perl $script_bin/sv2gene.pl -f $features -i $clustered_file -w all_samples_whitelist.txt
    else
      echo "perl $script_bin/sv2gene.pl -f $features -i $clustered_file"
      perl $script_bin/sv2gene.pl -f $features -i $clustered_file
    fi
    rm $clustered_file
  done

fi

if [[ $replace -eq 1 ]]
then
    echo "Adding any new CNV calls to data/cnv'"
    for annofile in *_annotated_SVs.txt
    do
      python $script_bin/getCNVs.py -f $annofile
    done

  echo "Removing calls marked as false positives in 'all_samples_false_calls.txt'"
  for annofile in *_annotated_SVs.txt
  do
    perl $script_bin/clean_files.pl $annofile
  done

  if [ -f "all_genes_filtered.txt" ] && [ -f "all_bps_filtered.txt" ]
  then
    rm "all_genes_filtered.txt"
    rm "all_bps_filtered.txt"
  fi

  for clean_file in *cleaned_SVs.txt
  do

    # Delete file if empty
    if [[ ! -s $clean_file ]]
    then
      rm $clean_file
    else
      # Annotate un-annotated (manually added) calls
      # Append any new hit genes to 'all_genes.txt'
      perl $script_bin/sv2gene.pl -r -f $features -i $clean_file
      rm $clean_file
    fi
  done

  echo "Writing bp info for cleaned, reannotated SV calls to 'all_bps_cleaned.txt')"

  # for reanno_file in *reannotated_SVs.txt
  # do
  #   # Grab some of the fields from these newly annotated files, and write them to 'all_bps_cleaned.txt'
  #   python $script_bin/getbps.py -f $reanno_file
  # done
  #
  # # This shouldn't be neccessary. All calls in this file are taken from 'reannotated' files, which should have FP removed already...
  # echo "Removing false positives from bp file 'all_bps_cleaned.txt', writing new bp file to 'all_bps_filtered.txt'"
  # python $script_bin/update_bps.py
  # rm 'all_bps_cleaned.txt'

  # Merge all samples
  echo "Merging all samples into single file..."
  perl $script_bin/merge_samples.pl *reannotated_SVs.txt

fi

exit 0
