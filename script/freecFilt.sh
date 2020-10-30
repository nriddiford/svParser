#!/usr/bin/env bash

## usage
usage() {
    echo -e "
options:
  -e    exclude bed file
  -d    directory containing cnv calls
  -h    show this message
"
}

exclude_file='/Users/Nick_curie/Documents/Curie/Data/Genomes/Dmel_v6.12/Mappability/dmel6_unmappable_100.bed'
dir=$(pwd)

while getopts 'hd:e:' flag; do
  case "${flag}" in
    d)  dir="$OPTARG" ;;
    e)  exclude_file="$OPTARG" ;;
    h)  usage
        exit 0 ;;
  esac
done


function getBase(){
  stem=$(basename "$1" )
  name=$(echo $stem | cut -d '_' -f 1)
  echo $name
}

function argv {
    for f in $dir/*_sig_cnvs.txt
    do
      s=$(getBase $f)
      echo "$s"
      bedtools subtract -a $f -b $exclude_file -f 0.25 -N > ${s}_filt_cnvs.txt
    done
    echo
}

argv
