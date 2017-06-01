#!/usr/bin/bash

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

cd filtered
pwd

echo "perl ../script/merge_vcf.pl"
perl ../script/merge_vcf.pl


cd summary
samples+=( $(ls -1 *.txt | cut -d '.' -f 1 | sort -u ) )

for ((i=0;i<${#samples[@]};++i))
do
  echo "perl ../../script/svMerger.pl -f ${samples[i]}*.txt"
  perl ../../script/svMerger.pl -f ${samples[i]}*.txt
done
