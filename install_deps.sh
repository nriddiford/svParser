#/usr/bin/bash

while read lib; do
  echo "cpanm $lib"
  cpanm $lib
done <requirements.txt
