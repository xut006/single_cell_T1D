#!/bin/sh

rename s/\ /_/ *.csv
rename s/\ /_/ *.csv

for f in *.csv 
do
  sed '1,11d' $f > $f.trimmed
  echo "Trimmed ... $f.trimmed"
done
