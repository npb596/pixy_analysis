#!/bin/bash

vcffolder=$1

vcf_folder_slug=$(echo $vcffolder | sed 's/..\/01_simulating-test-data\/data\///g' | sed 's/\//_/g')

echo $vcf_folder_slug

find $vcffolder -name "*.vcf*" > tmp/${vcf_folder_slug}_vcflist.txt

cd tmp

split -l 100 ${vcf_folder_slug}_vcflist.txt vcf_files_

cd ..

ls tmp/vcf_files_* > tmp/vcf_lists.tmp

while read vcffolder
do

./scripts/scikit_folder_of_vcfs.sh $vcffolder

done < tmp/vcf_lists.tmp

