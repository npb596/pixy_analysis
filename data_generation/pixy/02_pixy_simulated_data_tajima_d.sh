#!/usr/bin/bash
#SBATCH --mem=20GB


mkdir -p tmp
mkdir -p data

find ../01_simulating-test-data/data/simulated_var_only -type f -name ".vcf.gz" > tmp/vcf_var_only.txt
find ../01_simulating-test-data/data/simulated_invar -type f -name ".vcf.gz" > tmp/vcf_invar.txt
find ../01_simulating-test-data/data/simulated_missing_sites -type f -name ".vcf.gz" > tmp/vcf_missing_sites.txt
find ../01_simulating-test-data/data/simulated_missing_genos -type f -name ".vcf.gz" > tmp/vcf_missing_genos.txt
find ../01_simulating-test-data/data/accuracy_invar -type f -name ".vcf.gz" > tmp/vcf_accuracy.txt

##rm -r data/var_only
##mkdir -p data/var_only
#
#while read vcf
#do
#
#mkdir -p tmp/1
#rm -r tmp/1
#vcfslug=$(echo $vcf | sed 's/.*\///g')
#echo $vcf
#
#python pixy.py --stats tajima_d --vcf $vcf --zarr_path tmp/1 --window_size 10000 --populations populations_pi.txt --bypass_filtration yes --outfile_prefix data/var_only/$vcfslug
#
#rm -r tmp/1
#
#done < tmp/vcf_var_only.txt 
#
# Parallelize instead of loop as above, assumes files have already been bgzipped and indexed
#vcfs=(`grep -v "tbi" tmp/vcf_var_only.txt`)
# gnu-parallel must be available on machine to use this
#parallel '
#vcfslug=$(basename {})
# Check if output files already exist, simply comment out every line but pixy run if you want to overwrite files
#if [ ! -s data/var_only/${vcfslug}_tajima_d.txt ];
#then
#python pixy.py --stats tajima_d --vcf {} --window_size 10000 --populations populations_pi.txt --bypass_invariant_check yes --output_folder data/var_only/ --output_prefix ${vcfslug}
#fi' ::: ${vcfs[@]}
#
##rm -r data/invar
##mkdir -p data/invar
#
#while read vcf
#do
#
#mkdir -p tmp/1
#rm -r tmp/1
#vcfslug=$(echo $vcf | sed 's/.*\///g')
#echo $vcf
#
#python pixy.py --stats tajima_d --vcf $vcf --zarr_path tmp/1 --window_size 10000 --populations populations_pi.txt --bypass_filtration yes --outfile_prefix data/invar/$vcfslug
#
#rm -r tmp/1  
#
#done < tmp/vcf_invar.txt
#
# Parallelize instead of loop as above, assumes files have already been bgzipped and indexed
#vcfs=(`grep -v "tbi" tmp/vcf_invar.txt`)
# gnu-parallel must be available on machine to use this
#parallel '
#vcfslug=$(basename {})
# Check if output files already exist, simply comment out every line but pixy run if you want to overwrite files
#if [ ! -s data/invar/${vcfslug}_tajima_d.txt ];
#then
#python pixy.py --stats tajima_d --vcf {} --window_size 10000 --populations populations_pi.txt --bypass_invariant_check yes --output_folder data/invar/ --output_prefix ${vcfslug}
#fi' ::: ${vcfs[@]}
#
##rm -r data/missing_sites
##mkdir -p data/missing_sites
#
#while read vcf
#do
#
#mkdir -p tmp/1
#rm -r tmp/1
#vcfslug=$(echo $vcf | sed 's/.*\///g')
#echo $vcf
#
#python pixy.py --stats tajima_d --vcf $vcf --zarr_path tmp/1 --window_size 10000 --populations populations_pi.txt --bypass_filtration yes --outfile_prefix data/missing_sites/$vcfslug
#
#rm -r tmp/1  
#
#done < tmp/vcf_missing_sites.txt
#
# Parallelize instead of loop as above, assumes files have already been bgzipped and indexed
#vcfs=(`grep -v "tbi" tmp/vcf_missing_sites.txt`)
# gnu-parallel must be available on machine to use this
#parallel '
#vcfslug=$(basename {})
# Check if output files already exist, simply comment out every line but pixy run if you want to overwrite files
#if [ ! -s data/missing_sites/${vcfslug}_tajima_d.txt ];
#then
#python pixy.py --stats tajima_d --vcf {} --window_size 10000 --populations populations_pi.txt --bypass_invariant_check yes --output_folder data/missing_sites/ --output_prefix ${vcfslug}
#fi' ::: ${vcfs[@]}
#
##rm -r data/missing_genos
##mkdir -p data/missing_genos
#
#while read vcf
#do
#
#mkdir -p tmp/1
#rm -r tmp/1
#vcfslug=$(echo $vcf | sed 's/.*\///g')
#echo $vcf
#
#python pixy.py --stats tajima_d --vcf $vcf --zarr_path tmp/1 --window_size 10000 --populations populations_pi.txt --bypass_filtration yes --outfile_prefix data/missing_genos/$vcfslug
#
#rm -r tmp/1   
#
#done < tmp/vcf_missing_genos.txt
#
# Parallelize instead of loop as above, assumes files have already been bgzipped and indexed
#vcfs=(`grep -v "tbi" tmp/vcf_missing_genos.txt`)
# gnu-parallel must be available on machine to use this
#parallel '
#vcfslug=$(basename {})
# Check if output files already exist, simply comment out every line but pixy run if you want to overwrite files
#if [ ! -s data/missing_genos/${vcfslug}_tajima_d.txt ];
#then
#python pixy.py --stats tajima_d --vcf {} --window_size 10000 --populations populations_pi.txt --bypass_invariant_check yes --output_folder data/missing_genos/ --output_prefix ${vcfslug}
#fi' ::: ${vcfs[@]}
#


mkdir -p data/accuracy_invar
while read vcf
do

mkdir -p tmp/1
rm -r tmp/1
vcfslug=$(echo $vcf | sed 's/.*\///g')
echo $vcf

python pixy.py --stats tajima_d --vcf $vcf --zarr_path tmp/1 --window_size 10000 --populations populations_pi.txt --bypass_filtration yes --outfile_prefix data/accuracy_invar/$vcfslug

rm -r tmp/1   

done < tmp/vcf_accuracy.txt

