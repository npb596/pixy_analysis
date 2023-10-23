#!/usr/bin/bash
#SBATCH -J tajima_d_pooled
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=2G
#SBATCH -t 00:30:00
#SBATCH -p jro0014_amd
#SBATCH --mail-type=ALL
#SBATCH --mail-user=npb0015@auburn.edu

module load python/anaconda/3.8.6
module load htslib/1.17

#mkdir -p tmp
#mkdir -p data

#find ../01_simulating-test-data/data/simulated_var_only -type f > tmp/vcf_var_only.txt
#find ../01_simulating-test-data/data/simulated_invar -type f > tmp/vcf_invar.txt
#find ../01_simulating-test-data/data/simulated_missing_sites -type f > tmp/vcf_missing_sites.txt
find ../01_simulating-test-data/data/simulated_missing_genos -type f > tmp/vcf_missing_genos.txt
#find ../01_simulating-test-data/data/accuracy_invar -type f > tmp/vcf_accuracy.txt

#rm -r data/var_only
#mkdir -p data/var_only

#while read vcf
#do

#mkdir -p tmp/1
#rm -r tmp/1
#vcfslug=$(echo $vcf | sed 's/.*\///g')
#echo $vcf
#bgzip $vcf
#tabix -p vcf ${vcf}.gz
#python ../../../../pixy/pixy/__main__.py --stats tajima_d --vcf ${vcf}.gz --window_size 10000 --populations populations_pi.txt --bypass_invariant_check yes --output_folder data/var_only/ --output_prefix $vcfslug

#rm -r tmp/1

#done < tmp/vcf_var_only.txt 

##rm -r data/invar
#mkdir -p data/invar

#while read vcf
#do

#mkdir -p tmp/1
#rm -r tmp/1
#vcfslug=$(echo $vcf | sed 's/.*\///g')
#echo $vcf
#bgzip $vcf
#tabix -p vcf ${vcf}.gz
#python ../../../../pixy/pixy/__main__.py --stats tajima_d --vcf ${vcf}.gz --window_size 10000 --populations populations_pi.txt --bypass_invariant_check yes --output_folder data/invar/ --output_prefix $vcfslug

#rm -r tmp/1  

#done < tmp/vcf_invar.txt

##rm -r data/missing_sites
#mkdir -p data/missing_sites

#while read vcf
#do

#mkdir -p tmp/1
#rm -r tmp/1
#vcfslug=$(echo $vcf | sed 's/.*\///g')
#echo $vcf
#bgzip $vcf
#tabix -p vcf ${vcf}.gz
#python ../../../../pixy/pixy/__main__.py --stats tajima_d --vcf ${vcf}.gz --window_size 10000 --populations populations_pi.txt --bypass_invariant_check yes --output_folder data/missing_sites/ --output_prefix $vcfslug

#rm -r tmp/1  

#done < tmp/vcf_missing_sites.txt


##rm -r data/missing_genos
#mkdir -p data/missing_genos

while read vcf
do

mkdir -p tmp/1
#rm -r tmp/1
vcfslug=$(echo $vcf | sed 's/.*\///g')
#echo $vcf
if ( file $vcf | grep -v "compressed" )
then 
bgzip $vcf
tabix -p vcf ${vcf}.gz
python ../../../../pixy/pixy/__main__.py --stats tajima_d --vcf ${vcf}.gz --window_size 10000 --populations populations_pi.txt --bypass_invariant_check yes --output_folder data/missing_genos/ --output_prefix $vcfslug
fi

#rm -r tmp/1   

done < tmp/vcf_missing_genos.txt



#mkdir -p data/accuracy_invar
#while read vcf
#do

#mkdir -p tmp/1
#rm -r tmp/1
#vcfslug=$(echo $vcf | sed 's/.*\///g')
#echo $vcf

#python pixy.py --stats tajima_d --vcf $vcf --zarr_path tmp/1 --window_size 10000 --populations populations_pi.txt --bypass_invariant_check yes --output_prefix data/accuracy_invar/$vcfslug

#rm -r tmp/1   

#done < tmp/vcf_accuracy.txt

