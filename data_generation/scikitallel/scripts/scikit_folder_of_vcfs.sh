#!/bin/bash

while read vcf
do

echo $vcf

#python scripts/compute_dxy_scikit-allel.py --vcf $vcf
#python scripts/compute_pi_scikit-allel.py --vcf $vcf
python scripts/compute_watterson_theta_scikit-allel.py --vcf $vcf
python scripts/compute_tajima_d_scikit-allel.py --vcf $vcf

done < $1

