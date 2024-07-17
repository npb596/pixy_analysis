#!/bin/bash

#SBATCH -J download_mosquito_bams
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=10G
#SBATCH -t 3-00:00:00
#SBATCH -p normal
#SBATCH -e download_mosquito_bams.err
#SBATCH -o download_mosquito_bams.out

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1695535/AB0283-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1695536/AB0284-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1695498/AB0244-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1695526/AB0274-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1695530/AB0278-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1695517/AB0265-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1695520/AB0268-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1695492/AB0238-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1695460/AB0205-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1695461/AB0206-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1695462/AB0207-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1695463/AB0208-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1695466/AB0211-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1695472/AB0217-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1695453/AB0198-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1695372/AB0103-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1695409/AB0147-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1697155/AK0065-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1697157/AK0067-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1697182/AK0095-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1697175/AK0088-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1697174/AK0087-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1697173/AK0086-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1697166/AK0077-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1697196/AK0110-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1697190/AK0104-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1697194/AK0108-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1697184/AK0098-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1697183/AK0096-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1697180/AK0093-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1697163/AK0073-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1697200/AK0127-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1697162/AK0072-C.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ169/ERZ1697181/AK0094-C.bam
