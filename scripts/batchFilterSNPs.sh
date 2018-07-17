#!/bin/bash
#$ -S /bin/bash
#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]"
#$ -l vf=4G
#$ -l h_vmem=4G
#$ -pe threads 1

source /mnt/wigclust19/data/safe/deschene/software/lmod/Lmod-7.7/init_wigclust
module load RBio/3.5.0 

export PATH_DATA=/mnt/wigclust17/data/unsafe/belleau/TCGAExample
export PATH_REF=/mnt/wigclust19/data/unsafe/belleau
cd $PATH_DATA
Rscript filterSNPs.r \
>filterSnps.stdout \
2>filterSnps.stderr
done


