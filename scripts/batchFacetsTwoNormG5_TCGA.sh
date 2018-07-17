#!/bin/bash
#$ -S /bin/bash
#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]"
#$ -l vf=64G
#$ -l h_vmem=64G
#$ -pe threads 1

source /mnt/wigclust19/data/safe/deschene/software/lmod/Lmod-7.7/init_wigclust
module load RBio/3.5.0 
module load facets/0.5.14

export PATH_DATA=/mnt/wigclust20/data/unsafe/belleau/TCGAExample
export PATH_REF=/mnt/wigclust19/data/unsafe/belleau

cd $PATH_DATA
>facetsTwoNormadat.log
for base in `seq 1 7`
do
date >>facetsTwoNormadat.log
Rscript facetsTwo_TCGA.r ${base}
>facetsTest_${base}.stdout \
2>facetsTest_${base}.stderr
done
