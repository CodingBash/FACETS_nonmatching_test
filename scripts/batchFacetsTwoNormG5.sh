#!/bin/bash
#$ -S /bin/bash
#$ -t 1-11
#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]"
#$ -l vf=64G
#$ -l h_vmem=64G
#$ -pe threads 1

source /mnt/wigclust19/data/safe/deschene/software/lmod/Lmod-7.7/init_wigclust
module load RBio/3.5.0 

export PATH_DATA=/mnt/wigclust17/data/unsafe/belleau/targetSeq
export PATH_REF=/mnt/wigclust19/data/unsafe/belleau
VALR=$(echo $(($SGE_TASK_ID - 1)))
cd $PATH_DATA
# base between 0 and 10
>facetsTwoNormadat.log
refCur=$(echo $(($VALR + 40))) 
for base in `seq 1 10`
do
CUR=$(echo $((1 + 4 * (40 + (($base + $VALR) % 11)) )))
REF=$(echo $((1 + 4 * (40 + $VALR))))
sampleCur=$(echo $(((($base + $VALR) % 11) + 40)) )
paste -d$',' <(cut -d$',' -f1,2,3,4 snpOrgTarget20G5.csv) \
<(cut -d$',' -f$REF,$(( $REF + 1 )),$(( $REF + 2 )),$(( $REF + 3 )) snpOrgTarget20G5.csv) \
<(cut -d$',' -f$CUR,$(( $CUR + 1 )),$(( $CUR + 2 )),$(( $CUR + 3 )) snpOrgTarget20G5.csv) \
>snpTmpG5_${refCur}_${sampleCur}.csv
date >>facetsTwoNormadat.log
Rscript facetsTwo.r ${refCur} ${sampleCur} \
>facetsTest_${refCur}_${sampleCur}.stdout \
2>facetsTest_${refCur}_${sampleCur}.stderr
done


