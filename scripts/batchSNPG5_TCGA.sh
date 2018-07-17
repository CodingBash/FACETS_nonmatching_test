#!/bin/bash
#$ -S /bin/bash
#$ -t 1
#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]"
#$ -l vf=64G
#$ -l h_vmem=64G
#$ -pe threads 1

source /mnt/wigclust19/data/safe/deschene/software/lmod/Lmod-7.7/init_wigclust
module load zlib/1.2.11
module load bzip2/1.0.6
module load xz/5.2.4
module load pcre/8.42
module load curl/7.60.0
module load htslib/1.8
module load facets/0.5.14


export PATH_DATA=/mnt/wigclust20/data/unsafe/belleau/TCGAExample
export PATH_REF=/mnt/wigclust19/data/unsafe/belleau

cd ${PATH_DATA}

NBCOV=""
LISTNAME=""
COVMIN=20

for LISTFILES in `cat ${PATH_DATA}/listBam.txt `
do
# LISTFILES=${PATH_DATA}/listFastqR1.txt
FILECUR=$(echo $LISTFILES)
sample=$(echo $LISTFILES |perl -pe 's/^C546.TCGA-([A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+)-.+$/$1/g')
ORGID=$(echo $FILECUR|perl -n -e '/^([A-Z0-9]+-[A-Z0-9]+)/;print $1')
TYPES=$(echo $sample|perl -n -e '/([A-Z0-9]+)$/;print $1')
if ! [ "$NBCOV" = "" ];then
NBCOV=${NBCOV},
LISTNAME=${LISTNAME}" "
fi
if [ "$TYPES" = "10A" ] || [ "$TYPES" = "11A" ]; then
NBCOV=${NBCOV}${COVMIN}
else
NBCOV=${NBCOV}0
fi
echo $TYPES $sample
LISTNAME=${LISTNAME}${PATH_DATA}/$FILECUR
done

snp-pileup -d 1000 -r${NBCOV}  /mnt/wigclust19/data/unsafe/belleau/dbSNP/20180418/dbSnp_20180418Hg38_G5.vcf snpTCGATest20G5.csv ${LISTNAME}

