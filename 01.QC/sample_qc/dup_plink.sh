#!/bin/bash
#$ -l h_vmem=1.7G
#$ -pe smp 30
#$ -binding linear:30
#$ -l h_rt=10:00:00

source /broad/software/scripts/useuse
PLINK2="/stanley/huang_lab/home/myu/software_tools/plink2"
WORKING_DIR="/stanley/huang_lab/home/myu/IBD/preprocess"

# # Perform LD pruning: autosomal SNPs
# $PLINK2 \
# --pfile ${WORKING_DIR}/geno_data/1.plink_format/gnomAD.v4.ibd_subset \
# --memory 30000 \
# --autosome \
# --geno 0.02 \
# --maf 0.05 \
# --snps-only just-acgt \
# --indep-pairwise 200 100 0.1 \
# --make-bed \
# --out ${WORKING_DIR}/QC/dup_check/pruned 

use PLINK
# if [ !-f ${pop}_pbk_btqc_mgqc-ibd-min0.125.genome.gz ]; then
# Estimate IBD
plink \
--bfile ${WORKING_DIR}/QC/dup_check/pruned \
--memory 510000 \
--genome gz \
--min 0.25 \
--out ${WORKING_DIR}/QC/dup_check/ibd-min0.25