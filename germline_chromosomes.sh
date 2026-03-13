#!/bin/bash
#SBATCH --job-name=germline_chromosomes
#SBATCH --export=ALL
#SBATCH --mem=16gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=D.E.Joia@sms.ed.ac.uk
#SBATCH --output=/mnt/loki/hartfield/caenorhabditis/scripts/output/%x_%a.out
#SBATCH --error=/mnt/loki/hartfield/caenorhabditis/scripts/error/%x_%a.err
#SBATCH --array=1-6

set -e

cd /mnt/loki/hartfield/caenorhabditis/analyses/Darcey/Cbriggsae

# Create a VCF with just Australian and Indian strains using bcftools.
# Use the SNPs-only file (there should be no structural variants in GERMLINE)
bcftools view -S aus_india_keep.txt -Oz -o aus_india.vcf.gz WI.biallelic.snps.isotypes.vcf.gz
bcftools index aus_india.vcf.gz

# Determine chromosome for this array task
# This assigns each task ID to the relevant chromosome from my text file
CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" chrom_list.txt)

# CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" chrom_list.txt)Subset chromosome from aus_india.vcf.gz
bcftools view -r $CHR -Oz -o aus_india_${CHR}.vcf.gz aus_india.vcf.gz
bcftools index aus_india_${CHR}.vcf.gz

# Convert to PED/MAP for GERMLINE
plink --vcf aus_india_${CHR}.vcf.gz --recode --allow-extra-chr --out aus_india_${CHR}_germline
