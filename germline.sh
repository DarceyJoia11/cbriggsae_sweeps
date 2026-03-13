#!/bin/bash
#SBATCH --job-name=germline_run
#SBATCH --export=ALL
#SBATCH --mem=16gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=D.E.Joia@sms.ed.ac.uk
#SBATCH --output=/mnt/loki/hartfield/caenorhabditis/scripts/output/germline_%a.out
#SBATCH --error=/mnt/loki/hartfield/caenorhabditis/scripts/error/germline_%a.err
#SBATCH --array=1-6

set -e

cd /mnt/loki/hartfield/caenorhabditis/analyses/Darcey/Cbriggsae

# Pick the chromosome for this job
CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" chrom_list.txt)

# GERMLINE input files
MAP_FILE="aus_india_${CHR}_germline.map"
PED_FILE="aus_india_${CHR}_germline.ped"
OUT_PREFIX="germline_${CHR}_result"

# Run GERMLINE
# minimum length of 2cM/Mb, no more than 2 SNPs between strains
# treat strains as haploid
printf "1\n${MAP_FILE}\n${PED_FILE}\n${OUT_PREFIX}\n" | \
/mnt/loki/hartfield/caenorhabditis/scripts/Darcey/germline \
-min_m 2 \
-err_hom 2 \
-err_het 2 \
-bits 128 \
-haploid
