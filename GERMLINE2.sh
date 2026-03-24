#!/bin/bash
#SBATCH --job-name=cbriggsae_ibd
#SBATCH --export=ALL
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=D.E.Joia@sms.ed.ac.uk
#SBATCH --output=/mnt/loki/hartfield/caenorhabditis/scripts/output/%x_%j.out
#SBATCH --error=/mnt/loki/hartfield/caenorhabditis/scripts/error/%x_%j.err

# Stop if any command fails
set -e

# MAKE SURE TO ACTIVATE AN ENVIRONMENT WITH VCFTOOLS INSTALLED BEFORE RUNNING SCRIPT

# Move to working directory
cd /mnt/loki/hartfield/caenorhabditis/analyses/Darcey/Cbriggsae

# Define tool paths
VCFTOOLS=vcftools
SHAPEIT=/mnt/loki/hartfield/caenorhabditis/scripts/Darcey/shapeit
G2=/mnt/loki/hartfield/caenorhabditis/scripts/Darcey/g2

# Loop over chromosomes
for chr in I II III IV V X
do

# to keep track of progress:
echo "Processing chromosome $chr"

  # Filter VCF
  $VCFTOOLS --gzvcf australia_${chr}.vcf.gz \
    --max-missing 0.9 \
    --recode --recode-INFO-all \
    --out australia_${chr}_filtered

  # Phase with SHAPEIT
  $SHAPEIT \
    --input-vcf australia_${chr}_filtered.recode.vcf \
    -M maps/chr_${chr}.map \
    -O australia_${chr}_phased \
    --thread ${SLURM_CPUS_PER_TASK}

  # Run GERMLINE2
$G2 \
  -m 2 \
  australia_${chr}_phased.haps \
  australia_${chr}_phased.sample \
  maps/chr_${chr}.map \
  australia_${chr}_ibd

done

