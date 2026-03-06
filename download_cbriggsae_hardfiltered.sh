#!/bin/sh

# Script to download nematode data

#SBATCH --job-name=nematode_download
#SBATCH --export=ALL
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=D.E.Joia@sms.ed.ac.uk
#SBATCH --output=/mnt/loki/hartfield/caenorhabditis/scripts/output/%x_%a.out
#SBATCH --error=/mnt/loki/hartfield/caenorhabditis/scripts/error/%x_%a.err

# Download hard-filtered C. briggsae variants
wget -P /mnt/loki/hartfield/caenorhabditis/raw/Cbriggsae_hardfiltered2026/ \
https://caendr-open-access-data-bucket.s3.us-east-2.amazonaws.com/dataset_release/c_briggsae/20250626/variation/WI.20250626.hard-filter.vcf.gz
