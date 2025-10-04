#!/bin/bash
#SBATCH --job-name=hic_compare_test
#SBATCH --time=72:00:00
#SBATCH --output=./logs/test_hic_compare_R.out
#SBATCH --error=./logs/test_hic_compare_R.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

#Usage: sbatch 04_compare_hic_compare.sh T6182 SDO416 500000

SAMPLE1=$1
SAMPLE2=$2
RESOLUTION=$3

module load R
Rscript /gpfs/commons/groups/gursoy_lab/sconner/ALS_DeepSeq/scripts/comparative_analyses_table.R "$SAMPLE1" "$SAMPLE2" "$RESOLUTION"
