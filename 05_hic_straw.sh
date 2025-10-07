#!/bin/bash
#SBATCH --job-name=straw_hic
#SBATCH --output=./logs/05_straw_hic_%j.out
#SBATCH --error=./logs/05_straw_hic_%j.err
#SBATCH --time=60:00:00
#SBATCH --mem=240G
#SBATCH --cpus-per-task=24

SAMPLE1=$1
RESOLUTION=$2

for i in {1..22}
    do
        ./straw NONE $SAMPLE1 $i $i BP $RESOLUTION > $SAMPLE1.NONE.chr$i.$RESOLUTION.txt
    done
    #./straw NONE GSM2795535_Rao-2017-HIC001_30.hic X X BP 500000 > HIC001/HIC001.NONE.chrX.500000.txt

#module load R
#Rscript /gpfs/commons/groups/gursoy_lab/sconner/ALS_DeepSeq/scripts/05_hic_straw.R "$SAMPLE1" "$RESOLUTION"