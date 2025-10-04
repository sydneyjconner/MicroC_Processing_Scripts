#!/bin/bash
#SBATCH --job-name=juicer_hic
#SBATCH --output=./logs/04_juicer_hic_%j.out
#SBATCH --error=./logs/04_juicer_hic_%j.err
#SBATCH --time=60:00:00
#SBATCH --mem=240G
#SBATCH --cpus-per-task=24


module load java

SAMPLE="$1"
CHROM_SIZES="/gpfs/commons/groups/gursoy_lab/sliu/QC_MicroC/archived.scripts/ALS_hic/sorted.hg38.genome"


INPUT_DIR="/gpfs/commons/home/sliu/ALS_MicroC_analyses/output/03_mapped_pairs_per_sample/Sample_${SAMPLE}_mapped_pairs"
#BASE_DIR="/gpfs/commons/groups/gursoy_lab/sliu/ALS_hic/Sample_${SAMPLE_NAME}_files/"
INPUT_FILE="${INPUT_DIR}/mapped.pairs"

OUTPUT_DIR="/gpfs/commons/home/sliu/ALS_MicroC_analyses/output/04_hic_per_sample/Sample_${SAMPLE}"
mkdir -p $OUTPUT_DIR
OUTPUT_FILE="${OUTPUT_DIR}/Sample_${SAMPLE}.hic"

java -Xmx240000m  -Djava.awt.headless=true -jar /nfs/sw/easybuild/software/juicer/1.22.01/juicer_tools_1.22.01.jar pre --threads 24 "$INPUT_FILE" "$OUTPUT_FILE" "$CHROM_SIZES"
