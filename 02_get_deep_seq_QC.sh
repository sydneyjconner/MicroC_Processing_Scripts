#!/bin/bash
#SBATCH --job-name=get_QC
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --output=logs/get_QC_%A_%a.log
#SBATCH --error=logs/get_QC_%A_%a.err

# 100g and 24 cpus 
module load BWA
module load SAMtools
module load Python/3.12.3-GCCcore-13.3.0
module load bcl2fastq2

# Define paths
REF_GENOME="/gpfs/commons/home/cangel/g2lab/resources/GRCh38/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"
CHROM_SIZES="/gpfs/commons/home/sliu/resources/hg38.genome"
INPUT_DIR="/gpfs/commons/groups/gursoy_lab/sconner/ALS_DeepSeq/output/01_combined_lanes_data"
OUTPUT_DIR="/gpfs/commons/groups/gursoy_lab/sconner/ALS_DeepSeq/output/02_quality_check_tables"

# Get sample directories
SAMPLE_LIST=($(ls -d "$INPUT_DIR"/Sample_*))
SAMPLE_DIR=${SAMPLE_LIST[$SLURM_ARRAY_TASK_ID]}
SAMPLE_NAME=$(basename "$SAMPLE_DIR")

# Define input FASTQ files
R1_FILE="$SAMPLE_DIR/${SAMPLE_NAME}_combined_R1.fastq.gz"
R2_FILE="$SAMPLE_DIR/${SAMPLE_NAME}_combined_R2.fastq.gz"

# Define downstream output directory
DOWNSTREAM_DIR="$OUTPUT_DIR/$SAMPLE_NAME"
mkdir -p "$DOWNSTREAM_DIR"
mkdir -p "$DOWNSTREAM_DIR/other_files"
mkdir -p "$DOWNSTREAM_DIR/temp"

# Intermediate file paths
MAPPED_PAIRS="$DOWNSTREAM_DIR/other_files/${SAMPLE_NAME}_mapped.pairs"
MAPPED_PT_BAM="$DOWNSTREAM_DIR/other_files/${SAMPLE_NAME}_mapped.PT.bam"
DEDUP_STATS="$DOWNSTREAM_DIR/${SAMPLE_NAME}_stats.txt"
QC_OUTPUT="$DOWNSTREAM_DIR/${SAMPLE_NAME}_output.csv"

# 1. BWA-MEM, pairtools parse + sort + split 
echo  "Running full pipeline with cleanup..."

bwa mem -5SP -T0 -t16 "$REF_GENOME" "$R1_FILE" "$R2_FILE" | \
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 \
    --nproc-in 8 --nproc-out 8 --chroms-path "$CHROM_SIZES" 2> "$DOWNSTREAM_DIR/other_files/parse.log" | \
pairtools sort --nproc 16 --tmpdir="$DOWNSTREAM_DIR/temp" 2> "$DOWNSTREAM_DIR/other_files/sort.log" | \
pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups \
    --output-stats "$DEDUP_STATS" 2> "$DOWNSTREAM_DIR/other_files/dedup.log" | \
tee "$MAPPED_PAIRS" | \
pairtools split --nproc-in 8 --nproc-out 8 --output-pairs /dev/null --output-sam - 2> "$DOWNSTREAM_DIR/other_files/split.log" | \
samtools view -bS -@16 - 2> "$DOWNSTREAM_DIR/other_files/samtools_view.log" | \
samtools sort -@16 -m 4G -o "$MAPPED_PT_BAM" 2> "$DOWNSTREAM_DIR/other_files/samtools_sort.log"

# Verify BAM created
if [ ! -s "$MAPPED_PT_BAM" ]; then
    echo "BAM file not created for $SAMPLE_NAME"
    exit 1
fi

# Index BAM
samtools index "$MAPPED_PT_BAM"

# Run QC
if [ -f "$DEDUP_STATS" ]; then
    echo "Running QC script..."
    python /gpfs/commons/groups/gursoy_lab/sliu/QC_MicroC/Micro-C/modified_get_qc.py -p "$DEDUP_STATS" -o "$QC_OUTPUT"
else
    echo "Skipping QC — stats file missing"
fi

echo "Finished processing $SAMPLE_NAME"
