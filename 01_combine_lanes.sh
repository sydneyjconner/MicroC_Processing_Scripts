#!/bin/bash
#SBATCH --job-name=combine_lanes
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8  # Set reasonable CPUs per job
#SBATCH --mem=48G  # Reduce memory per job since more jobs are running
#SBATCH --array=1-10  
#SBATCH --output=logs/01_combine_lanes_%A_%a.log
#SBATCH --error=logs/01_combine_lanes_%A_%a.err

# Define paths 
RAW_FASTQ_DIR_NAME="$1"

INPUT_DIR="/gpfs/commons/groups/gursoy_lab/sconner/ALS_DeepSeq/data/${RAW_FASTQ_DIR_NAME}"

OUTPUT_ROOT_DIR="/gpfs/commons/groups/gursoy_lab/sconner/ALS_DeepSeq/output"
OUTPUT_DIR="${OUTPUT_ROOT_DIR}/01_combined_lanes_data"

mkdir -p "$OUTPUT_DIR"

# Get list of sample directories
SAMPLE_LIST=($(ls -d "$INPUT_DIR"/Sample_*))

# Calculate which job corresponds to which sample and read type
TOTAL_SAMPLES=${#SAMPLE_LIST[@]}  # Get total number of samples
SAMPLE_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) / 2 ))  # Determine sample index (0-based)
READ_TYPE=$(( (SLURM_ARRAY_TASK_ID - 1) % 2 ))  # Determine if processing R1 or R2

# Get the sample directory
SAMPLE_DIR=${SAMPLE_LIST[$SAMPLE_INDEX]}    
SAMPLE_NAME=$(basename "$SAMPLE_DIR")
FASTQ_DIR="$SAMPLE_DIR/fastq"

# Set read type
if [[ $READ_TYPE -eq 0 ]]; then
    READ_SUFFIX="R1"
else
    READ_SUFFIX="R2"
fi


# create the output folder for a sample
mkdir -p "$OUTPUT_DIR/$SAMPLE_NAME"

# generating the combined fastq at destinated output folder
zcat "$FASTQ_DIR"/*_L001_001.${READ_SUFFIX}.fastq.gz "$FASTQ_DIR"/*_L002_001.${READ_SUFFIX}.fastq.gz | gzip > "$OUTPUT_DIR/$SAMPLE_NAME/${SAMPLE_NAME}_combined_${READ_SUFFIX}.fastq.gz"






