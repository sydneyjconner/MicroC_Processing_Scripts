#!/bin/bash

# Define the sample input directory
INPUT_DIR="/gpfs/commons/groups/gursoy_lab/sconner/ALS_DeepSeq/output/01_combined_lanes_data/ALS5"

# Count the number of Sample_* directories
NUM_SAMPLES=$(ls -d "$INPUT_DIR"/Sample_* 2>/dev/null | wc -l)

if [ "$NUM_SAMPLES" -eq 0 ]; then
    echo "❌ No Sample_* directories found in $INPUT_DIR"
    exit 1
fi

# Submit the array job with correct range
echo "✅ Found $NUM_SAMPLES samples. Submitting array job..."
sbatch --array=0-$(($NUM_SAMPLES - 1)) /gpfs/commons/groups/gursoy_lab/sconner/ALS_DeepSeq/scripts/02_get_deep_seq_QC.sh
#sbatch /gpfs/commons/groups/gursoy_lab/sconner/ALS_DeepSeq/scripts/02_get_deep_seq_QC.sh
