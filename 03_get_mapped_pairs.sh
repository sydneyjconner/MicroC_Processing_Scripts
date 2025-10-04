#!/bin/bash
#SBATCH --job-name=mapped_pairs
#SBATCH --time=120:00:00
#SBATCH --output=./logs/03_get_mapped_pairs_%j.out
#SBATCH --error=./logs/03_get_mapped_pairs_%j.err
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G


module load BWA
module load SAMtools
module load Python/3.12.3-GCCcore-13.3.0

# ----------- INPUT ARGUMENT ----------- #
SAMPLE="$1"
if [[ -z "$SAMPLE" ]]; then
  echo "Error: No sample ID provided."
  exit 1
fi

# ----------------- PATHS ----------------#


SAMPLE="$1"

REF_GENOME="/gpfs/commons/home/cangel/g2lab/resources/GRCh38/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"

CHROM_SIZES="/gpfs/commons/home/sliu/resources/hg38.genome"

INPUT_DIR="/gpfs/commons/home/sliu/ALS_MicroC_analyses/output/01_combined_lanes_data"

#test
#INPUT_DIR="/gpfs/commons/home/sliu/ALS_MicroC_analyses/SXL_070925_test_output/01_combined_lanes_data"

OUTPUT_DIR="/gpfs/commons/home/sliu/ALS_MicroC_analyses/output/03_mapped_pairs_per_sample/Sample_${SAMPLE}_mapped_pairs"

mkdir -p "$OUTPUT_DIR"

# first library
#L1_R1_FILE="$FASTQ_DIR/$SAMPLE/${SAMPLE}E2_combined_R1.fastq.gz"
#L1_R2_FILE="$FASTQ_DIR/$SAMPLE/${SAMPLE}E2_combined_R2.fastq.gz"

# second library
#L2_R1_FILE="$FASTQ_DIR/$SAMPLE/${SAMPLE}E3_combined_R1.fastq.gz"
#L2_R2_FILE="$FASTQ_DIR/$SAMPLE/${SAMPLE}E3_combined_R2.fastq.gz"

# ----------- FIND MATCHING LIBRARIES ----------- #
# Find all R1 and R2 files for this sample across multiple libraries
R1_FILES=()
R2_FILES=()

for fq_dir in "$INPUT_DIR"/Sample_${SAMPLE}*/; do
    if [[ -d "$fq_dir" ]]; then
        for fq_file in "$fq_dir"/*_combined_R1.fastq.gz; do
            [[ -f "$fq_file" ]] && R1_FILES+=("$fq_file")
        done
        for fq_file in "$fq_dir"/*_combined_R2.fastq.gz; do
            [[ -f "$fq_file" ]] && R2_FILES+=("$fq_file")
        done
    fi
done

# --------------  Sanity check -------------------------------
if [[ ${#R1_FILES[@]} -eq 0 || ${#R2_FILES[@]} -eq 0 ]]; then
    echo "Error: No R1 or R2 files found for sample $SAMPLE."
    exit 1
fi
echo "Found ${#R1_FILES[@]} R1 files and ${#R2_FILES[@]} R2 files for $SAMPLE."


# Get the directory this script is in
#SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# Use it to define a TMP_DIR inside this script's directory
#TMP_DIR="${SCRIPT_DIR}/tmp"

# Later in the pipeline:
# pairtools sort --tmpdir "$TMP_DIR"

# Pipeline without intermediate files written
bwa mem -5SP -T0 -t 24 "$REF_GENOME" \
    <(zcat "${R1_FILES[@]}") \
    <(zcat "${R2_FILES[@]}") |
pairtools parse \
    --min-mapq 40 \
    --walks-policy 5unique \
    --max-inter-align-gap 30 \
    --nproc-in 8 --nproc-out 8 \
    --chroms-path "$CHROM_SIZES" |
pairtools sort \
    --nproc 24 --tmpdir=./ |
pairtools dedup \
    --nproc-in 8 --nproc-out 8 \
    --mark-dups \
    --output-stats stats.txt |pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats "${OUTPUT_DIR}/stats.txt"|pairtools split --nproc-in 8 --nproc-out 8 --output-pairs "${OUTPUT_DIR}/mapped.pairs" --output-sam -|samtools view -bS -@24 | samtools sort -@24 -o "${OUTPUT_DIR}/mapped.PT.bam"

# Index the final BAM
samtools index "${OUTPUT_DIR}/mapped.PT.bam"



# pipeline step by step

# bwa mem -5SP -T0 -t48 "$REF_GENOME" <(zcat "$L1_R1_FILE" "$L2_R1_FILE") <(zcat "$L1_R2_FILE" "$L2_R2_FILE") -o aligned.sam

# pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 24 --nproc-out 24 --chroms-path "$CHROM_SIZES" aligned.sam >  parsed.pairsam

# pairtools sort --nproc 24 --tmpdir=/gpfs/commons/groups/gursoy_lab/sliu/ALS_hic  parsed.pairsam > sorted.pairsam

# pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt --output dedup.pairsam sorted.pairsam

# pairtools split --nproc-in 8 --nproc-out 8 --output-pairs mapped.pairs --output-sam unsorted.bam dedup.pairsam

# samtools sort -@24 -T /gpfs/commons/groups/gursoy_lab/sliu/ALS_hic -o mapped.PT.bam unsorted.bam

# samtools index mapped.PT.bam
#pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path "$CHROM_SIZES"  | pairtools sort --tmpdir=/gpfs/commons/groups/gursoy_lab/sliu/ALS_hic --nproc 24|pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt|pairtools split --nproc-in 8 --nproc-out 8 --output-pairs mapped.pairs --output-sam -|samtools view -bS -@24 | samtools sort -@24 -o mapped.PT.bam
#samtools index mapped.PT.bam




