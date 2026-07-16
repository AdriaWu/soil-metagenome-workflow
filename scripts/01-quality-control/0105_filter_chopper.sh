#!/bin/bash
#PBS -P <PROJECT>
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=32GB
#PBS -l jobfs=400GB
#PBS -l ncpus=8
#PBS -l storage=scratch/<PROJECT>+gdata/<PROJECT>
#PBS -l wd

set -euo pipefail
set -x

# ============================================================================
# Chopper Read Filtering
# ============================================================================
# Description: Quality and length filtering for Oxford Nanopore reads
# Input: FASTQ files (gzipped or uncompressed)
# Output: Filtered FASTQ files (gzipped)
# ============================================================================

# Activate conda environment
source <CONDA_PATH>/etc/profile.d/conda.sh
conda activate <QC_ENV_PATH>

# Configuration
NTHREADS=${PBS_NCPUS:-8}

# Filtering parameters (MODIFY AS NEEDED)
MIN_LENGTH=2000      # Minimum read length (bp)
MIN_QUALITY=18       # Minimum average quality score (Phred)

# Input/Output paths (MODIFY THESE)
INPUT_DIR="<INPUT_FASTQ_DIR>"
OUT_DIR="<OUTPUT_FILTERED_DIR>"

mkdir -p "$OUT_DIR"

# Local temporary directory on jobfs
LOCAL_DIR="${PBS_JOBFS}/chopper_work"
mkdir -p "$LOCAL_DIR"

# Process each FASTQ file
for input_file in "$INPUT_DIR"/*.fastq.gz "$INPUT_DIR"/*.fq.gz; do
    [ -f "$input_file" ] || continue
    
    # Extract sample name
    sample=$(basename "$input_file" .fastq.gz)
    sample=${sample%.fq.gz}
    
    output_file="${OUT_DIR}/${sample}_filtered.fastq.gz"
    
    echo "Filtering $sample..."
    echo "  Input: $input_file"
    echo "  Output: $output_file"
    echo "  Parameters: length >=${MIN_LENGTH}bp, quality >=${MIN_QUALITY}"
    
    # Run Chopper with original options
    chopper \
        --input "$input_file" \
        --quality "$MIN_QUALITY" \
        --minlength "$MIN_LENGTH" \
        --threads "$NTHREADS" \
        | gzip > "$output_file"
    
    # Report filtering statistics
    input_reads=$(zcat "$input_file" | awk 'NR%4==1' | wc -l)
    output_reads=$(zcat "$output_file" | awk 'NR%4==1' | wc -l)
    retention=$((100 * output_reads / input_reads))
    
    echo "  Input reads: $input_reads"
    echo "  Output reads: $output_reads"
    echo "  Retention: ${retention}%"
    echo ""
done

echo "Chopper filtering completed at $(date) on host $(hostname)"