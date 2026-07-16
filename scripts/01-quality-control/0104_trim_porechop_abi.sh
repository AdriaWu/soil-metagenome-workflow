#!/bin/bash
#PBS -P <PROJECT>
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=48GB
#PBS -l jobfs=200GB
#PBS -l ncpus=12
#PBS -l storage=scratch/<PROJECT>+gdata/<PROJECT>
#PBS -l wd

set -euo pipefail
set -x

# ============================================================================
# Porechop_ABI - Adapter Trimming & Chimera Splitting
# ============================================================================
# Description: Remove Oxford Nanopore adapters and SPLIT reads at internal
#              (mid-read) adapters. Internal adapters mark chimeric /
#              concatenated molecules (two library fragments read through as
#              one); left in place they are assembled straight through and
#              show up as internal-adapter contamination at NCBI submission.
#
# Why this step exists:
#   Dorado --trim (used at basecalling) removes adapters only from read ENDS.
#   It does NOT split reads at internal adapters, and its chimera splitter is
#   not fully sensitive. Chopper (next step) filters on quality/length only
#   and never touches adapters. Porechop_ABI fills that gap by finding
#   adapters anywhere in the read and, by default, splitting the read where an
#   internal adapter is found.
#
# Placement: AFTER Dorado basecalling/demux (00-basecalling), BEFORE Chopper
#            quality/length filtering (0105_filter_chopper.sh). Keep Dorado
#            --trim all as-is; the two steps are complementary.
#
# Input:  Per-barcode FASTQ (gzipped), one file per sample
# Output: Adapter-trimmed, chimera-split FASTQ (gzipped)
# ============================================================================

# Activate conda environment (Porechop_ABI installed into the QC env)
source <CONDA_PATH>/etc/profile.d/conda.sh
conda activate <QC_ENV_PATH>

# Configuration
NTHREADS=${PBS_NCPUS:-12}

# Post-split fragment floor (bp). Pieces smaller than this after an internal
# split are discarded. Kept below the downstream Chopper length floor (2000 bp),
# so Chopper remains the effective minimum-length gate.
MIN_SPLIT_READ_SIZE=1000

# Input/Output paths (MODIFY THESE)
INPUT_DIR="<INPUT_FASTQ_DIR>"      # demuxed per-barcode FASTQ (.fastq.gz)
OUT_DIR="<OUTPUT_ADAPTERTRIM_DIR>" # adapter-trimmed output

mkdir -p "$OUT_DIR"

# ---------------------------------------------------------------------------
# NOTE ON SCALE
# The loop below processes every FASTQ in INPUT_DIR serially in one job.
# That is fine for a handful of samples / testing. For a full cohort
# (hundreds of samples, ~100s of Gbp) Porechop_ABI is the slow step, so run
# this as a PBS ARRAY instead — one sample per task, in parallel:
#
#   #PBS -J 0-<N-1>              # N = number of FASTQ files
#   mapfile -t FILES < <(ls "$INPUT_DIR"/*.fastq.gz)
#   input_file="${FILES[$PBS_ARRAY_INDEX]}"
#   ... (run the single-sample block below on "$input_file") ...
#
# Submit the array from this directory so .o/.e logs land beside the script.
# ---------------------------------------------------------------------------

for input_file in "$INPUT_DIR"/*.fastq.gz "$INPUT_DIR"/*.fq.gz; do
    [ -f "$input_file" ] || continue

    sample=$(basename "$input_file" .fastq.gz)
    sample=${sample%.fq.gz}

    output_file="${OUT_DIR}/${sample}_adaptertrim.fastq.gz"

    # Skip logic: don't redo a sample that already finished
    if [ -s "$output_file" ]; then
        echo "SKIP $sample — output already exists: $output_file"
        continue
    fi

    echo "Adapter-trimming $sample..."
    echo "  Input:  $input_file"
    echo "  Output: $output_file"

    # Flags:
    #   -abi                  ab-initio adapter discovery (in addition to the
    #                         built-in adapter set); also reports if a file is
    #                         already trimmed. Drop it for speed if ab-initio
    #                         consistently reports "already trimmed".
    #   (default)             reads with an internal adapter are SPLIT. Do NOT
    #                         pass --no_split (keeps chimeras) or --discard_middle
    #                         (throws the whole read away instead of splitting).
    porechop_abi \
        -abi \
        -i "$input_file" \
        -o "$output_file" \
        -t "$NTHREADS" \
        --min_split_read_size "$MIN_SPLIT_READ_SIZE"

    # Report read counts (out > in => chimeras were split; out < in => sub-floor
    # fragments dropped). Both are expected.
    in_reads=$(zcat "$input_file" | awk 'NR%4==1' | wc -l)
    out_reads=$(zcat "$output_file" | awk 'NR%4==1' | wc -l)
    echo "  Input reads:  $in_reads"
    echo "  Output reads: $out_reads"
    echo ""
done

echo "Porechop_ABI adapter trimming completed at $(date) on host $(hostname)"
