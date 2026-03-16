#!/bin/bash
#PBS -P <PROJECT>
#PBS -q hugemem
#PBS -l walltime=30:00:00
#PBS -l mem=1470GB
#PBS -l jobfs=400GB
#PBS -l ncpus=48
#PBS -l storage=scratch/<PROJECT>+gdata/<PROJECT>
#PBS -l wd

set -euo pipefail
set -x

# ============================================================================
# Map Assembled Contigs to SMAG + Plant Reference for Taxonomic Assignment
# ============================================================================
# Description: Map assembled contigs against the combined SMAG + plant genome
#              reference to assign each contig to a MAG. This produces the
#              mapped contigs CSV required by annotate_contigs.py (step 6).
#
# Input:  Assembled/polished contigs (FASTA, optionally gzipped)
# Output: Sorted BAM, idxstats, and mapped_contigs.csv (contig_id, MAG_ID)
#
# Note:   This uses the asm5 preset (assembly-to-reference alignment), NOT
#         the map-ont preset used in 02-alignment for read-level mapping.
#         --secondary=no retains only the best hit per contig.
# ============================================================================

# ---- Tool Setup ----
export PATH=<MINIMAP2_PATH>:$PATH
module load samtools/1.22

# ---- I/O Paths ----
# Assembled contigs (FASTA or FASTA.gz) from assembly/polishing step
CONTIG_FILE="<INPUT_CONTIGS_FASTA>"

# SMAG + plant genome reference (FASTA/FASTA.gz or prebuilt .mmi)
REF_FILE="<COMBINED_SMAG_PLANT_REFERENCE>"

# Output directory
OUTPUT_DIR="<OUTPUT_CONTIG_MAPPING_DIR>"

# ---- Parameters ----
THREADS=$PBS_NCPUS
TMPDIR="${TMPDIR:-/tmp}"

# ---- Setup ----
mkdir -p "$OUTPUT_DIR"

echo "============================================"
echo "Contig-to-SMAG Mapping"
echo "============================================"
echo "Contigs:   $CONTIG_FILE"
echo "Reference: $REF_FILE"
echo "Output:    $OUTPUT_DIR"
echo "Threads:   $THREADS"
echo "Started:   $(date)"
echo "============================================"

# ---- Decompress contigs if gzipped ----
if [[ "$CONTIG_FILE" == *.gz ]]; then
    echo "Decompressing contigs..."
    CONTIG_DECOMPRESSED="${OUTPUT_DIR}/$(basename "${CONTIG_FILE%.gz}")"
    gunzip -c "$CONTIG_FILE" > "$CONTIG_DECOMPRESSED"
    CONTIG_FILE="$CONTIG_DECOMPRESSED"
fi

# ---- Map contigs to reference ----
BASE_NAME=$(basename "$CONTIG_FILE" .fasta)
BASE_NAME=$(basename "$BASE_NAME" .fa)
OUTPUT_BAM="${OUTPUT_DIR}/${BASE_NAME}.sorted.bam"

echo "Mapping contigs with minimap2 (asm5 preset)..."

# asm5: assembly-to-reference alignment (≤5% divergence)
# --secondary=no: retain only the best hit per contig
# --split-prefix: required for large references to manage memory
minimap2 \
    -ax asm5 \
    -t "$THREADS" \
    --secondary=no \
    --split-prefix="${TMPDIR}/${BASE_NAME}.split" \
    "$REF_FILE" \
    "$CONTIG_FILE" \
| samtools sort -@ "$THREADS" -m 3G \
    -T "${TMPDIR}/${BASE_NAME}.sorttmp" \
    -o "$OUTPUT_BAM" -

echo "Indexing BAM..."
samtools index -c "$OUTPUT_BAM"

# ---- Generate mapping statistics ----
echo "Generating mapping statistics..."
samtools flagstat "$OUTPUT_BAM" > "${OUTPUT_BAM%.bam}.flagstat.txt"
samtools idxstats "$OUTPUT_BAM" > "${OUTPUT_BAM%.bam}.idxstats.txt"

# ---- Extract mapped contigs CSV ----
# Create the contig_id → MAG_ID mapping table required by annotate_contigs.py
# (step 7).
#
# Filters:
#   -q MIN_MAPQ: Minimum mapping quality threshold
#   -F 0x904:    Exclude unmapped (0x4), secondary (0x100), and
#                supplementary (0x800) alignments — primary only
#
# SMAG contig name parsing:
#   SMAG reference contig names encode the MAG ID before the first "k"
#   character (e.g., "12345k141_0" → MAG_ID "12345"). The awk block
#   strips the suffix from the first "k" onward, then removes any
#   remaining "k" characters to recover the numeric MAG ID.
#   >>> Adjust or remove this parsing if your reference uses a different
#   >>> naming convention. <<<

MIN_MAPQ=20
MAPPED_CSV="${OUTPUT_DIR}/mapped_contigs.csv"

echo "Extracting mapped contigs table (MAPQ ≥ ${MIN_MAPQ}, primary only)..."
printf "contig_id,MAG_ID,MAPQ\n" > "$MAPPED_CSV"

samtools view -@ "$THREADS" -q "$MIN_MAPQ" -F 0x904 "$OUTPUT_BAM" \
| awk '
BEGIN { FS="\t"; OFS="," }
{
    contig_id = $1
    mag_id    = $3
    mapq      = $5

    # Parse MAG ID from SMAG reference contig name:
    # Strip everything from the first "k" onward, then remove "k"
    sub(/k.*/, "k", mag_id)
    gsub("k", "", mag_id)

    printf "%s,%s,%s\n", contig_id, mag_id, mapq
}' >> "$MAPPED_CSV"

N_MAPPED=$(tail -n +2 "$MAPPED_CSV" | wc -l)
echo "Mapped contigs: $N_MAPPED"
echo "Output CSV: $MAPPED_CSV"

# ---- Cleanup ----
if [[ -v CONTIG_DECOMPRESSED ]]; then
    rm -f "$CONTIG_DECOMPRESSED"
fi

echo "============================================"
echo "Contig mapping complete"
echo "BAM:  $OUTPUT_BAM"
echo "CSV:  $MAPPED_CSV"
echo "Finished: $(date) on $(hostname)"
echo "============================================"