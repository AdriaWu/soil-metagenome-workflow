# Abundance Analysis Scripts

Scripts for mapping reads to assembly, differential abundance analysis, and taxonomic annotation.

## Prerequisites

- minimap2 (v2.30+)
- samtools (v1.19+)
- R (v4.0+) with packages: edgeR, tidyverse, ggplot2
- Python 3.7+ with pandas

## Workflow Overview
```
Assembled contigs + Filtered reads
    ↓
[1] map_reads_to_assembly.sh → BAM files (per sample)
    ↓
[2] build_count_matrix.sh → Count matrix (samples × contigs)
    ↓
[3] qc_and_filter_counts.R → QC plots + Filtered count matrix
    ↓
[4] run_edger_analysis.R → Differential abundance results
    ↓
[5] empirical_fdr_calibration.R → Calibrated p-value thresholds
    ↓
[6] map_contigs_to_smag.sh → Contig-to-MAG mapping (mapped_contigs.csv)
    ↓
[7] annotate_contigs.py → Taxonomic annotation
```

---

## Configuration

Replace placeholders in scripts:

- `<PROJECT>`: NCI project code
- `<MINIMAP2_PATH>`: minimap2 binary directory
- `<INPUT_CONTIGS_FASTA>`: Assembled/polished contigs
- `<INPUT_FILTERED_FASTQ_DIR>`: Quality-filtered reads from step 01
- `<OUTPUT_BAM_DIR>`: Directory for BAM files
- `<INPUT_METADATA_CSV>`: Sample metadata (barcode, treatment, etc.)
- `<COMBINED_SMAG_PLANT_REFERENCE>`: Combined SMAG + plant genome reference (FASTA/FASTA.gz or .mmi)
- `<OUTPUT_CONTIG_MAPPING_DIR>`: Directory for contig-to-SMAG mapping output
- `<SMAG_TAXONOMY_CSV>`: SMAG taxonomy database
- `<MIN_LIBRARY_SIZE>`: Minimum library size cutoff (e.g., 10000)
- `<MIN_PREVALENCE>`: Minimum samples a contig must appear in (e.g., 3)
- `<MIN_TOTAL_READS>`: Minimum total reads for a contig (e.g., 100)

---

## Step 1: Map Reads to Assembly

Map filtered reads back to assembled contigs using minimap2.
```bash
# Edit script
nano map_reads_to_assembly.sh

# Key parameters:
# CONTIGS_FILE: Assembled contigs (from step 03)
# READS_DIR: Filtered reads (from step 01)
# OUTPUT_DIR: Output directory for BAM files

# Submit job
qsub map_reads_to_assembly.sh
```

**Output**: Sorted and indexed BAM files (one per sample)

**Resource requirements**:
- Queue: hugemem
- Memory: 500 GB
- Walltime: 24 hours
- Uses minimap2 preset `-ax map-ont` for Nanopore reads

---

## Step 2: Build Count Matrix

Extract read mappings from BAM files and generate count matrix.
```bash
# Edit script
nano build_count_matrix.sh

# Set paths:
# BAM_DIR: Output from step 1
# OUTPUT_CSV: Mapped reads table
# MATRIX_OUTPUT: Count matrix (samples × contigs)

# Key parameter:
# MIN_MAPQ: 60 (high-quality mappings only)

# Submit job
qsub build_count_matrix.sh
```

**Outputs**:
- `mapped_reads_table.csv`: Individual read mappings
- `count_matrix.csv`: Aggregated counts (samples × contigs)

**Filtering**:
- MAPQ ≥ 60 (high-quality alignments)
- Excludes unmapped, secondary, and supplementary alignments

---

## Step 3: QC and Filter Count Matrix

Clean sample names, perform quality control, and filter low-quality samples and contigs.
```bash
# Edit script
nano qc_and_filter_counts.R

# Set paths:
# count_matrix_file: Output from step 2
# metadata_file: Sample metadata
# output_dir: Directory for QC plots and filtered outputs

# Key parameters:
# cutoff: Minimum library size (total counts per sample)
# min_samples: Minimum number of samples a contig must appear in
# min_total_reads: Minimum total reads for a contig across all samples

# Run
Rscript qc_and_filter_counts.R
```

**Outputs**:
1. **QC Plots** (in `<OUTPUT_DIR>/plots/`):
   - `library_size_distribution.png`: Distribution of total counts per sample
   - `contig_prevalence_distribution.png`: Distribution of contig prevalence
   - `contig_abundance_distribution.png`: Distribution of total counts per contig

2. **Filtered Data**:
   - `count_matrix_filtered.csv`: Filtered count matrix (samples × contigs)
   - `metadata_filtered.csv`: Metadata for retained samples

3. **Console Output**:
   - Summary statistics before and after filtering
   - Number of samples and contigs retained

**Filtering Strategy**:
1. **Sample filtering**: Remove samples with library size < cutoff
2. **Contig filtering**: Keep contigs that meet EITHER:
   - Present in ≥ min_samples samples, OR
   - Total reads ≥ min_total_reads across all samples

**Name Cleaning**:
- Standardizes sample names to format: `{ID}_{condition}_{rep1}_{rep2}`
- Ensures consistency between metadata and count matrix
- Matches samples between both files

---

## Step 4: Differential Abundance Analysis

Identify differentially abundant contigs using edgeR.
```bash
# Edit script
nano run_edger_analysis.R

# Set paths:
# count_matrix_file: Output from step 2
# metadata_file: Sample metadata
# output_dir: Results directory

# Run
Rscript run_edger_analysis.R
```

**Outputs**:
- `edger_results_categorical.csv`: Results for categorical comparison
- `edger_results_continuous.csv`: Results for continuous covariate
- `edger_*_fit.rds`: Saved edgeR objects for downstream analysis

**Analysis modes**:
1. **Categorical**: Compare groups (e.g., small vs. large biomass)
2. **Continuous**: Test correlation with continuous variable (e.g., shoot mass)

---

## Step 5: Empirical FDR Calibration

Calibrate FDR thresholds using permutation testing and generate comprehensive visualizations.
```bash
# Edit script
nano empirical_fdr_calibration.R

# IMPORTANT: Set paths to plotting functions
# These should point to the 05-visualization scripts
# PLOT_VOLCANO_R: ../05-visualization/plot_volcano.R
# PLOT_HEATMAP_R: ../05-visualization/plot_heatmap.R

# Set paths (same as step 3)

# Key parameter:
# B: 500 (number of permutations; increase to 1000 for final analysis)

# Run
Rscript empirical_fdr_calibration.R
```

**Outputs**:
1. **FDR Calibration**:
   - `empirical_fdr_curve.png`: FDR vs. p-value threshold
   - `observed_vs_null_discoveries.png`: Observed vs. null discoveries
   - `empirical_tau.csv`: Calibrated p-value threshold

2. **Volcano Plots**:
   - `volcano_empirical.png`: Using empirical τ threshold
   - `volcano_bh.png`: Using Benjamini-Hochberg FDR

3. **Heatmaps**:
   - `heatmap_empirical.png`: Top 50 contigs (empirical threshold)
   - `heatmap_bh.png`: Top 50 contigs (BH FDR)

4. **Results Tables**:
   - `DA_contigs_empirical_fdr.csv`: Significant contigs
   - `DA_contigs_list.txt`: List of significant contig IDs
   - `all_contigs_results.csv`: Complete results table

**Dependencies**:
- Requires `plot_volcano.R` and `plot_heatmap.R` from `05-visualization/`
- Make sure these are sourced correctly at the top of the script

**Why multiple visualization modes?**
- **Empirical (τ)**: Conservative, permutation-based, controls false discoveries empirically
- **BH FDR**: Standard, assumes theoretical null distribution, commonly used in publications
- **Comparison**: Shows difference between theoretical and empirical approaches

---

## Step 6: Map Contigs to SMAG + Plant Reference

Map assembled contigs against the combined SMAG + plant genome reference to assign each contig to a MAG. This generates the `mapped_contigs.csv` required by step 7.

> **Note**: This step uses the `asm5` preset (assembly-to-reference alignment), which is different from the `map-ont` preset used in `02-alignment` for read-level competitive mapping. Here we are aligning contigs (not reads) to identify which MAG each contig belongs to.

```bash
# Edit script
nano map_contigs_to_smag.sh

# Set paths:
# CONTIG_FILE: Assembled/polished contigs (from assembly step)
# REF_FILE: Combined SMAG + plant genome reference (FASTA/FASTA.gz or .mmi)
# OUTPUT_DIR: Output directory for mapping results

# Submit job
qsub map_contigs_to_smag.sh
```

**Outputs**:
- `<base_name>.sorted.bam`: Sorted BAM of contig-to-reference alignments
- `<base_name>.sorted.flagstat.txt`: Mapping statistics
- `<base_name>.sorted.idxstats.txt`: Per-reference alignment counts
- `mapped_contigs.csv`: Contig-to-MAG mapping table (input for step 7)

**Key parameters**:
- Preset: `asm5` (assembly-to-reference, ≤5% divergence)
- `--secondary=no`: Retain only the best hit per contig
- `--split-prefix`: Required for large combined references

**Comparison of alignment steps in this pipeline**:

This pipeline uses minimap2 in three distinct contexts. Each serves a different purpose and requires different presets and parameters:

| | Reads → SMAG + Plant (02-alignment) | Reads → Contigs (step 1) | Contigs → SMAG + Plant (this step) |
|---|---|---|---|
| **Input** | Raw/filtered reads | Filtered reads | Assembled contigs |
| **Reference** | Combined SMAG + plant genomes | Assembled contigs | Combined SMAG + plant genomes |
| **Preset** | `-ax map-ont` | `-ax map-ont` | `-ax asm5` |
| **Why this preset** | Long noisy reads to reference | Long noisy reads to reference | Assembly to reference (≤5% divergence) |
| **Secondary alignments** | Yes (competitive filtering) | No (count matrix needs unique assignments) | No (best hit only for annotation) |
| **MAPQ filter** | Not applied (all hits retained) | MAPQ ≥ 60 (step 2) | Not applied (best hit retained) |
| **Purpose** | Identify and remove host contamination from reads | Quantify per-sample contig abundance | Assign each contig to a MAG for taxonomy |
| **Output** | BAMs for contamination removal | Count matrix (samples × contigs) | `mapped_contigs.csv` (contig_id, MAG_ID) |
| **Downstream** | Clean reads → assembly | Differential abundance (steps 3–5) | Taxonomic annotation (step 7) |

---

## Step 7: Annotate Contigs

Assign taxonomy to contigs based on MAG mapping.
```bash
# Edit script
nano annotate_contigs.py

# Set paths:
# MAPPED_CONTIGS_CSV: Contig-to-MAG mapping
# SMAG_TAXONOMY_CSV: SMAG taxonomy database
# PLANT_MAG_FILES: Optional plant genome MAG lists

# Run
python annotate_contigs.py
```

**Outputs**:
- `annotated_contigs.csv`: Contigs with full taxonomy
- `taxonomy_summary_phylum.csv`: Counts by phylum
- `taxonomy_summary_family.csv`: Counts by family

**Taxonomy assignment**:
1. Match contigs to SMAGs → assign microbial taxonomy
2. Match contigs to plant genomes → assign host taxonomy
3. Clean and standardize taxonomy labels

---

## Troubleshooting

### Low Mapping Rate

**Symptoms**: <50% of reads map to assembly

**Possible causes**:
- Assembly incomplete
- Different samples used for assembly vs. mapping
- Contamination

**Solutions**:
1. Check assembly completeness (N50, total length)
2. Verify sample metadata
3. Run QC on both assembly and reads

---

### No Significant Contigs

**Symptoms**: All FDR > 0.05

**Possible causes**:
- Small treatment effect
- High biological variability
- Insufficient sample size
- Low coverage

**Solutions**:
1. Check power analysis (need n ≥ 5 per group)
2. Filter for high-coverage contigs (≥10× coverage)
3. Consider alternative analysis methods
4. Increase sample size if possible

---

### Permutations Too Slow

**Symptoms**: Empirical FDR calibration takes >24 hours

**Solutions**:
1. Reduce number of permutations (B = 100 for testing)
2. Filter contigs before analysis (keep high coverage only)
3. Use parallel processing (if available)

---

## Integration with Downstream Analysis

After abundance analysis:
```bash
# Abundance complete
qsub scripts/04-abundance-analysis/map_reads_to_assembly.sh
qsub scripts/04-abundance-analysis/build_count_matrix.sh
Rscript scripts/04-abundance-analysis/qc_and_filter_counts.R
Rscript scripts/04-abundance-analysis/run_edger_analysis.R
Rscript scripts/04-abundance-analysis/empirical_fdr_calibration.R
qsub scripts/04-abundance-analysis/map_contigs_to_smag.sh
python scripts/04-abundance-analysis/annotate_contigs.py

# Next: Binning and functional annotation (optional)
cd ../05-binning
# Or: Functional analysis of DA contigs
cd ../06-functional-annotation
```

---

## Parameter Tuning

### Library Size Cutoff

It depends on the overall situation of your dataset.

| Cutoff | Interpretation | Use Case |
|--------|---------------|----------|
| 1,000  | Very lenient  | Exploratory, high sample loss tolerance |
| 10,000 | Standard      | Recommended for most datasets |
| 50,000 | Conservative  | High-quality samples only |

### Contig Prevalence

| min_samples | Interpretation | Use Case |
|-------------|---------------|----------|
| 2           | Lenient       | Small sample sizes (n < 10) |
| 3           | Standard      | Recommended (n ≥ 10) |
| 5           | Conservative  | Large sample sizes (n ≥ 20) |

### Mapping Quality Threshold (MAPQ)

| MAPQ | Interpretation | Use Case |
|------|---------------|----------|
| 0    | Any alignment | Exploratory |
| 20   | ≥99% correct  | Standard |
| 60   | ≥99.9999% correct | Conservative (recommended) |

### edgeR Filtering
```R
# Default: filterByExpr (adaptive)
# Manual:
keep <- rowSums(cpm(y) > 1) >= 3  # CPM > 1 in ≥3 samples
```

### Log Fold Change Threshold
```R
# Conservative (recommended)
abs(logFC) >= 0.5  # 1.4-fold change

# Moderate
abs(logFC) >= 1.0  # 2-fold change

# Liberal
abs(logFC) >= 0    # Any change
```

---

## References

- edgeR: https://bioconductor.org/packages/edgeR/
- minimap2: https://github.com/lh3/minimap2
- Permutation FDR: https://doi.org/10.1093/bioinformatics/bti685