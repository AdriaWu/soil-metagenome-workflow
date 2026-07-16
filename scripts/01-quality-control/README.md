# Quality Control Scripts

PBS scripts for Oxford Nanopore read quality assessment and filtering.

## Prerequisites

- Conda environment with: `nanoplot`, `fastqc`, `multiqc`, `chopper`, `porechop_abi`
- PBS/Torque cluster with jobfs support

## Configuration

Before running, edit each script to replace placeholders:

- `<PROJECT>`: Your NCI project code (e.g., `xe2`)
- `<CONDA_PATH>`: Path to your conda installation
- `<QC_ENV_PATH>`: Path to your QC conda environment
- `<INPUT_FASTQ_DIR>`: Directory containing input FASTQ files
- `<OUTPUT_BASE_DIR>`: Base directory for results
- `<OUTPUT_ADAPTERTRIM_DIR>`: Directory for adapter-trimmed reads (porechop_abi only)
- `<OUTPUT_FILTERED_DIR>`: Directory for filtered reads (chopper only)

## Preserving BAM Tags (MM/ML) During FASTQ Conversion

If your input data is in unaligned BAM format with modified basecalling tags (e.g. MM/ML from Dorado's `--modified-bases 5mC_5hmC`), converting to FASTQ with a standard `samtools fastq` command will **drop these tags**. This matters if you need to preserve methylation calls for downstream analysis.

### Converting BAM to FASTQ While Preserving Tags

Use `samtools fastq -T` to carry specified tags into the FASTQ header comments:

```bash
# Preserve MM and ML tags
samtools fastq -T MM,ML unaligned.bam > reads_with_tags.fastq

# Preserve all tags
samtools fastq -T '*' unaligned.bam > reads_with_all_tags.fastq
```

The tags are written into the FASTQ header comment field (after the first whitespace in the `@` line). For example:
```
@read_id  MM:Z:C+h?,19,4;C+m?,19,4;  ML:B:C,234,51,20,1
ACGTACGT...
+
IIIIIIII...
```

### Aligning Tag-Preserved FASTQ with Minimap2

When aligning FASTQ files that contain tags in the header comments, use minimap2's `-y` flag to copy those comments into the SAM output:

```bash
samtools fastq -T MM,ML unaligned.bam \
  | minimap2 -y -ax map-ont reference.fa - \
  | samtools sort -o aligned.bam
```

- `-T MM,ML`: writes the MM and ML tags into the FASTQ header comment field
- `-y`: copies FASTQ header comments into the SAM auxiliary tags

**Note**: The `-y` flag must be a separate argument, not grouped with `-ax` (e.g. `-ax map-ont -y`, not `-ax -y map-ont`).

### Alternative: `dorado aligner`

Dorado's built-in aligner accepts unaligned BAM directly and preserves MM/ML tags natively, bypassing FASTQ conversion:

```bash
dorado aligner reference.fa unaligned.bam > aligned.bam
```

This uses the `lr:hq` minimap2 preset by default (suited for sup-basecalled reads). Override with `--mm2-opts "-x map-ont"` if needed.

### Verifying Tags Are Preserved

After alignment, confirm MM/ML tags are present:
```bash
samtools view aligned.bam | head -1 | tr '\t' '\n' | grep "^MM\|^ML"
```

### Relevance to This QC Workflow

The QC tools in this pipeline (NanoPlot, FastQC, Chopper, Porechop_ABI) operate on FASTQ files and do not use or require MM/ML tags. However, if you are converting BAMs to FASTQ as input for this QC workflow **and** you need to use those same FASTQ files for downstream alignment with tag preservation, make sure to use `samtools fastq -T MM,ML` during the conversion step.

**Extra caveat for Porechop_ABI:** unlike the other QC tools, Porechop_ABI *rewrites* reads — it trims adapter bases and renames/splits reads at internal adapters — so it does **not** preserve MM/ML tag comments, and the read coordinates no longer match the original BAM. The adapter-trimmed FASTQ is the correct input for the **assembly** branch (metaMDBG needs no methylation tags). Do **not** feed it into a methylation/tag-preserving alignment branch; that branch should start from the original unaligned BAM (or a `samtools fastq -T MM,ML` conversion), not from Porechop_ABI output.

If QC and alignment are separate workflows using different FASTQ files, this is not a concern — just ensure the alignment workflow starts from the original unaligned BAMs or from FASTQ files generated with `-T MM,ML`.

---

## Workflow

### Step 1: NanoPlot (Read Statistics)
```bash
# Edit script to set paths
nano run_nanoplot.sh

# Submit job
qsub run_nanoplot.sh
```

**Output**: `<OUTPUT_BASE_DIR>/nanoplot/<sample>.txt` and `<sample>.html`

---

### Step 2: FastQC (Per-base Quality)
```bash
# Edit script to set paths
nano run_fastqc.sh

# Submit job
qsub run_fastqc.sh
```

**Output**: `<OUTPUT_BASE_DIR>/fastqc/<sample>_fastqc.html` and `.zip`

---

### Step 3: Porechop_ABI (Adapter Trimming & Chimera Splitting)

Runs **before** Chopper. Removes Oxford Nanopore adapters and, by default,
**splits reads at internal (mid-read) adapters**. An internal adapter marks a
chimeric/concatenated molecule; if left in place it is assembled straight
through and surfaces as internal-adapter contamination at NCBI submission.
Dorado `--trim` (at basecalling) only trims adapters from read *ends*, so this
step is what catches the rest.

```bash
# Edit script to set paths
nano 0104_trim_porechop_abi.sh

# Submit job
qsub 0104_trim_porechop_abi.sh
```

- Uses `-abi` (ab-initio adapter discovery, on top of the built-in adapter set);
  it also reports when a file is already trimmed.
- Splitting on internal adapters is the **default** — do not pass `--no_split`
  (keeps chimeras) or `--discard_middle` (discards the whole read).
- Porechop_ABI is the slow step. For a full cohort, run it as a PBS **array**
  (one sample per task) — see the header comments in the script.

**Output**: `<OUTPUT_ADAPTERTRIM_DIR>/<sample>_adaptertrim.fastq.gz`

---

### Step 4: Chopper (Quality Filtering)
```bash
# Edit script to set paths and filtering parameters
nano filter_chopper.sh

# Adjust parameters:
MIN_LENGTH=2000    # Minimum read length (bp)
MIN_QUALITY=18     # Minimum average quality score

# Submit job (run on the adapter-trimmed reads from Step 3)
qsub filter_chopper.sh
```

**Output**: `<OUTPUT_FILTERED_DIR>/<sample>_filtered.fastq.gz`

---

### Step 5: MultiQC (Aggregate Report)

After running NanoPlot and/or FastQC:
```bash
# Edit script to set paths
nano multiqc_report.sh

# Submit job
qsub multiqc_report.sh
```

**Output**: `<OUTPUT_BASE_DIR>/multiqc/multiqc_report.html`

## Example: Complete QC Workflow
```bash
# 1. Run NanoPlot on raw reads
qsub run_nanoplot.sh

# 2. Run FastQC on raw reads
qsub run_fastqc.sh

# 3. Wait for jobs to complete, then aggregate
qsub multiqc_report.sh

# 4. Adapter-trim & split chimeras on the raw demuxed reads
qsub 0104_trim_porechop_abi.sh

# 5. Quality/length filter the adapter-trimmed reads
qsub filter_chopper.sh

# 6. Run QC on the final filtered reads
# (Edit scripts to point to filtered directory, then rerun)
qsub run_nanoplot.sh
qsub run_fastqc.sh
qsub multiqc_report.sh
```

## Resource Requirements

| Script             | Memory | CPUs | Walltime | Jobfs |
|--------------------|--------|------|----------|-------|
| run_nanoplot       | 64 GB  | 8    | 5 hours  | 400 GB|
| run_fastqc         | 64 GB  | 8    | 5 hours  | 400 GB|
| trim_porechop_abi  | 48 GB  | 12   | 10 hours | 200 GB|
| filter_chopper     | 32 GB  | 8    | 10 hours | 400 GB|
| multiqc_report     | 16 GB  | 1    | 1 hour   | N/A   |

> Porechop_ABI uses little memory (<1 GB peak in testing) but is CPU-bound and
> slow. The row above is per-sample; for a full cohort run it as a PBS array so
> samples process in parallel rather than one 10-hour job per sample serially.

## Troubleshooting

### FastQC Java Heap Space Error

If FastQC crashes with `OutOfMemoryError`:
- Increase `_JAVA_OPTIONS` in `run_fastqc.sh`
- Currently set to `-Xmx48G` (48 GB)

### NanoPlot Empty File Error

If NanoPlot fails with "no reads found":
- Check input file format (gzipped FASTQ expected)
- Verify reads start with `@` character
- Try uncompressed FASTQ

### MultiQC No Data

If MultiQC report is empty:
- Verify FastQC/NanoPlot directories exist
- Check file permissions
- Ensure FastQC/NanoPlot jobs completed successfully

## Output Structure
```
<OUTPUT_BASE_DIR>/
├── nanoplot/
│   ├── sample1.txt
│   ├── sample1.html
│   ├── sample2.txt
│   └── sample2.html
├── fastqc/
│   ├── sample1_fastqc.html
│   ├── sample1_fastqc.zip
│   ├── sample2_fastqc.html
│   └── sample2_fastqc.zip
└── multiqc/
    ├── multiqc_report.html
    └── multiqc_data/

<OUTPUT_ADAPTERTRIM_DIR>/
├── sample1_adaptertrim.fastq.gz
└── sample2_adaptertrim.fastq.gz

<OUTPUT_FILTERED_DIR>/
├── sample1_filtered.fastq.gz
└── sample2_filtered.fastq.gz
```

## Notes

- Scripts use jobfs for temporary storage (faster I/O)
- Input files are copied to jobfs before processing
- Results are copied back to final output directory
- Temporary files are cleaned up automatically
- All scripts include error handling (`set -euo pipefail`)