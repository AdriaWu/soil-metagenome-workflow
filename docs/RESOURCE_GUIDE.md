# PBS Resource Allocation Guide

Guidance for adjusting PBS resource requests across the pipeline. The default values in all scripts were calibrated for a specific dataset (~193 Gbp ONT R10 data, 358 samples, combined SMAG + plant reference of ~21,000 genomes). **Your data will differ — adjust resources accordingly.**

## General Principles

### Right-size your jobs

Over-requesting wastes KSU (service units) and may increase queue wait times. Under-requesting causes job failures. The goal is to request ~20% more than actual peak usage.

### Profile before scaling

Run a small test (e.g., 1–2 samples) and check actual usage before submitting batch jobs:

```bash
# After job completes, check actual resource usage
qstat -f <JOB_ID>         # Walltime, memory used
nci_jobstats <JOB_ID>     # Detailed resource breakdown

# For running jobs
qstat -f <JOB_ID> | grep "resources_used"
```

### Common mistakes to avoid

- **Requesting hugemem for CPU-only tasks**: QC scripts (NanoPlot, FastQC, Chopper) rarely need >64 GB RAM. Use `normal` queue.
- **Requesting GPUs for non-GPU tools**: Only Dorado basecalling uses GPUs. Minimap2, samtools, R, and Python are CPU-only.
- **Over-requesting CPUs**: Many tools don't scale linearly. Minimap2 scales well to 48 threads; R scripts are typically single-threaded unless explicitly parallelised.
- **Matching walltime to worst case by default**: Start conservative, then reduce based on profiled runtimes.

---

## Resource Requirements by Step

### What drives resource needs

| Factor | Affects | Example |
|--------|---------|---------|
| **Total data volume** (Gbp) | RAM, walltime for mapping/assembly | 10 Gbp vs. 200 Gbp |
| **Number of samples** | Walltime for per-sample loops | 24 vs. 358 barcodes |
| **Reference size** | RAM for minimap2 indexing | Single genome vs. 21,000 MAGs |
| **Read length / N50** | Assembly RAM and time | 5 kb vs. 20 kb N50 |
| **Community complexity** | Assembly RAM and contig count | Low-diversity enrichment vs. soil |

### 00-basecalling

| Script | Queue | Key resource | Scales with |
|--------|-------|-------------|-------------|
| `run_dorado_basecaller.sh` | gpuvolta / dgxa100 | GPU time | Total POD5 data volume |
| `run_dorado_demux.sh` | normal (or GPU queue) | CPU time | Number of reads × barcodes |

- GPU type matters more than count: A100 is ~2–4× faster than V100 for SUP basecalling.
- Memory is rarely the bottleneck; 96–384 GB is typically sufficient.
- Walltime scales roughly linearly with data volume (~500 kb/s per V100, ~1–2 Mb/s per A100).

### 01-quality-control

| Script | Queue | RAM | Notes |
|--------|-------|-----|-------|
| `run_nanoplot.sh` | normal | 32–64 GB | Scales with number of reads per sample |
| `run_fastqc.sh` | normal | 32–64 GB | Java heap (`-Xmx`) must fit in allocated RAM |
| `filter_chopper.sh` | normal | 16–32 GB | I/O-bound, not RAM-bound |
| `multiqc_report.sh` | normal | 8–16 GB | Aggregation only; very light |

These scripts **do not require hugemem**. The `normal` queue is sufficient for virtually all QC tasks.

### 02-contamination-removal

| Script | Queue | RAM | Scales with |
|--------|-------|-----|-------------|
| `competitive_mapping.sh` | hugemem | 500–1000 GB | **Reference size** (dominant factor) + data volume |
| `filter_contaminants.sh` | hugemem | 500–1000 GB | Number of samples × BAM size |

- The combined SMAG + plant reference (~21,000 genomes) is the primary reason for hugemem. A smaller reference (e.g., single plant genome + SMAG subset) would fit in the `normal` queue.
- If your reference is <50 GB uncompressed, try `normal` queue with 190 GB RAM first.

### 03-assembly

| Script | Queue | RAM | Scales with |
|--------|-------|-----|-------------|
| `run_nanomdbg.sh` | hugemem | 1000–1470 GB | Data volume × community complexity |
| `polish_nanomdbg.sh` (v1.2) | hugemem | 1000–1470 GB | Assembly size × read volume |
| `assembly_qc.sh` | normal | 64–128 GB | Assembly size |

- metaMDBG is memory-intensive for large, complex metagenomes. For smaller datasets (<20 Gbp), `normal` queue with 190 GB may suffice.
- MetaQUAST does not require hugemem.

### 04-abundance-analysis

| Script | Queue | RAM | Scales with |
|--------|-------|-----|-------------|
| `map_reads_to_assembly.sh` | normal–hugemem | 64–500 GB | Assembly size (reference) + data volume |
| `build_count_matrix.sh` | normal | 32–64 GB | Number of BAM files × contigs |
| `qc_and_filter_counts.R` | normal | 16–32 GB | Count matrix dimensions |
| `run_edger_analysis.R` | normal | 16–32 GB | Count matrix dimensions |
| `empirical_fdr_calibration.R` | normal | 32–64 GB | Permutations × matrix size |
| `map_contigs_to_smag.sh` | hugemem | 500–1470 GB | **Reference size** (same as 02) |
| `annotate_contigs.py` | normal | 8–16 GB | Number of contigs |

- Step 1 (read mapping) may or may not need hugemem depending on assembly size. A 7.69 Gbp assembly as reference needs ~50–100 GB for minimap2 indexing; the combined SMAG reference needs much more.
- R scripts are typically single-threaded — requesting 48 CPUs for them wastes resources. Use 1–4 CPUs.
- The FDR permutation step benefits from increased walltime more than increased RAM.

---

## Quick Reference: Matching Queue to Task

| Queue | Max RAM | Use when |
|-------|---------|----------|
| **normal** | 190 GB | QC, R/Python analysis, small references, MetaQUAST |
| **hugemem** | 1470 GB | Large reference mapping, metagenome assembly, large BAM operations |
| **gpuvolta** | 384 GB | Dorado basecalling (V100) |
| **dgxa100** | 384 GB | Dorado basecalling (A100, faster) |

---

## Estimating RAM for minimap2

Minimap2 RAM usage depends primarily on reference size:

```bash
# Check reference size
du -sh reference.fasta
# or for gzipped
zcat reference.fasta.gz | wc -c  # uncompressed bytes
```

**Rule of thumb**: minimap2 index requires ~5–8× the uncompressed reference size in RAM, plus overhead for alignment buffers.

| Reference size (uncompressed) | Approximate RAM needed | Suggested queue |
|-------------------------------|----------------------|-----------------|
| < 1 GB | 16–32 GB | normal |
| 1–10 GB | 64–128 GB | normal |
| 10–50 GB | 190 GB | normal (max) |
| > 50 GB | 500–1470 GB | hugemem |

---

## Walltime Estimation

| Task type | Rough rate | Example |
|-----------|-----------|---------|
| Basecalling (V100 SUP) | ~500 kb/s | 100 Gbp ≈ 56 hours |
| Basecalling (A100 SUP) | ~1.5 Mb/s | 100 Gbp ≈ 19 hours |
| Minimap2 read mapping | ~5–10 min/Gbp | 200 Gbp ≈ 17–33 hours |
| metaMDBG assembly | Highly variable | 200 Gbp complex soil ≈ 24–48 hours |
| edgeR + permutation FDR | ~1–4 hours | 500 permutations, 500k contigs |

These are rough estimates. Always profile with a subset first.

---

## Reducing KSU Usage

1. **Profile first**: Run 1–2 samples, check `nci_jobstats`, then scale.
2. **Use normal queue when possible**: hugemem nodes cost more KSU per hour.
3. **Match CPUs to tool parallelism**: Don't request 48 CPUs for single-threaded R scripts.
4. **Set realistic walltime**: Unused walltime isn't charged, but shorter requests may schedule faster.
5. **Use jobfs for I/O-heavy steps**: Copying data to local SSD (jobfs) before processing avoids network bottlenecks and reduces walltime.
