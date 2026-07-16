# Long-Read Soil Metagenomics Pipeline

[![DOI](https://zenodo.org/badge/1089933240.svg)](https://doi.org/10.5281/zenodo.18974276)

Workflow for high-quality genome assembly and differential abundance analysis 
from Oxford Nanopore soil metagenomic data.

## Overview

This repository contains the complete bioinformatic pipeline used in:

> Wu, H. (2025). Long-Read Soil Metagenomic Assembly and Differential 
> Abundance Analysis Across Plant Shoot Biomass Gradients. *Master's Thesis*, 
> Australian National University.

**📄 [Read the full thesis (PDF)](https://github.com/AdriaWu/soil-metagenome-workflow/blob/main/docs/HeyueWu_masters_thesis.pdf)**

## Data Availability

All raw sequencing data from this study are publicly available through NCBI:

**🔗 [NCBI BioProject PRJNA1400213](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1400213)**

- **358 BioSamples** across three experimental datasets
- **193 Gbp** of Oxford Nanopore R10 sequencing data
- Quality-filtered reads (≥300 bp, mean Q≥15)
- Comprehensive sample metadata included

### Datasets

The BioProject encompasses three complementary experiments:

1. **Boorowa agricultural soil time series** — Depth-stratified soil cores (surface and subsurface) sampled in 2019 and 2023, with paired soil chemistry data
2. **Plant productivity gradient** - Rhizosphere microbiomes across shoot biomass gradients in multiple crop hosts
3. **Forest vs. river soil comparison** - ACT soils under different fertilizer treatments

### Workflow 

### Pipeline Summary

```
Raw ONT Reads
    ↓
Adapter Trimming & Chimera Splitting (Porechop_ABI)
    ↓
Quality Control & Filtering (Chopper, NanoPlot)
    ↓
Host/Contamination Removal (minimap2 + competitive classification)
    ↓
Co-assembly (nanoMDBG)
    ↓
Assembly QC (MetaQUAST)
    ↓
Read Mapping & Count Matrix (minimap2, samtools)
    ↓
Compositional Transformation (GBM zero imputation, centred log-ratio)
    ↓
Differential Abundance (edgeR + empirical FDR)
    ↓
Community Analysis (PERMANOVA, betadisper, PCA)
    ↓
Taxonomic Annotation & Visualisation
```

### Key Features

- **High-contiguity assembly** using metaMDBG optimised for ONT R10 chemistry
- **Adapter-clean reads before assembly** — Porechop_ABI removes residual ONT adapters and splits chimeric reads at internal adapters, preventing the internal-adapter contigs that NCBI's Foreign Contamination Screen rejects
- **Rigorous contamination removal** through competitive mapping against host references
- **Compositionally-aware analysis** — GBM zero imputation and centred log-ratio transformation for sparse count data
- **Empirical FDR-controlled** differential abundance testing
- **Genome-resolved ecological interpretation** linking taxa to functional traits

## Key Results

- **556,967 contigs** totaling 7.69 Gbp
- **86 contigs >1 Mb** (longest: 4.52 Mb)
- **923 circular replicons** (complete genetic elements)
- **20% higher N50** vs. global SMAG catalogue
- **194 differentially abundant contigs** across biomass gradient

## Quick Start

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/)
- 32+ GB RAM recommended
- ~500 GB disk space for intermediate files

### Installation

```bash
# Clone repository
git clone https://github.com/AdriaWu/soil-metagenome-workflow.git
cd soil-metagenome-workflow

# Create conda environment
conda env create -f environment.yml
conda activate soil-meta
```

## Repository Structure

```
soil-metagenome-workflow/
├── scripts/
│   ├── 00-basecalling/          # Dorado basecalling and demultiplexing
│   ├── 01-quality-control/      # NanoPlot, FastQC, Porechop_ABI adapter trim, Chopper filtering
│   ├── 02-contamination-removal/# Competitive mapping and host read removal
│   ├── 03-assembly/             # metaMDBG assembly and QC
│   ├── 04-abundance-analysis/   # Read mapping, edgeR, FDR, annotation
│   └── 05-visualization/        # Volcano plots, heatmaps
├── docs/
│   ├── HeyueWu_masters_thesis.pdf
│   └── RESOURCE_GUIDE.md        # PBS resource allocation guidance
├── environment.yml
├── LICENSE
└── README.md
```

Each `scripts/` subdirectory contains its own README with step-by-step instructions, parameter tuning advice, and troubleshooting.

### PBS Resource Allocation

The default PBS resource requests in all scripts were calibrated for this study (~193 Gbp, 358 samples, combined reference of ~21,000 genomes). **Your dataset will have different characteristics — adjust resources accordingly.** See **[`docs/RESOURCE_GUIDE.md`](docs/RESOURCE_GUIDE.md)** for detailed guidance on right-sizing jobs, profiling resource usage, and reducing KSU costs.

## Citation

If you use this workflow, please cite both the thesis and the code repository:

```bibtex
@mastersthesis{wu2025longread,
  author  = {Wu, Heyue},
  title   = {Long-Read Soil Metagenomic Assembly and Differential
             Abundance Analysis Across Plant Shoot Biomass Gradients},
  school  = {Australian National University},
  year    = {2025},
  type    = {Master's Thesis}
}

@misc{wu2025soilmeta,
  author    = {Wu, Heyue},
  title     = {soil-metagenome-workflow},
  year      = {2025},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.18974276},
  url       = {https://doi.org/10.5281/zenodo.18974276},
  note      = {Software}
}
```

## License

[MIT License](LICENSE)

## Contact

- **Author**: Heyue Wu
- **Supervisor**: Professor Justin Borevitz
- **Institution**: Research School of Biology, Australian National University

## Acknowledgments

Special thanks to the Borevitz Lab, particularly Viraj Kolhapuri, Alek Meade, 
Saige Waugh, Ashley Jones, Christopher Bradley for their contributions.