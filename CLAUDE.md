# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**MethylMiner** is a bioinformatics pipeline for discovering rare methylation events in DNA methylation arrays. It detects Differentially Methylated Regions (DMRs) that cause disease, particularly focusing on cases like Fragile X syndrome where FMR1 hypermethylation leads to developmental disorders.

The codebase consists of two main components:
1. **methFlow** - NextFlow DSL2 pipeline for data processing
2. **MethVis** - Python Dash web application for interactive visualization

## Architecture

### Two-Component Design

**methFlow Pipeline** (`methFlow/` directory):
- NextFlow DSL2 orchestration (267 lines in `map.nf`)
- Sequential processing stages: QC → Normalization → Sorting → DMR Detection → Annotation
- Supports 3 workflows: QC_WF (QC only), DMR_WF (DMR only), QUANTILE_WF (complete)
- LSF cluster integration configured in `nextflow.config`

**MethVis Dashboard** (`MethVis/` directory):
- Jupyter Dash web application (`app.py` - 943 lines, largest single file)
- Interactive exploration of pipeline results
- PCA analysis, distribution plots, DMR browser

### Processing Pipeline (6 Stages)

1. **QC & Normalization** (`1_metharray_QC_norm.R`) - R/minfi processing
2. **Beta Value Sorting** (Bash) - Sort by genomic position for windowing
3. **BigWig Generation** (`getBigWig.py`) - Convert to genome browser format
4. **Rare Event Detection** (`2_findEpivariation.pl`) - Core algorithm (255 lines Perl)
5. **DMR Consolidation** (`3_getDMRlist.R`) - Merge adjacent significant windows
6. **Genomic Annotation** (`5_annotateDMRs.R`) - Annotate with genes/features

### Data Flow

```
IDAT files + metadata CSV → QC/Normalization → Sorted beta values →
BigWig files → Sliding window analysis → Significant regions →
Consolidated DMRs → Annotated Excel output
```

## Common Development Commands

### Running the Pipeline

**Load required modules (HPCF cluster):**
```bash
module load nextflow/21.10.5
module load R/4.1.0
module load perl/5.10.1
module load python/3.7.0
```

**Execute complete analysis:**
```bash
cd methFlow
nextflow run map.nf -entry QUANTILE_WF \
  --workdir /path/to/idat/files \
  --runName myanalysis \
  --outdir /path/to/results \
  --windowSize 1000 \
  --qCutMin 0.25 \
  --qCutMax 0.75 \
  --email your@email.com
```

**Run QC validation only:**
```bash
nextflow run map.nf -entry QC_WF \
  --workdir /path/to/idat/files \
  --runName myanalysis \
  --outdir /path/to/results
```

**Run DMR detection on existing QC output:**
```bash
nextflow run map.nf -entry DMR_WF \
  --runName myanalysis \
  --outdir /path/to/results \
  --windowSize 2000 \
  --qCutMin 0.1 \
  --qCutMax 0.9
```

### Testing with Example Data

**Use provided Fragile X test samples:**
```bash
cd methFlow
nextflow run map.nf -entry QUANTILE_WF \
  --workdir ../data/FragileX_exampleData \
  --runName fragileX_test \
  --outdir ./test_output
```

### Visualization Dashboard

**Setup Conda environment:**
```bash
conda create -n dash
conda activate dash
conda install jupyterLab numpy scipy matplotlib seaborn pandas \
  matplotlib-venn ipykernel plotly dash dash-bio \
  dash-core-components dash-bootstrap-components \
  notebook jupyter-dash -c conda-forge -c anaconda -c plotly
```

**Launch visualization:**
```bash
cp ./MethVis/app.py /path/to/runName_preprocessIllumina/
cd /path/to/runName_preprocessIllumina/
conda activate dash
python app.py
```

## Key Technologies

### Programming Languages
- **NextFlow DSL2** (267 lines) - Pipeline orchestration
- **R** (~500 lines) - Bioinformatics processing using Bioconductor/minfi
- **Python** (~1033 lines) - Web dashboard (Jupyter Dash)
- **Perl** (255 lines) - Core sliding-window algorithm for rare event detection
- **Bash** - File sorting and system operations

### Critical Dependencies

**R Bioconductor packages:**
- `minfi` - Methylation array processing
- `IlluminaHumanMethylationEPICmanifest` - EPIC probe definitions
- `IlluminaHumanMethylationEPICanno.ilm10b4.hg19` - EPIC annotations
- `GenomicRanges`, `BiocGenerics` - Genomic data structures

**Python packages:**
- `dash`, `plotly` - Web framework and interactive plots
- `numpy`, `pandas`, `scipy` - Scientific computing
- `sklearn` - PCA analysis

**System requirements:**
- LSF job scheduler (configured in `nextflow.config`)
- Internet access for manifest download
- 8GB memory per process (configurable)

## Input Requirements

### Required Files
- **IDAT files**: Paired Red/Green files per sample (`*_Red.idat`, `*_Grn.idat`)
- **Metadata CSV** with exact column names:
  - `Sample_Name` - Sample identifier
  - `Sentrix_ID` - 12-digit array ID
  - `Sentrix_Position` - Array position (R01C01 format)
  - `Sample_ID` - Concatenated Sentrix_ID_Position
  - `Reported_Sex` - Male/Female/Unknown
  - `Sample_Group` - Case/Control

### Current Platform Support
- **Array platform**: EPIC only (not 450K, EPICv2)
- **Genome build**: hg19 only (not hg38)
- **Sample types**: Any tissue type compatible with EPIC arrays

## Output Structure

Pipeline creates organized directory structure:
```
<outdir>/<runName>/<runName>_preprocessIllumina/
├── normalized_data/
│   ├── betaValues/           # Sorted beta value matrices
│   ├── bigWig/              # Genome browser tracks (per sample)
│   ├── dmr_<params>/        # DMR results and annotations
│   └── qcReports/           # QC visualizations (PDFs)
├── raw_data/                # Initial QC outputs
├── sex_prediction/          # Sex validation results
└── app.py                   # MethVis dashboard (copied here)
```

## Key Algorithm Details

### Core DMR Detection (`2_findEpivariation.pl`)
- **Sliding window approach**: Default 1kb windows across genome
- **Quantile-based outlier detection**: Identifies samples with extreme methylation
- **Consecutive probe requirement**: 3+ adjacent outlier probes needed
- **Directional consistency**: All outliers in same direction (hyper/hypo)
- **Statistical scoring**: Calculates significance for each window

### Parameters for DMR Sensitivity
- `--windowSize`: Window size in bp (default: 1000)
- `--qCutMin`: Lower quantile threshold (default: 0.25)
- `--qCutMax`: Upper quantile threshold (default: 0.75)

Lower quantile cutoffs = more sensitive detection but more false positives

## Test Data

**Fragile X validation samples** (`data/FragileX_exampleData/`):
- 3 samples with known FMR1 hypermethylation
- NA09145, NA09237 (males), NA07063 (female)
- Positive controls for pipeline validation
- Should detect hypermethylated DMR at Xq27.3 (FMR1 locus)

## Known Limitations

- EPIC platform only (no 450K, EPICv2 support)
- hg19 genome only
- No cell type deconvolution
- No CNV calling
- No Docker/containerization
- Cluster-dependent (LSF configuration)