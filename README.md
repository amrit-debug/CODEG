# WGCNA and Differential Gene Expression (DGE) Analysis Pipeline

## Overview

This pipeline performs:

1. **Differential Gene Expression (DGE) Analysis**
2. **Weighted Gene Co-expression Network Analysis (WGCNA)**

The workflow accepts a gene expression matrix and sample metadata (trait information) and generates differential expression results, co-expression modules, module-trait relationships, hub genes, and network files suitable for downstream biological interpretation.

---

# Table of Contents

- Overview
- Requirements
- Installation
- Input Files
- WGCNA Command Usage
- DGE Command Usage
- References

---

# Requirements

## Operating System

The pipeline has been tested on:

- Linux (recommended)
- macOS
- Windows

---

## Software Requirements

### R

Recommended:

```text
R >= 4.2
```

Check version:

```bash
R --version
```

---

# Installation

## Install Required CRAN Packages

Open R and run:

```r
install.packages(c(
  "WGCNA",
  "readr",
  "tidyverse",
  "dplyr",
  "pheatmap",
  "dynamicTreeCut",
  "cluster",
  "matrixStats",
  "gridExtra",
  "heatmaply",
  "sva",
  "GEOquery",
  "CorLevelPlot"
))
```

---

## Install Bioconductor

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```

---

## Install DESeq2

```r
BiocManager::install("DESeq2")
```

---

## Bioconductor Packages

```r
BiocManager::install(c(
  "limma",
  "edgeR",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi"
))
```

---

# Input Files

Two files are required:

## 1. Gene Expression Matrix

### Format

- CSV file
- Rows = genes
- Columns = samples
- First column = gene identifiers



### Requirements

- Gene IDs must be unique.
- Expression values must be numeric.
- Sample names must be unique.
- Missing values should be removed whenever possible.

---

## 2. Sample Metadata / Trait File

### Format

- CSV file
- One row per sample
- Sample IDs must exactly match expression matrix column names



### Requirements

- Sample IDs must match expression matrix columns.
- No duplicated sample names.
- Trait variables may be categorical or numeric.

---


# WGCNA Command Usage
## For AD_CA3
Rscript WGCNA.R -i total_gene_counts_WGCNA.csv -p phenotype_WGCNA.csv -r 36 -c 8 -s 8 -b 17000 -t 12

## For AD_hpNPCs
Rscript WGCNA.R -i total_gene_counts_WGCNA.csv -p phenotype_WGCNA.csv -r 30 -c 10 -s 14 -b 17000 -t 12

## For PD_neural_progenitors
Rscript WGCNA.R -i total_gene_counts_WGCNA.csv -p phenotype_WGCNA.csv -r 24 -c 12 -s 12 -b 15000 -t 12

## For PD_terminally_diff_neurons
Rscript WGCNA.R -i total_gene_counts_WGCNA.csv -p phenotype_WGCNA.csv -r 20 -c 12 -s 8 -b 24000 -t 12

## For ALS_frontal
Rscript WGCNA.R -i total_gene_counts_WGCNA.csv -p phenotype_WGCNA.csv -r 30 -c 50 -s 18 -b 19000 -t 12

## For ALS_motor
Rscript WGCNA.R -i total_gene_counts_WGCNA.csv -p phenotype_WGCNA.csv -r 30 -c 50 -s 18 -b 19000 -t 12


---

# DGE Command Usage
## For AD_CA3 and AD_hpNPCs
Rscript DGE.R -i total_gene_counts_FAD_DGE.csv -p phenotype_FAD_DGE.csv -r 30
Rscript DGE.R -i total_gene_counts_SAD_DGE.csv -p phenotype_SAD_DGE.csv -r 30

## For PD and ALS (all tissue types)
Rscript DGE.R -i total_gene_counts_WGCNA.csv -p phenotype_WGCNA.csv -r 30


---

# References

## WGCNA

Langfelder P, Horvath S. (2008)

WGCNA: an R package for weighted correlation network analysis.

BMC Bioinformatics 9:559.

---

## DESeq2

Love MI, Huber W, Anders S. (2014)

Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.

Genome Biology 15(12):550.

---

## clusterProfiler

Wu T et al. (2021)

clusterProfiler 4.0: A universal enrichment tool.

The Innovation 2(3):100141.

---

