# CODEG
This repository contains source codes (WGCNA.R and DGE.R) for performing the Weighted Gene Co-Expression, Differential Gene Expression and Functional Enrichment Analysis, specifically designed to run on the input files. 
WGCNA analysis generates co-expressed gene pairs from gene modules

The code WGCNA.R performs Weigted Gene Co-expression Network Analysis using the inputs total_gene_counts1.csv and trait.csv

The file total_gene_counts1.csv contains raw read gene counts of the Control, FAD and SAD patient samples.

The file trait.csv contains the sample metadata including the traits.

Following are the inputs of WGCNA.R

  make_option(c("-i", "--input"), type="character", help="read count"),
  make_option(c("-p", "--trait"), type="character", help="trait file"),
  make_option(c("-r", "--min_read_count"), type="integer", default=30, help="minimum read counts to remove outliers"),
  make_option(c("-c", "--min_sample_count"), type="integer", default=2, help="minimum sample counts to remove outliers"),
  make_option(c("-b", "--maxblocksize"), type="integer", default=2000, help="set blocksize"),
  make_option(c("-s", "--soft_threshold"), type="integer", default=8, help="soft_threshold_power"),
  make_option(c("-t", "--threads"), type="integer", default=1, help="Number of threads")

1. User must define the design type of experiment based on the trait.csv for the variable dds.
2. User must define the read count threshold and sample counts based on the total_gene_counts1.csv for variable dds75.
3. User must define the soft threshold power for their adjacency matrix based on the soft threshold indices.
4. User must define the maxblocksize in the blockwisemodules function based on the number of genes retained after outlier removal and also based on the computational complexity.

Pre-requisites

Install the following packages:

WGCNA
DESeq2
readr
tidyverse
heatmaply
gridExtra
dplyr
sva
pheatmap
readr
org.Hs.eg.db
cluster
limma

Differential gene expression analysis, gene ontology enrichment analysis and gene set enrichment analysis

Following are the inputs for DGE.R

  make_option(c("-i", "--input"), type = "character", help = "read_counts"),
  make_option(c("-p", "--phenotype"), type = "character", help = "phenotype file"),
  make_option(c("-r", "--read_count"), type = "integer", help = "minimum read counts to remove outliers")

1. User must define the design type of experiment based on the trait_SAD.csv for the variable dds_SAD.
2. User must define the read count threshold and sample counts based on the total_gene_counts_SAD.csv for variable dds_SAD.


Pre-requisites

Install the following packages:


DelayedArray
DESeq2
tidyverse
tools
clusterProfiler
org.Hs.eg.db
conflicted
tibble
AnnotationDbi
optparse
