#This code performs differential gene expression analysis, gene ontology enrichment analysis and gene set enrichment analysis
  # 1. user must define the design type for the variable dds_SAD
  # 2. user must define the threshold read count for genes


library("DelayedArray")
library(DESeq2)
library(tidyverse)
library(tools)
library(clusterProfiler)
library(org.Hs.eg.db)
library(conflicted)
library(tibble)
library(AnnotationDbi)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "read_counts"),
  make_option(c("-p", "--phenotype"), type = "character", help = "phenotype file"),
  make_option(c("-r", "--read_count"), type = "integer", help = "minimum read counts to remove outliers")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


counts_data_SAD <- read.csv(opt$input, header = T, sep = ',', row.names = NULL)

counts_data_SAD = na.omit(counts_data_SAD)
counts_data_SAD <- counts_data_SAD[!duplicated(rownames(counts_data_SAD)), ]
samp_SAD <- counts_data_SAD %>%
  distinct(gene_id, .keep_all = TRUE) %>%   # remove duplicates
  as.data.frame() %>%                     # convert tibble → data.frame
  column_to_rownames(var = "gene_id")
samp_SAD <- samp_SAD[ , !(names(samp_SAD) == "symbol")]

coldat_SAD <- read.csv(opt$phenotype, header = T, sep = ',', row.names = NULL)
coldat1_SAD <- coldat_SAD %>% remove_rownames %>% tibble::column_to_rownames(var = "sample_id")
coldat1_SAD
all(colnames(samp_SAD) %in% rownames(coldat1_SAD))
all(colnames(samp_SAD) == rownames(coldat1_SAD))
dds_SAD <- DESeqDataSetFromMatrix(countData = samp_SAD, colData = coldat1_SAD, design = ~ Type)
dds_SAD
keep <- rowSums(counts(dds_SAD)) >= opt$read_count
keep
dds_SAD <- dds_SAD[keep,]
nrow(dds_SAD)
assay(dds_SAD)
dds_SAD$Type <- relevel(dds_SAD$Type, ref = "Ctrl")
dds_SAD <- DESeq(dds_SAD)
res0.01_SAD <- results(dds_SAD, contrast = c("Type", "SAD", "Ctrl"), alpha = 0.05)
res0.01_SAD
summary(res0.01_SAD)

res0.01_SAD <- res0.01_SAD[!is.na(res0.01_SAD$baseMean) & !is.na(res0.01_SAD$log2FoldChange), ]

png("MAplot.png", width = 2300, height = 1200, res=300)
plotMA(res0.01_SAD)
dev.off()

res0.01_SAD <- as.data.frame(res0.01_SAD)
res0.01_SAD_dge <- res0.01_SAD[!is.na(res0.01_SAD$padj) & res0.01_SAD$padj <= 0.05, ]
head(res0.01_SAD, n=5)
write.csv(res0.01_SAD, "exp_profile.csv")
write.csv(res0.01_SAD_dge, "differentially_expressed_genes.csv")



# List of genes of interest
organism <- "org.Hs.eg.db" # Database for human genes (use appropriate one for your organism)

# Run enrichGO

gene_list_SAD <- rownames(res0.01_SAD)
ego_SAD <- enrichGO(
  gene          = gene_list_SAD,      # Input gene list
  OrgDb         = org.Hs.eg.db,       # Organism database (e.g., org.Hs.eg.db for human)
  keyType       = "ENSEMBL",       # Type of gene identifiers (e.g., "SYMBOL", "ENSEMBL")
  ont           = "All",           
  pAdjustMethod = "BH",          
  pvalueCutoff  = 0.05,          
  qvalueCutoff  = 0.2,           
  universe      = NULL            
)

# View results

ego_SAD_df <- as.data.frame(ego_SAD)
write.csv(ego_SAD_df, "GO.csv")
cnetplot(ego_SAD, foldChange = res0.01_SAD$log2FoldChange, showCategory = 5)

specific_term <- "GO:0034454"  # Replace with your GO ID

filtered_result <- ego_SAD@result[ego_SAD@result$ID == specific_term, "geneID"]
filtered_result <- unlist(strsplit(filtered_result, "/"))

filtered_fc_genes <- res0.01_SAD[rownames(res0.01_SAD) %in% filtered_result, ]


# Visualize results
png("dotplot.png", width = 4000, height = 3000, res=600)
dotplot(ego_SAD, showCategory = 10)
dev.off()

sorted_gGO_gene_rank_SAD <- res0.01_SAD[order(res0.01_SAD$log2FoldChange, decreasing = TRUE), ]
gGO_list_SAD <- sorted_gGO_gene_rank_SAD$log2FoldChange
names(gGO_list_SAD) <- rownames(sorted_gGO_gene_rank_SAD)



gGO_list_SAD <- na.omit(gGO_list_SAD)
gGO_list_SAD <- gGO_list_SAD[!duplicated(names(gGO_list_SAD))]
gGO_list_SAD <- sort(gGO_list_SAD, decreasing = TRUE)

set.seed(1234)
gse_result_list_SAD <- gseGO(
  geneList     = gGO_list_SAD,      # Ranked list of genes
  OrgDb        = org.Hs.eg.db,  # Organism annotation database
  keyType      = "ENSEMBL",             # Type of gene identifiers (e.g., SYMBOL, ENSEMBL)
  ont          = "All",          
  pAdjustMethod = "BH",           
  pvalueCutoff  = 0.05,           
  by = "fgsea"                  
)

# View Results
head(gse_result_list_SAD)
gse_result_list_SAD_df <- as.data.frame(gse_result_list_SAD)
write.csv(gse_result_list_SAD_df, "GO_module.csv")
# Visualization

gse_result_list_SAD_df <- gse_result_list_SAD_df[gse_result_list_SAD_df$pvalue < 0.05, ]

png("gseGO_dotplot.png", width = 4000, height = 3000, res=600)
dotplot(gse_result_list_SAD, showCategory = 10, color = "pvalue") # Show top 10 enriched terms# Show top 10 enriched terms
dev.off()
png("gseGO_ridgeplot.png", width = 2300, height = 2000, res=300)
ridgeplot(gse_result_list_SAD)
dev.off()
