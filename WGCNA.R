#This code performs WGCNA analysis and generates co-expressed gene pairs from gene modules
#Following conditions are necessary to run this program as per user:
# 1. user must define the design type of experiment for variable dds
# 2. user must define the read count threshold and sample counts for variable dds75
# 3. user must define the soft threshold power as per their topology
# 4. user must define the maxblocksize as stated for the variable bwnet

library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="read count"),
  make_option(c("-p", "--trait"), type="character", help="trait file"),
  make_option(c("-r", "--min_read_count"), type="integer", default=30, help="minimum read counts to remove outliers"),
  make_option(c("-c", "--min_sample_count"), type="integer", default=2, help="minimum sample counts to remove outliers"),
  make_option(c("-b", "--maxblocksize"), type="integer", default=2000, help="set blocksize"),
  make_option(c("-s", "--soft_threshold"), type="integer", default=8, help="soft_threshold_power"),
  make_option(c("-t", "--threads"), type="integer", default=1, help="Number of threads")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

library(WGCNA)
library(DESeq2)
library(readr)
library(tidyverse)
library(heatmaply)
library(gridExtra)
library(dplyr)
library(sva)
library(pheatmap)
library(readr)
library("org.Hs.eg.db")
library(cluster)
library(limma)
options(mc.cores = opt$threads)

# 1. Fetch Data ------------------------------------------------
data <- read.csv(opt$input, header = T, sep = ',', row.names = NULL)
head(data, n=5)
phenoData <- read.csv(opt$trait, header = T, sep = ',', row.names = NULL)
phenoData 
class(phenoData)
class(phenoData)
dat <- as.data.frame(data)
head (dat, n=5)
dat1 <- dat %>%
  column_to_rownames(var = 'gene_id')
head(dat1, n=5)


# detect outlier genes
dat1 <- na.omit(dat1)
dat1[] <- lapply(dat1, function(x) as.numeric(as.character(x)))
gsg <- goodSamplesGenes(t(dat1))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
dat1 <- dat1[gsg$goodGenes == TRUE, ]
dat1 <- as.data.frame(dat1)

# convert NaN to NA
dat1[is.nan(as.matrix(dat1))] <- NA

# remove columns with NA
dat1 <- dat1[, colSums(is.na(dat1)) == 0]

dat1
htree <- hclust(dist(t(dat1)), method = "average")
png("tree.png", units = "in", width = 16, height = 10, res = 600)
pltree <- plot(htree)
dev.off()

# exclude outlier samples
samples.to.be.excluded <- ""
data1.subset <- dat1[,!(colnames(dat1) %in% samples.to.be.excluded)]
data1.subset
ncol(data1.subset)
# 3. Normalization -------------------------------------------------------------------
# create a deseq2 dataset

# exclude outlier samples
phenoData
rownames(phenoData) <- phenoData[[1]]
phenoData <- phenoData[, -1]
colData <- phenoData %>% 
  filter(!row.names(.) %in% samples.to.be.excluded)
colData
nrow(colData)

# making the rownames and column names identical
all(rownames(colData) %in% colnames(data1.subset))
all(rownames(colData) == colnames(data1.subset))


# create dds
head(data1.subset)
dds <- DESeqDataSetFromMatrix(countData = data1.subset,
                              colData = colData,
                              design = ~ Type) # not spcifying model




dds75 <- dds[rowSums(counts(dds) >= opt$min_read_count) >= opt$min_sample_count,] #36/8
dds75
nrow(dds75)


# perform variance stabilization
dds_norm <- vst(dds75)
vst_mat <- assay(dds_norm)



norm.counts_corrected <- vst_mat %>% 
  t()

write.csv(norm.counts_corrected, "normalized_corrected_counts.csv", row.names = TRUE)

power <- c(c(1:16), seq(from = 18, to = 50, by = 2))
power

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts_corrected,
                         powerVector = power,
                         networkType = "signed",
                         corFnc = "bicor",
                         verbose = 5)

head(sft, n = 1)
sft.data1 <- sft$fitIndices
head(sft.data1, n = 28)
# visualization to pick power


a1 <- ggplot(sft.data1, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()



a2 <- ggplot(sft.data1, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

jpeg("softpowerplotmean_ac.jpeg")
grid.arrange(a1, a2, nrow = 2)
dev.off()

k=softConnectivity(norm.counts_corrected, 
                   corFnc = "bicor", 
                   type = "signed",
                   blockSize = opt$maxblocksize,
                   power=8,
                   verbose = 2)
#k = abs(k)
sizeGrWindow(10,5)
par(mfrow=c(1,2))
jpeg("histo_ac.jpeg")
hist(k)
dev.off()

jpeg("sft_ac.jpeg")
scaleFreePlot(k, main="Checksft\n")
dev.off()

# convert matrix to numeric
norm.counts_corrected[] <- sapply(norm.counts_corrected, as.numeric)

soft_power <- opt$soft_threshold  #8
temp_cor <- cor
cor <- WGCNA::cor




bwnet <- blockwiseModules(norm.counts_corrected,
                          maxBlockSize = opt$maxblocksize,
                          TOMType = "signed",
                          mergeCutHeight = 0.25,
                          power = soft_power,
                          saveTOMs = TRUE,
                          corType = "bicor",
                          networkType = "signed",
                          minKMEtoStay = 0,
                          saveTOMFileBase = "TOMfilebase",
                          pamRespectsDendro = FALSE,
                          numericLabels = FALSE,
                          randomSeed = 54321,
                          nThreads = opt$threads,
                          verbose = 3)

head(bwnet$MEs, n = 5)
#bwnet$TOMsimilarityFromExpr

#TOM = TOMsimilarity
cor <- temp_cor


# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs                         


# Print out a preview
head(module_eigengenes)

# get number of genes for each module
table(bwnet$colors)
corMatrix <- cor(bwnet$MEs)
setwd("/home/amrit/final_neurodegeneration_datasets/AD_NGS_data_new/CA3neurons/new/")
#library(WGCNA)
pMatrix <- corPvalueStudent(corMatrix, nrow(bwnet$MEs))
moduleColors <- gsub("ME", "", colnames(bwnet$MEs))  # Remove "ME" prefix to get colors

# Load pheatmap for heatmap visualization


# Plot heatmap
png("merge_module_heatmap1.png", width = 2800, height = 2000, res = 300)
pheatmap(corMatrix,
         labels_row = moduleColors,
         labels_col = moduleColors,
         color = colorRampPalette(c("white", "lightyellow", "yellow", "orange","red"))(50),
         main = "Correlation Heatmap of Module Eigengenes",
         display_numbers = TRUE)  # Remove "ME" prefix to get colors
dev.off()
# Load pheatmap for heatmap visualization


moduleTree <- hclust(as.dist(1 - corMatrix), method = "average")

png("merge_module_heatmap2.png", width = 2800, height = 2000, res = 300)
pheatmap(corMatrix,
         labels_row = moduleColors,
         labels_col = moduleColors,
         cluster_rows = moduleTree,
         cluster_cols = moduleTree,
         color = colorRampPalette(c("white", "lightyellow", "yellow", "orange","red"))(50),
         main = "Module Eigengene Dendrogram Heatmap",
         display_numbers = FALSE)
dev.off()
png("colordendogram_soft13_new.jpeg", units = "in", width = 6, height = 5, res = 600)
plt <- plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                           c("unmerged", "merged"),
                           dendroLabels = FALSE,
                           addGuide = TRUE,
                           hang= 0.03,
                           guideHang = 0.05)
dev.off()





head(colData, n=2)
traits <- colData %>% 
  mutate(disease_state_bin = ifelse(grepl('Alzheimer', disease), 1, 0))%>%
  dplyr::select(4)
head(traits,n=2)

# binarize categorical variables

colData$sample_group <- factor(colData$Type, levels = c("Ctrl", "FAD", "SAD"))

genotype.out <- binarizeCategoricalColumns(colData$sample_group,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)
head(genotype.out)


traits <- cbind(traits, genotype.out)
head(traits)
traits$Control <- 1- traits$disease_state_bin
traits <- traits[, !colnames(traits) %in% "disease_state_bin"]
colnames(traits) <- c('FAD', 'SAD', 'Control')
traits <- traits[, c("Control", "FAD", "SAD")]

# Define numbers of genes and samples
nSamples <- nrow(norm.counts_corrected)
nGenes <- ncol(norm.counts_corrected)

fad_only_idx <- colData$sample_group %in% c("FAD", "Ctrl")  # 60 samples
fad_eigengenes <- module_eigengenes[fad_only_idx, ]
fad_traits <- data.frame(FAD = colData$sample_group[fad_only_idx] == "FAD")
fad_module.trait.corr <- cor(fad_eigengenes, fad_traits, use='p')
fad_module.trait.corr.pvals <- corPvalueStudent(fad_module.trait.corr, nSamples)


fad_module.trait.fdr <- matrix(
  p.adjust(fad_module.trait.corr.pvals, method = "fdr"),
  nrow = nrow(fad_module.trait.corr),
  ncol = ncol(fad_module.trait.corr),
  dimnames = dimnames(fad_module.trait.corr)
)

fad_textMatrix <- matrix(
  paste(signif(fad_module.trait.corr[,1], 2), 
        "\n(", signif(fad_module.trait.corr.pvals[,1], 2), ")", sep=""),
  nrow = nrow(fad_module.trait.corr),
  ncol = 1,
  dimnames = list(rownames(fad_module.trait.corr), "FAD")
)


jpeg("fad_heatmap.jpeg", width = 6000, height = 15000, res = 400)

# Big bottom margin for angled labels
par(mar = c(30, 20, 8, 10))

labeledHeatmap(
  Matrix = fad_module.trait.corr,
  xLabels = "Control vs FAD",
  yLabels = rownames(fad_module.trait.corr),
  ySymbols = substr(rownames(fad_module.trait.corr), 3, 15),
  colorLabels = TRUE,
  colors = blueWhiteRed(50),
  textMatrix = fad_textMatrix,
  setStdMargins = FALSE,
  cex.text = 1.5,
  cex.lab = 2.0,
  zlim = c(-1, 1),
  xLabelsAngle = 45   # 🔑 keep diagonal labels but with more room
)

dev.off()

sad_only_idx <- colData$sample_group %in% c("SAD", "Ctrl")  # 60 samples
sad_eigengenes <- module_eigengenes[sad_only_idx, ]
sad_traits <- data.frame(SAD = colData$sample_group[sad_only_idx] == "SAD")
sad_module.trait.corr <- cor(sad_eigengenes, sad_traits, use='p')
sad_module.trait.corr.pvals <- corPvalueStudent(sad_module.trait.corr, nSamples)


sad_module.trait.fdr <- matrix(
  p.adjust(sad_module.trait.corr.pvals, method = "fdr"),
  nrow = nrow(sad_module.trait.corr),
  ncol = ncol(sad_module.trait.corr),
  dimnames = dimnames(sad_module.trait.corr)
)

sad_textMatrix <- matrix(
  paste(signif(sad_module.trait.corr[,1], 2), 
        "\n(", signif(sad_module.trait.corr.pvals[,1], 2), ")", sep=""),
  nrow = nrow(sad_module.trait.corr),
  ncol = 1,
  dimnames = list(rownames(sad_module.trait.corr), "SAD"))

jpeg("sad_heatmap.jpeg", width = 6000, height = 15000, res = 400)

# Big bottom margin for angled labels
par(mar = c(30, 20, 8, 10))

labeledHeatmap(
  Matrix = sad_module.trait.corr,
  xLabels = "Control vs SAD",
  yLabels = rownames(sad_module.trait.corr),
  ySymbols = substr(rownames(sad_module.trait.corr), 3, 15),
  colorLabels = TRUE,
  colors = blueWhiteRed(50),
  textMatrix = sad_textMatrix,
  setStdMargins = FALSE,
  cex.text = 1.5,
  cex.lab = 2.0,
  zlim = c(-1, 1),
  xLabelsAngle = 45   # 🔑 keep diagonal labels but with more room
)

dev.off()

module.gene.mapping <- as.data.frame(bwnet$colors)
head(module.gene.mapping, n=50)
write.csv(module.gene.mapping, "genes_found_in_various_modules.csv")
module.membership.measure <- cor(module_eigengenes, norm.counts_corrected, use = 'p')
#module.membership.measure
module.membership.measure.pvals <- corPvalueStudent(as.matrix(module.membership.measure), nSamples)
#module.membership.measure.pvals
write.csv(module.membership.measure.pvals, "MM_genes.csv")
module.membership.measure.pvals[1:6,1:10]
# Calculate the gene significance and associated p-values
head(traits)
gene.signf.corr1 <- cor(norm.counts_corrected, traits$data.FAD.vs.all, use = 'p')
colnames(gene.signf.corr1)[1] <- "Y"
gene.signf.corr.pvals1 <- corPvalueStudent(gene.signf.corr1, nSamples)
corp_pval_FAD_vs_Ctrl <- gene.signf.corr.pvals1 %>% 
  as.data.frame() %>% 
  arrange(Y)
write.csv(corp_pval_FAD_vs_Ctrl, "corp_pval_FAD_vs_Ctrl.csv", row.names = TRUE)



gene.signf.corr2 <- cor(norm.counts_corrected, traits$data.SAD.vs.all, use = 'p')
colnames(gene.signf.corr2)[1] <- "Z"
gene.signf.corr.pvals2 <- corPvalueStudent(gene.signf.corr2, nSamples)
corp_pval_SAD_vs_Ctrl <- gene.signf.corr.pvals2 %>% 
  as.data.frame() %>% 
  arrange(Z)
write.csv(corp_pval_SAD_vs_Ctrl, "corp_pval_SAD_vs_Ctrl.csv", row.names = TRUE)




TOM = TOMsimilarityFromExpr(norm.counts_corrected, power=8)
annot = read.csv(file = "ensemble_to_gene_ids_v113.csv")
setwd("/home/amrit/final_neurodegeneration_datasets/AD_NGS_data_new/CA3neurons/new/cytoscape")
module_FAD = c('green', 'blue', 'red')
module_SAD = c('skyblue', 'royalblue', 'floralwhite', 'lightcyan')
ens_mods = colnames(norm.counts_corrected)
ens_mods
inModule = is.finite(match(bwnet$colors, module_SAD))
modens <- ens_mods[inModule]
modGenes <- annot$gene_name[match(modens, annot$ensembl_id)]
modTOM <- TOM[inModule, inModule]
dimnames(modTOM) <- list(modens, modens)
cyt <- exportNetworkToCytoscape(modTOM,
                                edgeFile = paste("Cytoscape-edges15-SAD", paste(module_SAD, collapse = "-"), ".txt", sep=""),
                                nodeFile = paste("Cytoscape-nodes15-SAD-", paste(module_SAD, collapse = "-"), ".txt", sep=""),
                                weighted = TRUE,
                                threshold = 0.01,
                                nodeNames = modens,
                                altNodeNames = modGenes,
                                nodeAttr = bwnet$colors[inModule])
cyt


