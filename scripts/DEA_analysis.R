#Load required libraries
library(DESeq2)
library(tidyverse)
library(dplyr)
library(tibble)
#Loading required data files
args <- commandArgs(trailingOnly=TRUE)
count_file <- args[1]
raw_counts <- read.csv(count_file, header = TRUE, row.names = "Geneid", stringsAsFactors = FALSE)
raw_counts <- raw_counts[,sort(colnames(raw_counts))]
colSums(raw_counts)

condition <- c(rep("LNCAP_Hypoxia", 2), rep("LNCAP_Normoxia", 2), rep("PC3_Hypoxia", 2), rep("PC3_Normoxia", 2))

my_colData <- as.data.frame(condition)
rownames(my_colData) <- colnames(raw_counts)

#creating dds object
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = my_colData,
                              design = ~condition)

count_matrix <- counts(dds)

#Annotation file can be downloaded from Ensembl BioMart (http://uswest.ensembl.org/biomart/martview/)
#Loading annotation file and filtering for specific biotypes
annotation_file <- args[2]
annotation <- read.csv(annotation_file, header=TRUE, stringsAsFactors = FALSE)

counts_gse <- read.csv(count_file,
                     header = TRUE,
                     stringsAsFactors = FALSE)
#remove version numbers on gene from both dfs
counts_gse$Geneid <- sub("\\..*$", "", counts_gse$Geneid)
annotation$Geneid <- sub("\\..*$", "", annotation$Geneid)
annotated_counts <- left_join(counts_gse, annotation, by = "Geneid") %>%
  select(Geneid, Genesymbol, Genetype, 
         LNCAP_Hypoxia_S1, LNCAP_Hypoxia_S2, LNCAP_Normoxia_S1, LNCAP_Normoxia_S2, 
         PC3_Hypoxia_S1, PC3_Hypoxia_S2, PC3_Normoxia_S1, PC3_Normoxia_S2)

biotypes_to_keep <- c("protein_coding", "IG_J_gene", "IG_V_gene", "IG_C_gene", "IG_D_gene", "TR_D_gene", "TR_C_gene", "TR_V_gene", "TR_J_gene")

filtered_counts <- annotated_counts %>% filter(Genetype %in% biotypes_to_keep)


zero_counts1 <- rowSums(filtered_counts[, 4:11] == 0)
zero_summary2 <- table(zero_counts1)

#Filtering genes based on zero counts
#keep genes that have non zero counts in 2 samples 
keep_genes <- zero_counts1 < 7
filtered_counts_nozero <- filtered_counts[keep_genes, ]

filtered_counts_output <- args[3]
write.csv(filtered_counts_nozero, file = filtered_counts_output, sep = ",", row.names = FALSE)

#remove version numbers from dds
rownames(dds) <- sub("\\..*$", "", rownames(dds))
#filter dds object to keep only genes that we have saved after filtering zero counts
dds_filtered <- dds[rownames(dds) %in% filtered_counts_nozero$Geneid, ]

#run DESeq analysis and normalization 
dds <- DESeq(dds_filtered)
dds
normalized_counts <- counts(dds, normalized = T)
normalized_counts_df <- as.data.frame(normalized_counts)
normalized_counts_output <- args[4]
write.csv(normalized_counts_df, file = normalized_counts_output, row.names = TRUE)

#Visualizations for sample variability

#Variance stabilizing transformation
vsd <- vst(dds_filtered, blind = TRUE)  # blind=TRUE for exploratory PCA
plot_PCA = function (vsd.obj) {
  pcaData <- plotPCA(vsd.obj,  intgroup = c("condition"), returnData = T)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    labs(x = paste0("PC1: ",percentVar[1],"% variance"),
         y = paste0("PC2: ",percentVar[2],"% variance"),
         title = "PCA Plot colored by condition") +
    ggrepel::geom_text_repel(aes(label = name), color = "black")
}

png(filename = "pcab.png", 
    width = 2000, height = 2000, res = 300)  # adjust width/height as needed
plot_PCA(vsd)
dev.off()

plotDists = function (vsd.obj) {
  sampleDists <- dist(t(assay(vsd.obj)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd.obj$condition)
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(55)
  pheatmap::pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,  clustering_distance_cols = sampleDists, col = colors, fontsize_row = 4, fontsize_col = 4, fontsize_legend = 4, fontsize = 4)
}
png(filename = "distance_plot.png", width = 1000, height = 900, res = 300)  # adjust width/height as needed
plotDists(vsd)
dev.off()

#heatmap for top 500 most variable genes across all samples
variable_gene_heatmap <- function (vsd.obj, num_genes = 500, annotation, title = "Heatmap") {
  brewer_palette <- "RdBu"
  ramp <- colorRampPalette( RColorBrewer::brewer.pal(11, brewer_palette))
  mr <- ramp(256)[256:1]
  stabilized_counts <- assay(vsd.obj)
  row_variances <- rowVars(stabilized_counts)
  top_variable_genes <- stabilized_counts[order(row_variances, decreasing=T)[1:num_genes],]
  top_variable_genes <- top_variable_genes - rowMeans(top_variable_genes, na.rm=T)
  gene_names <- annotation$Genesymbol[match(rownames(top_variable_genes), annotation$Geneid)]
  rownames(top_variable_genes) <- gene_names
  coldata <- as.data.frame(vsd.obj@colData)
  coldata$sizeFactor <- NULL
  pheatmap::pheatmap(top_variable_genes, color = mr, annotation_col = coldata, fontsize_col = 8, fontsize_row = 250/num_genes, border_color = NA, main = title)
}

png(filename = "variable_gene_heatmap.png", 
    width = 1000, height = 1700, res = 200)  # adjust width/height as needed
variable_gene_heatmap(vsd, num_genes = 40, annotation = annotation)
dev.off()

#comparison of density plots for raw and vst normalized counts for all samples
raw_counts <- assay(dds)
vst_counts <- assay(vsd)

png("density_plots_raw_vst.png",
    width = 4000, height = 4000, res = 300)  # Adjust width, height (pixels), and resolution (dpi)

par(mfrow = c(4, 4), mar = c(3, 3, 2, 1))  # mar adjusts margins (bottom, left, top, right)

for (i in 1:8) {
  # Raw counts density
  plot(density(raw_counts[, i]),
       main = paste("Raw - Sample", colnames(raw_counts)[i]),
       xlab = "Expression",
       col = "red",
       lwd = 2,
       ylim = c(0, max(sapply(1:8, function(j) max(density(raw_counts[, j])$y, na.rm = TRUE))))  # Uniform y-axis
       )
  
  # VST counts density (next panel)
  plot(density(vst_counts[, i]),
       main = paste("VST - Sample", colnames(vst_counts)[i]),
       xlab = "Expression",
       col = "blue",
       lwd = 2,
       ylim = c(0, max(sapply(1:8, function(j) max(density(vst_counts[, j])$y, na.rm = TRUE))))  # Uniform y-axis
       )
}
dev.off()

#Differential expression analysis for each cell line separately
#Based on exploratory analysis, we can see that samples cluster by cell line (PC3 vs LNCAP) and then by condition (Hypoxia vs Normoxia) so we run DEA separately for each cell line
#for lncap cell line
dds_lncap <- dds[, grepl("LNCAP", colnames(dds))]
dds_lncap$condition <- droplevels(dds_lncap$condition)
dds_lncap$condition <- relevel(dds_lncap$condition, ref = "LNCAP_Normoxia")
dds_lncap <- DESeq(dds_lncap)
#extract results for lncap DEA
res_lncap <- results(dds_lncap, contrast = c("condition", "LNCAP_Hypoxia", "LNCAP_Normoxia"))

reslncapOrdered <- res_lncap[order(res_lncap$padj), ]

sum(reslncapOrdered$padj < 0.05, na.rm = TRUE)
write.csv(as.data.frame(reslncapOrdered), file = "DEGs_lncap.csv")

#Repeat for PC3
dds_pc3 <- dds[, grepl("PC3", colnames(dds))]
dds_pc3$condition <- droplevels(dds_pc3$condition)
dds_pc3$condition <- relevel(dds_pc3$condition, ref = "PC3_Normoxia")
dds_pc3 <- DESeq(dds_pc3)
#extract results for pc3 DEA
res_pc3 <- results(dds_pc3, contrast = c("condition", "PC3_Hypoxia", "PC3_Normoxia"))

respc3Ordered <- res_pc3[order(res_pc3$padj), ]

sum(respc3Ordered$padj < 0.05, na.rm = TRUE)
write.csv(as.data.frame(respc3Ordered), file = "DEGs_pc3.csv")

#volcano plot for lncap
res_df <- as.data.frame(reslncapOrdered)
res_df <- na.omit(res_df)
res_df$gene <- rownames(res_df)

res_df$regulation <- "Not Significant"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

qp <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "#FEA405", 
                                "Downregulated" = "purple", 
                                "Not Significant" = "gray")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  annotate("text", x = min(res_df$log2FoldChange), y = -log2(0.05) + 0.5,
           label = "padj = 0.05", hjust = 0, size = 3) +
  theme_minimal() +
  labs(title = "Volcano Plot", 
       x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-Value") +
  theme(plot.title = element_text(hjust = 0.5))

v_plot <- "vp_lncap.png"
ggsave(v_plot, plot = qp, width = 8, height = 6, dpi = 300)

#volcano plot for pc3
res_df <- as.data.frame(respc3Ordered)
res_df <- na.omit(res_df)
res_df$gene <- rownames(res_df)

res_df$regulation <- "Not Significant"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

qp <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "#FEA405", 
                                "Downregulated" = "purple", 
                                "Not Significant" = "gray")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  annotate("text", x = min(res_df$log2FoldChange), y = -log2(0.05) + 0.5,
           label = "padj = 0.05", hjust = 0, size = 3) +
  theme_minimal() +
  labs(title = "Volcano Plot", 
       x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-Value") +
  theme(plot.title = element_text(hjust = 0.5))

v_plot <- "vp_pc3.png"
ggsave(v_plot, plot = qp, width = 8, height = 6, dpi = 300)

#heatmap for top 20 DE genes in lncap
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(dplyr)

DE_gene_heatmap <- function(res_lncap, count_matrix, padj_cutoff = 0.0001, ngenes = 20) {
  # Generate the color palette
  brewer_palette <- "RdBu"
  ramp <- colorRampPalette(RColorBrewer::brewer.pal(11, brewer_palette))
  mr <- ramp(256)[256:1]  # Reversed palette (blue to red)

  # Convert DESeqResults to data frame and get significant genes
  significant_genes <- as.data.frame(res_lncap) %>%
    filter(padj < padj_cutoff) %>%
    arrange(desc(log2FoldChange)) %>%
    head(ngenes)

  # Extract count data for significant genes
  heatmap_values <- count_matrix[rownames(significant_genes), ]  # Use row names (Ensembl IDs)
  
  #map ensembl IDs in res_lncap to gene symbols
  gene_symbols <- annotation$Genesymbol[match(rownames(significant_genes), annotation$Geneid)]
  rownames(heatmap_values) <- gene_symbols
  
  # Scale rows for heatmap (z-score normalization)
  heatmap_values <- t(scale(t(heatmap_values)))

  # Plot the heatmap using pheatmap
  p <- pheatmap::pheatmap(heatmap_values,
                          color = mr,
                          scale = "none",  # Already scaled
                          cluster_rows = TRUE,
                          cluster_cols = TRUE,
                          fontsize_col = 10,
                          fontsize_row = max(6, 200/ngenes),  # Minimum font size of 6
                          border_color = NA,
                          main = paste("Top", ngenes, "DE Genes (padj <", padj_cutoff, ")"))

  # Return the pheatmap object
  return(invisible(p))
}

count_matrix <- assay(dds_lncap)  # Replace with your count matrix if different
png("de_gene_heatmap_lncap.png",
    width = 1000, height = 1700, res = 200)  # Adjust dimensions and resolution
d <- DE_gene_heatmap(res_lncap, count_matrix, padj_cutoff = 0.001, ngenes = 30)
dev.off()

#heatmap for top 20 DE genes in pc3
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(dplyr)

DE_gene_heatmap <- function(res_pc3, count_matrix, padj_cutoff = 0.0001, ngenes = 20) {
  # Generate the color palette
  brewer_palette <- "RdBu"
  ramp <- colorRampPalette(RColorBrewer::brewer.pal(11, brewer_palette))
  mr <- ramp(256)[256:1]  # Reversed palette (blue to red)

  # Convert DESeqResults to data frame and get significant genes
  significant_genes <- as.data.frame(res_pc3) %>%
    filter(padj < padj_cutoff) %>%
    arrange(desc(log2FoldChange)) %>%
    head(ngenes)

  # Extract count data for significant genes
  heatmap_values <- count_matrix[rownames(significant_genes), ]  # Use row names (Ensembl IDs)
  
  #map ensembl IDs in res_lncap to gene symbols
  gene_symbols <- annotation$Genesymbol[match(rownames(significant_genes), annotation$Geneid)]
  rownames(heatmap_values) <- gene_symbols
  
  # Scale rows for heatmap (z-score normalization)
  heatmap_values <- t(scale(t(heatmap_values)))

  # Plot the heatmap using pheatmap
  p <- pheatmap::pheatmap(heatmap_values,
                          color = mr,
                          scale = "none",  # Already scaled
                          cluster_rows = TRUE,
                          cluster_cols = TRUE,
                          fontsize_col = 10,
                          fontsize_row = max(6, 200/ngenes),  # Minimum font size of 6
                          border_color = NA,
                          main = paste("Top", ngenes, "DE Genes (padj <", padj_cutoff, ")"))

  # Return the pheatmap object
  return(invisible(p))
}

count_matrix <- assay(dds_pc3)  # Replace with your count matrix if different
png("de_gene_heatmap_pc3.png",
    width = 1000, height = 1700, res = 200)  # Adjust dimensions and resolution
d <- DE_gene_heatmap(res_pc3, count_matrix, padj_cutoff = 0.001, ngenes = 30)
dev.off()

#Gene set enrichment analysis for Lncap DEGs
#Downloaded hallmark gene sets from MSigDB (https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp)
library(fgsea)
pathway_file <- args[5]
hallmark_pathway <- gmtPathways(pathway_file)
library(org.Hs.eg.db)

res_lncap$symbol <- mapIds(org.Hs.eg.db,
                           keys = rownames(res_lncap),
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")
res_clean <- res_lncap[!is.na(res_lncap$symbol), ]        # remove NAs
res_clean <- res_clean[!duplicated(res_clean$symbol), ]   # remove duplicates

# named numeric vector: names = symbols, values = log2FC
lncap_ranked_list <- res_clean$log2FoldChange
names(lncap_ranked_list) <- res_clean$symbol

# sort decreasing (most up first)
lncap_ranked_list <- sort(lncap_ranked_list, decreasing = TRUE)

fgsea_results <- fgsea(pathways = hallmark_pathway,
                       stats = lncap_ranked_list,
                       minSize = 15,
                       maxSize = 500,
                       nperm = 1000)
fgsea_results_ordered <- fgsea_results[order(-NES)]
head(fgsea_results_ordered[, .(pathway, padj, NES)])
plotEnrichment(hallmark_pathway[["HALLMARK_HYPOXIA"]], lncap_ranked_list)

#waterfall plot for fgsea results
waterfall_plot <- function (fsgea_results, graph_title) {
  fgsea_results %>% 
    mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>% # removes 'HALLMARK_' from the pathway title 
    ggplot( aes(reorder(short_name,NES), NES)) +
      geom_bar(stat= "identity", aes(fill = padj<0.05))+
      coord_flip()+
      labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title)+
      theme(axis.text.y = element_text(size = 7), 
            plot.title = element_text(hjust = 1))
}
library(stringr)
p <- waterfall_plot(fgsea_results, "Hallmark pathways altered by hypoxia in LNCaP cells")
ggsave("GSEA_hallmark_lncap.png", p, width = 8, height = 6, dpi = 300)

#Gene set enrichment analysis for PC3 DEGs

library(fgsea)
hallmark_pathway <- gmtPathways(pathway_file)
library(org.Hs.eg.db)

res_pc3$symbol <- mapIds(org.Hs.eg.db,
                           keys = rownames(res_pc3),
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")
res_clean <- res_pc3[!is.na(res_pc3$symbol), ]        # remove NAs
res_clean <- res_clean[!duplicated(res_clean$symbol), ]   # remove duplicates

# named numeric vector: names = symbols, values = log2FC
pc3_ranked_list <- res_clean$log2FoldChange
names(pc3_ranked_list) <- res_clean$symbol

# sort decreasing (most up first)
pc3_ranked_list <- sort(pc3_ranked_list, decreasing = TRUE)

fgsea_results <- fgsea(pathways = hallmark_pathway,
                       stats = pc3_ranked_list,
                       minSize = 15,
                       maxSize = 500,
                       nperm = 1000)
fgsea_results_ordered <- fgsea_results[order(-NES)]
head(fgsea_results_ordered[, .(pathway, padj, NES)])
plotEnrichment(hallmark_pathway[["HALLMARK_HYPOXIA"]], pc3_ranked_list)

#waterfall plot for fgsea results
waterfall_plot <- function (fsgea_results, graph_title) {
  fgsea_results %>% 
    mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>% # removes 'HALLMARK_' from the pathway title 
    ggplot( aes(reorder(short_name,NES), NES)) +
      geom_bar(stat= "identity", aes(fill = padj<0.05))+
      coord_flip()+
      labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title)+
      theme(axis.text.y = element_text(size = 7), 
            plot.title = element_text(hjust = 1))
}
library(stringr)
p <- waterfall_plot(fgsea_results, "Hallmark pathways altered by hypoxia in PC3 cells")
ggsave("GSEA_hallmark_pc3.png", p, width = 8, height = 6, dpi = 300)

#GSEA using Reactome pathways for lncap DEGs
res_lncap <- read.csv("DEGs_lncap.csv", row.names = 1)

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(stats)

ncbi_list <- clusterProfiler::bitr(
  geneID = rownames(res_lncap),        # use Ensembl IDs from row names
  fromType = "ENSEMBL",          
  toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db
)

res_lncap$ENSEMBL <- rownames(res_lncap)

res_mapped <- res_lncap %>%
  left_join(ncbi_list, by = "ENSEMBL") %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(ENTREZID, .keep_all = TRUE)

ngenes <- res_mapped$log2FoldChange
names(ngenes) <- res_mapped$ENTREZID
ngenes <- sort(ngenes, decreasing = TRUE)

library(ReactomePA)
enp_gsea <- gsePathway(
  ngenes,
  organism = "human",
  #pvalueCutoff = 0.05,
  verbose = FALSE
)

pathways <- enp_gsea@result
pathways <- pathways[order(pathways$p.adjust), ]  # Sort by FDR (adjusted p-value)
top_pathways <- pathways[order(abs(pathways$NES), decreasing = TRUE), ]  # Sort by NES

library(dplyr)
library(forcats)

top20 <- top_pathways[1:20, ] %>%
  mutate(Description = fct_reorder(Description, NES))  # Reorder factor for y-axis

#ggplot for fgsea results
library(ggplot2)

r1 <- ggplot(top20, aes(x = NES,
                        y = Description,
                        color = p.adjust,
                        size = setSize)) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "#0072B2", high = "#D55E00", name = "FDR (p.adjust)") +
  scale_size(range = c(3, 10), name = "Gene Set Size") +
  labs(
    title = "Top 20 Enriched Pathways",
    subtitle = "Gene Set Enrichment Analysis (GSEA)",
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    caption = "Data source: clusterProfiler::gsePathway"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5),
    plot.title = element_text(face = "bold", size = 9),
    plot.subtitle = element_text(size = 7),
    legend.position = "right"
  )
ggsave("GSEA_reactome_lncap.png", r1, width = 8, height = 6, dpi = 300)

#GSEA using Reactome pathways for pc3 DEGs
res_pc3 <- read.csv("DEGs_pc3.csv", row.names = 1)
head(res_pc3)

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(stats)

ncbi_list <- clusterProfiler::bitr(
  geneID = rownames(res_pc3),        # use Ensembl IDs from row names
  fromType = "ENSEMBL",          
  toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db
)

res_pc3$ENSEMBL <- rownames(res_pc3)

res_mapped <- res_pc3 %>%
  left_join(ncbi_list, by = "ENSEMBL") %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(ENTREZID, .keep_all = TRUE)

ngenes <- res_mapped$log2FoldChange
names(ngenes) <- res_mapped$ENTREZID
ngenes <- sort(ngenes, decreasing = TRUE)

library(ReactomePA)
enp_gsea <- gsePathway(
  ngenes,
  organism = "human",
  #pvalueCutoff = 0.05,
  verbose = FALSE
)

pathways <- enp_gsea@result
pathways <- pathways[order(pathways$p.adjust), ]  # Sort by FDR (adjusted p-value)
top_pathways <- pathways[order(abs(pathways$NES), decreasing = TRUE), ]  # Sort by NES

library(dplyr)
library(forcats)

top20 <- top_pathways[1:20, ] %>%
  mutate(Description = fct_reorder(Description, NES))  # Reorder factor for y-axis

#ggplot for fgsea results
library(ggplot2)

r1 <- ggplot(top20, aes(x = NES,
                        y = Description,
                        color = p.adjust,
                        size = setSize)) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "#0072B2", high = "#D55E00", name = "FDR (p.adjust)") +
  scale_size(range = c(3, 10), name = "Gene Set Size") +
  labs(
    title = "Top 20 Enriched Pathways",
    subtitle = "Gene Set Enrichment Analysis (GSEA)",
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    caption = "Data source: clusterProfiler::gsePathway"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5),
    plot.title = element_text(face = "bold", size = 9),
    plot.subtitle = element_text(size = 7),
    legend.position = "right"
  )
ggsave("GSEA_reactome_pc3.png", r1, width = 8, height = 6, dpi = 300)