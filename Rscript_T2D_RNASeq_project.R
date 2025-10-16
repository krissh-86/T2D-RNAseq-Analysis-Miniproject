# Set working directory 
setwd("C:/Users/Dell/Desktop/Bversity/Mini project/Galaxy/Final project R script")  # Replace with your path

# Load libraries
library(dplyr)
library(DESeq2)
library(tidyverse)
library(GEOquery)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(biomaRt)

# Step 1: Load raw counts & metadata
counts <- read.table("Merged_Feature_Counts_Final.tabular", header=TRUE, row.names=1, sep="\t", check.names=FALSE)
metadata <- read.csv("metadata.csv", header=TRUE, row.names=1)
metadata$condition <- factor(metadata$condition, levels = c("Non_T2D", "T2D"))
all(colnames(counts) %in% rownames(metadata))  # Should return TRUE


# Step 2: Create DESeq2 data set & pre-filter
dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=metadata,
                              design=~condition)

# Step 3: Run DESeq2 analysis
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "T2D", "Non_T2D"))
res_df <- as.data.frame(res)
head(res)
write.csv(as.data.frame(res), "DESeq2_results.csv")


# Step 4: MA plot & volcano plot
plotMA(res)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Volcano Plot of DEGs')
# Order results by adjusted p-value (padj)
res_sorted <- res[order(res$padj), ]

# Step 5: Filter significant genes
sig_genes <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(sig_genes, "Significant_DEGs.csv", row.names = TRUE)

# Step 6: Separate upregulated and downregulated genes
upregulated_genes <- subset(res_df, padj < 0.05 & log2FoldChange > 0)
write.csv(upregulated_genes, "Upregulated_genes.csv", row.names = TRUE)

downregulated_genes <- subset(res_df, padj < 0.05 & log2FoldChange < 0)
write.csv(downregulated_genes, "Downregulated_genes.csv", row.names = TRUE)

# Step 7: Normalized counts heatmap for significant genes
vsd <- vst(dds, blind = FALSE)
norm_counts <- assay(vsd)
sig_gene_names <- rownames(sig_genes)
sig_norm_counts <- norm_counts[sig_gene_names,]

annotation_col <- data.frame(Type = factor(metadata$condition))
rownames(annotation_col) <- colnames(sig_norm_counts)

pheatmap(sig_norm_counts,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
         main = "Heatmap of Significant DEGs")

# Step 8: PCA plot for sample clustering
pca <- prcomp(t(norm_counts))
percent_var <- (pca$sdev^2) / sum(pca$sdev^2) * 100
pca_data <- as.data.frame(pca$x)
pca_data$Type <- meta_data$condition  # Use "condition" here, not "Type"

ggplot(pca_data, aes(PC1, PC2, color = Type)) +
  geom_point(size = 3) +
  labs(
    title = "PCA of RNA-seq samples",
    x = sprintf("PC1 (%.2f%%)", percent_var[1]),
    y = sprintf("PC2 (%.2f%%)", percent_var[2])
  ) +
  theme_minimal()

# Step 9: Sample distance heatmap
rld <- rlog(dds, blind = TRUE)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = "Sample-to-Sample Distance Heatmap")

# Step 10: Dispersion plot
plotDispEsts(dds, main = "Dispersion Plot")

# Step 11: Gene annotation using biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annot <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol', 'description'),
                              filters = 'entrezgene_id',
                              values = sig_gene_names,
                              mart = ensembl)
sig_genes_df <- as.data.frame(sig_genes)
sig_genes_df <- rownames_to_column(sig_genes_df, var = "entrezgene_id")
annot$entrezgene_id <- as.character(annot$entrezgene_id)
sig_genes_annot <- left_join(sig_genes_df, annot, by = "entrezgene_id")
head(sig_genes_annot)
write.csv(sig_genes_annot, "Annotated_Significant_DEGs.csv", row.names = FALSE)


# Step 12: Functional enrichment analysis (GO & KEGG)
entrez_ids <- sig_genes_annot$entrezgene_id %>% na.omit() %>% unique()

# Biological Process (bp)
ego_bp <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                   ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05, readable = TRUE)
# Molecular Function (mf)
ego_mf <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                   ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05, readable = TRUE)
# Cellular Component (cc)
ego_cc <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                   ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05, readable = TRUE)
#KEGG pathway enrichment
kegg_enrich <- enrichKEGG(gene = entrez_ids, # or filtered_entrez if you used that
                          organism = "hsa",  # Homo sapiens
                          pvalueCutoff = 0.05)

# Step 13: Visualization of enrichment results
dotplot(ego_bp, showCategory = 20) + ggtitle("GO: Biological Process")
dotplot(ego_mf, showCategory = 20) + ggtitle("GO: Molecular Function")
dotplot(ego_cc, showCategory = 20) + ggtitle("GO: Cellular Component")
dotplot(kegg_enrich, showCategory = 20) + ggtitle("KEGG Pathway Enrichment")

# Step 14: Load Enrichr clinical overlap results (manually exported CSV)
enrichr_table <- read.csv("Annotated_Significant_DEGs_WithEnrichrHits.csv")

# Replace NA with 0 and sort decreasing by overlap count
enrichr_table$Diabetes_Study_HitCount[is.na(enrichr_table$Diabetes_Study_HitCount)] <- 0
enrichr_table <- enrichr_table[order(-enrichr_table$Diabetes_Study_HitCount), ]

# Save sorted annotation
write.csv(enrichr_table, "Final_Annotated_DEGs_With_DiabetesHits_Sorted.csv", row.names=FALSE)

# Step 15: Biomarker prioritization - rank by overlap and fold change magnitude
enrichr_table$biomarker_rank_score <- enrichr_table$Diabetes_Study_HitCount * abs(enrichr_table$log2FoldChange)
enrichr_table <- enrichr_table[order(-enrichr_table$biomarker_rank_score), ]

# Select top 10 biomarker candidates
top_biomarkers <- head(enrichr_table, 10)
write.csv(top_biomarkers, "Top_10_Biomarker_Candidates.csv", row.names=FALSE)

# Step 16: Visualize top biomarker candidates (barplot)
library(ggplot2)
ggplot(top_biomarkers, aes(x = reorder(hgnc_symbol, Diabetes_Study_HitCount), y = Diabetes_Study_HitCount)) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() +
  labs(title="Top 10 Biomarker Candidates by Diabetes Study Hit Count",
       x="Gene Symbol",
       y="Number of Diabetes GEO Overlaps") +
  theme_minimal()

# Step 17: Top_10_Biomarker_Genes_Summary_table
# Load dplyr if not already loaded
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)

# Create summary table with selected columns
top10_genes <- head(enrichr_table, 10)

summary_table <- dplyr::select(top10_genes, hgnc_symbol, log2FoldChange, padj, Diabetes_Study_HitCount, description)

# Print summary table to console
print(summary_table)

# Export summary table to CSV for reporting/sharing
write.csv(summary_table, "Top_10_Biomarker_Genes_Summary.csv", row.names = FALSE)

