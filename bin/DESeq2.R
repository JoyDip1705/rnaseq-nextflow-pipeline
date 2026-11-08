#!/usr/bin/env Rscript

# ============================================================
# RNA-seq DIFFERENTIAL EXPRESSION ANALYSIS
# DESeq2
# Comparison: late_luad vs early_luad
# ============================================================

suppressPackageStartupMessages({
	    library(DESeq2)
	    library(ggplot2)
})

# ------------------------------------------------------------
# 1. INPUT FILES
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
	    stop(
		 "Usage: Rscript DESeq2.R ",
	         "raw_counts_matrix.csv sample_metadata.csv"
	    )
}

count_file <- args[1]
metadata_file <- args[2]

# ------------------------------------------------------------
# 2. READ COUNT MATRIX
# ------------------------------------------------------------

counts_data <- read.csv(
			count_file,
			header = TRUE,
			check.names = FALSE,
			stringsAsFactors = FALSE
			)

# First column contains gene IDs
gene_ids <- counts_data[[1]]

# Remaining columns contain sample counts
count_matrix <- as.matrix(counts_data[, -1, drop = FALSE])

rownames(count_matrix) <- gene_ids

# Ensure numeric counts
storage.mode(count_matrix) <- "numeric"

# ------------------------------------------------------------
# 3. READ SAMPLE METADATA
# ------------------------------------------------------------

metadata <- read.csv(
		     metadata_file,
		     header = TRUE,
		     stringsAsFactors = FALSE
)

# Check metadata columns
required_columns <- c("sample", "condition")

if (!all(required_columns %in% colnames(metadata))) {
	stop("Metadata must contain sample and condition columns")
}

# ------------------------------------------------------------
# 4. VALIDATE SAMPLE MATCHING
# ------------------------------------------------------------

if (anyDuplicated(metadata$sample)) {
	stop("Duplicate sample names found in metadata")
}

if (anyDuplicated(colnames(count_matrix))) {
	stop("Duplicate sample names found in count matrix")
}

missing_samples <- setdiff(
			   metadata$sample,
			   colnames(count_matrix)
)

if (length(missing_samples) > 0) {
	stop(
	     "Samples missing from count matrix: ",
	     paste(missing_samples, collapse = ", ")
	)
}

extra_samples <- setdiff(
			 colnames(count_matrix),
			 metadata$sample)

if (length(extra_samples) > 0) {
	stop(
	     "Samples in count matrix but absent from metadata: ",
	     paste(extra_samples, collapse = ", ")
	)
}

# Reorder metadata to match count matrix
metadata <- metadata[
		     match(colnames(count_matrix), metadata$sample),
		     ,
		     drop = FALSE
]

rownames(metadata) <- metadata$sample

# ------------------------------------------------------------
# 5. CHECK COUNT VALUES
# ------------------------------------------------------------

if (anyNA(count_matrix)) {
	stop("Missing values detected in count matrix")
}

if (any(count_matrix < 0)) {
	stop("Negative counts detected")
}

if (any(count_matrix != round(count_matrix))) {
	stop("DESeq2 requires raw integer counts")
}

count_matrix <- round(count_matrix)

# ------------------------------------------------------------
# 6. SET CONDITION LEVELS
# ------------------------------------------------------------

metadata$condition <- factor(
			     metadata$condition,
			     levels = c("early_luad", "late_luad")
)

if (anyNA(metadata$condition)) {
	stop("Unexpected condition found in metadata")
}

# ------------------------------------------------------------
# 7. CREATE DESEQ2 OBJECT
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
			      countData = count_matrix,
			      colData = metadata,
			      design = ~ condition
)

# ------------------------------------------------------------
# 8. FILTER LOW-COUNT GENES
# ------------------------------------------------------------

keep <- rowSums(counts(dds) >= 10) >= 2

dds <- dds[keep, ]

# ------------------------------------------------------------
# 9. RUN DESEQ2
# ------------------------------------------------------------

dds <- DESeq(dds)

saveRDS(dds, "DESeq2_object.rds")

# ------------------------------------------------------------
# 10. DIFFERENTIAL EXPRESSION
# ------------------------------------------------------------

res <- results(
	       dds,
	       contrast = c("condition", "late_luad", "early_luad"),
	       alpha = 0.05)

# ------------------------------------------------------------
# 11. SAVE COMPLETE RESULTS
# ------------------------------------------------------------

res_df <- as.data.frame(res)

res_df$gene_id <- rownames(res_df)

res_df <- res_df[
		 order(res_df$padj, na.last = TRUE),]

write.csv(res_df,"DESeq2_results.csv",row.names = FALSE)

# ------------------------------------------------------------
# 12. SAVE SIGNIFICANT DEGS
# ------------------------------------------------------------

sig_genes <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1,]

write.csv(sig_genes,"significant_DEGs.csv",row.names = FALSE)

# ------------------------------------------------------------
# 13. SAVE NORMALIZED COUNTS
# ------------------------------------------------------------

normalized_counts <- counts(dds, normalized = TRUE)

write.csv(normalized_counts, "normalized_counts.csv")

# ------------------------------------------------------------
# 14. VOLCANO PLOT
# ------------------------------------------------------------

res_df$significance <- "Not significant"

res_df$significance[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange >= 1] <- "Upregulated"

res_df$significance[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange <= -1] <- "Downregulated"

res_df$neg_log10_padj <- -log10(pmax(res_df$padj, 1e-300))

volcano <- ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj, color = significance)) +
	geom_point(alpha = 0.7, size = 1.5, na.rm = TRUE) +
	scale_color_manual(values = c("Downregulated" = "blue",
				      "Not significant" = "grey70",
				      "Upregulated" = "red")) +
	geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
	geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
	theme_bw() +
	labs(title = "late LUAD vs early LUAD",
	     x = "Log2 fold change",
	     y = "-Log10 adjusted P-value")
	 

ggsave("volcano_plot.png", volcano, width = 8, height = 6, dpi = 400)

# ------------------------------------------------------------
# 15. MA PLOT
# ------------------------------------------------------------

png("MA_plot.png", width = 1800, height = 1400, res = 200)

plotMA(res, ylim = c(-5, 5), main = "late LUAD vs early LUAD")
	       
dev.off()

# ------------------------------------------------------------
# 16. PCA PLOT
# ------------------------------------------------------------

vsd <- vst(dds, blind = TRUE)
	       
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
			    
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
	geom_point(size = 4) +
	geom_text(vjust = -1, size = 3) +
	theme_bw() +
	labs(title = "PCA of RNA-seq samples", x = paste0("PC1: ", percentVar[1], "% variance"), y = paste0("PC2: ", percentVar[2], "% variance"))

ggsave("PCA_plot.png", pca_plot, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------
# 17. HEATMAP (top 50 genes by adjusted p-value, ggplot2)
# ------------------------------------------------------------

top_genes <- head(res_df$gene_id[!is.na(res_df$padj)], 50)

heat_mat <- assay(vsd)[top_genes, , drop = FALSE]

# Z-score each gene across samples
heat_mat <- t(scale(t(heat_mat)))

heat_mat <- heat_mat[complete.cases(heat_mat), , drop = FALSE]  # drop zero-variance genes

# Hierarchical clustering to order rows and columns
row_order <- hclust(dist(heat_mat))$order

col_order <- hclust(dist(t(heat_mat)))$order

gene_levels   <- rownames(heat_mat)[row_order]
sample_levels <- colnames(heat_mat)[col_order]

# Long format without extra packages
heat_df <- data.frame(
		      	gene   = rep(rownames(heat_mat), times = ncol(heat_mat)),
				sample = rep(colnames(heat_mat), each = nrow(heat_mat)),
				zscore = as.vector(heat_mat),
					stringsAsFactors = FALSE
				)

heat_df$condition <- metadata$condition[match(heat_df$sample, metadata$sample)]
heat_df$sample_label <- paste0(heat_df$sample, "\n(", heat_df$condition, ")")

heat_df$gene <- factor(heat_df$gene, levels = rev(gene_levels))

heat_df$sample_label <- factor(
			       	heat_df$sample_label,
					levels = unique(heat_df$sample_label[order(match(heat_df$sample, sample_levels))])
				)

heatmap_plot <- ggplot(heat_df, aes(x = sample_label, y = gene, fill = zscore)) +
		geom_tile(color = "white", linewidth = 0.2) +
		scale_fill_gradient2(low = "blue", mid = "white", high = "red",
				     midpoint = 0, name = "Z-score") +
		theme_minimal() +
		theme(axis.text.y = element_text(size = 6),
					      	 axis.text.x = element_text(size = 8),
						 panel.grid = element_blank()) +
		labs(title = "Top 50 DE genes: late LUAD vs early LUAD",
						     	     x = NULL, y = NULL)

ggsave("heatmap_top50.png", heatmap_plot, width = 7, height = 10, dpi = 400)

# ------------------------------------------------------------
# 18. SUMMARY
# ------------------------------------------------------------

cat("\nDESeq2 analysis completed\n")

cat("Total genes tested:", nrow(res_df), "\n")

cat("Significant DEGs:", nrow(sig_genes), "\n")

cat("Upregulated:", sum(sig_genes$log2FoldChange > 0), "\n")

cat("Downregulated:", sum(sig_genes$log2FoldChange < 0), "\n")

cat("\nComparison: late_luad vs early_luad\n")

