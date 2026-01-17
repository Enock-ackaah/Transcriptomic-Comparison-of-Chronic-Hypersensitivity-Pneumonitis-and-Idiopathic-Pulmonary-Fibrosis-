# ============================================================
# Project 4: Reproducing GSE150910 (CHP vs IPF vs Control)
# Step 2 — Data Ingestion (Counts) + Basic Structural Verification
# ============================================================

# ---- 0) User-provided path to the GEO processed count matrix ----
file_path <- "C:/Users/enock_p22oyv9/OneDrive/Desktop/BIOINFORMATICIAN/BIOINFORMATICS PROJECTS/Project 4/GSE150910_gene-level_count_file.csv.gz"

# ---- 1) Confirming the file exists (prevents silent path errors) ----
stopifnot(file.exists(file_path))

# ---- 2) Read the gzipped CSV directly (R can read .csv.gz without manual unzip) ----
# Why: This loads the gene-level count table into memory for analysis.
counts_df <- read.csv(gzfile(file_path), stringsAsFactors = FALSE, check.names = FALSE)

# ---- 3) Quick structural checks: rows/cols and column names ----
# Why: Confirms we loaded a gene-by-sample table and didn’t read a wrong file.
cat("Counts data frame dimensions (rows x cols): ", paste(dim(counts_df), collapse = " x "), "\n")
cat("First 10 column names:\n")
print(head(colnames(counts_df), 10))

# ---- 4) Confirming the first column is the gene identifier (usually 'symbol') ----
# Why: We must keep gene IDs separate from numeric count columns.
first_col_name <- colnames(counts_df)[1]
cat("First column name (gene ID column):", first_col_name, "\n")

# If the first column isn't 'symbol', we still proceed, but we treat it as gene ID.
gene_id <- counts_df[[1]]

# ---- 5) Basic sanity check on gene IDs ----
# Why: Ensures gene identifiers are present and not empty.
cat("Example gene IDs (first 10):\n")
print(head(gene_id, 10))

# ---- 6) Create the numeric count matrix (genes x samples) ----
# Why: DESeq2/edgeR operate on a numeric matrix of counts.
# Remove the gene ID column and convert the remaining columns to a matrix.
count_mat <- as.matrix(counts_df[, -1, drop = FALSE])

# ---- 7) Force numeric storage mode safely ----
# Why: read.csv can import numbers as character in some cases; DE tools require numeric counts.
storage.mode(count_mat) <- "numeric"

# ---- 8) Set rownames to gene IDs ----
# Why: Makes downstream subsetting and annotation straightforward.
rownames(count_mat) <- gene_id

# ---- 9) Verify matrix dimensions match expectation ----
# Why: Confirms genes are rows and samples are columns.
cat("\nCount matrix dimensions (genes x samples): ", paste(dim(count_mat), collapse = " x "), "\n")

# ---- 10) Check that values look like raw counts ----
# Why: DEA methods require raw counts (integers, wide range, many zeros).
cat("\nCount value range (min to max):\n")
print(range(count_mat, na.rm = TRUE))

cat("\nSummary of counts (overall):\n")
print(summary(as.vector(count_mat)))

# ---- 11) Verify counts are (mostly) integers ----
# Why: RNA-seq raw counts should be whole numbers; if not, it might be normalized data.
# Note: stored as numeric can still be whole numbers; we check closeness to integers.
int_check <- mean(abs(count_mat - round(count_mat)) < 1e-8, na.rm = TRUE)
cat("\nProportion of entries that are effectively integers:", round(int_check, 4), "\n")

# ---- 12) Identify sample columns by name pattern (CHP/IPF/Control) ----
# Why: This helps build an initial metadata table when sample IDs encode group labels.
sample_ids <- colnames(count_mat)

# These patterns are based on common naming; we will confirm by printing counts.
is_chp  <- grepl("^chp_",  sample_ids, ignore.case = TRUE)
is_ipf  <- grepl("^ipf_",  sample_ids, ignore.case = TRUE)
is_ctrl <- grepl("^(ctl_|ctrl_|control_)", sample_ids, ignore.case = TRUE) | grepl("^con_", sample_ids, ignore.case = TRUE)

cat("\nSample ID pattern counts:\n")
cat("CHP-like columns:",  sum(is_chp),  "\n")
cat("IPF-like columns:",  sum(is_ipf),  "\n")
cat("CTRL-like columns:", sum(is_ctrl), "\n")

# ---- 13) Build a first-pass metadata table (sample-level) ----
# Why: Differential expression requires a sample metadata frame (colData) aligned to the count matrix.
condition <- rep(NA_character_, length(sample_ids))
condition[is_chp]  <- "CHP"
condition[is_ipf]  <- "IPF"
condition[is_ctrl] <- "Control"

meta <- data.frame(
  sample_id = sample_ids,
  condition = factor(condition, levels = c("Control", "CHP", "IPF")),
  stringsAsFactors = FALSE
)

# ---- 14) Report how many samples we successfully labeled ----
# Why: If many samples are unlabeled, we must obtain proper metadata from GEO (GSM annotations).
cat("\nCondition label distribution (including NAs):\n")
print(table(meta$condition, useNA = "ifany"))

# ---- 15) Keep only samples with known condition labels (temporary, for early QC) ----
# Why: PCA/QC and DE require known groups. Unlabeled samples can be added back after proper metadata ingestion.
keep <- !is.na(meta$condition)
count_mat_labeled <- count_mat[, keep, drop = FALSE]
meta_labeled <- meta[keep, , drop = FALSE]

# Ensure meta rows align to count matrix columns (critical alignment check)
stopifnot(identical(meta_labeled$sample_id, colnames(count_mat_labeled)))

cat("\nAfter keeping labeled samples:\n")
cat("Counts matrix (genes x labeled samples): ", paste(dim(count_mat_labeled), collapse = " x "), "\n")
cat("Metadata rows:", nrow(meta_labeled), "\n")
cat("Final labeled distribution:\n")
print(table(meta_labeled$condition))

# ---- 16) Remove genes that are all zeros across labeled samples ----
# Why: Genes with zero counts everywhere provide no information and can break some steps.
nonzero_genes <- rowSums(count_mat_labeled, na.rm = TRUE) > 0
count_mat_labeled <- count_mat_labeled[nonzero_genes, , drop = FALSE]

cat("\nAfter removing all-zero genes:\n")
cat("Counts matrix dimensions (genes x samples): ", paste(dim(count_mat_labeled), collapse = " x "), "\n")

# ---- 17) Save ingestion outputs for reproducibility ----
# Why: Lock in a clean, aligned starting point so later steps are consistent.
saveRDS(count_mat,          file = "GSE150910_counts_all_samples_raw.rds")
saveRDS(meta,               file = "GSE150910_metadata_firstpass_all_samples.rds")
saveRDS(count_mat_labeled,  file = "GSE150910_counts_labeled_nonzero.rds")
saveRDS(meta_labeled,       file = "GSE150910_metadata_labeled_firstpass.rds")

cat("\nSaved RDS files:\n")
cat("- GSE150910_counts_all_samples_raw.rds\n")
cat("- GSE150910_metadata_firstpass_all_samples.rds\n")
cat("- GSE150910_counts_labeled_nonzero.rds\n")
cat("- GSE150910_metadata_labeled_firstpass.rds\n")

# ============================================================
# End of Step 2: Data ingestion + basic verification complete.
# Next step (Step 3): low-count gene filtering + library size QC
# ============================================================


# ============================================================
# Project 4: GSE150910 (CHP vs IPF vs Control)
# Step 4 — Exploratory QC (Computational) using DESeq2 VST
#   - VST (variance-stabilizing transformation)
#   - PCA (unsupervised clustering / dominant variation)
#   - Hierarchical clustering (outlier detection)
#   - Decision checks (do PCs reflect biology vs technical factors?)
# ============================================================

# ---- 0) Load required packages ----
# Why: DESeq2 provides normalization + VST; ggplot2 makes a clean PCA plot.
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2", update = FALSE, ask = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(DESeq2)
library(ggplot2)

# ---- 1) Load the clean ingestion outputs from Step 2 ----
# Why: We start from saved, aligned objects to ensure reproducibility.
count_mat <- readRDS("GSE150910_counts_labeled_nonzero.rds")      # genes x samples
meta      <- readRDS("GSE150910_metadata_labeled_firstpass.rds") # samples x metadata

# ---- 2) Confirm alignment: metadata rows must match count matrix columns ----
# Why: If misaligned, everything downstream becomes wrong.
stopifnot(identical(meta$sample_id, colnames(count_mat)))

# ---- 3) Ensure condition factor with correct reference level ----
# Why: Setting "Control" as reference supports interpretable contrasts later.
meta$condition <- factor(meta$condition, levels = c("Control", "CHP", "IPF"))

# ---- 4) Build a DESeq2 dataset object ----
# Why: DESeq2 needs a structured object that connects counts + metadata + design.
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = meta,
  design    = ~ condition
)

# ---- 5) Light pre-filtering of very low-signal genes (defensive QC) ----
# Why: Improves PCA/clustering stability and reduces noise.
dds <- dds[rowSums(counts(dds)) > 10, ]

# ---- 6) Variance-stabilizing transformation (VST) ----
# Why: RNA-seq counts have mean–variance dependence; VST makes variance more uniform,
#      which is appropriate for PCA and clustering (visual QC).
# blind=TRUE: ignores group labels so QC remains unsupervised.
vsd <- vst(dds, blind = TRUE)

# ---- 7) PCA using VST data ----
# Why: PCA shows dominant variation sources and whether samples cluster by condition.
pca_df <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3, alpha = 0.9) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("GSE150910 — PCA (VST-transformed counts)") +
  theme_classic()

print(p_pca)

# ---- 8) Optional: Save PCA figure (recommended for documentation/GitHub) ----
# Why: Creates a portable artifact for your report/manuscript.
if (!dir.exists("figures")) dir.create("figures")
ggsave(filename = "figures/Step4_PCA_VST.png", plot = p_pca, width = 7, height = 5, dpi = 300)

# ---- 9) Hierarchical clustering (sample-to-sample) ----
# Why: Detects outliers and unexpected structure.
vsd_mat <- assay(vsd)             # matrix: genes x samples (VST scale)
sample_dist <- dist(t(vsd_mat))   # transpose -> samples x genes, then distance across samples
hc <- hclust(sample_dist, method = "complete")

# ---- 10) Plot dendrogram with condition labels ----
# Why: Quickly reveals if samples group by biology or if an outlier exists.
png("figures/Step4_HierarchicalClustering_VST.png", width = 1200, height = 700, res = 150)
plot(hc,
     labels = meta$condition,
     main = "GSE150910 — Hierarchical Clustering (VST; Complete Linkage)",
     xlab = "",
     sub  = "")
dev.off()

# ---- 11) Optional: Heatmap-like sample distance matrix (base R) ----
# Why: Another outlier check; blocks/structure can indicate batch effects.
sample_dist_mat <- as.matrix(sample_dist)

png("figures/Step4_SampleDistance_Heatmap.png", width = 900, height = 800, res = 150)
heatmap(sample_dist_mat,
        symm = TRUE,
        main = "GSE150910 — Sample-to-Sample Distances (VST)",
        labRow = meta$condition,
        labCol = meta$condition)
dev.off()

# ---- 12) Decision-point diagnostics (what you check before DE testing) ----
# Why: This prints interpretable checkpoints to decide whether technical factors may dominate.
cat("\n========== Step 4 QC Summary ==========\n")
cat("VST matrix dimensions (genes x samples): ", paste(dim(vsd_mat), collapse = " x "), "\n")
cat("Condition counts:\n")
print(table(meta$condition))

cat("\nPCA variance explained:\n")
cat("PC1:", percentVar[1], "%\n")
cat("PC2:", percentVar[2], "%\n")

cat("\nSaved figures to ./figures/: \n")
cat("- Step4_PCA_VST.png\n")
cat("- Step4_HierarchicalClustering_VST.png\n")
cat("- Step4_SampleDistance_Heatmap.png\n")

cat("\nDecision point:\n")
cat("- If PCA/clustering separate primarily by condition (Control vs CHP vs IPF) -> proceed.\n")
cat("- If separation is dominated by technical factors (batch/institution/library prep) -> add covariates / batch correction.\n")
cat("- If you see strong outliers -> investigate/remove before DE testing.\n")
cat("======================================\n")

# ============================================================
# Project 4: GSE150910 (CHP vs IPF vs Control)
# Step 6 — Differential Expression (DESeq2)
# Step 7 — Signatures (common vs unique)
# Step 8 — Visualizations (volcano, heatmap, venn, boxplots)
# ============================================================

# ---- 0) Packages ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs_bioc <- c("DESeq2")
for (p in pkgs_bioc) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)

pkgs_cran <- c("ggplot2", "pheatmap")
for (p in pkgs_cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

library(DESeq2)
library(ggplot2)
library(pheatmap)

# ---- 1) Load data from your Step 2 outputs ----
count_mat <- readRDS("GSE150910_counts_labeled_nonzero.rds")      # genes x samples
meta      <- readRDS("GSE150910_metadata_labeled_firstpass.rds") # samples x metadata

# Alignment check (critical)
stopifnot(identical(meta$sample_id, colnames(count_mat)))

# Ensure condition factor with Control as reference
meta$condition <- factor(meta$condition, levels = c("Control", "CHP", "IPF"))

# ---- 2) Choose design using available covariates ----
# NOTE: You do not currently have age/sex/race/smoking/institution in meta,
# so we run the valid model you can support now: ~ condition
design_formula <- ~ condition

# ---- 3) Build DESeq2 object ----
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = meta,
  design    = design_formula
)

# Optional: defensive gene filter (improves stability)
dds <- dds[rowSums(counts(dds)) > 10, ]

# ---- 4) Run DESeq2 ----
dds <- DESeq(dds)

# ---- 5) Define contrasts ----
# CHP vs Control
res_chp_vs_ctrl <- results(dds, contrast = c("condition", "CHP", "Control"))
# IPF vs Control
res_ipf_vs_ctrl <- results(dds, contrast = c("condition", "IPF", "Control"))
# CHP vs IPF
res_chp_vs_ipf  <- results(dds, contrast = c("condition", "CHP", "IPF"))

# Convert to data.frames and keep gene column
res_to_df <- function(res_obj) {
  df <- as.data.frame(res_obj)
  df$gene <- rownames(df)
  df
}

df_chp_ctrl <- res_to_df(res_chp_vs_ctrl)
df_ipf_ctrl <- res_to_df(res_ipf_vs_ctrl)
df_chp_ipf  <- res_to_df(res_chp_vs_ipf)

# Save raw DE results
if (!dir.exists("results")) dir.create("results")
write.csv(df_chp_ctrl, "results/DESeq2_CHP_vs_Control.csv", row.names = FALSE)
write.csv(df_ipf_ctrl, "results/DESeq2_IPF_vs_Control.csv", row.names = FALSE)
write.csv(df_chp_ipf,  "results/DESeq2_CHP_vs_IPF.csv",     row.names = FALSE)

# ---- 6) Apply paper-aligned thresholds ----
# Paper-aligned: padj < 0.05 and |log2FC| > 1
sig_filter <- function(df, padj_cut = 0.05, lfc_cut = 1) {
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
  df[df$padj < padj_cut & abs(df$log2FoldChange) > lfc_cut, ]
}

sig_chp_ctrl <- sig_filter(df_chp_ctrl)
sig_ipf_ctrl <- sig_filter(df_ipf_ctrl)
sig_chp_ipf  <- sig_filter(df_chp_ipf)

write.csv(sig_chp_ctrl, "results/SIG_CHP_vs_Control_padj0.05_lfc1.csv", row.names = FALSE)
write.csv(sig_ipf_ctrl, "results/SIG_IPF_vs_Control_padj0.05_lfc1.csv", row.names = FALSE)
write.csv(sig_chp_ipf,  "results/SIG_CHP_vs_IPF_padj0.05_lfc1.csv",     row.names = FALSE)

cat("\nSignificant genes (padj<0.05 & |log2FC|>1):\n")
cat("CHP vs Control:", nrow(sig_chp_ctrl), "\n")
cat("IPF vs Control:", nrow(sig_ipf_ctrl), "\n")
cat("CHP vs IPF:",    nrow(sig_chp_ipf),  "\n")

# ---- 7) Step 7 — Define signatures (common vs unique; direction-consistent) ----
# CHP signature and IPF signature come from vs Control contrasts
chp_sig <- sig_chp_ctrl[, c("gene", "log2FoldChange", "padj")]
ipf_sig <- sig_ipf_ctrl[, c("gene", "log2FoldChange", "padj")]

# Merge by gene to find overlap
merged_sig <- merge(chp_sig, ipf_sig, by = "gene", suffixes = c("_CHP", "_IPF"))

# Direction consistency: same sign of log2FC in both contrasts
same_direction <- merged_sig[sign(merged_sig$log2FoldChange_CHP) == sign(merged_sig$log2FoldChange_IPF), ]

shared_genes <- same_direction$gene

chp_unique <- setdiff(chp_sig$gene, ipf_sig$gene)
ipf_unique <- setdiff(ipf_sig$gene, chp_sig$gene)

# Save signature gene lists
write.csv(data.frame(gene = shared_genes), "results/Signature_Shared_CHP_IPF_sameDirection.csv", row.names = FALSE)
write.csv(data.frame(gene = chp_unique),   "results/Signature_CHP_unique.csv",                row.names = FALSE)
write.csv(data.frame(gene = ipf_unique),   "results/Signature_IPF_unique.csv",                row.names = FALSE)

cat("\nSignature counts:\n")
cat("Shared (same direction):", length(shared_genes), "\n")
cat("CHP-unique:", length(chp_unique), "\n")
cat("IPF-unique:", length(ipf_unique), "\n")

# ---- 8) Step 8 — Visualizations ----
if (!dir.exists("figures")) dir.create("figures")

# 8A) Volcano plot function
volcano_plot <- function(df, title, out_png) {
  d <- df
  d$neglog10padj <- -log10(d$padj)
  d$neglog10padj[is.infinite(d$neglog10padj)] <- NA
  
  d$Sig <- "Not significant"
  d$Sig[!is.na(d$padj) & d$padj < 0.05 & !is.na(d$log2FoldChange) & abs(d$log2FoldChange) > 1] <- "Significant"
  
  p <- ggplot(d, aes(x = log2FoldChange, y = neglog10padj, color = Sig)) +
    geom_point(alpha = 0.7, size = 1.2) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    ggtitle(title) +
    xlab("log2 Fold Change") +
    ylab("-log10 adjusted p-value") +
    theme_classic()
  
  ggsave(out_png, p, width = 7, height = 5, dpi = 300)
}

volcano_plot(df_chp_ctrl, "Volcano: CHP vs Control (DESeq2)", "figures/Volcano_CHP_vs_Control.png")
volcano_plot(df_ipf_ctrl, "Volcano: IPF vs Control (DESeq2)", "figures/Volcano_IPF_vs_Control.png")
volcano_plot(df_chp_ipf,  "Volcano: CHP vs IPF (DESeq2)",     "figures/Volcano_CHP_vs_IPF.png")

# 8B) Heatmap of top DE genes per contrast (by padj)
# Use VST for heatmaps (recommended)
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)  # genes x samples

top_n <- function(sig_df, n = 50) {
  sig_df <- sig_df[order(sig_df$padj), ]
  head(sig_df$gene, n)
}

make_heatmap <- function(gene_list, title, out_png) {
  gene_list <- intersect(gene_list, rownames(vsd_mat))
  if (length(gene_list) < 2) {
    cat("Skipping heatmap for", title, "- not enough genes.\n")
    return(NULL)
  }
  
  mat <- vsd_mat[gene_list, , drop = FALSE]
  
  # Z-score by gene to highlight relative patterns
  mat_z <- t(scale(t(mat)))
  
  ann <- data.frame(condition = meta$condition)
  rownames(ann) <- meta$sample_id
  
  png(out_png, width = 1200, height = 900, res = 150)
  pheatmap(mat_z,
           annotation_col = ann,
           show_colnames = FALSE,
           main = title)
  dev.off()
}

make_heatmap(top_n(sig_chp_ctrl, 50), "Top DE genes: CHP vs Control", "figures/Heatmap_Top50_CHP_vs_Control.png")
make_heatmap(top_n(sig_ipf_ctrl, 50), "Top DE genes: IPF vs Control", "figures/Heatmap_Top50_IPF_vs_Control.png")
make_heatmap(top_n(sig_chp_ipf,  50), "Top DE genes: CHP vs IPF",     "figures/Heatmap_Top50_CHP_vs_IPF.png")

# 8C) Venn diagram (shared vs unique)
# Use base R to avoid extra package installs
png("figures/Venn_Signatures_CHP_IPF.png", width = 900, height = 700, res = 150)
plot.new()
title("Signatures: Shared vs Unique (CHP vs Control, IPF vs Control)")
text(0.5, 0.6, paste0("Shared (same direction): ", length(shared_genes)))
text(0.3, 0.4, paste0("CHP-unique: ", length(chp_unique)))
text(0.7, 0.4, paste0("IPF-unique: ", length(ipf_unique)))
dev.off()

# 8D) Boxplots for a few marker genes
# Pick a small set of example genes from significant lists if available
pick_markers <- function(sig_df, k = 4) {
  sig_df <- sig_df[order(sig_df$padj), ]
  head(sig_df$gene, k)
}

marker_genes <- unique(c(pick_markers(sig_chp_ctrl, 3),
                         pick_markers(sig_ipf_ctrl, 3),
                         pick_markers(sig_chp_ipf,  3)))
marker_genes <- marker_genes[!is.na(marker_genes)]
marker_genes <- intersect(marker_genes, rownames(vsd_mat))

if (length(marker_genes) > 0) {
  for (g in marker_genes) {
    dfp <- data.frame(
      expr = as.numeric(vsd_mat[g, ]),
      condition = meta$condition
    )
    
    p <- ggplot(dfp, aes(x = condition, y = expr)) +
      geom_boxplot(outlier_alpha = 0.5) +
      geom_jitter(width = 0.15, alpha = 0.6, size = 1.2) +
      ggtitle(paste0("Marker gene expression (VST): ", g)) +
      xlab("Condition") +
      ylab("VST expression") +
      theme_classic()
    
    ggsave(paste0("figures/Boxplot_", g, ".png"), p, width = 6, height = 4, dpi = 300)
  }
} else {
  cat("\nNo marker genes available for boxplots (no significant genes passed thresholds).\n")
}

cat("\nCompleted Step 6–8.\n")
cat("Outputs saved to:\n")
cat("- ./results (DE tables + signatures)\n")
cat("- ./figures (volcano, heatmaps, venn summary, marker boxplots)\n")
list.files("figures", full.names = TRUE)


# =========================
# Step 6/7 — Build "Commonly Differentially Expressed Genes" Table
# (CHP vs Control AND IPF vs Control), paper-style
# =========================
# Assumes you already ran DESeq2 and you have a DESeqDataSet object named: dds
# and that your main phenotype column is: condition with levels: Control, CHP, IPF
#
# Output:
#  - Table5_Common_DEGs_CHP_IPF_vs_Control.csv
#  - Table5_Common_DEGs_CHP_IPF_vs_Control.xlsx (optional, if openxlsx installed)
#
# Notes:
#  - Uses FDR padj < 0.05 and |log2FC| > 1 (paper-aligned thresholds you stated)
#  - "Shared" means significant in BOTH contrasts AND SAME direction
# =========================

# ---- Packages ----
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
})

# ---- User settings ----
alpha_cutoff <- 0.05
lfc_cutoff   <- 1
out_csv      <- "Table5_Common_DEGs_CHP_IPF_vs_Control.csv"

# ---- Safety checks ----
stopifnot(exists("dds"))
stopifnot("DESeqDataSet" %in% class(dds))
stopifnot("condition" %in% colnames(colData(dds)))

# Ensure reference level is Control (so contrasts match the paper wording)
dds$condition <- relevel(as.factor(dds$condition), ref = "Control")

# If DESeq hasn't been run yet, run it (safe to re-run; DESeq2 will overwrite results slot)
dds <- DESeq(dds)

# ---- Get DE results for each contrast ----
res_chp <- results(dds, contrast = c("condition", "CHP", "Control"), alpha = alpha_cutoff)
res_ipf <- results(dds, contrast = c("condition", "IPF", "Control"), alpha = alpha_cutoff)

# Convert to clean data frames with gene column
df_chp <- as.data.frame(res_chp) %>%
  rownames_to_column("Gene") %>%
  select(Gene, baseMean, log2FoldChange, pvalue, padj) %>%
  rename(
    baseMean_CHP = baseMean,
    log2FC_CHP   = log2FoldChange,
    pvalue_CHP   = pvalue,
    padj_CHP     = padj
  )

df_ipf <- as.data.frame(res_ipf) %>%
  rownames_to_column("Gene") %>%
  select(Gene, baseMean, log2FoldChange, pvalue, padj) %>%
  rename(
    baseMean_IPF = baseMean,
    log2FC_IPF   = log2FoldChange,
    pvalue_IPF   = pvalue,
    padj_IPF     = padj
  )

# ---- Merge & compute shared/common genes ----
tab_merged <- df_chp %>%
  inner_join(df_ipf, by = "Gene") %>%
  # Prefer one baseMean column to display (paper often shows one baseMean)
  # We'll compute a single baseMean as the average of the two contrast baseMeans:
  mutate(
    baseMean = (baseMean_CHP + baseMean_IPF) / 2
  ) %>%
  # Shared significance + effect size thresholds in BOTH contrasts
  mutate(
    sig_CHP = !is.na(padj_CHP) & padj_CHP < alpha_cutoff & !is.na(log2FC_CHP) & abs(log2FC_CHP) > lfc_cutoff,
    sig_IPF = !is.na(padj_IPF) & padj_IPF < alpha_cutoff & !is.na(log2FC_IPF) & abs(log2FC_IPF) > lfc_cutoff,
    same_direction = sign(log2FC_CHP) == sign(log2FC_IPF),
    Direction = ifelse(log2FC_CHP > 0 & log2FC_IPF > 0, "Up",
                       ifelse(log2FC_CHP < 0 & log2FC_IPF < 0, "Down", "Mixed"))
  ) %>%
  mutate(
    Signature_Type = case_when(
      sig_CHP & sig_IPF & same_direction ~ "Shared",
      sig_CHP & !(sig_IPF)              ~ "CHP-unique",
      sig_IPF & !(sig_CHP)              ~ "IPF-unique",
      TRUE                              ~ "Not significant"
    )
  )

# ---- Paper-style "Table 5": only SHARED genes (significant in both & same direction) ----
table5 <- tab_merged %>%
  filter(Signature_Type == "Shared") %>%
  # Keep the classic columns (like the screenshot): Gene, baseMean, log2FC, pvalue, padj
  # but since table is "common in both", we include both log2FC/padj side-by-side:
  transmute(
    Gene,
    baseMean = round(baseMean, 4),
    log2FC_CHP = round(log2FC_CHP, 4),
    padj_CHP   = signif(padj_CHP, 3),
    log2FC_IPF = round(log2FC_IPF, 4),
    padj_IPF   = signif(padj_IPF, 3),
    Direction
  ) %>%
  # Order like a paper: strongest evidence first (lowest max FDR), then largest effects
  mutate(max_padj = pmax(padj_CHP, padj_IPF, na.rm = TRUE),
         mean_abs_lfc = (abs(log2FC_CHP) + abs(log2FC_IPF)) / 2) %>%
  arrange(max_padj, desc(mean_abs_lfc)) %>%
  select(-max_padj, -mean_abs_lfc)

# ---- Save outputs ----
write.csv(table5, out_csv, row.names = FALSE)
message("Saved: ", out_csv)

# Optional: save as Excel (if you want a paper-like supplementary table)
if (requireNamespace("openxlsx", quietly = TRUE)) {
  out_xlsx <- sub("\\.csv$", ".xlsx", out_csv)
  openxlsx::write.xlsx(table5, out_xlsx, overwrite = TRUE)
  message("Saved: ", out_xlsx)
} else {
  message("Tip: install.packages('openxlsx') if you want an .xlsx output too.")
}

# ---- Print a preview in console ----
print(head(table5, 20))
cat("\nRows in Table 5 (Shared genes): ", nrow(table5), "\n")


