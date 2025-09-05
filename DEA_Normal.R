# --- Limma DEA for NORMAL - FIXED + 5K Gene Selection ---

suppressPackageStartupMessages({
  library(limma)
})

# === Paths ===
expr_path <- "/Users/dyuth/Desktop/Data Science Project/Data/GTEX_Merged_ExprONLY_Filtered_ProteinCoding_MIR100HG_FINAL.csv"
output_path <- "/Users/dyuth/Desktop/Data Science Project/dea/limma_DEA_results_full.csv"

# === Load expression + metadata ===
expr_df <- read.csv(expr_path, check.names = FALSE)
stopifnot("sample" %in% colnames(expr_df))

# Make sample IDs unique before setting rownames (if needed)
expr_df$sample <- make.unique(as.character(expr_df$sample))

# Extract metadata columns by name pattern (case-insensitive match)
meta_cols <- grep("^sample$|^MIR100HG(_z)?$|^MIR100HG_group$", colnames(expr_df), ignore.case = TRUE, value = TRUE)
meta <- expr_df[, meta_cols, drop = FALSE]
meta$sample <- gsub("-", "_", meta$sample, fixed = TRUE)

# Subset metadata to only Q1 and Q4
meta_sub <- subset(meta, tolower(MIR100HG_group) %in% c("low (q1)", "high (q4)"))
meta_sub$sample <- gsub("-", "_", meta_sub$sample, fixed = TRUE)

# === Prepare expression matrix (exclude metadata columns) ===
expr_only <- expr_df[, !(colnames(expr_df) %in% meta_cols), drop = FALSE]
expr_only$sample <- expr_df$sample
expr_only$sample <- gsub("-", "_", expr_only$sample, fixed = TRUE)
rownames(expr_only) <- expr_only$sample
expr_only$sample <- NULL

# Find matching samples
common <- intersect(rownames(expr_only), meta_sub$sample)
if (length(common) < 4) stop("❌ Not enough overlapping samples between expression and metadata.")

# Subset data
expr_sub <- expr_only[common, , drop = FALSE]
meta_sub <- meta_sub[match(common, meta_sub$sample), , drop = FALSE]
rownames(meta_sub) <- meta_sub$sample

cat("✅ Samples:", nrow(expr_sub), " | Raw gene columns:", ncol(expr_sub), "\n")

# === Clean expression ===
clean_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- gsub(",", "", x, fixed = TRUE)
  x <- gsub("\\s+", "", x)
  suppressWarnings(as.numeric(x))
}

expr_sub <- as.data.frame(expr_sub, check.names = FALSE)
expr_num <- as.data.frame(lapply(expr_sub, clean_numeric), stringsAsFactors = FALSE)
rownames(expr_num) <- rownames(expr_sub)

# Drop low-quality genes
finite_counts <- vapply(expr_num, function(v) sum(is.finite(v)), integer(1))
keep_enough <- finite_counts >= pmax(3, ceiling(0.5 * nrow(expr_num)))
expr_num <- expr_num[, keep_enough, drop = FALSE]

# Drop zero-variance genes
sds <- apply(expr_num, 2, sd, na.rm = TRUE)
keep_var <- is.finite(sds) & sds > 0
expr_num <- expr_num[, keep_var, drop = FALSE]

# Ensure unique gene IDs
if (anyDuplicated(colnames(expr_num))) {
  colnames(expr_num) <- make.unique(colnames(expr_num))
}

cat("✅ After filtering: ", ncol(expr_num), " genes remain\n")

# === CRITICAL FIX: Ensure sample order matches between expression and metadata ===
# Sort both datasets by sample names to ensure proper alignment
sample_order <- sort(common)
expr_num <- expr_num[sample_order, , drop = FALSE]
meta_sub <- meta_sub[sample_order, , drop = FALSE]

# Double-check alignment
if (!identical(rownames(expr_num), rownames(meta_sub))) {
  stop("❌ Sample order mismatch between expression and metadata!")
}

cat("✅ Sample alignment verified\n")

# === Debugging: Check dimensions before DEA ===
cat("🔍 Debug info:\n")
cat("  Expression matrix: ", nrow(expr_num), "samples x", ncol(expr_num), "genes\n")
cat("  Metadata: ", nrow(meta_sub), "samples\n")
cat("  Group levels: ", levels(factor(meta_sub$MIR100HG_group)), "\n")
cat("  Group counts: ", table(meta_sub$MIR100HG_group), "\n")

# === DEA (limma) - FIXED ===
expr_mat <- as.matrix(expr_num)
storage.mode(expr_mat) <- "double"

# Transpose for limma (genes as rows, samples as columns)
eset <- t(expr_mat)

# Create design matrix and check for issues
# Fix: Use the actual case-sensitive group names from your data
group <- factor(meta_sub$MIR100HG_group, levels = c("low (q1)", "high (q4)"))

# Debug the group factor
cat("🔍 Group factor debug:\n")
cat("  Length of group factor: ", length(group), "\n")
cat("  Number of samples in meta_sub: ", nrow(meta_sub), "\n")
cat("  Group factor levels: ", levels(group), "\n")
cat("  First few group values: ", head(group), "\n")

# Remove any NA values from group and corresponding samples
valid_samples <- !is.na(group)
if (sum(valid_samples) != length(group)) {
  cat("⚠️  Found NA values in group factor, removing affected samples\n")
  group <- group[valid_samples]
  eset <- eset[, valid_samples, drop = FALSE]
  meta_sub <- meta_sub[valid_samples, , drop = FALSE]
  cat("  After removing NAs: ", length(group), "samples remain\n")
}

design <- model.matrix(~ group)

# Only assign rownames if dimensions match
if (nrow(design) == nrow(meta_sub)) {
  rownames(design) <- rownames(meta_sub)
} else {
  cat("⚠️  Design matrix rows (", nrow(design), ") != metadata rows (", nrow(meta_sub), ")\n")
  cat("  Using default rownames for design matrix\n")
}

# Final dimension check
cat("🔍 Final dimensions:\n")
cat("  Expression set (eset): ", nrow(eset), "genes x", ncol(eset), "samples\n")
cat("  Design matrix: ", nrow(design), "samples x", ncol(design), "coefficients\n")

# Verify dimensions match
if (ncol(eset) != nrow(design)) {
  stop("❌ Dimension mismatch: ncol(eset) = ", ncol(eset), ", nrow(design) = ", nrow(design))
}

# Run limma
fit <- lmFit(eset, design)
fit <- eBayes(fit)

cat("✅ limma analysis completed successfully!\n")

# === Collect results ===
tab_bh <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH", sort.by = "P")
raw_p <- fit$p.value[, 2]
bonf <- p.adjust(raw_p, method = "bonferroni")

tab <- tab_bh
tab$P.Value.Raw <- raw_p[rownames(tab)]
tab$adj.P.Val.Bonf <- bonf[rownames(tab)]
tab$adj.P.Val.BH <- tab$adj.P.Val
tab$adj.P.Val <- NULL
tab$SYMBOL <- rownames(tab)

# Significance flags
for (alpha in c(0.05, 0.01, 0.001, 0.0001, 0.00001)) {
  tab[[paste0("BH_sig_", alpha)]]   <- tab$adj.P.Val.BH   <= alpha
  tab[[paste0("Bonf_sig_", alpha)]] <- tab$adj.P.Val.Bonf <= alpha
}

# === Save main results ===
write.csv(tab, output_path, row.names = TRUE)
cat("\n✅ Saved main results to:", output_path, "\n")

# === Summary ===
cat("\n📊 Significant genes (counts):\n")
for (alpha in c(0.05, 0.01, 0.001, 0.0001, 0.00001)) {
  cat(sprintf("  BH  ≤ %-8s : %d\n", alpha, sum(tab[[paste0("BH_sig_", alpha)]], na.rm = TRUE)))
  cat(sprintf("  Bonf≤ %-8s : %d\n", alpha, sum(tab[[paste0("Bonf_sig_", alpha)]], na.rm = TRUE)))
}

# ========================================================================================
# === THREE APPROACHES FOR GETTING ~5000 GENES ===
# ========================================================================================

cat("\n🎯 === APPROACH 1: OPTIMAL P-VALUE THRESHOLD ===\n")

# Find the optimal Bonferroni p-value threshold that gets closest to 5000 genes
target_n <- 5000
best_thresh <- 0.001
best_diff <- Inf

# Test a range of thresholds
thresholds <- seq(0.00001, 0.01, by = 0.00001)
for (thresh in thresholds) {
  n_genes <- sum(tab$adj.P.Val.Bonf <= thresh, na.rm = TRUE)
  diff <- abs(n_genes - target_n)
  if (diff < best_diff) {
    best_diff <- diff
    best_thresh <- thresh
  }
}

optimal_n <- sum(tab$adj.P.Val.Bonf <= best_thresh, na.rm = TRUE)
cat(sprintf("🎯 Optimal Bonferroni threshold: %.6f gives %d genes (difference: %d from target)\n", 
            best_thresh, optimal_n, abs(optimal_n - target_n)))

# Save the optimized gene set
signif_optimal <- tab[tab$adj.P.Val.Bonf <= best_thresh, ]
optimal_filename <- sprintf("limma_DEA_Bonf_%.6f_%dgenes.csv", best_thresh, optimal_n)
write.csv(signif_optimal, optimal_filename, row.names = TRUE)
cat("✅ Saved optimal threshold gene set to:", optimal_filename, "\n")

cat("\n🎯 === APPROACH 2: LOG FOLD CHANGE + P-VALUE FILTERING ===\n")

# Test different combinations of p-value and logFC thresholds
p_thresh <- 0.001
logfc_thresholds <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)

cat("Testing combinations with Bonferroni p <", p_thresh, ":\n")
best_lfc_combo <- NULL
best_lfc_diff <- Inf

for (lfc_thresh in logfc_thresholds) {
  n_genes <- sum(tab$adj.P.Val.Bonf <= p_thresh & abs(tab$logFC) > lfc_thresh, na.rm = TRUE)
  diff <- abs(n_genes - target_n)
  cat(sprintf("  |logFC| > %.1f: %d genes (diff: %d)\n", lfc_thresh, n_genes, diff))
  
  if (diff < best_lfc_diff && n_genes > 0) {
    best_lfc_diff <- diff
    best_lfc_combo <- list(p_thresh = p_thresh, lfc_thresh = lfc_thresh, n_genes = n_genes)
  }
}

# Save the best logFC combination
if (!is.null(best_lfc_combo)) {
  signif_combined <- tab[tab$adj.P.Val.Bonf <= best_lfc_combo$p_thresh & 
                           abs(tab$logFC) > best_lfc_combo$lfc_thresh, ]
  combined_filename <- sprintf("limma_DEA_Bonf%.3f_logFC%.1f_%dgenes.csv", 
                               best_lfc_combo$p_thresh, best_lfc_combo$lfc_thresh, 
                               best_lfc_combo$n_genes)
  write.csv(signif_combined, combined_filename, row.names = TRUE)
  cat(sprintf("✅ Best combination: Bonf p < %.3f AND |logFC| > %.1f (%d genes)\n", 
              best_lfc_combo$p_thresh, best_lfc_combo$lfc_thresh, best_lfc_combo$n_genes))
  cat("✅ Saved combined filter gene set to:", combined_filename, "\n")
}

cat("\n🎯 === APPROACH 3: TOP 5000 BY SIGNIFICANCE ===\n")

# Simply take the top 5000 most significant genes by Bonferroni p-value
tab_sorted <- tab[order(tab$adj.P.Val.Bonf), ]
top_5000 <- tab_sorted[1:min(5000, nrow(tab_sorted)), ]

top5k_filename <- "limma_DEA_top5000_genes.csv"
write.csv(top_5000, top5k_filename, row.names = TRUE)

cat("✅ Top 5000 genes by Bonferroni p-value saved to:", top5k_filename, "\n")
cat("📊 P-value range in top 5000: ", sprintf("%.2e", min(top_5000$adj.P.Val.Bonf)), 
    " to ", sprintf("%.2e", max(top_5000$adj.P.Val.Bonf)), "\n")

# ========================================================================================
# === SUMMARY AND RECOMMENDATIONS ===
# ========================================================================================

cat("\n📋 === SUMMARY OF 5K GENE SELECTION METHODS ===\n")
cat("1. Optimal p-value threshold:", optimal_filename, "\n")
if (!is.null(best_lfc_combo)) {
  cat("2. Combined p-value + logFC:", combined_filename, "\n")
} else {
  cat("2. Combined p-value + logFC: No suitable combination found\n")
}
cat("3. Top 5000 by significance:", top5k_filename, "\n")

cat("\n🔬 === RECOMMENDATIONS ===\n")
cat("• For statistical rigor: Use the optimal p-value threshold approach (#1)\n")
cat("• For biological relevance: Use the combined p-value + logFC approach (#2)\n")
cat("• For simplicity: Use the top 5000 approach (#3)\n")

# Also save traditional thresholds for comparison
signif_bonf_001 <- tab[tab$adj.P.Val.Bonf <= 0.001, ]
write.csv(signif_bonf_001, "limma_DEA_Bonf001_traditional.csv", row.names = TRUE)
cat("\n✅ Traditional Bonferroni p < 0.001 set saved (", nrow(signif_bonf_001), " genes)\n")

# === Preview ===
cat("\n📋 Top 10 results preview:\n")
print(head(tab[, c("SYMBOL", "logFC", "P.Value.Raw", "adj.P.Val.BH", "adj.P.Val.Bonf")], 10))