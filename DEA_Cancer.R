# ============================================================================
# UNIVERSAL LIMMA DEA SCRIPT FOR ANY CANCER TYPE
# ============================================================================
# Change ONLY the CANCER_TYPE variable to run for different cancers

# === CONFIGURATION - CHANGE THIS LINE ONLY ===
CANCER_TYPE <- "PRAD"  # OPTIONS: "PAAD", "STAD", "PRAD"

# ============================================================================
# AUTOMATIC PATH CONFIGURATION
# ============================================================================

# Base directory (adjust if needed)
base_dir <- "/Users/dyuth/Desktop/Data Science Project"

# Generate cancer-specific file paths
expr_path <- file.path(base_dir, "Data", paste0(CANCER_TYPE, "_Merged_Filtered_ProteinCoding_MIR100HG_FIXED.csv"))
dea_dir <- file.path(base_dir, "dea")
output_path <- file.path(dea_dir, paste0("limma_DEA_results_full_", CANCER_TYPE, ".csv"))

# Create DEA directory if it doesn't exist
if (!dir.exists(dea_dir)) {
  dir.create(dea_dir, recursive = TRUE)
}

cat("=================================================================\n")
cat("UNIVERSAL DEA PIPELINE FOR:", CANCER_TYPE, "\n")
cat("=================================================================\n")
cat("Input file:", expr_path, "\n")
cat("Output directory:", dea_dir, "\n")

# Check if input file exists
if (!file.exists(expr_path)) {
  stop("ERROR: Input file not found: ", expr_path, 
       "\nPlease ensure you have run the Python preprocessing pipeline first.")
}

# ============================================================================
# LOAD LIBRARIES
# ============================================================================

suppressPackageStartupMessages({
  library(limma)
})

# ============================================================================
# LOAD AND PREPARE DATA
# ============================================================================

cat("\nLoading expression data for", CANCER_TYPE, "...\n")

# Load expression + metadata
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

cat("Total samples in dataset:", nrow(expr_df), "\n")
cat("Samples with Q1/Q4 stratification:", nrow(meta_sub), "\n")

# === Prepare expression matrix (exclude metadata columns) ===
expr_only <- expr_df[, !(colnames(expr_df) %in% meta_cols), drop = FALSE]
expr_only$sample <- expr_df$sample
expr_only$sample <- gsub("-", "_", expr_only$sample, fixed = TRUE)
rownames(expr_only) <- expr_only$sample
expr_only$sample <- NULL

# Find matching samples
common <- intersect(rownames(expr_only), meta_sub$sample)
if (length(common) < 4) {
  stop("ERROR: Not enough overlapping samples between expression and metadata (need ≥4, found ", 
       length(common), ")")
}

# Subset data
expr_sub <- expr_only[common, , drop = FALSE]
meta_sub <- meta_sub[match(common, meta_sub$sample), , drop = FALSE]
rownames(meta_sub) <- meta_sub$sample

cat("Samples used for analysis:", nrow(expr_sub), " | Raw gene columns:", ncol(expr_sub), "\n")

# ============================================================================
# CLEAN EXPRESSION DATA
# ============================================================================

cat("Cleaning expression data...\n")

# Function to clean numeric data
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

cat("After filtering:", ncol(expr_num), "genes remain\n")

# ============================================================================
# ENSURE SAMPLE ALIGNMENT
# ============================================================================

# Sort both datasets by sample names to ensure proper alignment
sample_order <- sort(common)
expr_num <- expr_num[sample_order, , drop = FALSE]
meta_sub <- meta_sub[sample_order, , drop = FALSE]

# Double-check alignment
if (!identical(rownames(expr_num), rownames(meta_sub))) {
  stop("ERROR: Sample order mismatch between expression and metadata!")
}

cat("Sample alignment verified\n")

# ============================================================================
# DIAGNOSTIC INFORMATION
# ============================================================================

cat("\nDiagnostic Information for", CANCER_TYPE, ":\n")
cat("  Expression matrix:", nrow(expr_num), "samples x", ncol(expr_num), "genes\n")
cat("  Metadata:", nrow(meta_sub), "samples\n")

group_table <- table(meta_sub$MIR100HG_group)
cat("  Group counts:", paste(names(group_table), "=", group_table, collapse=", "), "\n")

# ============================================================================
# DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================================

cat("\nRunning differential expression analysis...\n")

# Prepare data for limma
expr_mat <- as.matrix(expr_num)
storage.mode(expr_mat) <- "double"

# Transpose for limma (genes as rows, samples as columns)
eset <- t(expr_mat)

# Create group factor
group <- factor(meta_sub$MIR100HG_group, levels = c("Low (Q1)", "High (Q4)"))

# Remove any NA values from group and corresponding samples
valid_samples <- !is.na(group)
if (sum(valid_samples) != length(group)) {
  cat("WARNING: Found NA values in group factor, removing affected samples\n")
  group <- group[valid_samples]
  eset <- eset[, valid_samples, drop = FALSE]
  meta_sub <- meta_sub[valid_samples, , drop = FALSE]
  cat("  After removing NAs:", length(group), "samples remain\n")
}

# Create design matrix
design <- model.matrix(~ group)
if (nrow(design) == nrow(meta_sub)) {
  rownames(design) <- rownames(meta_sub)
}

# Final dimension check
cat("Final dimensions:\n")
cat("  Expression set (eset):", nrow(eset), "genes x", ncol(eset), "samples\n")
cat("  Design matrix:", nrow(design), "samples x", ncol(design), "coefficients\n")

# Verify dimensions match
if (ncol(eset) != nrow(design)) {
  stop("ERROR: Dimension mismatch: ncol(eset) = ", ncol(eset), ", nrow(design) = ", nrow(design))
}

# Run limma
fit <- lmFit(eset, design)
fit <- eBayes(fit)

cat("Limma analysis completed successfully for", CANCER_TYPE, "!\n")

# ============================================================================
# COLLECT AND PREPARE RESULTS
# ============================================================================

cat("\nProcessing results...\n")

# Get results
tab_bh <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH", sort.by = "P")
raw_p <- fit$p.value[, 2]
bonf <- p.adjust(raw_p, method = "bonferroni")

# Prepare final results table
tab <- tab_bh
tab$P.Value.Raw <- raw_p[rownames(tab)]
tab$adj.P.Val.Bonf <- bonf[rownames(tab)]
tab$adj.P.Val.BH <- tab$adj.P.Val
tab$adj.P.Val <- NULL
tab$SYMBOL <- rownames(tab)

# Add significance flags
for (alpha in c(0.05, 0.01, 0.001, 0.0001, 0.00001)) {
  tab[[paste0("BH_sig_", alpha)]]   <- tab$adj.P.Val.BH   <= alpha
  tab[[paste0("Bonf_sig_", alpha)]] <- tab$adj.P.Val.Bonf <= alpha
}

# Save main results
write.csv(tab, output_path, row.names = TRUE)
cat("Saved main results to:", output_path, "\n")

# ============================================================================
# MULTIPLE FILTERING APPROACHES
# ============================================================================

cat("\n", rep("=", 60), "\n", sep="")
cat("LOG FOLD CHANGE FILTERING FOR", CANCER_TYPE, "\n")
cat(rep("=", 60), "\n", sep="")

# === APPROACH 1: MATCH NORMAL DATASET (|logFC| > 2 + Bonf p < 0.001) ===
cat("\nAPPROACH 1: EXACT MATCH TO NORMAL DATASET\n")
p_thresh <- 0.001
lfc_thresh <- 2.0

signif_consistent <- tab[tab$adj.P.Val.Bonf <= p_thresh & abs(tab$logFC) > lfc_thresh, ]
consistent_n <- nrow(signif_consistent)

cat(sprintf("Consistent filtering (Bonf p < %.3f AND |logFC| > %.1f): %d genes\n", 
            p_thresh, lfc_thresh, consistent_n))

consistent_filename <- file.path(dea_dir, sprintf("limma_DEA_%s_Bonf%.3f_logFC%.1f_%dgenes.csv", 
                                                  CANCER_TYPE, p_thresh, lfc_thresh, consistent_n))
write.csv(signif_consistent, consistent_filename, row.names = TRUE)
cat("Saved consistent gene set to:", consistent_filename, "\n")

# === APPROACH 2: EXPLORE OTHER logFC THRESHOLDS ===
cat("\nAPPROACH 2: EXPLORE OTHER logFC THRESHOLDS\n")

p_thresh <- 0.001
logfc_thresholds <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)

cat("Gene counts with different |logFC| thresholds (Bonferroni p < 0.001):\n")
for (lfc_thresh in logfc_thresholds) {
  n_genes <- sum(tab$adj.P.Val.Bonf <= p_thresh & abs(tab$logFC) > lfc_thresh, na.rm = TRUE)
  cat(sprintf("  |logFC| > %.1f: %d genes\n", lfc_thresh, n_genes))
  
  # Save each threshold for comparison
  if (n_genes > 0) {
    signif_lfc <- tab[tab$adj.P.Val.Bonf <= p_thresh & abs(tab$logFC) > lfc_thresh, ]
    lfc_filename <- file.path(dea_dir, sprintf("limma_DEA_%s_Bonf%.3f_logFC%.1f_%dgenes.csv", 
                                               CANCER_TYPE, p_thresh, lfc_thresh, n_genes))
    write.csv(signif_lfc, lfc_filename, row.names = TRUE)
  }
}

# === APPROACH 3: TRADITIONAL BONFERRONI ONLY ===
cat("\nAPPROACH 3: TRADITIONAL BONFERRONI ONLY\n")
signif_bonf_001 <- tab[tab$adj.P.Val.Bonf <= 0.001, ]
traditional_n <- nrow(signif_bonf_001)
cat(sprintf("Traditional Bonferroni p < 0.001 (no logFC filter): %d genes\n", traditional_n))

traditional_filename <- file.path(dea_dir, paste0("limma_DEA_", CANCER_TYPE, "_Bonf001_traditional.csv"))
write.csv(signif_bonf_001, traditional_filename, row.names = TRUE)

# === APPROACH 4: OPTIMIZE FOR TARGET GENE COUNT ===
cat("\nAPPROACH 4: OPTIMIZE FOR SIMILAR GENE COUNT\n")

# Target number of genes (adjust based on your normal dataset)
target_n <- 1500  # ADJUST THIS BASED ON YOUR NORMAL DATASET RESULTS

best_lfc_target <- NULL
best_diff_target <- Inf
lfc_fine_thresholds <- seq(0.1, 3.0, by = 0.1)

for (lfc_thresh in lfc_fine_thresholds) {
  n_genes <- sum(tab$adj.P.Val.Bonf <= 0.001 & abs(tab$logFC) > lfc_thresh, na.rm = TRUE)
  diff <- abs(n_genes - target_n)
  
  if (diff < best_diff_target && n_genes > 0) {
    best_diff_target <- diff
    best_lfc_target <- list(lfc_thresh = lfc_thresh, n_genes = n_genes, diff = diff)
  }
}

if (!is.null(best_lfc_target)) {
  cat(sprintf("Optimal |logFC| > %.1f gives %d genes (target: %d, diff: %d)\n",
              best_lfc_target$lfc_thresh, best_lfc_target$n_genes, target_n, best_lfc_target$diff))
  
  signif_optimal <- tab[tab$adj.P.Val.Bonf <= 0.001 & abs(tab$logFC) > best_lfc_target$lfc_thresh, ]
  optimal_filename <- file.path(dea_dir, sprintf("limma_DEA_%s_Bonf001_logFC%.1f_%dgenes_OPTIMAL.csv",
                                                 CANCER_TYPE, best_lfc_target$lfc_thresh, best_lfc_target$n_genes))
  write.csv(signif_optimal, optimal_filename, row.names = TRUE)
  cat("Saved optimal gene set to:", optimal_filename, "\n")
}

# ============================================================================
# SUMMARY AND RECOMMENDATIONS
# ============================================================================

cat("\n", rep("=", 60), "\n", sep="")
cat("SUMMARY FOR", CANCER_TYPE, "DATASET\n")
cat(rep("=", 60), "\n", sep="")

cat("Significant genes (counts):\n")
for (alpha in c(0.05, 0.01, 0.001, 0.0001, 0.00001)) {
  cat(sprintf("  BH  <= %-8s : %d\n", alpha, sum(tab[[paste0("BH_sig_", alpha)]], na.rm = TRUE)))
  cat(sprintf("  Bonf<= %-8s : %d\n", alpha, sum(tab[[paste0("Bonf_sig_", alpha)]], na.rm = TRUE)))
}

cat("\nFILTERED GENE SETS CREATED:\n")
cat(sprintf("1. Consistent with normal dataset (|logFC| > 2.0): %d genes\n", consistent_n))
cat(sprintf("2. Traditional Bonferroni only: %d genes\n", traditional_n))
if (!is.null(best_lfc_target)) {
  cat(sprintf("3. Optimized for target count (|logFC| > %.1f): %d genes\n", 
              best_lfc_target$lfc_thresh, best_lfc_target$n_genes))
}

cat("\nRECOMMENDATIONS FOR", CANCER_TYPE, ":\n")
cat("• For pipeline consistency: Use", basename(consistent_filename), "\n")
cat("• For maximum genes: Use", basename(traditional_filename), "\n")
if (!is.null(best_lfc_target)) {
  cat("• For balanced approach: Use", basename(optimal_filename), "\n")
}

cat("\nIMPORTANT: Choose ONE approach and use it consistently across all datasets!\n")
cat("Using |logFC| > 2.0 + Bonf p < 0.001 matches your normal dataset methodology.\n")

# ============================================================================
# PREVIEW RESULTS
# ============================================================================

cat("\nTop 10 results preview:\n")
print(head(tab[, c("SYMBOL", "logFC", "P.Value.Raw", "adj.P.Val.BH", "adj.P.Val.Bonf")], 10))

cat("\nPreview of consistent filtering results (|logFC| > 2.0):\n")
if (nrow(signif_consistent) > 0) {
  print(head(signif_consistent[, c("SYMBOL", "logFC", "P.Value.Raw", "adj.P.Val.BH", "adj.P.Val.Bonf")], 10))
} else {
  cat("No genes meet the |logFC| > 2.0 + Bonf p < 0.001 criteria for", CANCER_TYPE, "\n")
}

cat("\n", rep("=", 60), "\n", sep="")
cat("DEA ANALYSIS COMPLETE FOR", CANCER_TYPE, "!\n")
cat("Check the", dea_dir, "directory for all output files.\n")
cat(rep("=", 60), "\n", sep="")