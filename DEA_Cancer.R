

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(gridExtra)
  library(RColorBrewer)
})


CANCER_TYPE <- "STAD"  


base_dir <- "/Users/dyuth/Desktop/Data Science Project"


expr_path <- file.path(base_dir, "Data", paste0(CANCER_TYPE, "_Merged_Filtered_ProteinCoding_MIR100HG_FIXED.csv"))
dea_dir <- file.path(base_dir, "dea")
validation_dir <- file.path(dea_dir, "validation")


if (!dir.exists(dea_dir)) dir.create(dea_dir, recursive = TRUE)
if (!dir.exists(validation_dir)) dir.create(validation_dir, recursive = TRUE)


output_path_complete <- file.path(dea_dir, paste0("limma_DEA_results_complete_", CANCER_TYPE, ".csv"))
output_path_filtered <- file.path(dea_dir, paste0("limma_DEA_results_filtered_", CANCER_TYPE, ".csv"))
validation_summary_path <- file.path(validation_dir, paste0("dea_validation_summary_", CANCER_TYPE, ".csv"))



validation_metrics <- list()
start_time <- Sys.time()


if (!file.exists(expr_path)) {
  stop(" Input file not found: ", expr_path, 
       "\n run the Python preprocessing pipeline first.")
}

cat("Input file:", expr_path, "\n")
cat("Output directory:", dea_dir, "\n")
cat("Complete results:", basename(output_path_complete), "\n")
cat("Filtered results:", basename(output_path_filtered), "\n")



cat("\nLoading and validating expression data for", CANCER_TYPE, "...\n")


expr_df <- read.csv(expr_path, check.names = FALSE)
validation_metrics$initial_samples <- nrow(expr_df)
validation_metrics$initial_genes <- ncol(expr_df) - 1  

cat("Initial data validation:\n")
cat("  Samples:", validation_metrics$initial_samples, "\n")
cat("  Genes:", validation_metrics$initial_genes, "\n")


missing_pct <- sum(is.na(expr_df)) / (nrow(expr_df) * ncol(expr_df)) * 100
validation_metrics$initial_missing_pct <- missing_pct
cat("  Missing data:", round(missing_pct, 2), "%\n")


stopifnot("sample" %in% colnames(expr_df))


mir_check <- grep("MIR100HG", colnames(expr_df), ignore.case = TRUE, value = TRUE)
validation_metrics$mir100hg_columns_detected <- length(mir_check)
cat("MIR100HG columns detected:", paste(mir_check, collapse = ", "), "\n")


expr_df$sample <- make.unique(as.character(expr_df$sample))




meta_cols <- grep("^sample$|^MIR100HG(_z)?$|^MIR100HG_group$", colnames(expr_df), ignore.case = TRUE, value = TRUE)
meta <- expr_df[, meta_cols, drop = FALSE]
meta$sample <- gsub("-", "_", meta$sample, fixed = TRUE)


if ("MIR100HG_group" %in% colnames(meta)) {
  group_dist <- table(meta$MIR100HG_group, useNA = "ifany")
  validation_metrics$group_low_q1 <- as.numeric(group_dist["Low (Q1)"])
  validation_metrics$group_high_q4 <- as.numeric(group_dist["High (Q4)"])
  validation_metrics$group_total <- sum(group_dist)
  cat("Group distribution:", paste(names(group_dist), "=", group_dist, collapse=", "), "\n")
} else {
  stop(" No MIR100HG_group column found for stratification")
}


meta_sub <- subset(meta, tolower(MIR100HG_group) %in% c("low (q1)", "high (q4)"))
meta_sub$sample <- gsub("-", "_", meta_sub$sample, fixed = TRUE)

validation_metrics$analysis_samples <- nrow(meta_sub)
cat("Total samples in dataset:", nrow(expr_df), "\n")
cat("Samples with Q1/Q4 stratification:", nrow(meta_sub), "\n")


expr_only <- expr_df[, !(colnames(expr_df) %in% meta_cols), drop = FALSE]
expr_only$sample <- expr_df$sample
expr_only$sample <- gsub("-", "_", expr_only$sample, fixed = TRUE)
rownames(expr_only) <- expr_only$sample
expr_only$sample <- NULL


mir_expr_cols <- grep("^MIR100HG(_x|_y|_z)?$", colnames(expr_only), ignore.case = TRUE, value = TRUE)
preserved_mir <- intersect(mir_expr_cols, colnames(expr_only))
validation_metrics$mir100hg_preserved <- length(preserved_mir)
cat("MIR100HG expression columns preserved:", paste(preserved_mir, collapse = ", "), "\n")


common <- intersect(rownames(expr_only), meta_sub$sample)
validation_metrics$matched_samples <- length(common)

if (length(common) < 4) {
  stop("Not enough overlapping samples between expression and metadata (need ≥4, found ", 
       length(common), ")")
}


expr_sub <- expr_only[common, , drop = FALSE]
meta_sub <- meta_sub[match(common, meta_sub$sample), , drop = FALSE]
rownames(meta_sub) <- meta_sub$sample

cat("Samples used for analysis:", nrow(expr_sub), " | Raw gene columns:", ncol(expr_sub), "\n")



cat("Cleaning expression data with validation\n")


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


protected_genes <- mir_expr_cols[mir_expr_cols %in% colnames(expr_num)]
validation_metrics$protected_genes <- length(protected_genes)


pre_filter_metrics <- list(
  total_genes = ncol(expr_num),
  genes_with_na = sum(sapply(expr_num, function(x) any(is.na(x)))),
  genes_all_zero = sum(sapply(expr_num, function(x) all(x == 0, na.rm = TRUE))),
  genes_no_variance = sum(sapply(expr_num, function(x) var(x, na.rm = TRUE) == 0))
)


finite_counts <- vapply(expr_num, function(v) sum(is.finite(v)), integer(1))
keep_enough <- finite_counts >= pmax(3, ceiling(0.5 * nrow(expr_num)))
keep_enough[protected_genes] <- TRUE  
expr_num <- expr_num[, keep_enough, drop = FALSE]


sds <- apply(expr_num, 2, sd, na.rm = TRUE)
keep_var <- is.finite(sds) & sds > 0
keep_var[protected_genes] <- TRUE  

expr_num <- expr_num[, keep_var, drop = FALSE]


validation_metrics$genes_after_filtering <- ncol(expr_num)
validation_metrics$filtering_retention_rate <- (ncol(expr_num) / pre_filter_metrics$total_genes) * 100

final_mir <- intersect(protected_genes, colnames(expr_num))
validation_metrics$mir100hg_final <- length(final_mir)


if (anyDuplicated(colnames(expr_num))) {
  colnames(expr_num) <- make.unique(colnames(expr_num))
}

cat("Filtering results:\n")
cat("Genes before filtering:", pre_filter_metrics$total_genes, "\n")
cat("Genes after filtering:", validation_metrics$genes_after_filtering, "\n")
cat("Retention rate:", round(validation_metrics$filtering_retention_rate, 1), "%\n")
cat("MIR100HG genes preserved:", paste(final_mir, collapse = ", "), "\n")


cat("\nEnsuring sample alignment and assessing statistical power\n")


sample_order <- sort(common)
expr_num <- expr_num[sample_order, , drop = FALSE]
meta_sub <- meta_sub[sample_order, , drop = FALSE]


if (!identical(rownames(expr_num), rownames(meta_sub))) {
  stop("Sample order mismatch between expression and metadata")
}

validation_metrics$sample_alignment_verified <- TRUE
cat("Sample alignment verified\n")


group_table <- table(meta_sub$MIR100HG_group)
min_group_size <- min(group_table)
max_group_size <- max(group_table)

validation_metrics$min_group_size <- min_group_size
validation_metrics$max_group_size <- max_group_size
validation_metrics$group_balance <- abs(max_group_size - min_group_size) / max_group_size

cat("Statistical Power Assessment:\n")
cat("Group sizes - Min:", min_group_size, "Max:", max_group_size, "\n")
cat("Group imbalance:", round(validation_metrics$group_balance * 100, 1), "%\n")


effect_sizes <- c(0.5, 1.0, 1.5, 2.0)
power_estimates <- sapply(effect_sizes, function(es) {

  pooled_n <- min_group_size
  df <- 2 * pooled_n - 2
  ncp <- es * sqrt(pooled_n / 2)
  power <- 1 - pt(qt(0.975, df), df, ncp) + pt(qt(0.025, df), df, ncp)
  return(power)
})


validation_metrics$power_effect_size_0.5 <- power_estimates[1]
validation_metrics$power_effect_size_1.0 <- power_estimates[2]
validation_metrics$power_effect_size_1.5 <- power_estimates[3]
validation_metrics$power_effect_size_2.0 <- power_estimates[4]

cat("Statistical power estimates:\n")
for (i in seq_along(effect_sizes)) {
  cat("Effect size", effect_sizes[i], ":", round(power_estimates[i] * 100, 1), "%\n")
}



cat("\nDiagnostic Information for", CANCER_TYPE, ":\n")
cat("Expression matrix:", nrow(expr_num), "samples x", ncol(expr_num), "genes\n")
cat("Metadata:", nrow(meta_sub), "samples\n")
cat("Group counts:", paste(names(group_table), "=", group_table, collapse=", "), "\n")



cat("\nRunning enhanced differential expression analysis\n")


expr_mat <- as.matrix(expr_num)
storage.mode(expr_mat) <- "double"


eset <- t(expr_mat)


group <- factor(meta_sub$MIR100HG_group, levels = c("Low (Q1)", "High (Q4)"))


valid_samples <- !is.na(group)
if (sum(valid_samples) != length(group)) {
  cat("found NA values in group factor, removing affected samples\n")
  group <- group[valid_samples]
  eset <- eset[, valid_samples, drop = FALSE]
  meta_sub <- meta_sub[valid_samples, , drop = FALSE]
  cat("  After removing NAs:", length(group), "samples remain\n")
}


design <- model.matrix(~ group)
if (nrow(design) == nrow(meta_sub)) {
  rownames(design) <- rownames(meta_sub)
}


validation_metrics$design_rank <- qr(design)$rank
validation_metrics$design_full_rank <- validation_metrics$design_rank == ncol(design)

cat("Design matrix validation:\n")
cat("Rank:", validation_metrics$design_rank, "/", ncol(design), "\n")
cat("Full rank:", validation_metrics$design_full_rank, "\n")


cat("Final dimensions:\n")
cat("Expression set (eset):", nrow(eset), "genes x", ncol(eset), "samples\n")
cat("Design matrix:", nrow(design), "samples x", ncol(design), "coefficients\n")

validation_metrics$final_samples <- ncol(eset)
validation_metrics$final_genes <- nrow(eset)


if (ncol(eset) != nrow(design)) {
  stop("dimension mismatch: ncol(eset) = ", ncol(eset), ", nrow(design) = ", nrow(design))
}


tryCatch({
  fit <- lmFit(eset, design)
  fit <- eBayes(fit)
  validation_metrics$limma_success <- TRUE
  validation_metrics$limma_error <- "None"
}, error = function(e) {
  validation_metrics$limma_success <- FALSE
  validation_metrics$limma_error <- as.character(e)
  stop("Limma analysis failed: ", e)
})




cat("\nProcessing complete results with validation\n")


tab_complete <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH", sort.by = "none")
raw_p <- fit$p.value[, 2]
bonf <- p.adjust(raw_p, method = "bonferroni")


tab_complete$P.Value.Raw <- raw_p[rownames(tab_complete)]
tab_complete$adj.P.Val.Bonf <- bonf[rownames(tab_complete)]
tab_complete$adj.P.Val.BH <- tab_complete$adj.P.Val
tab_complete$adj.P.Val <- NULL
tab_complete$SYMBOL <- rownames(tab_complete)


alpha_levels <- c(0.05, 0.01, 0.001, 0.0001, 0.00001)
for (alpha in alpha_levels) {
  tab_complete[[paste0("BH_sig_", alpha)]]   <- tab_complete$adj.P.Val.BH   <= alpha
  tab_complete[[paste0("Bonf_sig_", alpha)]] <- tab_complete$adj.P.Val.Bonf <= alpha
}



validation_metrics$total_tests <- nrow(tab_complete)
validation_metrics$significant_bh_0.05 <- sum(tab_complete$adj.P.Val.BH <= 0.05)
validation_metrics$significant_bonf_0.001 <- sum(tab_complete$adj.P.Val.Bonf <= 0.001)


validation_metrics$mean_abs_logfc <- mean(abs(tab_complete$logFC))
validation_metrics$median_abs_logfc <- median(abs(tab_complete$logFC))
validation_metrics$max_abs_logfc <- max(abs(tab_complete$logFC))


validation_metrics$p_value_dist_uniform_test <- ks.test(tab_complete$P.Value.Raw, "punif")$p.value


up_genes <- sum(tab_complete$logFC > 0)
down_genes <- sum(tab_complete$logFC < 0)
validation_metrics$up_regulated <- up_genes
validation_metrics$down_regulated <- down_genes
validation_metrics$direction_bias <- abs(up_genes - down_genes) / (up_genes + down_genes)

cat("Statistical validation results:\n")
cat("  Total genes tested:", validation_metrics$total_tests, "\n")
cat("  Significant (BH p<0.05):", validation_metrics$significant_bh_0.05, "\n")
cat("  Significant (Bonf p<0.001):", validation_metrics$significant_bonf_0.001, "\n")
cat("  Up-regulated genes:", validation_metrics$up_regulated, "\n")
cat("  Down-regulated genes:", validation_metrics$down_regulated, "\n")
cat("  Direction bias:", round(validation_metrics$direction_bias * 100, 1), "%\n")


mir_in_results <- tab_complete[grep("MIR100HG", tab_complete$SYMBOL, ignore.case = TRUE), ]
validation_metrics$mir100hg_in_results <- nrow(mir_in_results)

if (nrow(mir_in_results) > 0) {
  validation_metrics$mir100hg_logfc <- mir_in_results$logFC[1]
  validation_metrics$mir100hg_pvalue <- mir_in_results$P.Value.Raw[1]
  validation_metrics$mir100hg_bonf_pvalue <- mir_in_results$adj.P.Val.Bonf[1]
  validation_metrics$mir100hg_significant_bonf_0.001 <- mir_in_results$adj.P.Val.Bonf[1] <= 0.001
  
  cat("MIR100HG analysis results:\n")
  cat("LogFC:", round(validation_metrics$mir100hg_logfc, 4), "\n")
  cat("Raw p-value:", format(validation_metrics$mir100hg_pvalue, scientific = TRUE), "\n")
  cat("Bonferroni p-value:", format(validation_metrics$mir100hg_bonf_pvalue, scientific = TRUE), "\n")
  cat("Significant (Bonf p<0.001):", validation_metrics$mir100hg_significant_bonf_0.001, "\n")
} else {
  validation_metrics$mir100hg_logfc <- NA
  validation_metrics$mir100hg_pvalue <- NA
  validation_metrics$mir100hg_bonf_pvalue <- NA
  validation_metrics$mir100hg_significant_bonf_0.001 <- FALSE
  cat(" MIR100HG not found in results\n")
}



cat("\nGenerating plots\n")


pdf(file.path(validation_dir, paste0("dea_diagnostics_cancer_", CANCER_TYPE, ".pdf")), 
    width = 12, height = 8)

# P-value distribution
par(mfrow = c(2, 3))
hist(tab_complete$P.Value.Raw, breaks = 20, main = paste("P-value Distribution -", CANCER_TYPE), 
     xlab = "P-value", col = "lightblue", border = "black")
abline(h = nrow(tab_complete)/20, col = "red", lty = 2, lwd = 2)

# Volcano plot
plot(tab_complete$logFC, -log10(tab_complete$P.Value.Raw), 
     main = paste("Volcano Plot -", CANCER_TYPE), xlab = "Log2 Fold Change", ylab = "-log10(P-value)",
     pch = 20, col = ifelse(tab_complete$adj.P.Val.Bonf <= 0.001, "red", "gray"))
abline(v = c(-0.5, 0.5), col = "blue", lty = 2)
abline(h = -log10(0.001), col = "blue", lty = 2)

# Mean-variance relationship
plot(fit$Amean, sqrt(fit$sigma), main = "Mean-Variance Relationship",
     xlab = "Average Expression", ylab = "Standard Deviation", pch = 20)

# Effect size distribution
hist(tab_complete$logFC, breaks = 50, main = "Effect Size Distribution",
     xlab = "Log2 Fold Change", col = "lightgreen", border = "black")
abline(v = 0, col = "red", lwd = 2)

# Sample correlation heatmap
if (ncol(eset) <= 50) {  # Only if manageable number of samples
  sample_cors <- cor(eset)
  heatmap(sample_cors, main = "Sample Correlation Heatmap", 
          col = colorRampPalette(c("blue", "white", "red"))(50))
}

# Top genes heatmap
if (validation_metrics$significant_bonf_0.001 > 0 && validation_metrics$significant_bonf_0.001 <= 50) {
  top_genes <- rownames(tab_complete[tab_complete$adj.P.Val.Bonf <= 0.001, ])[1:min(20, validation_metrics$significant_bonf_0.001)]
  heatmap(eset[top_genes, ], main = "Top Significant Genes", 
          col = colorRampPalette(c("blue", "white", "red"))(50),
          scale = "row")
}

dev.off()



write.csv(tab_complete, output_path_complete, row.names = TRUE)
cat(" Saved COMPLETE results (", nrow(tab_complete), " genes) to:", basename(output_path_complete), "\n")



cat("\nCreating filtered results for pipeline compatibility\n")


p_threshold <- 0.001
logfc_threshold <- 1.5  


tab_filtered <- tab_complete[
  tab_complete$adj.P.Val.Bonf <= p_threshold & 
    abs(tab_complete$logFC) >= logfc_threshold, 
]


tab_filtered <- tab_filtered[order(tab_filtered$adj.P.Val.Bonf), ]


write.csv(tab_filtered, output_path_filtered, row.names = TRUE)
cat("Saved FILTERED results (", nrow(tab_filtered), " genes) to:", basename(output_path_filtered), "\n")


legacy_filename <- file.path(dea_dir, paste0("limma_DEA_", CANCER_TYPE, "_Bonf", p_threshold, "_logFC", logfc_threshold, "_", nrow(tab_filtered), "genes.csv"))
write.csv(tab_filtered, legacy_filename, row.names = TRUE)
cat("Saved LEGACY format to:", basename(legacy_filename), "\n")



cat("\n", rep("=", 60), "\n", sep="")
cat("ADDITIONAL FILTERING VARIATIONS FOR", CANCER_TYPE, "\n")
cat(rep("=", 60), "\n", sep="")


filtering_configs <- list(
  list(p = 0.001, lfc = 0.5, name = "relaxed"),
  list(p = 0.001, lfc = 1.0, name = "moderate"),
  list(p = 0.001, lfc = 2.0, name = "strict"),
  list(p = 0.01, lfc = 1.0, name = "less_stringent")
)

variant_results <- list()
cat("Creating additional filtered versions:\n")
for (config in filtering_configs) {
  filtered_genes <- tab_complete[
    tab_complete$adj.P.Val.Bonf <= config$p & 
      abs(tab_complete$logFC) >= config$lfc, 
  ]
  
  variant_results[[config$name]] <- nrow(filtered_genes)
  
  if (nrow(filtered_genes) > 0) {
    variant_filename <- file.path(dea_dir, paste0("limma_DEA_", CANCER_TYPE, "_", config$name, "_", nrow(filtered_genes), "genes.csv"))
    write.csv(filtered_genes, variant_filename, row.names = TRUE)
    cat(sprintf("  %s (p<%.3f, |logFC|>%.1f): %d genes -> %s\n", 
                config$name, config$p, config$lfc, nrow(filtered_genes), basename(variant_filename)))
  }
}


validation_metrics$filtered_relaxed <- variant_results$relaxed
validation_metrics$filtered_moderate <- variant_results$moderate
validation_metrics$filtered_strict <- variant_results$strict
validation_metrics$filtered_less_stringent <- variant_results$less_stringent


p001_results <- tab_complete[
  tab_complete$adj.P.Val.Bonf <= 0.001 & 
    abs(tab_complete$logFC) >= 0.5, 
]
p001_filename <- file.path(dea_dir, paste0("limma_DEA_cancer_", CANCER_TYPE, "_Bonf0.001_logFC0.5_genes.csv"))
write.csv(p001_results, p001_filename, row.names = TRUE)
cat("Created p<0.001, |logFC|>=0.5 file:", nrow(p001_results), "genes\n")

validation_metrics$p001_genes <- nrow(p001_results)




end_time <- Sys.time()
validation_metrics$runtime_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
validation_metrics$timestamp <- as.character(Sys.time())
validation_metrics$r_version <- R.version.string
validation_metrics$limma_version <- as.character(packageVersion("limma"))
validation_metrics$cancer_type <- CANCER_TYPE


validation_flags <- c()

if (validation_metrics$group_balance > 0.2) {
  validation_flags <- c(validation_flags, "HIGH_GROUP_IMBALANCE")
}

if (validation_metrics$direction_bias > 0.8) {
  validation_flags <- c(validation_flags, "EXTREME_DIRECTION_BIAS")
}

if (validation_metrics$p_value_dist_uniform_test < 0.01) {
  validation_flags <- c(validation_flags, "NON_UNIFORM_P_DISTRIBUTION")
}

if (validation_metrics$mir100hg_in_results == 0) {
  validation_flags <- c(validation_flags, "MIR100HG_MISSING")
}

if (validation_metrics$min_group_size < 10) {
  validation_flags <- c(validation_flags, "SMALL_SAMPLE_SIZE")
}

if (validation_metrics$significant_bonf_0.001 == 0) {
  validation_flags <- c(validation_flags, "NO_SIGNIFICANT_GENES")
}

validation_metrics$validation_flags <- if(length(validation_flags) > 0) paste(validation_flags, collapse = ";") else "None"
validation_metrics$validation_passed <- length(validation_flags) == 0


safe_convert <- function(x) {
  if (is.null(x)) return("NULL")
  if (is.logical(x)) return(as.character(x))
  if (is.numeric(x)) {
    if (length(x) == 1) return(as.character(x))
    else return(paste(x, collapse = ","))
  }
  return(as.character(x))
}

validation_df <- data.frame(
  metric = names(validation_metrics),
  value = sapply(validation_metrics, safe_convert, USE.NAMES = FALSE),
  stringsAsFactors = FALSE
)

write.csv(validation_df, validation_summary_path, row.names = FALSE)



cat("Complete dataset statistics:\n")
cat("Total genes tested:", nrow(tab_complete), "\n")
cat("Genes with p < 0.05:", sum(tab_complete$P.Value.Raw < 0.05, na.rm = TRUE), "\n")
cat("Genes with Bonf p < 0.001:", sum(tab_complete$adj.P.Val.Bonf < 0.001, na.rm = TRUE), "\n")
cat("Genes with |logFC| > 1:", sum(abs(tab_complete$logFC) > 1, na.rm = TRUE), "\n")
cat("Genes with |logFC| > 2:", sum(abs(tab_complete$logFC) > 2, na.rm = TRUE), "\n")

cat("\nFiltered dataset:\n")
cat("  Significant genes:", nrow(tab_filtered), "\n")
if (nrow(tab_filtered) > 0) {
  cat("Upregulated (logFC > 0):", sum(tab_filtered$logFC > 0, na.rm = TRUE), "\n")
  cat("Downregulated (logFC < 0):", sum(tab_filtered$logFC < 0, na.rm = TRUE), "\n")
}

cat("\nData Quality Metrics:\n")
cat("Initial retention rate:", round(validation_metrics$filtering_retention_rate, 1), "%\n")
cat("Group balance score:", round((1 - validation_metrics$group_balance) * 100, 1), "%\n")
cat("Direction bias:", round(validation_metrics$direction_bias * 100, 1), "%\n")


cat("Analysis Summary:\n")
cat("Cancer Type:", CANCER_TYPE, "\n")
cat("Runtime:", round(validation_metrics$runtime_minutes, 2), "minutes\n")
cat("Timestamp:", validation_metrics$timestamp, "\n")

cat("\nData Quality:\n")
cat("Initial samples:", validation_metrics$initial_samples, "\n")
cat("Analysis samples:", validation_metrics$final_samples, "\n")
cat("Initial genes:", validation_metrics$initial_genes, "\n")
cat("Final genes:", validation_metrics$final_genes, "\n")
cat("Gene retention:", round(validation_metrics$filtering_retention_rate, 1), "%\n")

cat("\nStatistical Results:\n")
cat("Total tests:", validation_metrics$total_tests, "\n")
cat("Significant (Bonf p<0.001):", validation_metrics$significant_bonf_0.001, "\n")
cat("Direction bias:", round(validation_metrics$direction_bias * 100, 1), "%\n")
cat("Min group size:", validation_metrics$min_group_size, "\n")

cat("\nMIR100HG Status:\n")
if (validation_metrics$mir100hg_in_results > 0) {
  cat("  Status: YES\n")
  cat("  LogFC:", round(validation_metrics$mir100hg_logfc, 4), "\n")
  cat("  Bonferroni p-value:", format(validation_metrics$mir100hg_bonf_pvalue, scientific = TRUE), "\n")
  cat("  Significant:", validation_metrics$mir100hg_significant_bonf_0.001, "\n")
} else {
  cat("  Status: NO\n")
}

cat("\nValidation Status:\n")
if (validation_metrics$validation_passed) {
  cat("  Overall: YES\n")
} else {
  cat("  Overall: NO\n")
  cat("  Flags:", validation_metrics$validation_flags, "\n")
}

cat("\nOutput Files Generated:\n")
cat("  Complete results:", basename(output_path_complete), "\n")
cat("  Filtered results:", basename(output_path_filtered), "\n")
cat("  Legacy format:", basename(legacy_filename), "\n")
cat("  Validation summary:", basename(validation_summary_path), "\n")
cat("  Diagnostic plots:", paste0("dea_diagnostics_cancer_", CANCER_TYPE, ".pdf"), "\n")


cat("Top 10 most significant genes (complete results):\n")
top_genes <- head(tab_complete[order(tab_complete$P.Value.Raw), ], 10)
print(top_genes[, c("SYMBOL", "logFC", "P.Value.Raw", "adj.P.Val.BH", "adj.P.Val.Bonf", "Bonf_sig_0.001")])

if (nrow(tab_filtered) > 0) {
  cat("\nTop 10 filtered results (existing pipeline):\n")
  print(head(tab_filtered[, c("SYMBOL", "logFC", "P.Value.Raw", "adj.P.Val.BH", "adj.P.Val.Bonf")], 10))
} else {
  cat("\nNo genes pass the filtered criteria (p <", p_threshold, ", |logFC| >", logfc_threshold, ")\n")
}

