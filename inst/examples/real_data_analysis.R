# ============================================================================
# SPHERE: Spatial Variable Gene Detection in Human Breast Cancer
# ============================================================================
# This script demonstrates a complete SPHERE analysis workflow using
# human breast cancer spatial transcriptomics data. It covers:
#   1. Loading and preparing the expression data
#   2. Handling missing pathway annotations
#   3. Filtering pathways by size
#   4. Fitting the SPHERE model
#   5. Inspecting and saving results
#
# NOTE: Running the full model fit may take several hours depending on
# the number of genes, spots, and chains.
# ============================================================================


# ----------------------------------------------------------------------------
# 0. Load package
# ----------------------------------------------------------------------------

library(SPHERE)


# ----------------------------------------------------------------------------
# 1. Load data
# ----------------------------------------------------------------------------

# Full human breast cancer spatial transcriptomics dataset
load("humanBC_dataset.RData")

# Pre-processed subset with selected genes and spots
load("BrC_New_set.RData")


# ----------------------------------------------------------------------------
# 2. Prepare expression data and spatial coordinates
# ----------------------------------------------------------------------------

# Expression count matrix (spots x genes)
# round() ensures integer counts for the Poisson likelihood
data_mat <- round(BrC_final)

# Spatial coordinates for each spot (x, y)
spot <- BrC_spot

# Quick sanity check
cat("Spots:", nrow(data_mat), "\n")
cat("Genes:", ncol(data_mat), "\n")
cat("Library size range:",
    paste(range(rowSums(data_mat)), collapse = " - "), "\n")


# ----------------------------------------------------------------------------
# 3. Prepare pathway annotations
# ----------------------------------------------------------------------------

# pathway_df is loaded from the .RData file and must contain a column
# named "Pathway" with one entry per gene (same order as columns of data_mat)

# 3a. Handle genes with missing pathway annotations
#     In the Stan CAR prior, a gene with no neighbors (n_j = 0) receives
#     an independent Normal(0, 10) prior on Beta[j], so these genes are
#     effectively treated as isolated nodes in the pathway graph.
na_idx <- which(is.na(pathway_df$Pathway))

if (length(na_idx) > 0) {
  # Give each unannotated gene a unique singleton pathway label
  pathway_df$Pathway[na_idx] <- paste0("isolated_", seq_along(na_idx))
  cat("Assigned", length(na_idx),
      "unannotated genes to individual singleton pathways.\n")
} else {
  cat("No missing pathway annotations found.\n")
}

# 3b. Filter pathways by size
#     Removes pathways that are too small or too large to be informative
#     for the CAR prior
gen_path           <- filter_pathways_by_limit(pathway_df)
gene_pathway_final <- gen_path$Pathway

cat("Genes retained    :", length(gene_pathway_final), "\n")
cat("Pathways retained :", length(unique(gene_pathway_final)), "\n")

# 3c. Subset data_mat to retained genes (if filter reduced gene set)
if (length(gene_pathway_final) < ncol(data_mat)) {
  data_mat <- data_mat[, names(gene_pathway_final)]
  cat("data_mat subsetted to", ncol(data_mat), "retained genes.\n")
}


# ----------------------------------------------------------------------------
# 4. Fit the SPHERE model
# ----------------------------------------------------------------------------

# NOTE: This step may take several hours for large datasets.
# Consider running overnight or on a computing cluster.

BrC_fit <- fit_sphere(
  data_mat          = data_mat,
  spot              = spot,
  gene_group        = gene_pathway_final,
  iter_sampling     = 5000,         # posterior draws per chain
  iter_warmup       = 2000,         # warmup (burn-in) per chain
  chains            = 3,            # number of independent MCMC chains
  seed              = 8,            # for reproducibility
  knots             = 30,           # inducing points for low-rank GP
  alpha             = c(10, 3),     # Dirichlet prior on mixture weights
  refresh           = 100,          # print progress every 100 iterations
  # Observation noise prior (half-normal)
  mu_noise          = 0,
  sd_noise          = 1,
  # GP lengthscale prior (log-normal)
  mu_gp_lengthscale = 0,
  sd_gp_lengthscale = 3,
  # GP amplitude prior (half-normal)
  mu_gp_amplitude   = 0,
  sd_gp_amplitude   = 12,
  # Global intercept prior (normal)
  mu_intercept      = 0,
  sd_intercept      = 1,
  # CAR spatial correlation prior (Beta)
  shape_beta_rho    = 5,
  rate_beta_rho     = 2,
  # CAR precision prior (half-normal)
  mu_beta_sig       = 0,
  sd_beta_sig       = 1
)

cat("Model fitting completed in", round(BrC_fit$runtime / 3600, 2),
    "hours.\n")


# ----------------------------------------------------------------------------
# 5. Inspect results
# ----------------------------------------------------------------------------

summary_df <- BrC_fit$summary

# 5a. Convergence diagnostics
cat("\nParameters with Rhat > 1.05:\n")
high_rhat <- summary_df[!is.na(summary_df$rhat) &
                          summary_df$rhat > 1.05, ]
if (nrow(high_rhat) == 0) {
  cat("None — all parameters converged!\n")
} else {
  print(high_rhat)
}

# 5b. Extract posterior gene classifications
#     Z[j] = 1 → gene j is non-spatially expressed (non-SE)
#     Z[j] = 2 → gene j is spatially expressed (SE)
z_rows    <- grep("^Z\\[", summary_df$variable)
z_summary <- summary_df[z_rows, c("variable", "mean", "median", "sd",
                                  "q2.5", "q97.5", "rhat")]

cat("\nPosterior Z summary (first 10 genes):\n")
print(head(z_summary, 10))

# 5c. Identify SE genes (posterior mean Z > 1.5 as threshold)
se_genes    <- z_summary$variable[z_summary$mean > 1.5]
nonse_genes <- z_summary$variable[z_summary$mean <= 1.5]

cat("\nDetected SE genes    :", length(se_genes), "\n")
cat("Detected non-SE genes:", length(nonse_genes), "\n")
cat("\nSE gene names:\n")
print(se_genes)

# 5d. GP lengthscale and amplitude summaries
ell_rows <- grep("^ell_gs\\[",     summary_df$variable)
sig_rows <- grep("^sig_eta_gs\\[", summary_df$variable)

cat("\nGP lengthscale summary (SE genes):\n")
print(summary_df[ell_rows, c("variable", "mean", "sd", "q2.5", "q97.5")])

cat("\nGP amplitude summary (SE genes):\n")
print(summary_df[sig_rows, c("variable", "mean", "sd", "q2.5", "q97.5")])


# ----------------------------------------------------------------------------
# 6. Save results
# ----------------------------------------------------------------------------

save(BrC_fit, z_summary, se_genes, nonse_genes,
     file = "SPHERE_BrC_results.RData")
cat("\nResults saved to SPHERE_BrC_results.RData\n")
