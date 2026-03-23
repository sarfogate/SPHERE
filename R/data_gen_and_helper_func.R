## ============================================================
## SPHERE: Data Generation and Helper Functions
## ============================================================


## ------------------------------------------------------------
## Function to generate fixed number of spots (coordinates)
## ------------------------------------------------------------

#' Generate Spatial Spot Coordinates
#'
#' Generates a data frame of random spatial coordinates for a given number
#' of spots within a specified bounding box.
#'
#' @param num_spots Integer. Number of spatial spots to generate.
#' @param x1 Numeric. Minimum x-coordinate (default: 0).
#' @param x2 Numeric. Maximum x-coordinate (default: 1).
#' @param y1 Numeric. Minimum y-coordinate (default: 0).
#' @param y2 Numeric. Maximum y-coordinate (default: 1).
#' @param seed Integer or NULL. Random seed for reproducibility (default: NULL).
#'
#' @return A data frame with columns \code{x} and \code{y} containing
#'   the spatial coordinates of each spot.
#'
#' @importFrom stats runif
#' @export
generate_spatial_spots <- function(num_spots, x1 = 0, x2 = 1,
                                   y1 = 0, y2 = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  data.frame(
    x = runif(num_spots, min = x1, max = x2),
    y = runif(num_spots, min = y1, max = y2)
  )
}


## ------------------------------------------------------------
## Function to simulate a named spatial pattern
## ------------------------------------------------------------

#' Simulate a Spatial Expression Pattern
#'
#' Assigns a spatial expression pattern (e.g. hotspot, streak, gradient)
#' to a set of spatial spots. Used internally to generate spatially
#' expressed gene patterns for simulation studies.
#'
#' @param spots A data frame with columns \code{x} and \code{y} giving
#'   spatial coordinates.
#' @param pat_type Character. Type of spatial pattern. One of
#'   \code{"hotspot"}, \code{"streak"}, \code{"gradient"}, \code{"ring"},
#'   or \code{"wave"}. Use \code{NULL} for no pattern (all low expression).
#' @param high_exp_prop Numeric. Proportion of spots assigned to the
#'   high-expression region (default: 0.3).
#' @param grad_percent Numeric. Gradient threshold in normalized coordinates
#'   for the \code{"gradient"} pattern (default: 0.3).
#' @param gene_id Integer or NULL. Gene index used to assign gene-specific
#'   random pattern centers. If NULL, uses the center of the tissue
#'   (default: NULL).
#' @param seed Integer or NULL. Random seed for reproducibility (default: NULL).
#'
#' @return The input \code{spots} data frame with two additional columns:
#'   \code{eta} (numeric spatial effect) and \code{mark} (character:
#'   \code{"high"}, \code{"medium"}, or \code{"low"}).
#'
#' @importFrom stats runif median
#' @export
simulate_expression <- function(spots, pat_type, high_exp_prop = 0.3,
                                grad_percent = 0.3, gene_id = NULL,
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  x1 <- min(spots$x); x2 <- max(spots$x)
  y1 <- min(spots$y); y2 <- max(spots$y)
  x_range <- x2 - x1; y_range <- y2 - y1

  # Normalize coordinates to [0,1] for scale-invariant patterns
  spots$x_norm <- (spots$x - x1) / x_range
  spots$y_norm <- (spots$y - y1) / y_range

  num_high  <- round(high_exp_prop * nrow(spots))
  spots$eta <- 0  # Default to low expression

  if (is.null(pat_type)) {
    spots$eta <- 0
  } else {
    center_x <- 0.5
    center_y <- 0.5

    if (!is.null(gene_id)) {
      center_x <- runif(1, 0.1, 0.9)
      center_y <- runif(1, 0.1, 0.9)
    }

    if (pat_type == "hotspot") {
      radius <- sqrt(sort((spots$x_norm - center_x)^2 +
                            (spots$y_norm - center_y)^2)[num_high])
      spots$eta <- as.integer((spots$x_norm - center_x)^2 +
                                (spots$y_norm - center_y)^2 <= radius^2)

    } else if (pat_type == "streak") {
      streak_x    <- if (is.null(gene_id)) center_x else runif(1, 0.1, 0.9)
      streak_spots <- order(abs(spots$x_norm - streak_x))[1:num_high]
      spots$eta[streak_spots] <- 1

    } else if (pat_type == "gradient") {
      threshold_x <- grad_percent
      spots$eta   <- pmax(0, pmin(1, (spots$x_norm - threshold_x) /
                                    (1 - threshold_x)))
      if (!is.null(gene_id)) {
        shift     <- runif(1, -0.2, 0.2)
        spots$eta <- pmax(0, pmin(1, spots$eta + shift))
      }

    } else if (pat_type == "ring") {
      distances   <- (spots$x_norm - center_x)^2 +
        (spots$y_norm - center_y)^2
      ring_spots  <- order(abs(distances - median(distances)))[1:num_high]
      spots$eta[ring_spots] <- 1

    } else if (pat_type == "wave") {
      phase <- if (!is.null(gene_id)) runif(1, -1, 1) else 0
      freq  <- if (!is.null(gene_id)) runif(1, 2.5, 5.5) else 4
      w_y   <- sin(freq * pi * (spots$x_norm + phase))
      y_scaled      <- (w_y - min(w_y)) / (max(w_y) - min(w_y))
      dist_to_curve <- abs(spots$y_norm - y_scaled)
      wave_spots    <- order(dist_to_curve)[1:num_high]
      spots$eta[wave_spots] <- 1
    }
  }

  # Assign expression marks
  if (is.null(pat_type)) {
    spots$mark <- "low"
  } else if (pat_type == "gradient") {
    spots$mark <- ifelse(spots$eta >= 0.75, "high",
                         ifelse(spots$eta <= 0.40, "low", "medium"))
  } else {
    spots$mark <- ifelse(spots$eta >= 0.75, "high", "low")
  }

  spots$x_norm <- NULL
  spots$y_norm <- NULL

  return(spots)
}


## ------------------------------------------------------------
## Function -- Pattern-based ST count data generator
## ------------------------------------------------------------

#' Generate Pattern-Based Spatial Transcriptomics Count Data
#'
#' Simulates spatial transcriptomics count data with a user-specified
#' spatial expression pattern. Genes are classified as spatially expressed
#' (SE) or non-spatially expressed (non-SE) according to \code{prop}.
#' SE genes follow a pattern defined by \code{pat_type}, while non-SE
#' genes have flat Poisson expression.
#'
#' @param spots A data frame or matrix of spatial coordinates
#'   (columns \code{x} and \code{y}).
#' @param num_genes Integer. Total number of genes to simulate.
#' @param prop Numeric vector of length 2 summing to 1. Proportions of
#'   non-SE and SE genes respectively (default: \code{c(0.8, 0.2)}).
#' @param pat_type Character. Spatial pattern type for SE genes. One of
#'   \code{"hotspot"}, \code{"streak"}, \code{"gradient"}, \code{"ring"},
#'   or \code{"wave"}.
#' @param seed Integer or NULL. Random seed for reproducibility
#'   (default: NULL).
#' @param G Integer. Number of biological pathway groups.
#' @param gene_grp Integer vector of length \code{num_genes} giving
#'   pathway membership for each gene.
#' @param cor_grp Numeric. CAR spatial correlation parameter \eqn{\rho}
#'   for pathway effects (default: see \code{create_car_beta2}).
#' @param tau_beta Numeric. CAR precision parameter for pathway effects.
#' @param high_exp_prop Numeric. Proportion of spots in the
#'   high-expression region (default: 0.3).
#' @param grad_percent Numeric. Gradient threshold for the
#'   \code{"gradient"} pattern (default: 0.2).
#' @param beta0 Numeric. Log-scale baseline expression for non-SE genes
#'   and low-expression spots (default: 2).
#' @param beta1 Numeric or NULL. Log-scale high expression for SE genes.
#'   If NULL, defaults to \code{beta0 + log(3)} (3x fold change).
#' @param sigma Numeric. Standard deviation of log-scale noise
#'   (default: 0.3).
#'
#' @return A named list with elements:
#' \describe{
#'   \item{Y}{Integer matrix (\eqn{n \times p}) of simulated counts.}
#'   \item{se_idx}{Integer vector of SE gene column indices.}
#'   \item{nonse_idx}{Integer vector of non-SE gene column indices.}
#'   \item{Beta}{Numeric vector of CAR pathway effects.}
#'   \item{epsilon}{Numeric matrix of log-scale noise.}
#'   \item{beta0}{Baseline log-mean used.}
#'   \item{beta1}{High-expression log-mean used.}
#' }
#'
#' @importFrom stats rnorm rpois quantile
#' @export
generate_genedata_pattern <- function(
    spots, num_genes, prop = c(0.8, 0.2), pat_type, seed = NULL,
    G, gene_grp, cor_grp, tau_beta,
    high_exp_prop = 0.3,
    grad_percent  = 0.2,
    beta0         = 2,
    beta1         = NULL,
    sigma         = 0.3) {

  if (!is.null(seed)) set.seed(seed)

  if (length(prop) != 2 || any(prop < 0) || abs(sum(prop) - 1) > 1e-6)
    stop("prop must be a vector of length 2 summing to 1")

  spots <- as.data.frame(spots)
  colnames(spots)[1:2] <- c("x", "y")

  n <- nrow(spots)
  p <- num_genes

  if (is.null(beta1)) beta1 <- beta0 + log(3)

  # Gene counts per category
  num_cat    <- floor(p * prop)
  num_cat[which.max(prop * p - num_cat)] <-
    num_cat[which.max(prop * p - num_cat)] + (p - sum(num_cat))
  n_nonse <- num_cat[1]
  n_se    <- num_cat[2]

  # Random gene arrangement
  se_idx    <- sort(sample(seq_len(p), n_se))
  nonse_idx <- setdiff(seq_len(p), se_idx)
  Z         <- rep(1L, p)
  Z[se_idx] <- 2L

  # CAR gene-level effects
  Beta <- create_car_beta2(p, gene_grp, cor_grp, tau_beta)

  # Noise
  epsilon <- matrix(rnorm(n * p, mean = 0, sd = sigma), nrow = n, ncol = p)

  # Generate counts
  result_mat <- matrix(0L, nrow = n, ncol = p)

  for (j in nonse_idx) {
    loglambda_j     <- beta0 + Beta[j] + epsilon[, j]
    result_mat[, j] <- rpois(n, lambda = exp(loglambda_j))
  }

  for (j in se_idx) {
    pattern_data <- simulate_expression(
      spots, pat_type = pat_type, high_exp_prop = high_exp_prop,
      grad_percent = grad_percent, gene_id = j)

    loglambda_j            <- beta0 + Beta[j] + epsilon[, j]
    high_mask              <- pattern_data$mark == "high"
    med_mask               <- pattern_data$mark == "medium"
    loglambda_j[high_mask] <- beta1 + Beta[j] + epsilon[high_mask, j]
    loglambda_j[med_mask]  <- (beta0 + beta1) / 2 + Beta[j] +
      epsilon[med_mask, j]

    result_mat[, j] <- rpois(n, lambda = exp(loglambda_j))
  }

  # Gene naming
  gene_names            <- character(p)
  gene_names[nonse_idx] <- paste0("Gene", nonse_idx, "_NonSE")
  gene_names[se_idx]    <- paste0("Gene", se_idx, "_SE_", pat_type)
  colnames(result_mat)  <- gene_names

  cat("---- Pattern Simulation Summary ----\n")
  cat("n spots   :", n, " | P genes:", p, "\n")
  cat("SE genes  :", n_se, " | Non-SE genes:", n_nonse, "\n")
  cat("SE indices:", paste(se_idx, collapse = ", "), "\n")
  cat("Pattern   :", pat_type, "\n")
  cat("beta0 =", round(beta0, 3), "| beta1 =", round(beta1, 3),
      "| fold change =", round(exp(beta1 - beta0), 2), "x\n")
  cat("sigma     :", sigma, "\n")
  cat("Observed counts: median =", median(result_mat),
      "| 95% =", paste(round(quantile(as.vector(result_mat),
                                      c(.025, .975))), collapse = ", "),
      "| max =", max(result_mat), "\n")
  cat("------------------------------------\n")

  return(list(
    Y         = result_mat,
    se_idx    = se_idx,
    nonse_idx = nonse_idx,
    Beta      = Beta,
    epsilon   = epsilon,
    beta0     = beta0,
    beta1     = beta1
  ))
}


## ------------------------------------------------------------
## Function -- SPHERE model data generator
## ------------------------------------------------------------

#' Generate SPHERE Model Simulation Data
#'
#' Simulates spatial transcriptomics count data consistent with the
#' SPHERE model. Includes library size variation, CAR pathway effects,
#' gene-specific Gaussian Process spatial effects, and auto-calibrated
#' baseline expression.
#'
#' @param spots A matrix or data frame of spatial coordinates
#'   (at least 2 columns: x, y).
#' @param num_genes Integer. Total number of genes to simulate.
#' @param prop Numeric vector of length 2 summing to 1. Proportions of
#'   non-SE and SE genes respectively (default: \code{c(0.8, 0.2)}).
#' @param G Integer. Number of biological pathway groups.
#' @param gene_grp Integer vector of length \code{num_genes} giving
#'   pathway membership for each gene.
#' @param rho Numeric. CAR spatial correlation for pathway effects
#'   (default: 0.85).
#' @param tau_beta Numeric. CAR precision for pathway effects
#'   (default: 10).
#' @param tau_gs Numeric vector or NULL. GP amplitude for each SE gene.
#'   If NULL, sampled from \code{Uniform(0.05, 0.30)}.
#' @param ell_gs Numeric vector or NULL. GP lengthscale for each SE gene.
#'   If NULL, sampled from \code{Uniform(0.50, 2.00)}.
#' @param nugget Numeric. Nugget added to GP covariance diagonal for
#'   numerical stability (default: 1e-6).
#' @param eps_sd Numeric. Standard deviation of iid log-scale noise
#'   (default: 0.20).
#' @param depth_model Character. Model for sequencing depth. One of
#'   \code{"poisson"}, \code{"negbin"}, or \code{"fixed"}
#'   (default: \code{"poisson"}).
#' @param depth_lambda Numeric. Mean for Poisson or NegBin depth model
#'   (default: 5).
#' @param depth_size Numeric. Overdispersion parameter for NegBin depth
#'   model — smaller values give more overdispersion (default: 5).
#' @param depth_fixed Numeric. Fixed library size when
#'   \code{depth_model = "fixed"} (default: 1).
#' @param target_mean Numeric. Desired median count used to auto-calibrate
#'   the global intercept \code{mu0} (default: 10).
#' @param seed Integer or NULL. Random seed for reproducibility
#'   (default: NULL).
#' @param verbose Logical. If TRUE, prints a simulation summary
#'   (default: TRUE).
#'
#' @return A named list with elements:
#' \describe{
#'   \item{Y}{Integer matrix (\eqn{n \times p}) of simulated counts.}
#'   \item{N_i}{Integer vector of per-spot library sizes.}
#'   \item{mu0}{Numeric. Auto-calibrated global log-expression intercept.}
#'   \item{Beta}{Numeric vector of CAR pathway effects.}
#'   \item{epsilon}{Numeric matrix of iid log-scale noise.}
#'   \item{eta}{Numeric matrix of GP spatial effects.}
#'   \item{Z}{Integer vector of true gene classifications (1=non-SE, 2=SE).}
#'   \item{se_idx}{Integer vector of SE gene indices.}
#'   \item{tau_gs}{Numeric vector of GP amplitudes for SE genes.}
#'   \item{ell_gs}{Numeric vector of GP lengthscales for SE genes.}
#'   \item{lambda}{Numeric matrix of latent Poisson rates.}
#'   \item{Y_mean}{Numeric matrix of expected counts.}
#' }
#'
#' @importFrom stats rnorm rpois rnbinom quantile median
#' @export
gen_genedata_model <- function(
    spots, num_genes, prop = c(0.8, 0.2), G, gene_grp,
    rho = 0.85, tau_beta = 10, tau_gs = NULL, ell_gs = NULL,
    nugget = 1e-6, eps_sd = 0.20,
    depth_model  = c("poisson", "negbin", "fixed"),
    depth_lambda = 5,
    depth_size   = 5,
    depth_fixed  = 1,
    target_mean  = 10,
    seed = NULL, verbose = TRUE) {

  depth_model <- match.arg(depth_model)

  stopifnot(is.matrix(spots) || is.data.frame(spots))
  spots <- as.matrix(spots)
  stopifnot(ncol(spots) >= 2)

  if (length(prop) != 2 || any(prop < 0) || abs(sum(prop) - 1) > 1e-8)
    stop("prop must be length-2 and sum to 1 (Non-SE, SE).")

  if (length(gene_grp) != num_genes)
    stop("gene_grp must have length num_genes.")

  if (!is.null(seed)) set.seed(seed)

  n <- nrow(spots)
  P <- num_genes

  # Sequencing depth
  if (depth_model == "poisson") {
    N_i <- rpois(n, lambda = depth_lambda) + 1L
  } else if (depth_model == "negbin") {
    N_i <- rnbinom(n, mu = depth_lambda, size = depth_size) + 1L
  } else {
    N_i <- rep.int(as.integer(depth_fixed), n)
  }

  # Gene-level state
  n_se      <- round(P * prop[2])
  se_idx    <- sort(sample(seq_len(P), n_se))
  nonse_idx <- setdiff(seq_len(P), se_idx)
  Z         <- rep(1L, P)
  Z[se_idx] <- 2L

  # CAR Beta
  Beta <- create_car_beta2(num_genes = P, gene_grp = gene_grp,
                           rho = rho, tau_beta = tau_beta)

  # iid noise
  epsilon <- matrix(rnorm(n * P, mean = 0, sd = eps_sd), nrow = n, ncol = P)

  # Spatial effects
  eta <- matrix(0, nrow = n, ncol = P)

  if (n_se > 0) {
    if (is.null(tau_gs)) tau_gs <- runif(n_se, 0.05, 0.30)
    if (is.null(ell_gs)) ell_gs <- runif(n_se, 0.50, 2.00)

    if (length(tau_gs) != n_se) stop("tau_gs must have length = number of SE genes.")
    if (length(ell_gs) != n_se) stop("ell_gs must have length = number of SE genes.")

    eta_se        <- sim_gs_eta2(spots = spots[, 1:2, drop = FALSE],
                                 num_genes = n_se,
                                 tau_gs = tau_gs, ell_gs = ell_gs,
                                 nugget = nugget)
    eta[, se_idx] <- eta_se
  } else {
    tau_gs <- numeric(0)
    ell_gs <- numeric(0)
  }

  # Auto-calibrate mu0
  log_part <- sweep(epsilon + eta, 2, Beta, `+`)
  mu0      <- log(target_mean) - log(median(N_i)) - median(log_part)

  # Generate counts
  lambda <- exp(mu0 + log_part)
  Y_mean <- outer(N_i, rep(1, P)) * lambda
  Y      <- matrix(rpois(n * P, lambda = as.vector(Y_mean)),
                   nrow = n, ncol = P)

  # Gene naming
  gene_names            <- character(P)
  gene_names[nonse_idx] <- paste0("Gene", nonse_idx, "_NonSE")
  gene_names[se_idx]    <- paste0("Gene", se_idx, "_SE_model")
  colnames(Y)           <- gene_names

  if (verbose) {
    cat("---- ST Simulation Summary ----\n")
    cat("n spots:", n, " | P genes:", P, "\n")
    cat("SE genes:", n_se, " (", round(100 * n_se / P, 1), "% )\n", sep = "")
    cat("Depth model:", depth_model,
        " | median N_i:", median(N_i),
        " | range:", paste(range(N_i), collapse = " - "), "\n")
    cat("mu0 (auto-calibrated):", signif(mu0, 4), "\n")
    cat("Expected mean counts : median =", signif(median(Y_mean), 4),
        " | 95% =", paste(signif(quantile(as.vector(Y_mean),
                                          c(.025, .975)), 4),
                          collapse = ", "),
        " | max =", signif(max(Y_mean), 4), "\n")
    cat("Observed counts      : median =", median(Y),
        " | 95% =", paste(quantile(as.vector(Y), c(.025, .975)),
                          collapse = ", "),
        " | max =", max(Y), "\n")
    cat("-------------------------------\n")
  }

  return(list(
    Y       = Y,
    N_i     = N_i,
    mu0     = mu0,
    Beta    = Beta,
    epsilon = epsilon,
    eta     = eta,
    Z       = Z,
    se_idx  = se_idx,
    tau_gs  = tau_gs,
    ell_gs  = ell_gs,
    lambda  = lambda,
    Y_mean  = Y_mean
  ))
}


## ------------------------------------------------------------
## Plot all Figures
## ------------------------------------------------------------

pal <- colorRampPalette(c('#00274c', '#00274c', "lightyellow2",
                          '#ffcb05', '#ffcb05'))

#' Plot Gene Expression Data
#'
#' Creates spatial expression plots for selected genes arranged in a grid.
#'
#' @param gene_df A matrix or data frame of gene expression values
#'   (spots x genes).
#' @param spot A matrix of spatial coordinates with two columns (x, y).
#' @param genes Character vector of gene names to plot. If NULL, plots
#'   all genes (default: NULL).
#' @param nrow Integer. Number of rows in the plot grid (default: 1).
#'
#' @return A grid of ggplot objects arranged by \pkg{gridExtra}.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradientn
#'   theme_minimal labs theme element_blank
#' @export
plot_gene_data <- function(gene_df, spot, genes = NULL, nrow = 1) {
  x <- y <- NULL  # fix R CMD CHECK note
  gene_df <- as.matrix(gene_df)
  if (is.null(genes)) genes <- colnames(gene_df)

  rel_expr <- apply(gene_df, 2, function(x) {
    if (max(x) == min(x)) return(rep(0, length(x)))
    (x - min(x)) / (max(x) - min(x))
  })
  rel_expr <- as.matrix(rel_expr)

  plots <- lapply(genes, function(gene_name) {
    df_plot <- data.frame(x = spot[, 1], y = spot[, 2],
                          expression = rel_expr[, gene_name])
    ggplot(df_plot, aes(x, y, color = expression)) +
      geom_point(size = 4) +
      scale_color_gradientn(colours = pal(5), limits = c(0, 1)) +
      theme_minimal() +
      labs(title = paste("Relative Expression for\n", gene_name),
           color = "Relative\nExpression") +
      theme(axis.title = element_blank(),
            axis.text  = element_blank(),
            axis.ticks = element_blank())
  })
  do.call(gridExtra::grid.arrange, c(plots, nrow = nrow))
}


## ------------------------------------------------------------
## Function to filter genes with less total counts
## ------------------------------------------------------------

#' Filter Genes with Low Total Counts
#'
#' Removes genes (columns) whose total count across all spots falls
#' below a specified threshold.
#'
#' @param df A matrix or data frame of gene expression counts.
#' @param threshold Numeric. Minimum total count required to keep a
#'   gene (default: 10).
#' @param keep_first Logical. If TRUE, always keeps the first column
#'   regardless of its sum (default: TRUE).
#'
#' @return A filtered \code{data.table} with low-count genes removed.
#'
#' @importFrom data.table as.data.table
#' @export
filter_columns_by_sum <- function(df, threshold = 10, keep_first = TRUE) {
  ..cols_keep <- NULL  # fix R CMD CHECK note
  dt <- as.data.table(df)
  if (keep_first) {
    numeric_part <- dt[, -1, with = FALSE]
    sums         <- colSums(numeric_part, na.rm = TRUE)
    cols_keep    <- c(names(dt)[1], names(numeric_part)[sums > threshold])
  } else {
    sums      <- colSums(dt, na.rm = TRUE)
    cols_keep <- names(dt)[sums > threshold]
  }
  return(dt[, ..cols_keep])
}


## ------------------------------------------------------------
## Function to make RBF basis matrix
## ------------------------------------------------------------

#' Construct Low-Rank RBF Gaussian Process Basis Matrix
#'
#' Computes a radial basis function (RBF) basis matrix using \code{r}
#' knot centers selected by k-means clustering. Used to construct the
#' low-rank GP approximation in SPHERE.
#'
#' @param coords A numeric matrix of spatial coordinates (\eqn{n \times 2}).
#' @param r Integer. Number of knots (inducing points) (default: 30).
#' @param lengthscale Numeric or NULL. RBF lengthscale. If NULL, uses
#'   the median pairwise distance between knots (default: 1.5).
#'
#' @return A numeric matrix of dimension \eqn{n \times r} containing
#'   the RBF basis function values.
#'
#' @importFrom stats kmeans dist median
#' @export
make_rbf_basis <- function(coords, r = 30, lengthscale = 1.5) {
  N  <- nrow(coords)
  km <- kmeans(coords, centers = r)
  knots <- km$centers
  d  <- as.matrix(dist(rbind(coords, knots)))
  D  <- d[1:N, (N + 1):(N + r)]
  if (is.null(lengthscale)) {
    lengthscale <- median(as.matrix(dist(knots)))
  }
  Phi <- exp(-(D^2) / (2 * lengthscale^2))
  return(Phi)
}


## ------------------------------------------------------------
## Function to make RBF distance matrix
## ------------------------------------------------------------

#' Compute Spot-to-Knot Distance Matrix for GP Approximation
#'
#' Computes the matrix of Euclidean distances from each spatial spot
#' to each of \code{r} knot centers selected by k-means clustering.
#'
#' @param coords A numeric matrix of spatial coordinates (\eqn{n \times 2}).
#' @param r Integer. Number of knots (inducing points) (default: 30).
#'
#' @return A numeric matrix of dimension \eqn{n \times r} containing
#'   Euclidean distances from spots to knots.
#'
#' @importFrom stats kmeans dist
#' @export
make_rbf_dist <- function(coords, r = 30) {
  coords <- as.matrix(coords)
  N      <- nrow(coords)
  km     <- kmeans(coords, centers = r)
  knots  <- km$centers
  all_d  <- as.matrix(dist(rbind(coords, knots)))
  D      <- all_d[1:N, (N + 1):(N + r)]
  D
}


## ------------------------------------------------------------
## Function to generate gene pathway group
## ------------------------------------------------------------

#' Generate Gene Pathway Group Membership
#'
#' Assigns genes to \code{G} pathway groups in a round-robin fashion.
#' Genes are assigned sequentially: gene 1 to group 1, gene 2 to group
#' 2, ..., gene G to group G, gene G+1 to group 1, and so on.
#'
#' @param p Integer. Total number of genes.
#' @param G Integer. Number of pathway groups (default: 3).
#'
#' @return An integer vector of length \code{p} giving pathway group
#'   membership for each gene.
#'
#' @export
generate_gene_grp <- function(p, G = 3) {
  reps <- ceiling(p / G)
  rep(1:G, each = reps)[1:p]
}


## ------------------------------------------------------------
## Function to get the tau and ell for SE genes
## ------------------------------------------------------------

#' Get Tau and Ell Parameters for Spatially Expressed Genes
#'
#' Extracts the GP amplitude (tau) and lengthscale (ell) parameters
#' for spatially expressed (SE) genes based on a given proportion.
#'
#' @param p Integer. Total number of genes.
#' @param prop Numeric vector of length 2 giving the proportion of
#'   non-SE and SE genes respectively (e.g. \code{c(0.7, 0.3)}).
#' @param t_gs Numeric vector of GP amplitude values for SE genes.
#' @param l_gs Numeric vector of GP lengthscale values for SE genes.
#'
#' @return A named list with elements \code{tau_gs} and \code{ell_gs},
#'   each a numeric vector of length \code{round(prop[2] * p)}.
#'
#' @export
get_tau_ell <- function(p, prop, t_gs, l_gs) {
  se_count <- round(prop[2] * p)
  list(
    tau_gs = t_gs[1:se_count],
    ell_gs = l_gs[1:se_count]
  )
}


## ------------------------------------------------------------
## Relative expression
## ------------------------------------------------------------

#' Compute Relative Gene Expression
#'
#' Applies a log1p transformation followed by min-max normalization
#' to convert raw counts to relative expression values in \eqn{[0, 1]}.
#'
#' @param raw_exp A numeric vector or matrix of raw expression counts.
#'
#' @return A numeric vector or matrix of relative expression values
#'   scaled to \eqn{[0, 1]}.
#'
#' @export
relative_expr <- function(raw_exp) {
  raw_exp_log <- log1p(raw_exp)
  (raw_exp_log - min(raw_exp_log)) / (max(raw_exp_log) - min(raw_exp_log))
}


## ------------------------------------------------------------
## Create pathway effects (internal)
## ------------------------------------------------------------

#' Simulate CAR Pathway Effects for Genes
#'
#' Simulates gene-level pathway effects \eqn{\beta_j} from a
#' Conditional Autoregressive (CAR) prior. Genes in the same pathway
#' group are treated as spatial neighbors and share information.
#'
#' @param num_genes Integer. Total number of genes.
#' @param gene_grp Integer vector of length \code{num_genes} giving
#'   pathway group membership for each gene.
#' @param rho Numeric. CAR spatial correlation parameter (default: 0.9).
#' @param tau_beta Numeric. CAR precision parameter (default: 1).
#'
#' @return A numeric vector of length \code{num_genes} containing
#'   simulated gene-level pathway effects.
#'
#' @importFrom stats rnorm
#' @keywords internal
create_car_beta2 <- function(num_genes, gene_grp, rho = 0.9, tau_beta = 1) {
  stopifnot(length(gene_grp) == num_genes)

  W <- matrix(0, num_genes, num_genes)
  for (i in 1:(num_genes - 1)) {
    for (j in (i + 1):num_genes) {
      if (gene_grp[i] == gene_grp[j]) {
        W[i, j] <- 1
        W[j, i] <- 1
      }
    }
  }

  deg  <- rowSums(W)
  D_mat <- diag(deg)
  Q0   <- D_mat - rho * W

  iso <- which(deg == 0)
  if (length(iso) > 0) Q0[iso, iso] <- 1

  Q  <- tau_beta * (Q0 + t(Q0)) / 2
  ev <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 0)
    stop(sprintf("Q is not positive definite; min eigenvalue = %.3e. Reduce rho.",
                 min(ev)))

  R    <- chol(Q)
  z    <- rnorm(num_genes)
  Beta <- backsolve(R, z, transpose = TRUE)
  Beta <- backsolve(R, Beta)
  Beta
}


## ------------------------------------------------------------
## Create gene-specific spatial variability (internal)
## ------------------------------------------------------------

#' Simulate Gene-Specific Gaussian Process Spatial Effects.
#'
#' Simulates spatial expression effects for SE genes using
#' gene-specific Gaussian Process covariance functions with an
#' RBF (squared exponential) kernel.
#'
#' @param spots A numeric matrix of spatial coordinates (\eqn{n \times 2}).
#' @param num_genes Integer. Number of SE genes to simulate.
#' @param tau_gs Numeric vector of length \code{num_genes}. GP amplitude
#'   (signal variance) for each gene.
#' @param ell_gs Numeric vector of length \code{num_genes}. GP lengthscale
#'   for each gene.
#' @param nugget Numeric. Nugget added to diagonal for numerical stability
#'   (default: 1e-6).
#'
#' @return A numeric matrix of dimension \eqn{n \times} \code{num_genes}
#'   containing simulated GP spatial effects.
#'
#' @importFrom stats dist rnorm
#' @keywords internal
sim_gs_eta2 <- function(spots, num_genes, tau_gs, ell_gs, nugget = 1e-6) {
  n <- nrow(spots)
  if (length(tau_gs) != num_genes || length(ell_gs) != num_genes)
    stop("tau_gs and ell_gs must have length num_genes")

  c_mat   <- matrix(0, n, num_genes)
  dist_sq <- as.matrix(dist(spots))^2

  for (j in 1:num_genes) {
    K_j        <- tau_gs[j] * exp(-dist_sq / (2 * ell_gs[j]^2))
    diag(K_j)  <- diag(K_j) + nugget
    R          <- chol(K_j)
    z          <- rnorm(n)
    c_mat[, j] <- t(R) %*% z
  }
  return(c_mat)
}
