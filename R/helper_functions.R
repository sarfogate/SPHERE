## ============================================================
## SPHERE: Helper Functions
## ============================================================


## ------------------------------------------------------------
## Color palette for spatial expression plots
## ------------------------------------------------------------

pal <- colorRampPalette(c('#00274c', '#00274c', "lightyellow2",
                          '#ffcb05', '#ffcb05'))


## ------------------------------------------------------------
## Plot spatial gene expression
## ------------------------------------------------------------

#' @keywords internal
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
    ggplot2::ggplot(df_plot, ggplot2::aes(x, y, color = expression)) +
      ggplot2::geom_point(size = 4) +
      ggplot2::scale_color_gradientn(colours = pal(5), limits = c(0, 1)) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = paste("Relative Expression for\n", gene_name),
                    color = "Relative\nExpression") +
      ggplot2::theme(axis.title = ggplot2::element_blank(),
                     axis.text  = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank())
  })
  do.call(gridExtra::grid.arrange, c(plots, nrow = nrow))
}


## ------------------------------------------------------------
## Filter genes with low total counts
## ------------------------------------------------------------

#' @keywords internal
filter_columns_by_sum <- function(df, threshold = 10, keep_first = TRUE) {
  ..cols_keep <- NULL  # fix R CMD CHECK note
  dt <- data.table::as.data.table(df)
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
## RBF basis matrix
## ------------------------------------------------------------

#' @keywords internal
make_rbf_basis <- function(coords, r = 30, lengthscale = 1.5) {
  N     <- nrow(coords)
  km    <- stats::kmeans(coords, centers = r)
  knots <- km$centers
  d     <- as.matrix(stats::dist(rbind(coords, knots)))
  D     <- d[1:N, (N + 1):(N + r)]
  if (is.null(lengthscale)) {
    lengthscale <- stats::median(as.matrix(stats::dist(knots)))
  }
  Phi <- exp(-(D^2) / (2 * lengthscale^2))
  return(Phi)
}


## ------------------------------------------------------------
## RBF distance matrix
## ------------------------------------------------------------

#' @keywords internal
make_rbf_dist <- function(coords, r = 30) {
  coords <- as.matrix(coords)
  N      <- nrow(coords)
  km     <- stats::kmeans(coords, centers = r)
  knots  <- km$centers
  all_d  <- as.matrix(stats::dist(rbind(coords, knots)))
  D      <- all_d[1:N, (N + 1):(N + r)]
  D
}


## ------------------------------------------------------------
## Generate gene pathway group membership
## ------------------------------------------------------------

#' @keywords internal
generate_gene_grp <- function(p, G = 3) {
  reps <- ceiling(p / G)
  rep(1:G, each = reps)[1:p]
}


## ------------------------------------------------------------
## Get tau and ell for SE genes
## ------------------------------------------------------------

#' @keywords internal
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

#' @keywords internal
relative_expr <- function(raw_exp) {
  raw_exp_log <- log1p(raw_exp)
  (raw_exp_log - min(raw_exp_log)) / (max(raw_exp_log) - min(raw_exp_log))
}


## ------------------------------------------------------------
## Simulate CAR pathway effects (internal)
## ------------------------------------------------------------

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

  deg   <- rowSums(W)
  D_mat <- diag(deg)
  Q0    <- D_mat - rho * W

  iso <- which(deg == 0)
  if (length(iso) > 0) Q0[iso, iso] <- 1

  Q  <- tau_beta * (Q0 + t(Q0)) / 2
  ev <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 0)
    stop(sprintf("Q is not positive definite; min eigenvalue = %.3e. Reduce rho.",
                 min(ev)))

  R    <- chol(Q)
  z    <- stats::rnorm(num_genes)
  Beta <- backsolve(R, z, transpose = TRUE)
  Beta <- backsolve(R, Beta)
  Beta
}


## ------------------------------------------------------------
## Simulate GP spatial effects (internal)
## ------------------------------------------------------------

#' @keywords internal
sim_gs_eta2 <- function(spots, num_genes, tau_gs, ell_gs, nugget = 1e-6) {
  n <- nrow(spots)
  if (length(tau_gs) != num_genes || length(ell_gs) != num_genes)
    stop("tau_gs and ell_gs must have length num_genes")

  c_mat   <- matrix(0, n, num_genes)
  dist_sq <- as.matrix(stats::dist(spots))^2

  for (j in 1:num_genes) {
    K_j       <- tau_gs[j] * exp(-dist_sq / (2 * ell_gs[j]^2))
    diag(K_j) <- diag(K_j) + nugget
    R         <- chol(K_j)
    z         <- stats::rnorm(n)
    c_mat[, j] <- t(R) %*% z
  }
  return(c_mat)
}
