
##---------------------------------------------------------
## Function to generate fixed number of spots (coordinates)
##---------------------------------------------------------

generate_spatial_spots <- function(num_spots, x1 = 0, x2 = 1, y1 = 0, y2 = 1, seed =NULL) {
  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  data.frame(
    x = runif(num_spots, min = x1, max = x2),
    y = runif(num_spots, min = y1, max = y2)
  )
}



filter_columns_by_sum <- function(df, threshold = 10, keep_first = TRUE) {
  # Convert to data.table for convenience
  library(data.table)
  dt <- as.data.table(df)
  
  # Identify columns to check (skip first if keep_first = TRUE)
  if (keep_first) {
    numeric_part <- dt[, -1, with = FALSE]
    sums <- colSums(numeric_part, na.rm = TRUE)
    cols_keep <- c(names(dt)[1], names(numeric_part)[sums > threshold])
  } else {
    sums <- colSums(dt, na.rm = TRUE)
    cols_keep <- names(dt)[sums > threshold]
  }
  # Return filtered table
  return(dt[, ..cols_keep])
}


make_rbf_basis <- function(coords, r = 30, lengthscale = 1.5) {
  N <- nrow(coords)
  
  # choose r knot centers using k-means or random subset
  km <- kmeans(coords, centers = r)
  knots <- km$centers
  
  # compute pairwise distance matrix between points and knots
  d <- as.matrix(dist(rbind(coords, knots)))
  D <- d[1:N, (N+1):(N+r)]  # N x r
  
  # set lengthscale = median distance between knots if not provided
  if (is.null(lengthscale)) {
    lengthscale <- median(as.matrix(dist(knots)))
  }
  
  # RBF basis
  Phi <- exp( - (D^2) / (2 * lengthscale^2) )
  
  return(Phi)   # N x r
}


make_rbf_dist <- function(coords, r = 30) {
  coords <- as.matrix(coords)
  N <- nrow(coords)
  
  km <- kmeans(coords, centers = r)
  knots <- km$centers
  
  all_d <- as.matrix(dist(rbind(coords, knots)))
  D <- all_d[1:N, (N+1):(N+r)]
  
  D
}


generate_gene_grp <- function(p, G = 3) {
  reps <- ceiling(p / G)
  rep(1:G, each = reps)[1:p]
}


get_tau_ell <- function(p, prop) {
  se_count <- round(prop[2] * p)
  list(
    tau_gs = t_gs[1:se_count],
    ell_gs = l_gs[1:se_count]
  )
}













