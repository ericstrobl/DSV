parallel_analysis_DSV <- function(Y, labels,
                                   B = 1000, alpha = 0.95, ridge = 1e-8){
  # Parallel analysis with the DSV algorithm
  # Inputs:
  #   Y      : (n x p) item matrix (e.g., ABC items), rows = subjects
  #   labels : length-n group labels used for within-group rank-standardization before PCA
  #   B      : number of permutations for the null distribution
  #   alpha  : null quantile used as cutoff (e.g., 0.95 keeps components above 95% null eigenvalue)
  #   ridge  : unused here (kept for compatibility / future stabilizing tweaks)
  
  # Within-group rank-standardize before building the permutation null
  Y1 <- rank_normalizeData_byLabels(Y, labels)
  Y1 <- as.matrix(Y1)
  n  <- nrow(Y1); p <- ncol(Y1)
  
  # Observed eigenvalues of Cov(Y) (PCA-based parallel analysis target)
  eig_obs <- sort(eigen(cov(Y1))$values, decreasing = TRUE)
  
  # Permutation null: independently permute rows within each column (breaks dependence structure)
  # then recompute eigenvalues; keep order-specific null distributions.
  null_eigs <- matrix(NA_real_, nrow = B, ncol = p)
  for (b in seq_len(B)) {
    cat("\r", b)
    
    Yb <- Y
    for (j in seq_len(p)) Yb[, j] <- Y[sample.int(n), j]
    Yb <- rank_normalizeData_byLabels(Yb, labels)
    
    null_eigs[b, ] <- sort(eigen(cov(Yb))$values, decreasing = TRUE)
  }
  
  # Retain components whose observed eigenvalues exceed the alpha-quantile of the null
  cutoff <- apply(null_eigs, 2, stats::quantile, probs = alpha, type = 7, na.rm = TRUE)
  k      <- sum(eig_obs > cutoff)
  
  # Upper-tail permutation p-values per component index (smoothed with +1)
  obs_mat       <- matrix(eig_obs, nrow = B, ncol = p, byrow = TRUE)
  p_value_upper <- (1 + colSums(null_eigs >= obs_mat)) / (B + 1)
  
  list(
    k = k,
    eig_obs = eig_obs,
    cutoff = cutoff,
    alpha = alpha,
    B = B,
    p_value_upper = p_value_upper
  )
}
