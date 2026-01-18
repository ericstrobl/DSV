permutation_testing_DSV <- function(X, Y, labels,
                                    nperms = 1e5, nc = 3) {
  # Permutation testing with the DSV algorithm
  # Inputs:
  #   X      : (n x pX) severity score matrix (e.g., CASI-5 domains), rows = subjects
  #   Y      : (n x pY) item matrix (e.g., ABC items), rows = subjects
  #   labels : length-n group labels to be permuted under exchangeability (e.g., ASD vs TD)
  #   nperms : number of label permutations for the null distribution (paper uses 100,000)
  #   nc     : number of components retained in DSV during each refit
  
  # helper: convert an empirical FDR curve into q-values at observed stats
  # q(s) = min_{t <= s} FDRhat(t), enforcing monotonicity via cummin.
  q_from_curve <- function(stat_vec, t_grid, fdr_grid) {
    o <- order(t_grid)
    t <- t_grid[o]
    f <- pmin(1, fdr_grid[o])
    
    f_mon <- cummin(f)  # monotone envelope over increasing thresholds
    
    idx <- findInterval(stat_vec, t)  # largest threshold <= stat
    q <- rep(1, length(stat_vec))
    ok <- idx > 0
    q[ok] <- f_mon[idx[ok]]
    q
  }
  
  # 0) Build a stable reference orientation, then refit observed with alignment
  base0  <- DSV(X, Y, labels, nc = nc)    # no reference
  MR_ref <- base0$MR
  
  # Fix the PCA basis (ee) and align to MR_ref so permutation refits are comparable
  mod <- DSV(X, Y, labels, nc = nc, ee = base0$eigen, MR_ref = MR_ref)
  
  n  <- nrow(Y)
  pX <- ncol(X)
  pY <- ncol(Y)
  
  # Observed test statistics (absolute values for two-sided / magnitude testing) 
  absMR_obs   <- abs(mod$MR)   # pX x nc : severity->component differentials
  absRtW_obs  <- abs(mod$RtW)  # nc x pY : component->item weights
  
  T0_obs       <- sum(absMR_obs)      # omnibus magnitude across all entries
  Tfac_obs     <- colSums(absMR_obs)  # component-wise magnitudes
  Tmaxfac_obs  <- Tfac_obs            # used for maxT-style FWER adjustment
  
  # Counters for smoothed permutation p-values (>= observed) 
  c0         <- 0L
  c_fac      <- integer(nc)
  c_fac_fwer <- integer(nc)
  
  c_MR  <- matrix(0L, nrow = pX, ncol = nc)
  c_RtW <- matrix(0L, nrow = nc, ncol = pY)
  
  # Method-B FDR: choose threshold grids + accumulate expected null exceedances 
  # Factors: small family, use observed Tfac values as the step grid.
  fac_grid <- sort(unique(as.numeric(Tfac_obs)))
  EV_fac   <- numeric(length(fac_grid))  # sum_b V_b(t) across perms
  
  # Entrywise: use quantile grids of observed absolute stats (keeps evaluation stable/cheap).
  MR_vec   <- as.vector(absMR_obs)
  RtW_vec  <- as.vector(absRtW_obs)
  
  MR_pos   <- MR_vec[is.finite(MR_vec)]
  RtW_pos  <- RtW_vec[is.finite(RtW_vec)]
  
  MR_grid  <- unique(as.numeric(quantile(MR_pos,  probs = seq(0.10, 0.99, length.out = 60), na.rm = TRUE)))
  RtW_grid <- unique(as.numeric(quantile(RtW_pos, probs = seq(0.10, 0.99, length.out = 60), na.rm = TRUE)))
  
  if (length(MR_grid)  < 2) MR_grid  <- sort(unique(c(0, MR_pos)))
  if (length(RtW_grid) < 2) RtW_grid <- sort(unique(c(0, RtW_pos)))
  
  EV_MR  <- numeric(length(MR_grid))
  EV_RtW <- numeric(length(RtW_grid))
  
  # Permutation loop: permute group labels, refit DSV with fixed basis, align to MR_ref 
  for (pp in seq_len(nperms)) {
    cat("\r", pp)
    
    perm <- sample.int(n)
    
    mod1 <- DSV(X, Y, labels[perm],
                nc = nc, ee = base0$eigen, MR_ref = MR_ref)
    
    absMR_1  <- abs(mod1$MR)
    absRtW_1 <- abs(mod1$RtW)
    
    # omnibus (global differential association strength)
    T0_1 <- sum(absMR_1)
    c0   <- c0 + (T0_1 >= T0_obs)
    
    # component-wise
    Tfac_1 <- colSums(absMR_1)
    c_fac  <- c_fac + (Tfac_1 >= Tfac_obs)
    
    # maxT for component-wise FWER (compare each observed T_j to permutation max over j)
    maxT_1      <- max(Tfac_1)
    c_fac_fwer  <- c_fac_fwer + (maxT_1 >= Tmaxfac_obs)
    
    # entrywise exceedance counts (for permutation p-values)
    c_MR  <- c_MR  + (absMR_1  >= absMR_obs)
    c_RtW <- c_RtW + (absRtW_1 >= absRtW_obs)
    
    # Method-B accumulators: expected number of null exceedances at each threshold 
    for (k in seq_along(fac_grid)) EV_fac[k] <- EV_fac[k] + sum(Tfac_1 >= fac_grid[k])
    
    MR_perm_vec <- as.vector(absMR_1)
    for (k in seq_along(MR_grid))  EV_MR[k]  <- EV_MR[k]  + sum(MR_perm_vec  >= MR_grid[k])
    
    RtW_perm_vec <- as.vector(absRtW_1)
    for (k in seq_along(RtW_grid)) EV_RtW[k] <- EV_RtW[k] + sum(RtW_perm_vec >= RtW_grid[k])
  }
  
  # Smoothed permutation p-values (add-one)
  p0          <- (c0 + 1) / (nperms + 1)
  p_fac       <- (c_fac + 1) / (nperms + 1)
  p_fac_fwer  <- (c_fac_fwer + 1) / (nperms + 1)
  
  p_MR        <- (c_MR + 1) / (nperms + 1)
  p_RtW       <- (c_RtW + 1) / (nperms + 1)
  
  dimnames(p_MR)  <- dimnames(absMR_obs)
  dimnames(p_RtW) <- dimnames(absRtW_obs)
  
  # Method B (empirical FDR via thresholds): build curves + q-values for each family
  
  # Factors (nc tests)
  R_fac      <- vapply(fac_grid, function(t) sum(Tfac_obs >= t), numeric(1))
  EV_fac_bar <- EV_fac / nperms
  FDRhat_fac <- EV_fac_bar / pmax(R_fac, 1)
  q_fac      <- q_from_curve(Tfac_obs, fac_grid, FDRhat_fac)
  
  # MR entries (pX * nc tests)
  R_MR      <- vapply(MR_grid, function(t) sum(MR_vec >= t), numeric(1))
  EV_MR_bar <- EV_MR / nperms
  FDRhat_MR <- EV_MR_bar / pmax(R_MR, 1)
  qMR_vec   <- q_from_curve(MR_vec, MR_grid, FDRhat_MR)
  
  q_MR <- matrix(qMR_vec, nrow = pX, ncol = nc)
  dimnames(q_MR) <- dimnames(absMR_obs)
  
  # RtW entries (nc * pY tests)
  R_RtW      <- vapply(RtW_grid, function(t) sum(RtW_vec >= t), numeric(1))
  EV_RtW_bar <- EV_RtW / nperms
  FDRhat_RtW <- EV_RtW_bar / pmax(R_RtW, 1)
  qRtW_vec   <- q_from_curve(RtW_vec, RtW_grid, FDRhat_RtW)
  
  q_RtW <- matrix(qRtW_vec, nrow = nc, ncol = pY)
  dimnames(q_RtW) <- dimnames(absRtW_obs)
  
  # Output (keeps both p-values and Method-B q-values + curves for transparency) 
  list(
    omnibus = list(stat = T0_obs, pval = p0),
    
    factors = list(
      stat = Tfac_obs,
      pval = p_fac,
      pFWER = p_fac_fwer,
      qB = q_fac,
      fdr_curve = list(t = fac_grid, R = R_fac, EV = EV_fac_bar, FDRhat = FDRhat_fac)
    ),
    
    MR_ref = MR_ref,
    mod = mod,
    
    MR_entrywise = list(
      stat = absMR_obs,
      pval = p_MR,
      qB   = q_MR,
      fdr_curve = list(t = MR_grid, R = R_MR, EV = EV_MR_bar, FDRhat = FDRhat_MR)
    ),
    
    RtW_entrywise = list(
      stat = absRtW_obs,
      pval = p_RtW,
      qB   = q_RtW,
      fdr_curve = list(t = RtW_grid, R = R_RtW, EV = EV_RtW_bar, FDRhat = FDRhat_RtW)
    )
  )
}
