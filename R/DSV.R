DSV <- function(X,Y,labels,MR_ref = NULL, nc=3, ee=NULL){
  # Differentially Supervised Variamax algorithm
  # Inputs:
  #   X      : (n x pX) severity score matrix (e.g., CASI-5 domain sums), rows = subjects
  #   Y      : (n x pY) item matrix (e.g., ABC items), rows = subjects
  #   labels : length-n vector with group labels (e.g., "ASD" vs "TD") used for within-group standardization
  #   MR_ref : optional (pX x nc) reference MR matrix for Procrustes alignment (stabilizes signs/order across refits)
  #   nc     : number of retained components / factors (m in the paper)
  #   ee     : optional eigen(cov(Y)) object; supply to reuse the same PCA basis across runs/permutations
  
  # 1) Within-group rank-standardize so covariances equal Spearman correlations
  #    (and ASD/TD differences are not driven by scale/level differences).
  Y = rank_normalizeData_byLabels(Y,labels)
  X = rank_normalizeData_byLabels(X,labels)
  
  # 2) Differential Spearman correlation matrix B = Cov_ASD(X,Y) - Cov_TD(X,Y)
  coef = cov(X[labels=="ASD",],Y[labels=="ASD",]) - cov(X[labels=="TD",],Y[labels=="TD",])
  
  # 3) PCA basis of Y (pooled): reuse ee across calls/permutations for stability/comparability
  if (is.null(ee)){
    ee = eigen(cov(Y))
  }
  
  # 4) Map differential associations into the m=nc orthonormal component space
  #    M: severity -> components (differential); W: components -> original Y items
  M = coef%*% ee$vectors[,1:nc,drop=FALSE] %*% diag(1/sqrt(ee$values[1:nc]))
  W = diag(sqrt(ee$values[1:nc]))%*% t(ee$vectors[,1:nc,drop=FALSE])
  
  # 5) Varimax rotation chosen to sparsify differential severity->component associations (M)
  R = my_varimax(M,normalize=FALSE)$rotmat
  
  # 6) Rotated component scores ("optimal components") for subjects
  factors = sweep(Y %*% ee$vectors[,1:nc,drop=FALSE],2,1/sqrt(ee$values[1:nc]),FUN="*") %*% R
  
  # 7) Loadings from X to rotated components (differential severity->component associations)
  MR = M %*% R
  rownames(MR) = colnames(X)
  
  # 8) Optional Procrustes alignment to a reference solution (fixes sign/permutation ambiguity)
  Q <- diag(nc)
  if (!is.null(MR_ref)) {
    sv <- svd(crossprod(MR, MR_ref))   # t(Mc) %*% Mrc
    Q  <- sv$u %*% t(sv$v)
  }
  
  # Apply alignment consistently to scores, loadings, and rotation
  factors <- factors %*% Q
  MR      <- MR %*% Q
  R <- R %*% Q
  
  # 9) Rotated item-weight map W* = R' W (named RtW here), so Y_hat = factors %*% RtW
  RtW <- t(R) %*% W
  colnames(RtW) <- colnames(Y)
  
  # Implied differential X->Y effects through components: (X->C*) (C*->Y)
  effects <- MR %*% RtW
  
  # 10) Variance accounting (trace of Cov(Y)) and loadings-style matrice
  total_var <- sum(diag(cov(Y)))   # trace; ~ q if perfectly standardized
  
  L_unrot <- t(W)                  # q x nc (unrotated loadings)
  L_rot   <- t(RtW)                # q x nc (rotated + aligned loadings)
  
  Vaccounted_unrot <- variance_accounting(L_unrot, total_var = total_var)
  Vaccounted_rot   <- variance_accounting(L_rot,   total_var = total_var)
  
  list(
    MR = MR, RtW = RtW, effects = effects,
    optimal_components = factors,
    R = R, Q = Q,
    eigen = ee, X = X, Y = Y,
    total_variance_Y = total_var,
    loadings_unrot = L_unrot,
    loadings_rot = L_rot,
    Vaccounted_unrot = Vaccounted_unrot,
    Vaccounted_rot = Vaccounted_rot
  )
}
