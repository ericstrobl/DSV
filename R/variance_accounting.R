variance_accounting <- function(loadings, total_var = NULL) {
  loadings <- as.matrix(loadings)            # q x m
  ss <- colSums(loadings^2)                  # SS loadings per factor (length m)
  
  if (is.null(total_var)) {
    # If variables are standardized and you used a correlation matrix, total_var = q.
    # Here we default to q, but in your pipeline trace(cov(Y)) is safer:
    total_var <- nrow(loadings)
  }
  
  prop_total <- ss / total_var               # proportion of TOTAL variance
  cum_total  <- cumsum(prop_total)
  
  prop_retained <- ss / sum(ss)              # within retained factor space
  cum_retained  <- cumsum(prop_retained)
  
  out <- data.frame(
    Factor = seq_along(ss),
    SS_Loadings = ss,
    PropVar_Total = prop_total,
    CumVar_Total  = cum_total,
    PropVar_Retained = prop_retained,
    CumVar_Retained  = cum_retained
  )
  
  # add percent versions for easy reporting
  out$PctVar_Total     <- 100 * out$PropVar_Total
  out$PctVar_Retained  <- 100 * out$PropVar_Retained
  out$CumPctVar_Total  <- 100 * out$CumVar_Total
  out$CumPctVar_Retained <- 100 * out$CumVar_Retained
  
  out
}
