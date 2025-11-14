#  KO_FDR(W, alpha_nominal)
#  --------------------------------------------------------------
#  Implements the original knockoff filter introduced here:
#. https://doi.org/10.1111/rssb.12265.
#  Controls the FDR at a nominal level (alpha_nominal) using a data-driven
#  threshold on the feature importance statistics (W).
#  - Input:
#     W : numeric vector of feature statistics.
#     alpha_nominal : desired FDR control level (e.g., 0.1).
#  - Output:
#     Vector of selected (rejected) feature indices.
KO_FDR <- function(W, alpha_nominal){
  ts = sort(c(0, abs(W)))
  ts = ts[ts>0]
  ratio = sapply(ts, function(t)
    (1 + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= alpha_nominal)
  threshold = ifelse(length(ok) > 0, ts[ok[1]], Inf)
  S = sort(which(W >= threshold))
  return(S)
}



#  KO_FDR_posthoc(W, alpha_nominal)
#  --------------------------------------------------------------
#  Implements our post-hoc adaptive version of the original knockoff filter.
#  This version computes a data-dependent alpha (alpha_dep), adjusting
#  the level based on the observed data configuration.
#  - Input:
#     W : numeric vector of feature statistics.
#     alpha_nominal : initial significance level.
#  - Output:
#      A list with:
#        • alpha_dep : data-dependent nominal type-I error
#        • R : indices of selected features.
KO_FDR_posthoc <- function(W, alpha_nominal){
  
  ts = sort(abs(W))
  ts = ts[ts>0]
  ratio = sapply(ts, function(t)
    (1 + sum(W <= -t)) / max(1, sum(W >= t)))
  sum_neg = sapply(ts, function(t)
    sum(W <= -t))
  
  ok = union(which(ratio <= alpha_nominal), which(sum_neg == 0))
  
  threshold = ifelse(length(ok) > 0, ts[ok[1]], Inf)
  
  R = which(W >= threshold)
  
  if(length(R)>0){
    U <- which(W<=-threshold)
    alpha_dep <- (1+length(U))/length(R)
  }else{
    alpha_dep <- alpha_nominal
  }
  
  # Returns the data-dependent alpha and the rejection set R
  return(list(alpha_dep = alpha_dep, R = R))
}



#  KO_FDR_derand(W, alpha_nominal, alpha_ebh)
#  --------------------------------------------------------------
#  Implements the derandomized knockoff procedure, introduced here:
#  https://doi.org/10.1093/jrsssb/qkad085.
#  Uses an e-value formulation and the Benjamini–Hochberg
#  (eBH) procedure for FDR control.
#  - Input:
#      W : matrix of W-statistics (rows = features, columns = knockoff runs)
#     alpha_nominal : nominal FDR control level for each replicate.
#        alpha_ebh : significance level for the eBH step.
#   - Output:
#      A list with:
#        • alpha_dep : alpha used in the eBH step.
#        • R : indices of selected features.
KO_FDR_derand <- function(W, alpha_nominal, alpha_ebh){
  k = ncol(W)
  p = nrow(W)
  thresholds = rep(0, k)
  for(j in 1:k){
    ts = sort(abs(W[, j]))
    ts = ts[ts>0]
    ratio = sapply(ts, function(t)
      (1 + sum(W[, j] <= -t)) / max(1, sum(W[, j] >= t)))
    sum_neg = sapply(ts, function(t)
      sum(W[, j] <= -t))
    
    ok = union(which(ratio <= alpha_nominal), which(sum_neg == 0))
    
    thresholds[j] = ifelse(length(ok) > 0, ts[ok[1]], Inf)
  }
  
  R_mat = (W >= matrix(thresholds, nrow = p, ncol = k, byrow = TRUE))
  U_vec = 1+colSums((W <= matrix(-thresholds, nrow = p, ncol = k, byrow = TRUE)))
  E_avg = rowSums(p * R_mat / matrix(U_vec, nrow = p, ncol = k, byrow = TRUE))/k
  
  E_avg_sort = sort(E_avg, decreasing = TRUE)
  
  idx_star = max(which((1:p)*E_avg_sort >= p/alpha_ebh), default = 0)
  
  R = which(E_avg >= p/(alpha_ebh*idx_star))
  
  # Returns the data-dependent alpha and the rejection set R
  return(list(alpha_dep = alpha_ebh, R = R))
}




#  KO_FDR_derand_posthoc(W, alpha_nominal, alpha_ebh)
#  --------------------------------------------------------------
#  Implements our posthoc derandomized knockoff procedure
#  - Input:
#      W : matrix of W-statistics (rows = features, columns = knockoff runs)
#     alpha_nominal : knockoff level alpha
#        alpha_ebh : initial significance level for the eBH step
#   - Output:
#      A list with:
#        • alpha_dep : alpha used in the eBH step.
#        • R : indices of selected features.

KO_FDR_derand_posthoc <- function(W, FDR_nominal, alpha_ebh = 0){
  k = ncol(W)
  p = nrow(W)
  thresholds = rep(0, k)
  for(j in 1:k){
    ts = sort(abs(W[, j]))
    ts = ts[ts>0]
    ratio = sapply(ts, function(t)
      (1 + sum(W[, j] <= -t)) / max(1, sum(W[, j] >= t)))
    sum_neg = sapply(ts, function(t)
      sum(W[, j] <= -t))
    
    ok = union(which(ratio <= FDR_nominal), which(sum_neg == 0))
    
    thresholds[j] = ifelse(length(ok) > 0, ts[ok[1]], Inf)
  }
  
  R_mat = (W >= matrix(thresholds, nrow = p, ncol = k, byrow = TRUE))
  U_vec = 1+colSums((W <= matrix(-thresholds, nrow = p, ncol = k, byrow = TRUE)))
  E_avg = rowSums(p * R_mat / matrix(U_vec, nrow = p, ncol = k, byrow = TRUE))/k
  
  E_avg_sort = sort(E_avg, decreasing = TRUE)
  
  idx_ebh = max(which((1:p)*E_avg_sort >= p/alpha_ebh), default = 0)
  
  R = which(E_avg >= p/(alpha_ebh*idx_ebh))
  
  if(length(R) > 0){
    alpha_dep = p/(idx_ebh * E_avg_sort[idx_ebh])
  }else{
    
    idx_star = max(which(((1:p)*E_avg_sort) == max((1:p)*E_avg_sort)))
    
    if(idx_star * E_avg_sort[idx_star] < p){
      alpha_dep = FDR_nominal
      R = c()
    }else{
      alpha_dep = p/(idx_star * E_avg_sort[idx_star])
      R = which(E_avg >= E_avg_sort[idx_star])
    }
  }
  
  # Returns the data-dependent alpha and the rejection set R
  return(list(alpha_dep = alpha_dep, R = R, R_mat = R_mat))
}
