# FDR standard
KO_FDR <- function(W, FDR_nominal){
  ts = sort(c(0, abs(W)) )
  ts = ts[ts>0]
  ratio = sapply(ts, function(t)
    (1 + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= FDR_nominal)
  threshold = ifelse(length(ok) > 0, ts[ok[1]], Inf)
  S = sort(which(W >= threshold))
  return(S)
}

KO_FDR_posthoc <- function(W, alpha_nominal, stopping_time="ph"){
  
  ts = sort(abs(W))
  ts = ts[ts>0]
  ratio = sapply(ts, function(t)
    (1 + sum(W <= -t)) / max(1, sum(W >= t)))
  sum_neg = sapply(ts, function(t)
    sum(W <= -t))
  sum_pos = sapply(ts, function(t)
    sum(W >= t))
  
  if(stopping_time=="RB"){
    ok = union(which(ratio <= alpha_nominal), which(sum_pos < 1/alpha_nominal))
  }else{
    ok = union(which(ratio <= alpha_nominal), which(sum_neg == 0))
  }
  
  threshold = ifelse(length(ok) > 0, ts[ok[1]], Inf)
  
  R = which(W >= threshold)
  
  if(length(R)>0){
    U <- which(W<=-threshold)
    alpha_dep <- (1+length(U))/length(R)
  }else{
    alpha_dep <- alpha_nominal
  }
  
  if(alpha_dep>1){
    R = c()
    alpha_dep = alpha_nominal
  }
  # Returns the data-dependent alpha and the rejection set R
  return(list(alpha_dep = alpha_dep, R = R))
}

# Function to control PFER
KO_select_PFER <- function(W, PFER_nominal){
  order_w <- order(abs(W),decreasing = TRUE)
  sorted_w <- W[order_w]
  negid <- which(sorted_w<0)
  S <- c()
  if(PFER_nominal> length(negid)){ 
    #if the total number of negitives is less than PFER_nominal, select all positives
    S <-  which(W>0)
  }else{
    TT <- negid[PFER_nominal]
    S <- which(sorted_w[1:TT]>0)
    S <- order_w[S]
  }
  return(S)
}


# Function for FDR control with derandomized knockoffs
# W is now a matrix (one column for every knockoff run)

KO_FDR_derand <- function(W, FDR_nominal, alpha_ebh){
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
  
  idx_star = max(which((1:p)*E_avg_sort >= p/alpha_ebh), default = 0)
  
  R = which(E_avg >= p/(alpha_ebh*idx_star))
  
  # Returns the data-dependent alpha and the rejection set R
  return(list(alpha_dep = alpha_ebh, R = R, R_mat = R_mat))
}



# Function for post-hoc FDR control with derandomized knockoffs
# W is now a matrix (one column for every knockoff run)

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
  E_threshold <- p/(alpha_ebh*idx_ebh)
  
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
      E_threshold <- E_avg_sort[idx_star]
    }
  }
  
  # Returns the data-dependent alpha and the rejection set R
  return(list(alpha_dep = alpha_dep, R = R, R_mat = R_mat, E_avg_sort = E_avg_sort, E_threshold = E_threshold))
}


# Function for post-hoc FDR control with derandomized knockoffs
# W is now a matrix (one column for every knockoff run)

KO_FDR_derand_posthoc_given_rejection_set <- function(W, FDR_nominal, rejection_set){
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
  
  level_to_return = p/(length(rejection_set)*min(E_avg[rejection_set]))

  # Returns the data-dependent alpha and the rejection set R
  return(list(alpha_dep = level_to_return, FDR_nominal = FDR_nominal, rejection_set = rejection_set))
}

# Function for PFER control with derandomized knockoffs                  
KO_PFER_derand <- function(W, PFER_nominal, eta = 0.5){
  k = ncol(W)
  p = nrow(W)
  thresholds = rep(0, k)
  for(j in 1:k){
    ts = sort(abs(W[, j]))
    ts = ts[ts>0]
    sum_neg = sapply(ts, function(t)
      sum(W[, j] <= -t))
    
    ok = which(sum_neg <= PFER_nominal-1)
    
    thresholds[j] = ifelse(length(ok) > 0, ts[ok[1]], Inf)
  }
  
  R_mat = (W >= matrix(thresholds, nrow = p, ncol = k, byrow = TRUE))
  
  R_sums = rowSums(R_mat)
  
  R = which(R_sums >= (eta * ncol(W)))
  
  # Returns the data-dependent eta and the rejection set R
  return(list(PFER_nominal = PFER_nominal, eta_dep = eta, R = R, R_mat = R_mat))
}

# Function for post-hoc PFER control with derandomized knockoffs
KO_PFER_derand_posthoc <- function(W, PFER_nominal){
  k = ncol(W)
  p = nrow(W)
  thresholds = rep(0, k)
  for(j in 1:k){
    ts = sort(abs(W[, j]))
    ts = ts[ts>0]
    sum_neg = sapply(ts, function(t)
      sum(W[, j] <= -t))
    
    ok = which(sum_neg <= PFER_nominal-1)
    
    thresholds[j] = ifelse(length(ok) > 0, ts[ok[1]], Inf)
  }
  
  R_mat = (W >= matrix(thresholds, nrow = p, ncol = k, byrow = TRUE))
  
  R_sums = rowSums(R_mat)
  
  R_times_eta = sapply(1:ncol(W), function(i) sum(R_sums >= i)) * (1:ncol(W)) / ncol(W)
  
  eta_dep = max(which(R_times_eta == max(R_times_eta))) / ncol(W)
  
  R = which(R_sums >= (eta_dep * ncol(W)))
  
  if(length(R) == 0){
    eta_dep = 0.5
  }
  
  
  # Returns the data-dependent eta and the rejection set R
  return(list(PFER_nominal = PFER_nominal, eta_dep = eta_dep, R = R, R_mat = R_mat))
}


#--------------------------------------------------------------
# Function: evenly_spaced_indices
# Purpose:
#   Return indices of p_relevant nonzero elements evenly spaced
#   among p total features, leaving roughly equal gaps of zeros.
#
# Example:
#   evenly_spaced_indices(600, 50)
#   -> returns indices spaced about 12 apart: 12, 24, 36, ...
#--------------------------------------------------------------

evenly_spaced_indices <- function(p, p_relevant) {
  if (p_relevant <= 0) stop("p_relevant must be positive")
  if (p_relevant > p) stop("p_relevant cannot exceed p")
  
  # Compute spacing
  step <- p / p_relevant
  
  # Start from the end of the first block (≈ step)
  indices <- round(seq(step, p, by = step))
  
  # Ensure exactly p_relevant indices (in case rounding changes count)
  indices <- indices[seq_len(p_relevant)]
  
  return(indices)
}

# Example:
evenly_spaced_indices(600, 50)
# [1]  12  24  36  48  60  72  84  96 108 120 132 ...





# generate design matrix randomly, taken for the conditional KO simulations
gene_X <- function(X_type = "IID_Normal", n, p, X_seed = NULL, rho = 0.50){
  if(!is.null(X_seed)){
    set.seed(X_seed)
  }
  
  X_type <- stringr::str_split(X_type, pattern = "_D_")[[1]][1]
  
  model_X <- F # set False if experiment with fixed-X
  if(model_X){
    basis <- matrix(rnorm(n*p), n)
  } else{
    basis <- qr.Q(qr(matrix(rnorm(n*p), n)))
  }
  
  cor_radius <- 5
  if(X_type == "IID_Normal"){
    cov_mat <- diag(p)
    X <- matrix(rnorm(n*p), n)
  } else if(X_type == "Coef_AR"){
   
    
    cov_mat <- solve(rho^(abs(outer(1:p, 1:p, "-"))))
    normalizer <- diag(1 / sqrt(diag(cov_mat)))
    cov_mat <- normalizer %*% cov_mat %*% normalizer
    
    R <- chol(cov_mat)
    X <- basis %*% R
  } else if(X_type == "X_AR"){
   
    
    cov_mat <- rho^(abs(outer(1:p, 1:p, "-")))
    normalizer <- diag(1 / sqrt(diag(cov_mat)))
    cov_mat <- normalizer %*% cov_mat %*% normalizer
    
    R <- chol(cov_mat)
    X <- basis %*% R
  } else if(X_type == "Homo_Block"){
   
    block_size <- 10
    
    blockSigma <- matrix(rho, block_size, block_size)
    diag(blockSigma) <- 1
    
    cov_mat <- as.matrix(diag(p / block_size) %x% blockSigma)
    normalizer <- diag(1 / sqrt(diag(cov_mat)))
    cov_mat <- normalizer %*% cov_mat %*% normalizer
    
    R <- chol(cov_mat)
    X <- basis %*% R
  } else if(X_type == "MCC"){
    if(n %% (p+1) == 0){
      X <- lapply(1:(n/(p+1)), function(i){
        rbind(diag(rep(1, p)), rep(0, p))
      })
      X <- do.call(rbind, X)
      X <- scale(X, center = T, scale = F)
      X <- scale(X, center = F, scale = sqrt(colSums(X^2)))
    } else{
      cov_mat <- matrix(-1/p, p, p)
      diag(cov_mat) <- 1
      
      R <- chol(cov_mat)
      X <- basis %*% R
    }
  } else if(X_type == "MCC_Block"){
    block_size <- 2
    
    blockSigma <- matrix(-1/block_size, block_size, block_size)
    diag(blockSigma) <- 1
    
    cov_mat <- as.matrix(diag(p / block_size) %x% blockSigma)
    
    R <- chol(cov_mat)
    X <- basis %*% R
  } else if(X_type == "Sparse"){
    sparsity <- 0.01
    X <- diag(1, nrow = n, ncol = p)
    lower_tri <- lower.tri(X)
    X[lower_tri] <- replicate(sum(lower_tri), rbinom(1, 1, sparsity))
  }
  # X <- scale(X, center = FALSE, scale = sqrt(colSums(X^2)))
  if(!exists("cov_mat")) cov_mat <- NA
  
  return(list(X = X, Xcov.true = cov_mat))
}


plot_selection_heatmap  <- function(R_mat, feature_names, title = "Selection Heatmap") {
  
  library(grid)
  library(pheatmap)
  
  # Assign names
  sel_matrix <- R_mat
  rownames(sel_matrix) <- feature_names
  colnames(sel_matrix) <- paste0("R", seq_len(ncol(R_mat)))
  
  # Convert logical → numeric
  sel_matrix_num <- apply(sel_matrix, 2, as.numeric)
  rownames(sel_matrix_num) <- rownames(sel_matrix)
  colnames(sel_matrix_num) <- colnames(sel_matrix)
  
  # Frequency counts
   row_freq <- rowSums(sel_matrix_num)
    col_freq <- colSums(sel_matrix_num)
  
  # Order by frequency
   sel_matrix_ordered <- sel_matrix_num[order(row_freq, decreasing = TRUE),
                                         order(col_freq, decreasing = TRUE)]
  
  # Plot heatmap
  p <- pheatmap::pheatmap(
    sel_matrix_ordered,
    color = c("white", "navy"),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    # fontsize_row = 20,
    #  fontsize_col = 20,
    #  cellheight = 15.5, 
    main = title,
    legend_breaks = c(0, 1),
    legend_labels = c("Not selected", "Selected"),
    legend = FALSE
  )
  
  # Return both matrix and plot
  return(p)
}
