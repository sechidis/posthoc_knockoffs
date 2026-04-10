#!/usr/bin/env Rscript

# ============================================================
# Parallel simulation study for knockoff methods (LSF array version)
# ============================================================

suppressPackageStartupMessages({
  library(knockoff)
  library(tibble)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(scales)
  library(digest)
})

source("../utils.R")
# ------------------------------------------------------------
# User-defined parameters
# ------------------------------------------------------------
M <- 100
rho <- 0.5
PFER_nominal <- 1


params = commandArgs(trailingOnly=TRUE) # take seeds index
iter <- params[1]
print(iter)

models=c("gaussian", "logistic")
for (model in models){ 
  p_relevant = c(5,10)#seq(1,15,1)
  p          = 50
  n          = 250
  # ------------------------------------------------------------
  # Define simulation grid
  # ------------------------------------------------------------
  if (model == "gaussian") {
    amplitude  = seq(2, 6, 0.50)
  } else {
    amplitude  = seq(4, 12, 1)
  }
  
  cases <- expand.grid(p_relevant = p_relevant,p = p,n = n,amplitude =amplitude)
  
  
  
  # 1) Initialize an empty tibble with the correct columns & types
  results <- tibble::tibble(
    n          = integer(),
    p          = integer(),
    p_relevant = integer(),
    amplitude  = numeric(),
    model      = character(),
    iter  = character(),
    seed       = integer(),
    method     = factor(levels = c("Original", "Posthoc")),
    PFER_nominal = numeric(),
    R          = integer(),
    TD         = integer(),
    FD         = integer(),
    FDP        = numeric(),
    power      = numeric(),
    eta_values = numeric()
  )
  
  for (cases_index in seq(1,dim(cases)[1],1)){
    p_relevant = cases[cases_index,'p_relevant']
    p  = cases[cases_index,'p']
    n  = cases[cases_index,'n']
    amplitude  = cases[cases_index,'amplitude']
    
    
    # ------------------------- Seed generation -------------------------
    seed_digest <- substr(digest::digest(
      paste(p_relevant, p, n, iter, amplitude, sep = "_"),
      algo = "xxhash32"
    ), 1, 7) # shorter to avoid NA
    
    seed_i <- strtoi(seed_digest, base = 16L)
    if (is.na(seed_i)) stop("❌ Seed generation failed for job digest=", seed_digest)
    
    set.seed(seed_i)
    cat("✅ Using seed:", seed_i, "\n")
    
    # ------------------------- Core simulation function -------------------------
    run_one <- function(amplitude, model, rho, n, p, M, p_relevant, PFER_nominal, seed = NULL, iter) {
      if (!is.null(seed)) set.seed(seed)
      model <- match.arg(model, choices = c("gaussian", "logistic"))
      
      # --- 1. Checks and setup ---
      stopifnot(p > 1, p_relevant >= 0, p_relevant <= p)
      if (p_relevant == 0) warning("p_relevant = 0: no true signals.")
      
      # --- 2. Covariance and data ---
      Sigma <- toeplitz(rho^(0:(p - 1)))
      X <- matrix(rnorm(n * p), n, p) %*% chol(Sigma)
      
      # --- 3. Coefficients ---
      beta_true <- rep(0, p)
      nonzero   <- evenly_spaced_indices(p, p_relevant)
      sign_loc  <- if (length(nonzero) > 1) nonzero[seq(2, length(nonzero), 2)] else integer(0)
      
      beta_true[nonzero] <- rnorm(p_relevant, mean = amplitude, sd = 1) / sqrt(n)
      beta_true[sign_loc] <- -beta_true[sign_loc]
      
      # --- 4. Response ---
      if (model == "gaussian") {
        y <- as.vector(X %*% beta_true + rnorm(n))
      } else {
        linpred <- X %*% beta_true
        prob <- 1 / (1 + exp(-linpred))
        y <- rbinom(n, 1, prob)
      }
      
      # --- 5. Knockoff statistics ---
      W <- matrix(NA, nrow = p, ncol = M)
      for (M_index in seq_len(M)) {
        Xk <- knockoff::create.second_order(X)
        W[, M_index] <- knockoff::stat.glmnet_coefdiff(X, Xk, y)
      }
      
      
      res_trad <- KO_PFER_derand(W, PFER_nominal = PFER_nominal, eta = 0.50)
      res_post <- KO_PFER_derand_posthoc(W, PFER_nominal = PFER_nominal)
      
      # --- 6. Traditional knockoff ---
      sel_trad   <- res_trad$R
      R_trad     <- length(sel_trad)
      FD_trad    <- sum(!(sel_trad %in% nonzero))
      TD_trad    <- sum(sel_trad %in% nonzero)
      FDP_trad   <- ifelse(R_trad > 0, FD_trad / R_trad, 0)
      power_trad <- ifelse(p_relevant > 0, TD_trad / p_relevant, 0)
      
      # --- 7. Posthoc knockoff ---
      sel_post   <- as.integer(res_post$R)
      alpha_post <- res_post$alpha_dep
      R_post     <- length(sel_post)
      FD_post    <- sum(!(sel_post %in% nonzero))
      TD_post    <- sum(sel_post %in% nonzero)
      FDP_post   <- ifelse(R_post > 0, FD_post / R_post, 0)
      power_post <- ifelse(p_relevant > 0, TD_post / p_relevant, 0)
      
      # --- 8. Results ---
      tibble::tibble(
        n = n,
        p = p,
        p_relevant = p_relevant,
        amplitude = amplitude,
        model = model,
        iter = iter,
        seed = seed,
        method = c("Original derandomised PFER", "Posthoc derandomised PFER"),
        PFER_nominal = c(res_trad$PFER_nominal, res_post$PFER_nominal),
        R = c(length(res_trad$R), length(res_post$R)),
        TD = c(sum(res_trad$R %in% nonzero), sum(res_post$R %in% nonzero)),
        FD = c(sum(!(res_trad$R %in% nonzero)), sum(!(res_post$R %in% nonzero))),
        FDP = c(ifelse(length(res_trad$R)>0, sum(!(res_trad$R %in% nonzero))/length(res_trad$R), 0),
                ifelse(length(res_post$R)>0, sum(!(res_post$R %in% nonzero))/length(res_post$R), 0)),
        power = c(sum(res_trad$R %in% nonzero)/length(nonzero), sum(res_post$R %in% nonzero)/length(nonzero)),
        eta_values = c(res_trad$eta_dep, res_post$eta_dep)
      )
    }
    
    # ------------------------- Run and save -------------------------
    pair <- run_one(
      model = model,
      rho = rho,
      n = n,
      p = p,
      amplitude = amplitude,
      M = M,
      p_relevant = p_relevant,
      PFER_nominal = PFER_nominal,
      seed = seed_i,
      iter = iter
    )
    
    
    # pair has two rows: Original, Posthoc
    results <- dplyr::bind_rows(results, pair)
    
  }
  outdir <- file.path("./results/", model)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  outfile <- file.path(outdir, paste0("result_", model,"_iter_",iter,".rds"))
  saveRDS(results, file = outfile)
  cat("✅ Saved:", outfile, "\n")
}
