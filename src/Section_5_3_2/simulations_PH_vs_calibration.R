# ============================================================
# Parallel simulation study for knockoff methods
# ============================================================

library(parallel)
library(knockoff)
library(cknockoff)
library(tibble)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(digest)

set.seed(123)

# ------------------------------------------------------------
# User-defined parameters
# ------------------------------------------------------------
X_types <- c("IID_Normal")
models=c("gaussian","logistic")
rho         <- 0.5
FDR_nominal <- 0.2


source("../utils.R")  # ensure KO_FDR and KO_FDR_posthoc are defined there


# ------------------------------------------------------------
# Define simulation grid
# ------------------------------------------------------------
for (X_type in X_types){ 
  iter       = 1:2000
  p_relevant = c(5,10,15)
  p          = 100
  n          = 500

  for (model in models){ 
    
    # ------------------------------------------------------------
    # Define simulation grid
    # ------------------------------------------------------------
    if (model == "gaussian") {
        amplitude  = seq(2, 6, 0.50)
    } else {
      amplitude  = seq(4, 12, 1)
    }
    
    jobs <- expand.grid(iter = iter,p_relevant = p_relevant,p = p,n = n,amplitude =amplitude)
    
    # ------------------------------------------------------------
    # Simulation function
    # ------------------------------------------------------------
    run_one <- function(p, p_relevant, amplitude,
                        model = c("gaussian", "logistic"),
                        rho = 0.5, n,
                        FDR_nominal = 0.1, seed = NULL) {
      
      if (!is.null(seed)) set.seed(seed)
      model <- match.arg(model)
      
      # --- 1. Checks and setup ---
      stopifnot(p > 1, p_relevant >= 0, p_relevant <= p)
      if (p_relevant == 0) warning("p_relevant = 0: no true signals.")
      
      # --- 2. Covariance and data ---
      #Sigma <- toeplitz(rho^(0:(p - 1)))
      #X <- matrix(rnorm(n * p), n, p) %*% chol(Sigma)
      
      X_Sigma <- gene_X(X_type = X_type, n, p, rho = rho)
      X <- X_Sigma$X    
      Sigma <- X_Sigma$Xcov.true
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
      kn.result <- knockoff.filter(X, y,
                                   knockoffs = cknockoff::ckn.create.fixed,
                                   statistic = cknockoff::stat.glmnet_coefdiff_tiebreak, # must specify this argument explicitly
                                   fdr = FDR_nominal # must specify this argument explicitly
      )
      
      
      
      # --- 6. Calibration knockoff ---
      sel_calibration  <- cknockoff::cknockoff(prelim_result = kn.result)$selected
      R_calibration     <- length(sel_calibration)
      FD_calibration    <- sum(!(sel_calibration %in% nonzero))
      TD_calibration    <- sum(sel_calibration %in% nonzero)
      FDP_calibration   <- ifelse(R_calibration > 0, FD_calibration / R_calibration, 0)
      power_calibration <- ifelse(p_relevant > 0, TD_calibration / p_relevant, 0)
      
      # --- 7. Posthoc knockoff ---
      res_post   <- KO_FDR_posthoc(kn.result$statistic, FDR_nominal)
      sel_post   <- as.integer(res_post$R)
      alpha_post <- res_post$alpha_dep
      R_post     <- length(sel_post)
      FD_post    <- sum(!(sel_post %in% nonzero))
      TD_post    <- sum(sel_post %in% nonzero)
      FDP_post   <- ifelse(R_post > 0, FD_post / R_post, 0)
      power_post <- ifelse(p_relevant > 0, TD_post / p_relevant, 0)
      
      # --- 9. Results ---
      results <- tibble(
        n, p, p_relevant, amplitude, model,seed,
        method = c("Calibration", "Posthoc"),
        R       = c(R_calibration, R_post),
        TD      = c(TD_calibration, TD_post),
        FD      = c(FD_calibration, FD_post),
        FDP     = c(FDP_calibration, FDP_post),
        power   = c(power_calibration, power_post),
        alpha_values = c(FDR_nominal, alpha_post)
      )
      
      list(results = results)
    }
    
    # ------------------------------------------------------------
    # Parallel execution
    # ------------------------------------------------------------
    cl <- makeCluster(detectCores() - 1)
    
    environment(run_one) <- .GlobalEnv
    
    base_seed <- sample.int(.Machine$integer.max, 1)
    
    clusterExport(cl, varlist = c(
      "jobs", "run_one", "FDR_nominal", "evenly_spaced_indices",
      "KO_FDR", "KO_FDR_posthoc", "model", "rho", "base_seed", "gene_X", "X_type"
    ))
    clusterEvalQ(cl, {
      library(knockoff)
      library(tibble)
      library(digest)
    })
    
    
    results_list <- parLapply(cl, seq_len(nrow(jobs)), function(i) {
      p_relevant_i <- jobs$p_relevant[i]
     p_i          <- jobs$p[i]
      n_i          <- jobs$n[i]
      iter_i       <- jobs$iter[i]
      amplitude_i  <- jobs$amplitude[i]
      
      seed_i <- as.integer(
        (strtoi(substr(
          digest::digest(
            paste(p_relevant_i, p_i, n_i, iter_i, amplitude_i, sep = "_"),
            algo = "xxhash64"
          ), 1, 8), base = 16L
        ) + base_seed + i) %% .Machine$integer.max
      )
      
      if (is.na(seed_i)) seed_i <- sample.int(.Machine$integer.max, 1)
      
      out <- run_one(
        p = p_i, p_relevant = p_relevant_i, amplitude = amplitude_i,
        model = model, rho = rho, n = n_i, FDR_nominal = FDR_nominal, seed = seed_i
      )
      out$results
    })
    
    stopCluster(cl)
    
    # ------------------------------------------------------------
    # Combine and save
    # ------------------------------------------------------------
    results <- bind_rows(results_list)
    
    if (!dir.exists("./results")) dir.create("./results", recursive = TRUE)
    
    file_path <- file.path("./results", paste0("results_", X_type, "_", model, ".rds"))
    saveRDS(results, file = file_path)
    
    cat("✅ Results saved to:", file_path, "\n")
  }
}
