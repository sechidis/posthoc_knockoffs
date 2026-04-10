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

# ------------------------- User parameters -------------------------
M <- 100
rho <- 0.5
FDR_nominal <- 0.1


models=c("gaussian", "logistic")
for (model in models){ 
  # ------------------------------------------------------------
  # Define simulation grid
  # ------------------------------------------------------------
  if (model == "gaussian") {
    jobs <- expand.grid(
      iter = 1:200,
      p_relevant = c(20, 40, 80),
      n = c(1000),
      p = c(800),
      amplitude = seq(2, 8, 1)
    )
  } else {
    jobs <- expand.grid(
      iter = 1:200,
      p_relevant = c(15, 30, 60),
      n = c(1000),
      p = c(600),
      amplitude = seq(4, 20, 2)
    )
    
  }
  # ------------------------- Command-line argument -------------------------
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop("Usage: Rscript simulations_derandomised_gaussian_lsf.R <job_index>")
  
  job_index <- as.integer(args[1])
  if (is.na(job_index) || job_index < 1 || job_index > nrow(jobs)) {
    stop("Invalid job index: ", job_index)
  }
  
  job <- jobs[job_index, ]
  
  cat("====================================================\n")
  cat("Running job index:", job_index, "\n")
  cat("Job parameters:\n")
  print(job)
  cat("====================================================\n")
  
  # ------------------------- Seed generation -------------------------
  seed_digest <- substr(digest::digest(
    paste(job$p_relevant, job$p, job$n, job$iter, job$amplitude, sep = "_"),
    algo = "xxhash32"
  ), 1, 7) # shorter to avoid NA
  
  seed_i <- strtoi(seed_digest, base = 16L)
  if (is.na(seed_i)) stop("❌ Seed generation failed for job ", job_index, " digest=", seed_digest)
  
  set.seed(seed_i)
  cat("✅ Using seed:", seed_i, "\n")
  
  # ------------------------- Core simulation function -------------------------
  run_one <- function(amplitude, model, rho, n, p, M, p_relevant, FDR_nominal, seed = NULL, iter) {
    if (!is.null(seed)) set.seed(seed)
    model <- match.arg(model, choices = c("gaussian", "logistic"))
    
    Sigma <- toeplitz(rho^(0:(p-1)))
    diags <- knockoff::create.solve_asdp(Sigma)
    X <- matrix(rnorm(n * p), n, p) %*% chol(Sigma)
    
    beta_true <- rep(0, p)
    nonzero <- evenly_spaced_indices(p, p_relevant)
    sign_loc <- nonzero[seq(2, length(nonzero), 2)]
    beta_true[nonzero] <- rnorm(p_relevant, mean = amplitude, sd = 1) / sqrt(n)
    beta_true[sign_loc] <- -beta_true[sign_loc]
    
    y <- if(model == "gaussian") X %*% beta_true + rnorm(n) else rbinom(n, 1, 1 / (1 + exp(-X %*% beta_true)))
    
    W <- matrix(NA, nrow = p, ncol = M)
    for (M_index in seq_len(M)) {
      Xk <- knockoff::create.gaussian(X, mu = rep(0,p), Sigma, diag_s = diags)
      W[, M_index] <- knockoff::stat.glmnet_coefdiff(X, Xk, y)
    }
    
    alpha_kn <- FDR_nominal / 2
    alpha_ebh <- FDR_nominal
    
    res_trad <- KO_FDR_derand(W, FDR_nominal = alpha_kn, alpha_ebh = alpha_ebh)
    res_post <- KO_FDR_derand_posthoc(W, FDR_nominal = alpha_kn, alpha_ebh = alpha_ebh)
    
    tibble(
      n = n,
      p = p,
      p_relevant = p_relevant,
      amplitude = amplitude,
      model = model,
      iter = iter,
      seed = seed,
      method = c("Original derandomised FDR", "Posthoc derandomised FDR"),
      R = c(length(res_trad$R), length(res_post$R)),
      TD = c(sum(res_trad$R %in% nonzero), sum(res_post$R %in% nonzero)),
      FD = c(sum(!(res_trad$R %in% nonzero)), sum(!(res_post$R %in% nonzero))),
      FDP = c(ifelse(length(res_trad$R)>0, sum(!(res_trad$R %in% nonzero))/length(res_trad$R), 0),
              ifelse(length(res_post$R)>0, sum(!(res_post$R %in% nonzero))/length(res_post$R), 0)),
      power = c(sum(res_trad$R %in% nonzero)/length(nonzero), sum(res_post$R %in% nonzero)/length(nonzero)),
      alpha_values = c(res_trad$alpha_dep, res_post$alpha_dep)
    )
  }
  
  # ------------------------- Run and save -------------------------
  res <- run_one(
    model = model,
    rho = rho,
    M = M,
    n = job$n,
    p = job$p,
    amplitude = job$amplitude,
    p_relevant = job$p_relevant,
    FDR_nominal = FDR_nominal,
    seed = seed_i,
    iter = job$iter
  )
  
  outdir <- file.path("./results/", model)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  outfile <- file.path(outdir, sprintf("result_%s_job%04d.rds", model, job_index))
  saveRDS(res, file = outfile)
  cat("✅ Saved:", outfile, "\n")
  
}