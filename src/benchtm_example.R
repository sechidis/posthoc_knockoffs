library(dplyr)
library(pheatmap)
library(grid)
library(parallel)

# We use benchtm to generate data that mimic clinical trials.
# Details on installing this package can be found here https://github.com/Sophie-Sun/benchtm
library(benchtm) 

# We used knockofftools to generate knockoffs when we have mixed variable types.
# Details on installing this package can be found here https://github.com/Novartis/knockofftools/tree/main
library(knockofftools)

source("./utils.R", echo=TRUE)
set.seed(123)
FDR_nominal = 0.20
### 1. Load scenario and filter case
data(scen_param)

cases <- scen_param %>% 
  filter(type == "continuous",
         pred == "pnorm(20*(X11-0.5))",
         b1_rel == 1)

R=10 # number of runs

W_list <- vector("list", R)   # preallocate

for (i in 1:R) {
  
  dat_1 <- generate_scen_data(scen = cases, include_truth = FALSE) %>%
    mutate_if(is.character, as.factor) %>%
    mutate(trt = as.factor(trt))
  dat_2 <- generate_scen_data(scen = cases, include_truth = FALSE) %>%
    mutate_if(is.character, as.factor) %>%
    mutate(trt = as.factor(trt))
  
  dat <- rbind(dat_1, dat_2)
  
  X <- dat %>% select(starts_with("X"))
  T <- dat$trt
  Y <- dat$Y
  X_all <- cbind(X, T)
  
  X_all_k <- knockofftools::knockoffs_seq(X_all)
  W <- knockofftools::stat_glmnet(y = Y, X = X_all, X_k = X_all_k)
  
  W_list[[i]] <- W
}

# Convert list to matrix
W <- do.call(cbind, W_list)

# Run the original and post-hoc KO
R_original <- matrix(0, nrow = nrow(W), ncol = ncol(W))
R_ph <- matrix(0, nrow = nrow(W), ncol = ncol(W))
alpha_ph  <- matrix(0, nrow = 1, ncol = ncol(W))
for (m in seq(1,M,1)){
  R_original[KO_FDR(W[,m], FDR_nominal),m] = 1 
  R_ph[KO_FDR_posthoc(W[,m], FDR_nominal)$R,m] = 1
  alpha_ph[m] = KO_FDR_posthoc(W[,m], FDR_nominal)$alpha_dep  
}

# Plot the selections
# Assign row names
rownames(R_ph) <- names(X_all)
# Create a color palette: white -> blue
my_colors <- colorRampPalette(c("white", "blue"))(256)

# Plot selections for original
heatmap(R_original, 
        Rowv = NA,          # do not reorder rows
        Colv = NA,          # do not reorder columns
        col = my_colors,    # custom color palette
        scale = "none",     # don't scale the values
        margins = c(8, 5))  # increase left margin for row names

# Plot selections for posthoc
heatmap(R_ph, 
        Rowv = NA,          # do not reorder rows
        Colv = NA,          # do not reorder columns
        col = my_colors,    # custom color palette
        scale = "none",     # don't scale the values
        margins = c(8, 5))  # increase left margin for row names
