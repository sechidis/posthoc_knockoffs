library(benchtm)
library(dplyr)
library(pheatmap)
library(grid)
library(parallel)
library(knockofftools)

source("../utils.R", echo=TRUE)


set.seed(12)
FDR_nominal = 0.20
### 1. Load scenario and filter case
data(scen_param)

cases <- scen_param %>% 
  filter(type == "continuous",
         pred == "pnorm(20*(X11-0.5))",
         b1_rel == 1)

dat_1 <- generate_scen_data(scen = cases, include_truth = FALSE) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(trt = as.factor(trt))
dat_2 <- generate_scen_data(scen = cases, include_truth = FALSE) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(trt = as.factor(trt))

dat <- rbind(dat_1,dat_2)

X <- dat %>% select(starts_with("X"))
T <- dat$trt
Y <- dat$Y
X_all <- cbind(T,X)


M=10
# Number of cores to use
ncores <- detectCores() - 1  # leave 1 core free
print(paste("Using", ncores, "cores"))

W_list <- mclapply(1:M, function(i) {
  set.seed(2+i)
  dat_1 <- generate_scen_data(scen = cases, include_truth = FALSE) %>%
    mutate_if(is.character, as.factor) %>%
    mutate(trt = as.factor(trt))
  dat_2 <- generate_scen_data(scen = cases, include_truth = FALSE) %>%
    mutate_if(is.character, as.factor) %>%
    mutate(trt = as.factor(trt))
  dat <- rbind(dat_1,dat_2)
  
  X <- dat %>% select(starts_with("X"))
  T <- dat$trt
  Y <- dat$Y
  X_all <- cbind(T,X)
  
  X_all_k <- knockofftools::knockoffs_seq(X_all)
  W <- knockofftools::stat_glmnet(y = Y, X = X_all, X_k = X_all_k)
  return(W)
}, mc.cores = ncores)
# Convert list to matrix
W <- do.call(cbind, W_list)  # each row = one iteration W

# Run the original and post-hoc KO
R_original <- matrix(FALSE, nrow = nrow(W), ncol = ncol(W))
R_ph <- matrix(0, nrow = nrow(W), ncol = ncol(W))
alpha_ph  <- matrix(0, nrow = 1, ncol = ncol(W))
for (m in seq(1,M,1)){
  R_original[KO_FDR(W[,m], FDR_nominal),m] = TRUE 
  R_ph[KO_FDR_posthoc(W[,m], FDR_nominal)$R,m] = TRUE
  alpha_ph[m] = KO_FDR_posthoc(W[,m], FDR_nominal)$alpha_dep  
}


pdf("./figures/benchtm_ph.pdf", width = 6.3, height = 8.1)
mean_val <- round(mean(alpha_ph), 2)
plot_selection_heatmap(R_ph, names(X_all), title = "") 
grid.text(label = bquote("Selection Heatmap with ph-KO and average " * alpha^ph * " = " * .(mean_val)),
          x = 0.5, y = 0.975, gp = gpar(fontsize = 18, fontface = "bold" ))
dev.off()


pdf("./figures/benchtm_original.pdf", width = 6.3, height = 8.1)
p_original <- plot_selection_heatmap(R_original, names(X_all), title = "")
grid.text(label = bquote("Selection Heatmap with original KO and " * alpha^kn * " = " * .(FDR_nominal)),
          x = 0.5, y = 0.975, gp = gpar(fontsize = 18, fontface = "bold" ))
dev.off()



plot_selection_heatmap <- function(R_mat, feature_names, title = "Selection Heatmap") {
  
  library(grid)
  library(pheatmap)
  
  # Assign names
  sel_matrix <- R_mat
  rownames(sel_matrix) <- feature_names
  colnames(sel_matrix) <- paste0("RCT ", seq_len(ncol(R_mat)))
  
  # Convert logical → numeric
  sel_matrix_num <- apply(sel_matrix, 2, as.numeric)
  rownames(sel_matrix_num) <- rownames(sel_matrix)
  colnames(sel_matrix_num) <- colnames(sel_matrix)
  
  # Frequency counts
  # row_freq <- rowSums(sel_matrix_num)
  #  col_freq <- colSums(sel_matrix_num)
  
  # Order by frequency
   # sel_matrix_ordered <- sel_matrix_num[order(row_freq, decreasing = TRUE),
  #                                       order(col_freq, decreasing = TRUE)]
  
  # Plot heatmap
  p <- pheatmap::pheatmap(
    sel_matrix_num,
    color = c("white", "navy"),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 14,
  #  fontsize_col = 20,
    cellheight = 16, 
  cellwidth = 40,
    main = title,
    legend_breaks = c(0, 1),
    legend_labels = c("Not selected", "Selected"),
    legend = FALSE
  )
  
  # Return both matrix and plot
  return(p)
}
