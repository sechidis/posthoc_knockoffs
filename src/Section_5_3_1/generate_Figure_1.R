# Before running this script, you must first generate the simulated outcomes.
# To do so, run the following bash file in the terminal:
#   simulations_PH_vs_original_lsf.sh

library(ggplot2)
library(scales)
library(cowplot)
library(dplyr)

model="gaussian"
p_relevant_fixed = 6
amplitude_fixed = 8

outdir <- file.path("./results/", model)
files <- list.files(outdir, pattern = "result_.*\\.rds$", full.names = TRUE)
results <- dplyr::bind_rows(lapply(files, readRDS))
saveRDS(results, paste0("./results/results_PH_vs_original_",model,".rds"))

# ------------------------------------------------------------
# Summaries across iterations
# ------------------------------------------------------------
summary_df <- results |>
  dplyr::group_by(p_relevant, p, n, amplitude, method, model) |>
  dplyr::summarize(
    avg_R      = mean(R),
    avg_TD     = mean(TD),
    avg_FD     = mean(FD),
    avg_FDP    = mean(FDP),
    avg_power  = mean(power),
    avg_alpha  = mean(alpha_values, na.rm = TRUE),
    exp_ratio  = mean(ifelse(alpha_values > 0, FDP / alpha_values, NA_real_), na.rm = TRUE),
    .groups = "drop"
  )

print(summary_df)


nice_theme <- theme(
  legend.position = "top",
  axis.title = element_text(size = 24, face = "bold"),
  axis.title.x = element_text(size = 24),
  axis.title.y = element_text(size = 24),
  axis.text = element_text(size = 24),
  legend.title = element_blank(),
  legend.text = element_text(size = 24),
  panel.grid.major = element_line(size = 0.3),
  panel.grid.minor = element_blank(),
  text = element_text(family = "Times")
)

nice_guides <- guides(
  linetype = guide_legend(keywidth = 3, keyheight = 1.5, ncol = 1), # stacked vertically
  color = guide_legend(keywidth = 3, keyheight = 1.2),
  shape = guide_legend(keywidth = 3, keyheight = 1.2)
)


################################
### Plots against p_relevant ###
################################
variable_to_plot_x = 'p_relevant'
# Filter results to plut 
summary_df_filtered <- summary_df  |>
  dplyr::filter(amplitude == amplitude_fixed &
                  p_relevant>=3 & p_relevant<=8  )


p_relevant_power <- ggplot(
  summary_df_filtered ,
  aes(x = !!sym(variable_to_plot_x), y = avg_power,
      color = method, shape = method)
) +
  geom_line(size = 2) +
  geom_point(size = 6) +
  scale_color_manual(
    values = c("Original" = "#E41A1C", "Posthoc" = "#377EB8"),
    labels = c(expression("Original KO"), expression("Posthoc KO"))
  ) +
  scale_shape_manual(
    values = c("Original" = 16, "Posthoc" = 17),
    labels = c(expression("Original KO"), expression("Posthoc KO"))
  ) +
  labs(
    x = expression(p["relevant"]),
    y = expression("Power"),
    color = expression("Method"),
    shape = expression("Method")
  ) +
  ylim(0, 1) +
  scale_x_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)) +
  theme_bw() +
  nice_theme +
  nice_guides 

print(p_relevant_power)




p_relevant_ratio <- ggplot(
  summary_df_filtered ,
  aes(x = !!sym(variable_to_plot_x), y = exp_ratio,
      color = method, shape = method)
) +
  geom_line(size = 2) +
  geom_point(size = 6) +
  scale_color_manual(
    values = c("Original" = "#E41A1C", "Posthoc" = "#377EB8"),
    labels = c(expression("Original KO"), expression("Posthoc KO"))
  ) +
  scale_shape_manual(
    values = c("Original" = 16, "Posthoc" = 17),
    labels = c(expression("Original KO"), expression("Posthoc KO"))
  ) +
  labs(
    x = expression(p["relevant"]),
    y = expression("Average ratio FDP/" * alpha),
    color = expression("Method"),
    shape = expression("Method")
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),  # ticks every 0.2
    limits = c(0, 1.1),              # ensure axis includes 1
    expand = expansion(mult = 0)   # avoid extra padding
  ) + 
  scale_x_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)) +
  theme_bw() +
  nice_theme + 
  nice_guides 


print(p_relevant_ratio)

# Combined FDP and alpha plot
plot_df <- summary_df_filtered  |>
  tidyr::pivot_longer(
    cols = c(avg_FDP, avg_alpha),
    names_to = "metric",
    values_to = "value"
  ) |>
  dplyr::mutate(metric = dplyr::recode(metric, avg_FDP = "FDR", avg_alpha = "alpha"))

# Base plot with all aesthetics
base_plot <- ggplot(plot_df, aes(
  x = !!sym(variable_to_plot_x),
  y = value,
  color = method,
  shape = method,
  linetype = metric
)) +
  geom_line(size = 2) +
  geom_point(size = 6) +
  scale_color_manual(values = c("Original" = "#E41A1C", "Posthoc" = "#377EB8"),
                     labels = list(expression("Original KO"), expression("Posthoc KO"))) +
  scale_shape_manual(values = c("Original" = 16, "Posthoc" = 17),
                     labels = list(expression("Original KO"), expression("Posthoc KO"))) +
  scale_linetype_manual(values = c("alpha" = "solid", "FDR" = "dotted"),
                        labels = list(expression(alpha), expression("Average FDP"))) +
  labs( x = expression(p["relevant"]),,
        y = expression("FDR"),
        color = expression("Method"),
        shape = expression("Method"),
        linetype = expression("Metric")) +
  theme_bw(base_size = 16) + 
  scale_x_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)) +
  nice_theme + nice_guides

# Extract legends
method_legend <- get_legend(base_plot + guides(linetype = "none"))
metric_legend <- get_legend(base_plot + guides(color = "none", shape = "none"))

# Main plot without legends
main_plot <- base_plot + theme(legend.position = "none")

# Combine plot and legends using cowplot
p_relevant_FDR <- plot_grid(
  method_legend,
  ggdraw() + draw_plot(main_plot) + draw_plot(metric_legend, x = 0.65, y = 0.7, width = 0.2, height = 0.2),
  ncol = 1, rel_heights = c(0.1, 1)
)

print(p_relevant_FDR)


# ------------------------------------------------------------
# Save plots
# ------------------------------------------------------------
dir.create("./figures", showWarnings = FALSE)

ggsave(
  filename = paste0("intro_power.pdf") ,
  plot = p_relevant_power,
  device = "pdf",
  path = "./figures_intro/",
  width = 6.5, height = 5.0, units = "in",  scale = 1
)
ggsave(
  filename = paste0("intro_FDR.pdf") ,
  plot = p_relevant_FDR,
  device = "pdf",
  path = "./figures_intro/",
  width = 6.5, height = 5.0, units = "in",  scale = 1
)

table(results[results$p_relevant == '3' & results$amplitude ==amplitude_fixed & results$method =="Original",'R'])
table(results[results$p_relevant == '3' & results$amplitude ==amplitude_fixed & results$method =="Posthoc",'R'])
