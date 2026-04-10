# Before running this script, you must first generate the simulated outcomes.
# To do so, run the following R file: simulations_PH_vs_calibration.sh


library(ggplot2)
library(scales)
library(cowplot)
library(dplyr)

X_types <- c("IID_Normal")
for (X_type in X_types){ 
  for (model in models){ 
    
    file_name <- paste0("results_",  X_type, "_",  model, ".rds")
    
    file_path <- file.path("./results/", file_name)
    results <- readRDS( file = file_path) |>
      mutate(method = factor(method,
                             levels = c("Calibration", "Posthoc")))
    
    
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
      axis.title = element_text(size = 30, face = "bold"),
      axis.title.x = element_text(size = 30),
      axis.title.y = element_text(size = 30),
      axis.text = element_text(size = 22),
      legend.title = element_blank(),
      legend.text = element_text(size = 22),
      panel.grid.major = element_line(size = 0.3),
      panel.grid.minor = element_blank(),
      text = element_text(family = "Times"),
      strip.text = element_text(face = "bold", size = 30),
      strip.background = element_rect(fill = "grey90", color = "grey70", linewidth = 0.8),
      panel.spacing = unit(2, "lines")
    )
    nice_guides <- guides(
      linetype = guide_legend(keywidth = 3, keyheight = 1.5, ncol = 1), # stacked vertically
      color = guide_legend(keywidth = 3, keyheight = 1.2),
      shape = guide_legend(keywidth = 3, keyheight = 1.2)
    )
    
    
    ###############################
    ### Plots against amplitude ###
    ###############################
    variable_to_plot_x = 'amplitude'
    # Filter results to plut 
    summary_df_filtered <- summary_df  |>
      dplyr::filter(p_relevant %in% c(10, 15)) |>
      mutate(p_relevant_lab = factor(p_relevant,
                                     levels = c(10, 15),
                                     labels = c("p[relevant]==10", "p[relevant]==15")))
    
    
    amplitude_power <- ggplot(
      summary_df_filtered ,
      aes(x = !!sym(variable_to_plot_x), y = avg_power,
          color = method, shape = method)
    ) +
      geom_line(size = 2) +
      geom_point(size = 6) +
      scale_color_manual(
        values = c("Calibration" = "#E41A1C",  "Posthoc" = "#377EB8"),
        labels = c(expression("Calibration KO"),  expression("Posthoc KO"))
      ) +
      scale_shape_manual(
        values = c("Calibration" = 16,  "Posthoc" = 17),
        labels = c(expression("Calibration KO"),  expression("Posthoc KO"))
      ) +
      labs(
        x = expression("Signal amplitude " * A),
        y = expression("Power"),
        color = expression("Method"),
        shape = expression("Method")
      ) +
      ylim(0, 1) +
      theme_bw() +
      facet_grid(~ p_relevant_lab, labeller = label_parsed) +
      nice_theme +
      nice_guides 
    
    print(amplitude_power)
    
    
    
    
    amplitude_ratio <- ggplot(
      summary_df_filtered ,
      aes(x = !!sym(variable_to_plot_x), y = exp_ratio,
          color = method, shape = method)
    ) +
      geom_line(size = 2) +
      geom_point(size = 6) +
      scale_color_manual(
        values = c("Calibration" = "#E41A1C",  "Posthoc" = "#377EB8"),
        labels = c(expression("Calibration KO"),  expression("Posthoc KO"))
      ) +
      scale_shape_manual(
        values = c("Calibration" = 16,  "Posthoc" = 17),
        labels = c(expression("Calibration KO"),  expression("Posthoc KO"))
      ) +
      labs(
        x = expression("Signal amplitude " * A),
        y = expression("Average ratio FDP/" * alpha),
        color = expression("Method"),
        shape = expression("Method")
      ) +
      scale_y_continuous(
        breaks = seq(0, 1, by = 0.2),  # ticks every 0.2
        limits = c(-.05, 1.1),              # ensure axis includes 1
        expand = expansion(mult = 0)   # avoid extra padding
      ) + 
      theme_bw() +
      facet_grid(~ p_relevant_lab, labeller = label_parsed) +
      nice_theme + 
      nice_guides 
    
    
    print(amplitude_ratio)
    
    # Combined FDP and alpha plot
    plot_df <- summary_df_filtered  %>%
      tidyr::pivot_longer(
        cols = c(avg_FDP, avg_alpha),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(metric = recode(metric, avg_FDP = "FDP", avg_alpha = "alpha"))
    
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
      scale_color_manual(values = c("Calibration" = "#E41A1C",  "Posthoc" = "#377EB8"),
                         labels = list(expression("Calibration KO"),  expression("Posthoc KO"))) +
      scale_shape_manual(values = c("Calibration" = 16,  "Posthoc" = 17),
                         labels = list(expression("Calibration KO"),   expression("Posthoc KO"))) +
      scale_linetype_manual(values = c("alpha" = "solid", "FDP" = "dotted"),
                            labels = list(expression(alpha), expression("Average " * FDP))) +
      labs(x = expression("Signal amplitude " * A),
           y = expression("FDR"),
           color = expression("Method"),
           shape = expression("Method"),
           linetype = expression("Metric")) +
      scale_y_continuous(
        breaks = seq(0, 0.40, by = 0.2),  # ticks every 0.2
        limits = c(0.1, 0.42),              # ensure axis includes 1
        expand = expansion(mult = 0)   # avoid extra padding
      ) + 
      theme_bw(base_size = 16) + 
      facet_grid(~ p_relevant_lab, labeller = label_parsed) +
      nice_theme + nice_guides
    
    # Extract legends
    method_legend <- get_legend(base_plot + guides(linetype = "none"))
    
    metric_legend <- get_legend(base_plot + guides(color = "none", shape = "none"))
    
    # Main plot without legends
    main_plot <- base_plot + theme(legend.position = "none")
    
    # Combine plot and legends using cowplot
    amplitude_FDR <- plot_grid(
      method_legend,
      ggdraw() + draw_plot(main_plot) + 
        draw_plot(metric_legend, x = 0.75, y = 0.57, width = 0.2, height = 0.2)+
        draw_plot(metric_legend, x = 0.28, y = 0.57, width = 0.2, height = 0.2),
      ncol = 1, rel_heights = c(0.1, 1)
    )
    
    print(amplitude_FDR)
    
    
    # ------------------------------------------------------------
    # Save plots
    # ------------------------------------------------------------
    dir.create("./figures", showWarnings = FALSE)
    
    ggsave(
    #  filename = paste0("results_",  X_type, "_",model,"_calibration_power.pdf") ,
      filename = paste0( model,"_calibration_power.pdf") ,
      plot = amplitude_power,
      device = "pdf",
      path = "./figures/",
      width = 12.5, height = 5.5, units = "in",  scale = 0.8
    )
    ggsave(
    #  filename = paste0("results_",  X_type, "_",model,"_calibration_FDR.pdf") ,
      filename = paste0(model,"_calibration_FDR.pdf") ,
      plot = amplitude_FDR,
      device = "pdf",
      path = "./figures/",
      width = 12.5, height = 5.5, units = "in",  scale = 0.8
    )
    
    ggsave(
    #  filename = paste0("results_",  X_type, "_",model,"_calibration_ratio.pdf") ,
      filename = paste0(model,"_calibration_ratio.pdf") ,
      plot = amplitude_ratio,
      device = "pdf",
      path = "./figures/",
      width = 12.5, height = 5.5, units = "in",  scale = 0.8
    )
  }
  
}
