##########################################
### pathreat.analysis.sfa.fig.R #########
### SFA figures (all panels) #############
##########################################

library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(sfaR)

source("code/pathreat.analysis.config.R")

cat("=== SFA Figures ===\n")

###############################
### load data #################
###############################

d_pa      <- readRDS(paths$est_sfa)
d_country <- readRDS(paths$est_sfa_country)
models    <- readRDS(paste0(paths$results_dir, "pathreat.sfa.models.est.Rds"))
g_global  <- readRDS(paths$est_global)

cat(sprintf("  PA-level: %s PAs\n", format(nrow(d_pa), big.mark = ",")))
cat(sprintf("  Country-level: %d countries\n", nrow(d_country)))

# total_pa_area already in d_country from estimation script

# global control mean for normalization
g_tc <- g_global[g_global$variable == "threat_composite", ]
ctrl_mean <- g_tc$control_mean

# global aggregates (pixel-weighted)
global_delta         <- weighted.mean(d_pa$delta, d_pa$n_pixels)
global_frontier      <- weighted.mean(d_pa$frontier, d_pa$n_pixels)
global_ctrl_mean     <- weighted.mean(d_pa$tc_control, d_pa$n_pixels)
global_delta_norm    <- global_delta / global_ctrl_mean * 100
global_frontier_norm <- global_frontier / global_ctrl_mean * 100
global_effect_norm   <- g_tc$coef / ctrl_mean * 100

# country-level normalized columns (each country's own control mean)
d_scatter_norm <- d_country %>%
  mutate(
    delta_norm    = delta_mean / tc_control_mean * 100,
    frontier_norm = frontier_mean / tc_control_mean * 100
  )


###############################
### shared theme & helpers ####
###############################

fig_theme <- theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  )

sfa_comparison_theme <- theme(
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 13),
  plot.margin = margin(5.5, 28, 5.5, 5.5, "pt")
)

size_scale <- scale_size_continuous(
  range = c(1.5, 8),
  name = expression("Total PA area (km"^2*")"),
  trans = scales::trans_new("log2", log2, function(x) 2^x),
  breaks = c(100, 1000, 10000, 100000),
  labels = scales::comma
)

# y-axis limits shared across country scatter panels (exclude Libya outlier)
d_scatter_ylim <- d_scatter_norm %>% filter(country_name != "Libya")
y_lo <- min(d_scatter_ylim$delta_norm, na.rm = TRUE)
y_hi <- max(d_scatter_ylim$delta_norm, na.rm = TRUE)
y_pad <- (y_hi - y_lo) * 0.03
y_lim <- c(y_lo - y_pad, y_hi + y_pad)

# helper: select outskirt points for labeling
label_outskirts <- function(df, x, y, top_n = 25) {
  xv <- df[[x]]
  yv <- df[[y]]
  xz <- (xv - median(xv, na.rm = TRUE)) / IQR(xv, na.rm = TRUE)
  yz <- (yv - median(yv, na.rm = TRUE)) / IQR(yv, na.rm = TRUE)
  df$outskirt_dist <- sqrt(xz^2 + yz^2)
  df[rank(-df$outskirt_dist, ties.method = "first") <= top_n, ]
}

# IUCN class labels and colors (for PA scatter)
iucn_labels <- c(
  "Strict protection" = "Strict (Ia\u2013IV)",
  "Less strict"       = "Multi-use (V\u2013VI)",
  "Not Reported"      = "Not Reported"
)
iucn_colors <- c(
  "Strict (Ia\u2013IV)"   = "dodgerblue3",
  "Multi-use (V\u2013VI)" = "darkorange2",
  "Not Reported"           = "grey55"
)


###########################################
### Fig 1: PA scatter + marginals #########
###########################################
{
  cat("  Fig 1: PA scatter + marginals\n")

  d_pa$iucn_label <- factor(iucn_labels[d_pa$iucn_class], levels = iucn_labels)

  pars_main <- coef(models$uhet, extraPar = TRUE)
  gamma_val <- pars_main["gamma"]

  ax_lo <- -0.22
  ax_hi <- 0.12

  p_main <- ggplot(d_pa, aes(x = frontier, y = delta)) +
    annotate("polygon",
             x = c(ax_lo, ax_hi, ax_hi, ax_lo),
             y = c(ax_lo, ax_hi, ax_hi + 0.5, ax_lo + 0.5),
             fill = "coral", alpha = 0.04) +
    annotate("polygon",
             x = c(ax_lo, ax_hi, ax_hi, ax_lo),
             y = c(ax_lo, ax_hi, ax_hi - 0.5, ax_lo - 0.5),
             fill = "steelblue", alpha = 0.04) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid",
                color = "grey30", linewidth = 0.5) +
    geom_point(aes(fill = iucn_label, size = log(n_pixels)),
               shape = 21, alpha = 0.35, stroke = 0.1, color = "grey40") +
    scale_fill_manual(values = iucn_colors, name = NULL) +
    scale_size_continuous(range = c(0.4, 3.5), guide = "none") +
    annotate("text", x = ax_lo + 0.01, y = ax_hi - 0.005,
             label = "Underperformance\n(actual > frontier)",
             hjust = 0, vjust = 1, size = 2.8, fontface = "italic",
             color = "coral4") +
    annotate("text", x = ax_hi - 0.01, y = ax_lo + 0.005,
             label = "Overperformance\n(actual < frontier)",
             hjust = 1, vjust = 0, size = 2.8, fontface = "italic",
             color = "steelblue4") +
    annotate("text", x = ax_hi - 0.005, y = ax_hi - 0.005,
             label = sprintf("Mean eff. = %.3f\n\u03B3 = %.3f\nN = %s",
                             mean(d_pa$te_jlms), gamma_val,
                             format(nrow(d_pa), big.mark = ",")),
             hjust = 1, vjust = 1, size = 2.8, color = "grey20") +
    labs(x = expression("Frontier " * hat(delta)[frontier] *
                        " (maximum achievable threat reduction)"),
         y = expression("Observed " * delta *
                        " (actual treatment effect)")) +
    coord_fixed(ratio = 1, xlim = c(ax_lo, ax_hi), ylim = c(ax_lo, ax_hi)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
    ) +
    guides(fill = guide_legend(override.aes = list(size = 3, alpha = 0.7)))

  p_top <- ggplot(d_pa, aes(x = frontier)) +
    geom_density(aes(fill = iucn_label), alpha = 0.4, linewidth = 0.3,
                 color = "grey40") +
    scale_fill_manual(values = iucn_colors) +
    scale_x_continuous(limits = c(ax_lo, ax_hi)) +
    theme_void() +
    theme(legend.position = "none")

  p_right <- ggplot(d_pa, aes(y = delta)) +
    geom_density(aes(fill = iucn_label), alpha = 0.4, linewidth = 0.3,
                 color = "grey40") +
    scale_fill_manual(values = iucn_colors) +
    scale_y_continuous(limits = c(ax_lo, ax_hi)) +
    theme_void() +
    theme(legend.position = "none")

  p_combined <- p_top + plot_spacer() +
    p_main + p_right +
    plot_layout(ncol = 2, nrow = 2,
                widths = c(4, 1), heights = c(1, 4))

  ggsave(paste0(paths$figures_dir, "fig.sfa.actual_vs_frontier.jpg"),
         plot = p_combined, width = 20, height = 20, units = "cm", dpi = 300)
  cat("    Saved: fig.sfa.actual_vs_frontier.jpg\n")
}


############################################################
### Fig 1b: Levels frontier curve (uhet, appendix) #########
############################################################
{
  cat("  Fig 1b: Levels frontier curve (uhet, appendix)\n")

  sf_main <- models[["uhet"]]
  frontier_coefs <- coef(sf_main)
  get_coef <- function(nm) {
    if (nm %in% names(frontier_coefs)) frontier_coefs[nm] else 0
  }

  tc_grid <- seq(0, max(d_pa$tc_control, na.rm = TRUE) * 1.02, length.out = 400)
  d_frontier_levels <- data.frame(
    tc_control   = tc_grid,
    tc_protected = get_coef("tc_control") * tc_grid +
      get_coef("tc_control_sq") * tc_grid^2 +
      get_coef("(Intercept)")
  )
  frontier_label_x <- 0.205
  frontier_label_y <- get_coef("tc_control") * frontier_label_x +
    get_coef("tc_control_sq") * frontier_label_x^2 +
    get_coef("(Intercept)")

  d_country_levels <- d_country %>%
    select(country_name, tc_control_mean, tc_protected_mean, total_pa_area)

  d_pa_levels <- d_pa %>%
    select(wdpaid, tc_control, tc_protected, size)

  select_countries <- c("Poland", "Cameroon", "Guinea", "Malaysia")
  d_select_levels <- d_country_levels %>%
    filter(country_name %in% select_countries)

  ax_hi <- max(
    quantile(d_country_levels$tc_control_mean,   0.99, na.rm = TRUE),
    quantile(d_country_levels$tc_protected_mean, 0.99, na.rm = TRUE)
  ) * 1.05
  ax_lim <- c(-0.005, ax_hi)

  col_country  <- "#5E3C99"
  col_frontier <- "#E69F00"
  col_ols      <- "#0072B2"

  line_cols <- scale_color_manual(
    values = c("SFA frontier" = col_frontier,
               "OLS fitted line" = col_ols),
    name = NULL
  )
  line_types <- scale_linetype_manual(
    values = c("SFA frontier" = "solid",
               "OLS fitted line" = "longdash"),
    name = NULL
  )

  p_country_levels <- ggplot(d_country_levels,
                             aes(x = tc_control_mean, y = tc_protected_mean)) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "grey55") +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "grey55") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "grey60", linewidth = 0.4) +
    geom_point(aes(size = total_pa_area), shape = 21,
               fill = "grey82", alpha = 0.65, stroke = 0.2,
               color = "grey55") +
    geom_smooth(aes(color = "OLS fitted line", linetype = "OLS fitted line"),
                method = "lm", formula = y ~ x, se = FALSE,
                linewidth = 0.8) +
    geom_line(data = d_frontier_levels,
              aes(x = tc_control, y = tc_protected,
                  color = "SFA frontier", linetype = "SFA frontier"),
              linewidth = 0.9, inherit.aes = FALSE) +
    geom_point(data = d_select_levels, aes(size = total_pa_area), shape = 21,
               fill = col_country, alpha = 0.9, stroke = 0.4,
               color = "grey20") +
    ggrepel::geom_text_repel(
      data = d_select_levels,
      aes(label = country_name),
      size = 2.8, color = "grey15", fontface = "bold",
      max.overlaps = Inf, segment.color = "grey70", segment.size = 0.2,
      min.segment.length = 0.1, box.padding = 0.25, seed = 42
    ) +
    scale_size_continuous(
      range = c(1.5, 7),
      trans = scales::trans_new("log2", log2, function(x) 2^x),
      breaks = c(100, 1000, 10000, 100000),
      guide = "none"
    ) +
    line_cols +
    line_types +
    labs(x = "Matched control threat level (0-1 scale)",
         y = "Threat level inside PAs (0-1 scale)") +
    coord_cartesian(xlim = ax_lim, ylim = ax_lim) +
    fig_theme +
    theme(legend.position = "none")

  p_pa_levels <- ggplot(d_pa_levels, aes(x = tc_control, y = tc_protected)) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "grey55") +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "grey55") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "grey30", linewidth = 0.65) +
    geom_point(aes(size = size), color = "grey65", alpha = 0.14) +
    geom_line(data = d_frontier_levels,
              aes(x = tc_control, y = tc_protected),
              color = col_frontier, linewidth = 0.9,
              inherit.aes = FALSE) +
    scale_size_continuous(
      name = expression("PA area (km"^2*")"),
      range = c(0.2, 2.4),
      trans = "log10",
      breaks = c(10, 100, 1000, 10000),
      labels = scales::comma
    ) +
    labs(x = "Matched control threat level (0-1 scale)",
         y = NULL) +
    coord_fixed(ratio = 1, xlim = ax_lim, ylim = ax_lim) +
    fig_theme +
    theme(legend.position = "none")

  p_pa_levels_standalone <- p_pa_levels +
    labs(y = "Actual threat level inside PAs (0–1 scale)") +
    sfa_comparison_theme

  p_levels_frontier <- p_country_levels + p_pa_levels +
    plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")

  ggsave(paste0(paths$figures_dir, "fig.sfa.levels.frontier-curve.uhet.country.jpg"),
         plot = p_country_levels, width = 12, height = 10,
         units = "cm", dpi = 300)
  cat("    Saved: fig.sfa.levels.frontier-curve.uhet.country.jpg\n")

  ggsave(paste0(paths$figures_dir, "fig.sfa.levels.frontier-curve.uhet.pa.jpg"),
         plot = p_pa_levels_standalone, width = 14, height = 14,
         units = "cm", dpi = 300)
  cat("    Saved: fig.sfa.levels.frontier-curve.uhet.pa.jpg\n")

  ggsave(paste0(paths$figures_dir, "fig.sfa.levels.frontier-curve.uhet.jpg"),
         plot = p_levels_frontier, width = 24, height = 12,
         units = "cm", dpi = 300)
  cat("    Saved: fig.sfa.levels.frontier-curve.uhet.jpg\n")
}


#################################################
### Fig 2: Frontier curve 2D-binned #############
#################################################
{
  cat("  Fig 2: Frontier curve (2D-binned, normalized)\n")

  sf_main <- models[["uhet"]]
  size_med <- median(d_pa$log_size)
  frontier_label <- sprintf("Frontier (median PA: %.0f km\u00B2)", exp(size_med) - 1)
  mean_eff <- mean(d_pa$te_jlms)
  pars <- coef(sf_main, extraPar = TRUE)
  gamma_val_2 <- pars["gamma"]

  # Build the delta-equivalent frontier curve from the level-space SFA fit.
  # sf_main fits tc_protected ~ 0 + tc_control + tc_control_sq (uhet on
  # log_size); convert to delta-space via fr_delta = fr_level - tc_control.
  frontier_coefs     <- coef(sf_main)
  frontier_intercept <- ifelse("(Intercept)" %in% names(frontier_coefs),
                               frontier_coefs["(Intercept)"], 0)
  get_coef <- function(nm) if (nm %in% names(frontier_coefs)) frontier_coefs[nm] else 0
  tc_seq    <- seq(min(d_pa$tc_control), max(d_pa$tc_control), length.out = 200)
  fr_level  <- frontier_intercept +
    get_coef("tc_control")    * tc_seq +
    get_coef("tc_control_sq") * tc_seq^2 +
    get_coef("log_size")      * size_med
  frontier_line   <- fr_level - tc_seq
  d_frontier_norm <- data.frame(
    tc_control    = tc_seq,
    frontier_norm = frontier_line / ctrl_mean * 100
  )

  d_pa$delta_norm <- d_pa$delta / ctrl_mean * 100

  # 2D quantile binning
  nx <- 25
  ny <- 6
  d_pa$xbin <- ntile(d_pa$tc_control, nx)
  d_pa <- d_pa %>%
    group_by(xbin) %>%
    mutate(ybin = ntile(delta_norm, ny)) %>%
    ungroup()

  d_binned_2d <- d_pa %>%
    group_by(xbin, ybin) %>%
    summarise(
      tc_mid     = mean(tc_control),
      delta_norm = mean(delta_norm),
      size       = mean(size),
      n          = n(),
      .groups = "drop"
    )

  p_2d <- ggplot(d_binned_2d, aes(x = tc_mid, y = delta_norm)) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 11, face = "bold")
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = global_effect_norm, linetype = "dashed",
               color = "steelblue", linewidth = 0.5) +
    geom_line(data = d_frontier_norm, aes(x = tc_control, y = frontier_norm),
              color = "#B2182B", linewidth = 1.2) +
    geom_point(aes(size = size), alpha = 0.5, shape = 16, color = "grey30") +
    scale_size_continuous(name = expression("PA area (km"^2*")"), range = c(0.5, 4),
                          trans = "log10", breaks = c(10, 100, 1000, 10000)) +
    geom_smooth(data = d_pa, aes(x = tc_control, y = delta_norm),
                method = "loess", span = 0.75, linewidth = 1,
                se = TRUE, alpha = 0.15, color = "#2166AC", fill = "#2166AC") +
    annotate("text", x = -Inf, y = global_effect_norm, label = "Global average",
             hjust = -0.05, vjust = -0.5, size = 3, color = "steelblue4",
             fontface = "italic") +
    annotate("text", x = max(d_pa$tc_control) * 0.65,
             y = min(d_frontier_norm$frontier_norm[
               d_frontier_norm$tc_control <= max(d_binned_2d$tc_mid)]) * 0.9,
             label = frontier_label, color = "#B2182B", fontface = "bold",
             size = 3, hjust = 0) +
    annotate("text", x = max(d_pa$tc_control) * 0.55,
             y = min(d_binned_2d$delta_norm) * 0.65,
             label = sprintf("\u03B3 = %.3f\nMean eff. = %.3f", gamma_val_2, mean_eff),
             size = 3, color = "grey30", hjust = 0) +
    labs(
      x = "Control threat level (counterfactual threat composite)",
      y = "Effect on overall threat index (% of control mean)"
    )

  ggsave(paste0(paths$figures_dir, "fig.sfa.frontier_curve.norm.2dbin.jpg"),
         plot = p_2d, width = 18, height = 14, units = "cm", dpi = 300)
  cat("    Saved: fig.sfa.frontier_curve.norm.2dbin.jpg\n")
}


####################################################
### Fig 3: Country scatter — all countries #########
####################################################
{
  cat("  Fig 3: Country scatter (all countries, normalized)\n")

  d_scatter_fig3 <- d_scatter_norm %>%
    filter(country_name != "Libya")

  d_label_all <- label_outskirts(d_scatter_fig3, "frontier_norm", "delta_norm", top_n = 25)

  p_country_scatter <- ggplot(d_scatter_fig3, aes(x = frontier_norm, y = delta_norm)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "black") +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "grey30", linewidth = 0.5) +
    geom_vline(xintercept = global_frontier_norm, linetype = "dashed",
               color = "firebrick", linewidth = 0.4, alpha = 0.7) +
    geom_hline(yintercept = global_delta_norm, linetype = "dashed",
               color = "steelblue", linewidth = 0.4, alpha = 0.7) +
    annotate("text", x = global_frontier_norm, y = y_lim[1] + y_pad,
             label = sprintf("Global potential: %.1f%%", global_frontier_norm),
             hjust = -0.05, vjust = 0, size = 2.5, color = "firebrick4") +
    annotate("text", x = min(d_scatter_fig3$frontier_norm),
             y = global_delta_norm,
             label = sprintf("Global actual: %.1f%%", global_delta_norm),
             hjust = 0, vjust = -0.5, size = 2.5, color = "steelblue4") +
    geom_point(aes(size = total_pa_area), shape = 21,
               fill = "steelblue", alpha = 0.5, stroke = 0.3, color = "grey30") +
    size_scale +
    ggrepel::geom_text_repel(
      data = d_label_all,
      aes(label = country_name), size = 2.2, color = "grey30",
      max.overlaps = 25, segment.color = "grey70", segment.size = 0.2,
      min.segment.length = 0.1, box.padding = 0.25, seed = 42) +
    annotate("text", x = min(d_scatter_fig3$frontier_norm) + 1,
             y = y_lim[2] - y_pad,
             label = "Underperformance",
             hjust = 0, vjust = 1, size = 2.8, fontface = "italic", color = "coral4") +
    annotate("text", x = max(d_scatter_fig3$frontier_norm) - 1,
             y = y_lim[1] + (y_lim[2] - y_lim[1]) * 0.12,
             label = "Overperformance",
             hjust = 1, vjust = 0, size = 2.8, fontface = "italic", color = "steelblue4") +
    labs(x = "Potential effect (% of control mean)",
         y = "Effect on threat index (% of control mean)") +
    coord_cartesian(ylim = y_lim) +
    fig_theme

  ggsave(paste0(paths$figures_dir, "fig.sfa.country_scatter.norm.jpg"),
         plot = p_country_scatter, width = 16, height = 14, units = "cm", dpi = 300)
  cat("    Saved: fig.sfa.country_scatter.norm.jpg\n")
}


####################################################
### Fig 4: Country scatter — 3 selected ############
####################################################
{
  cat("  Fig 4: Country scatter (3-country select)\n")

  # use raw (non-normalized) effect sizes
  d_scatter_raw <- d_country %>%
    select(country_name, delta_mean, frontier_mean, total_pa_area)

  select_countries <- c("Poland", "Cameroon", "Guinea", "Malaysia")
  d_select <- d_scatter_raw %>%
    filter(country_name %in% select_countries)

  # Poland for efficiency gap annotation
  d_gap <- d_select %>%
    filter(country_name == "Poland")

  # efficiency gap ] bracket on Poland (from dot to 45-degree line)
  y_top <- d_gap$delta_mean
  y_bot <- d_gap$frontier_mean
  y_mid <- (y_top + y_bot) / 2
  x_base <- d_gap$frontier_mean + 0.004  # shift right of the dot
  cap_w <- 0.002  # cap length in x-axis units

  # y-axis limits for raw scale
  y_lo_raw <- min(d_scatter_raw$delta_mean, na.rm = TRUE)
  y_hi_raw <- max(d_scatter_raw$delta_mean, na.rm = TRUE)
  y_pad_raw <- (y_hi_raw - y_lo_raw) * 0.03
  y_lim_raw <- c(y_lo_raw - y_pad_raw, y_hi_raw + y_pad_raw)

  p_select <- ggplot(d_scatter_raw, aes(x = frontier_mean, y = delta_mean)) +
    # zero reference lines
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey50") +
    # 45-degree frontier line
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "grey30", linewidth = 0.5) +
    # all countries as faint background
    geom_point(aes(size = total_pa_area), shape = 21,
               fill = "grey80", alpha = 0.25, stroke = 0.2, color = "grey60") +
    # measured effect ] bracket on Poland (from zero to dot, with gap)
    annotate("segment", x = x_base + cap_w, xend = x_base + cap_w,
             y = 0, yend = y_top + 0.001,
             color = "steelblue3", linewidth = 0.5) +
    annotate("segment", x = x_base, xend = x_base + cap_w,
             y = 0, yend = 0,
             color = "steelblue3", linewidth = 0.5) +
    annotate("segment", x = x_base, xend = x_base + cap_w,
             y = y_top + 0.001, yend = y_top + 0.001,
             color = "steelblue3", linewidth = 0.5) +
    annotate("text",
             x = x_base + cap_w + 0.001,
             y = (0 + y_top + 0.001) / 2,
             label = "Measured\neffect",
             hjust = 0, size = 2.8, color = "steelblue4", fontface = "italic") +
    # efficiency gap ] bracket on Poland (from dot to 45-degree line, with gap)
    annotate("segment", x = x_base + cap_w, xend = x_base + cap_w,
             y = y_top - 0.001, yend = y_bot,
             color = "coral3", linewidth = 0.5) +
    annotate("segment", x = x_base, xend = x_base + cap_w,
             y = y_top - 0.001, yend = y_top - 0.001,
             color = "coral3", linewidth = 0.5) +
    annotate("segment", x = x_base, xend = x_base + cap_w,
             y = y_bot, yend = y_bot,
             color = "coral3", linewidth = 0.5) +
    annotate("text",
             x = x_base + cap_w + 0.001,
             y = (y_top - 0.001 + y_bot) / 2,
             label = "Efficiency\ngap",
             hjust = 0, size = 2.8, color = "coral4", fontface = "italic") +
    # selected country dots
    geom_point(data = d_select, aes(size = total_pa_area), shape = 21,
               fill = "steelblue", alpha = 0.8, stroke = 0.4, color = "grey20") +
    size_scale +
    # labels (no connector lines)
    geom_text(data = d_select,
              aes(label = country_name),
              size = 3.2, color = "grey10", fontface = "bold",
              vjust = -1.2) +
    labs(x = expression("Potential effect (frontier " * hat(delta)[frontier] * ")"),
         y = expression("Actual effect (" * delta * ")")) +
    coord_cartesian(ylim = y_lim_raw) +
    fig_theme +
    guides(size = "none")

  ggsave(paste0(paths$figures_dir, "fig.sfa.country_scatter.norm-select.jpg"),
         plot = p_select, width = 16, height = 14, units = "cm", dpi = 300)
  cat("    Saved: fig.sfa.country_scatter.norm-select.jpg\n")
}


###########################################
### Fig 5: Country lollipop ###############
###########################################
{
  cat("  Fig 5: Country lollipop\n")

  global_delta_lp    <- weighted.mean(d_country$delta_mean, d_country$n_pixels)
  global_frontier_lp <- weighted.mean(d_country$frontier_mean, d_country$n_pixels)

  d_lollipop <- d_country %>%
    mutate(
      country_label = ifelse(is.na(country_name),
                             as.character(country_rast), country_name)
    ) %>%
    arrange(gap_mean)

  d_lollipop$country_label <- factor(d_lollipop$country_label,
                                      levels = d_lollipop$country_label)

  p_lollipop <- ggplot(d_lollipop, aes(y = country_label)) +
    geom_vline(xintercept = global_delta_lp, linewidth = 0.5, linetype = "dashed",
               color = "steelblue") +
    geom_vline(xintercept = global_frontier_lp, linewidth = 0.5, linetype = "dashed",
               color = "firebrick") +
    annotate("text", x = global_delta_lp, y = Inf,
             label = sprintf("Global actual: %.3f", global_delta_lp),
             vjust = -0.3, hjust = 0.5, size = 2.5, color = "steelblue4") +
    annotate("text", x = global_frontier_lp, y = Inf,
             label = sprintf("Global frontier: %.3f", global_frontier_lp),
             vjust = -1.5, hjust = 0.5, size = 2.5, color = "firebrick4") +
    geom_segment(aes(x = delta_mean, xend = frontier_mean,
                     yend = country_label),
                 color = "grey60", linewidth = 0.4) +
    geom_point(aes(x = delta_mean, color = "Actual effect", shape = "Actual effect"),
               size = 1.8) +
    geom_point(aes(x = frontier_mean, color = "Frontier", shape = "Frontier"),
               size = 1.8) +
    scale_color_manual(values = c("Actual effect" = "steelblue",
                                  "Frontier" = "firebrick"),
                       name = NULL) +
    scale_shape_manual(values = c("Actual effect" = 16,
                                  "Frontier" = 1),
                       name = NULL) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
    labs(x = expression("Treatment effect " * delta),
         y = NULL) +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      axis.title.x = element_text(size = 10),
      axis.text.y = element_text(size = 5),
      axis.text.x = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.2),
      plot.margin = margin(20, 5, 5, 5, "pt")
    )

  ggsave(paste0(paths$figures_dir, "fig.sfa.country_lollipop.jpg"),
         plot = p_lollipop, width = 18, height = 30, units = "cm", dpi = 300)
  cat("    Saved: fig.sfa.country_lollipop.jpg\n")
}


###################################################
### Fig 5b: Country lollipop — threat levels ######
###################################################
{
  cat("  Fig 5b: Country lollipop (threat levels)\n")

  # Whole-PA country aggregates have their own estimand. A separate country
  # regression's p-value cannot label these plotted levels.
  d_lollipop_lv <- d_country %>%
    mutate(
      country_label = ifelse(is.na(country_name), as.character(country_rast), country_name),
      actual_level = tc_protected_mean,
      potential_level = tc_control_mean + frontier_mean
    ) %>%
    arrange(actual_level)

  d_lollipop_lv$country_label <- factor(d_lollipop_lv$country_label,
                                         levels = d_lollipop_lv$country_label)

  p_lollipop_lv <- ggplot(d_lollipop_lv, aes(y = country_label)) +
    geom_segment(aes(x = potential_level, xend = tc_control_mean,
                     yend = country_label),
                 color = "grey75", linewidth = 0.4) +
    geom_point(aes(x = tc_control_mean),
               color = "grey40", shape = 16, size = 1.8) +
    geom_point(aes(x = actual_level),
               color = "steelblue", shape = 16, size = 1.8) +
    geom_point(aes(x = potential_level),
               color = "firebrick", shape = 1, size = 1.8) +
    # manual legend via dummy points (off-plot, hidden by coord_cartesian)
    geom_point(aes(x = -Inf, color = "Control"), shape = 16, size = 0, alpha = 0) +
    geom_point(aes(x = -Inf, color = "Protected"), shape = 16, size = 0, alpha = 0) +
    geom_point(aes(x = -Inf, color = "Potential"), shape = 1, size = 0, alpha = 0) +
    scale_color_manual(
      breaks = c("Control", "Protected", "Potential"),
      values = c("Control"            = "grey40",
                 "Protected" = "steelblue",
                 "Potential"           = "firebrick"),
      guide = guide_legend(
        override.aes = list(
          alpha = c(1, 1, 1),
          shape = c(16, 16, 1),
          size  = c(2.5, 2.5, 2.5)
        )
      ),
      name = NULL
    ) +
    labs(x = "Threat composite level",
         y = NULL) +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      axis.title.x = element_text(size = 10),
      axis.text.y = element_text(size = 5),
      axis.text.x = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.2),
      plot.margin = margin(20, 5, 5, 5, "pt")
    )

  ggsave(paste0(paths$figures_dir, "fig.sfa.country_lollipop-level.jpg"),
         plot = p_lollipop_lv, width = 18, height = 30, units = "cm", dpi = 300)
  cat("    Saved: fig.sfa.country_lollipop-level.jpg\n")
}


###################################################
### Fig 6: Efficiency gap — split sub-panels ######
### A) concept scatter   B) country leaderboard ###
### Saved as two JPGs for LaTeX subfigure layout ##
###################################################
{
  cat("  Fig 6: Efficiency gap (concept + leaderboard, split sub-panels)\n")

  # hero palette: teal = achieved, amber = closable gap
  col_achieved <- "#2C7DA0"
  col_gap      <- "#E69F00"

  # closable_total (gap_mean x area) is used to rank Panel B
  # delta_mean and frontier_mean are plotted directly in raw threat-composite units
  d_eg <- d_country %>%
    mutate(closable_total = gap_mean * total_pa_area)

  ## Panel A: C.8a-style closable-gap explainer in LEVEL space ##
  # x = matched-control threat level, y = actual within-PA threat level.
  # The dashed diagonal is the no-reduction line; the frontier curve is the
  # SFA-predicted minimum achievable within-PA threat level. For the annotated
  # country, the teal bracket shows the observed threat reduction and the amber
  # bracket shows the remaining closable gap to the frontier.
  col_country <- "#5E3C99"   # violet, for selected-country dots
  col_frontier <- col_gap

  sf_main <- models[["uhet"]]
  frontier_coefs <- coef(sf_main)
  get_coef <- function(nm) {
    if (nm %in% names(frontier_coefs)) frontier_coefs[nm] else 0
  }

  tc_grid <- seq(0, max(d_pa$tc_control, na.rm = TRUE) * 1.02, length.out = 400)
  d_frontier_levels <- data.frame(
    tc_control   = tc_grid,
    tc_protected = get_coef("tc_control") * tc_grid +
      get_coef("tc_control_sq") * tc_grid^2 +
      get_coef("(Intercept)")
  )

  d_scatter <- d_eg %>%
    select(country_name, tc_control_mean, tc_protected_mean,
           frontier_level_mean, total_pa_area) %>%
    rename(
      control_level   = tc_control_mean,
      actual_level    = tc_protected_mean,
      frontier_level  = frontier_level_mean
    )

  select_countries <- c("Poland", "Cameroon", "Guinea", "Malaysia")
  d_select <- d_scatter %>%
    filter(country_name %in% select_countries)

  # Closable-gap annotation on Poland.
  d_gap_pol <- d_select %>% filter(country_name == "Poland")

  # axis range: zoom around the bulk of countries while preserving a true
  # 45-degree no-reduction line via coord_fixed().
  ax_lim <- c(-0.005, 0.300)

  if (nrow(d_gap_pol) == 1) {
    x_pol      <- d_gap_pol$control_level
    y_control  <- d_gap_pol$control_level
    y_actual   <- d_gap_pol$actual_level
    y_frontier <- d_gap_pol$frontier_level
    x_off    <- ax_lim[2] * 0.020
    cap_w    <- ax_lim[2] * 0.008
    bracket_layers <- list(
      annotate("segment",
               x = x_pol + x_off, xend = x_pol + x_off,
               y = y_actual + 0.001, yend = y_control - 0.001,
               color = col_achieved, linewidth = 0.5),
      annotate("segment",
               x = x_pol + x_off - cap_w, xend = x_pol + x_off,
               y = y_control - 0.001, yend = y_control - 0.001,
               color = col_achieved, linewidth = 0.5),
      annotate("segment",
               x = x_pol + x_off - cap_w, xend = x_pol + x_off,
               y = y_actual + 0.001, yend = y_actual + 0.001,
               color = col_achieved, linewidth = 0.5),
      annotate("text",
               x = x_pol + x_off + cap_w * 1.5,
               y = (y_actual + y_control) / 2,
               label = "Estimated\nthreat reduction",
               hjust = 0, size = 3.6,
               color = col_achieved, fontface = "italic"),
      annotate("segment",
               x = x_pol + x_off, xend = x_pol + x_off,
               y = y_frontier + 0.001, yend = y_actual - 0.001,
               color = col_gap, linewidth = 0.5),
      annotate("segment",
               x = x_pol + x_off - cap_w, xend = x_pol + x_off,
               y = y_frontier + 0.001, yend = y_frontier + 0.001,
               color = col_gap, linewidth = 0.5),
      annotate("segment",
               x = x_pol + x_off - cap_w, xend = x_pol + x_off,
               y = y_actual - 0.001, yend = y_actual - 0.001,
               color = col_gap, linewidth = 0.5),
      annotate("text",
               x = x_pol + x_off + cap_w * 1.5,
               y = (y_frontier + y_actual) / 2,
               label = "Closable\ngap",
               hjust = 0, size = 3.8,
               color = col_gap, fontface = "italic")
    )
  } else {
    bracket_layers <- list()
  }

  pA <- ggplot(d_scatter, aes(x = control_level, y = actual_level)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey50") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "grey30", linewidth = 0.5) +
    geom_line(data = d_frontier_levels,
              aes(x = tc_control, y = tc_protected),
              color = col_frontier, linewidth = 0.9,
              inherit.aes = FALSE) +
    annotate("text",
             x = frontier_label_x, y = frontier_label_y + 0.018,
             label = "Minimum-threat\nfrontier",
             hjust = 0.5, size = 3.6,
             color = col_frontier, fontface = "bold") +
    geom_point(aes(size = total_pa_area), shape = 21,
               fill = "grey80", alpha = 0.30, stroke = 0.2,
               color = "grey60") +
    geom_point(data = d_select, aes(size = total_pa_area), shape = 21,
               fill = col_country, alpha = 0.85, stroke = 0.4,
               color = "grey20") +
    geom_text(data = d_select,
              aes(label = country_name),
              size = 4.2, color = "grey10", fontface = "bold",
              vjust = -1.2) +
    scale_size_continuous(
      range  = c(1.5, 8),
      trans  = scales::trans_new("log2", log2, function(x) 2^x),
      breaks = c(100, 1000, 10000, 100000),
      guide  = "none"
    ) +
    labs(
      x = "Matched control threat level (0-1 scale)",
      y = "Actual threat level inside PAs (0–1 scale)"
    ) +
    coord_fixed(ratio = 1, xlim = ax_lim, ylim = ax_lim, clip = "off") +
    fig_theme +
    sfa_comparison_theme

  pA <- pA + bracket_layers

  ## Panel B: top-25 leaderboard with Global row on top; raw delta units; no cap ##
  d_top25 <- d_eg %>%
    arrange(desc(closable_total)) %>%
    slice_head(n = 25) %>%
    select(country_name, delta_mean, frontier_mean)

  # Global aggregate row (pixel-weighted mean of country-level estimates)
  d_global <- tibble(
    country_name  = "Global",
    delta_mean    = weighted.mean(d_country$delta_mean,    d_country$n_pixels),
    frontier_mean = weighted.mean(d_country$frontier_mean, d_country$n_pixels)
  )

  d_bars <- bind_rows(d_global, d_top25) %>%
    mutate(country_name = factor(country_name, levels = rev(country_name)))
  n_rows <- nrow(d_bars)   # 26 = 25 countries + Global on top

  d_long <- bind_rows(
    d_bars %>%
      mutate(component = "Measured effect",
             x_start = 0,
             x_end   = delta_mean),
    d_bars %>%
      mutate(component = "Closable gap",
             x_start = delta_mean,
             x_end   = frontier_mean)
  ) %>%
    mutate(component = factor(component,
                              levels = c("Measured effect", "Closable gap")))

  pB <- ggplot(d_long, aes(y = country_name)) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
    geom_segment(aes(x = x_start, xend = x_end,
                     yend = country_name, color = component),
                 linewidth = 4) +
    geom_point(data = d_bars,
               aes(x = delta_mean),
               color = col_achieved, size = 1.4) +
    geom_point(data = d_bars,
               aes(x = frontier_mean),
               color = col_gap, size = 1.6) +
    # separator line between Global row and the country block
    geom_hline(yintercept = n_rows - 0.5, linewidth = 0.3,
               color = "grey50", linetype = "dashed") +
    scale_color_manual(
      values = c(`Measured effect` = col_achieved,
                 `Closable gap`    = col_gap),
      name = NULL
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.02))) +
    labs(
      x = expression("Treatment effect on threat composite (" * delta * ")"),
      y = NULL
    ) +
    fig_theme +
    theme(
      axis.text.y = element_text(size = 6.5),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom"
    )

  ## Save as two standalone sub-panels for LaTeX subfigure layout ##
  ggsave(paste0(paths$figures_dir, "fig.sfa.closable-gap.jpg"),
         plot = pA, width = 14, height = 14, units = "cm", dpi = 300)
  cat("    Saved: fig.sfa.closable-gap.jpg\n")

  # Keep the legacy filename current for any drafts that still include it.
  ggsave(paste0(paths$figures_dir, "fig.sfa.efficiency-gap.concept.jpg"),
         plot = pA, width = 14, height = 14, units = "cm", dpi = 300)
  cat("    Saved: fig.sfa.efficiency-gap.concept.jpg\n")

  ggsave(paste0(paths$figures_dir, "fig.sfa.efficiency-gap.leaderboard.jpg"),
         plot = pB, width = 14, height = 14, units = "cm", dpi = 300)
  cat("    Saved: fig.sfa.efficiency-gap.leaderboard.jpg\n")
}


##################################################################
### Fig 7: Share of potential reduction achieved #################
### Sub-panel B of manuscript Figure 5 (efficiency gap) ##########
##################################################################
{
  cat("  Fig 7: Share of potential reduction achieved\n")

  # Share closed is the fraction of the frontier-implied potential threat
  # reduction that is already achieved by actual PA outcomes:
  # (control - actual) / (control - frontier).
  # A positive finite denominator is needed to interpret this as a share of
  # potential reduction. The 2026-09-14 sensitivity check found little change
  # in the central distribution after removing the PA-count and 0.005 cutoffs.
  # See ideas/20260914_sfa_share_filter_sensitivity.md.
  d_share <- d_country %>%
    mutate(
      potential_reduction = tc_control_mean - frontier_level_mean,
      observed_reduction  = tc_control_mean - tc_protected_mean,
      share_closed        = observed_reduction / potential_reduction
    ) %>%
    filter(
      is.finite(potential_reduction),
      potential_reduction > 0,
      is.finite(share_closed)
    )

  med_share <- median(d_share$share_closed, na.rm = TRUE)
  share_plot_range <- range(d_share$share_closed) +
    c(-1, 1) * 3 * stats::bw.nrd0(d_share$share_closed)

  p_share_closed <- ggplot(d_share, aes(x = share_closed)) +
    geom_density(fill = "#2C7DA0", color = "#2C7DA0",
                 alpha = 0.45, linewidth = 0.5) +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "grey45") +
    geom_vline(xintercept = 1, linewidth = 0.35, color = "grey45") +
    geom_vline(xintercept = med_share, linetype = "dashed",
               linewidth = 0.4, color = "#2C7DA0") +
    annotate("text", x = share_plot_range[1], y = Inf,
             label = sprintf("Median %.1f%%", med_share * 100),
             vjust = 1.4, hjust = 0, size = 2.7,
             color = "#2C7DA0", fontface = "bold") +
    annotate("text", x = 0, y = Inf, label = "0%",
             vjust = 1.4, hjust = -0.1, size = 2.8, color = "grey35") +
    annotate("text", x = 1, y = Inf, label = "100%",
             vjust = 1.4, hjust = 1.1, size = 2.8, color = "grey35") +
    scale_x_continuous(
      labels = scales::label_percent(accuracy = 1),
      breaks = scales::breaks_width(0.5),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    # Expand the display to show the full density without dropping observations
    # before smoothing, as fixed scale limits would do.
    expand_limits(x = share_plot_range) +
    labs(
      x = "Share of potential threat reduction achieved",
      y = "Density"
    ) +
    fig_theme +
    theme(legend.position = "none")

  ggsave(paste0(paths$figures_dir, "fig.sfa.share-closed.density.jpg"),
         plot = p_share_closed, width = 14, height = 6.4, units = "cm", dpi = 300)
  cat("    Saved: fig.sfa.share-closed.density.jpg\n")
}


##################################################################
### Fig 7b: Global frontier budget ###############################
### Sub-panel C of manuscript Figure 5 (efficiency gap) ##########
##################################################################
{
  cat("  Fig 7b: Global frontier budget\n")

  # Aggregate the budget at the PA level so larger protected areas carry
  # proportionally more weight. The achieved segment is minus the PA-level
  # treatment effect times PA area. The estimated-underperformance segment is
  # the SFA conditional mean u_hat times PA area; it is not the positive part
  # of the raw deterministic-frontier residual.
  d_budget_pa <- d_pa %>%
    mutate(
      achieved_reduction = -delta,
      estimated_underperformance = u_hat
    )

  display_scale <- 1000
  achieved_total <- sum(
    d_budget_pa$achieved_reduction * d_budget_pa$size,
    na.rm = TRUE
  ) /
    display_scale
  underperformance_total <- sum(
    d_budget_pa$estimated_underperformance * d_budget_pa$size,
    na.rm = TRUE
  ) /
    display_scale
  budget_total <- achieved_total + underperformance_total

  d_budget <- data.frame(
    component = factor(
      c("Achieved", "Estimated underperformance"),
      levels = c("Achieved", "Estimated underperformance")
    ),
    value = c(achieved_total, underperformance_total)
  ) %>%
    mutate(
      xmin = cumsum(lag(value, default = 0)),
      xmax = cumsum(value),
      xmid = (xmin + xmax) / 2,
      share = value / budget_total,
      label = sprintf("%s\n%.1f\n%.1f%%", component, value, share * 100)
    )

  budget_cols <- c(
    "Achieved"                    = "#2C7DA0",
    "Estimated underperformance" = "#E69F00"
  )

  p_budget <- ggplot(d_budget) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1,
                  fill = component),
              color = "white", linewidth = 0.7) +
    geom_text(aes(x = xmid, y = 0.5, label = label,
                  color = component),
              size = 3.2, fontface = "bold", lineheight = 0.95) +
    annotate("segment", x = budget_total, xend = budget_total,
             y = -0.08, yend = 1.08,
             linewidth = 0.45, color = "grey30") +
    scale_fill_manual(values = budget_cols, guide = "none") +
    scale_color_manual(
      values = c("Achieved" = "white",
                 "Estimated underperformance" = "white"),
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = c(0, achieved_total,
                 budget_total),
      labels = scales::number_format(accuracy = 0.1),
      expand = expansion(mult = c(0.025, 0.01))
    ) +
    labs(x = expression("Threat-index area (" * 10^3 * " km"^2 * ")"),
         y = NULL) +
    coord_cartesian(ylim = c(-0.25, 1.05), clip = "off") +
    fig_theme +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5.5, 12, 20, 14, "pt")
    )

  ggsave(paste0(paths$figures_dir, "fig.sfa.global-frontier-budget.jpg"),
         plot = p_budget, width = 14, height = 5.6, units = "cm", dpi = 300)
  cat("    Saved: fig.sfa.global-frontier-budget.jpg\n")
}


##################################################################
### Fig C.8c: Efficiency density (country-level) #################
### Former sub-panel B of manuscript Figure 5 ####################
##################################################################
{
  cat("  Fig C.8c: Efficiency density (country-level)\n")

  # Countries reuse Panel A's teal (col_achieved = "#2C7DA0"); PAs in a deep violet
  # to stay distinguishable from both the teal and the amber closable-gap color.
  pal_lvl <- c("Protected areas" = "#5E3C99", "Countries" = "#2C7DA0")

  d_eff <- bind_rows(
    # archived: PA-level density (manuscript shows country-level only).
    # data.frame(efficiency = d_pa$te_jlms,      level = "Protected areas"),
    data.frame(efficiency = d_country$te_mean, level = "Countries")
  ) %>%
    filter(!is.na(efficiency)) %>%
    mutate(level = factor(level, levels = c("Protected areas", "Countries")))

  med_lines <- d_eff %>%
    group_by(level) %>%
    summarise(med_eff = median(efficiency, na.rm = TRUE), .groups = "drop")

  p_density <- ggplot(d_eff, aes(x = efficiency, fill = level, color = level)) +
    geom_density(alpha = 0.45, linewidth = 0.5) +
    geom_vline(
      data = med_lines,
      aes(xintercept = med_eff, color = level),
      linetype = "dashed", linewidth = 0.4, show.legend = FALSE
    ) +
    scale_fill_manual(values = pal_lvl, name = NULL) +
    scale_color_manual(values = pal_lvl, guide = "none") +
    coord_cartesian(xlim = c(0.88, 1.00)) +
    labs(
      x = "Efficiency (1 = on frontier)",
      y = "Density"
    ) +
    fig_theme +
    # archived: legend hidden when only the country-level density is plotted.
    # theme(legend.position = "bottom")
    theme(legend.position = "none")

  ggsave(paste0(paths$figures_dir, "fig.sfa.efficiency-gap.density.jpg"),
         plot = p_density, width = 14, height = 9, units = "cm", dpi = 300)
  cat("    Saved: fig.sfa.efficiency-gap.density.jpg\n")
}


######################################################################
### Fig 8: Efficiency density by taxon — ARCHIVED ####################
### Per-species efficiency = overlap-area-weighted mean of te_jlms ###
### across PAs that intersect the species' geographic range. #########
### Four overlaid taxon densities (mammal/bird/amphibian/reptile). ###
### Wrapped in `if (FALSE)` so it does not run by default; remove ####
### the wrapper to regenerate fig.sfa.efficiency-gap.bytaxa.density ##
### .jpg. Not used in the current manuscript. ########################
######################################################################
if (FALSE) {
  cat("  Fig 8: Efficiency density by taxon (species overlap-weighted)\n")

  sp_pa <- readRDS("data/store/pathreat.data.species-pa.Rds")

  d_eff_tax <- sp_pa %>%
    filter(!is.na(efficiency_s)) %>%
    mutate(taxon = tools::toTitleCase(taxon))

  tax_order <- c("Mammal", "Bird", "Amphibian", "Reptile")
  d_eff_tax$taxon <- factor(d_eff_tax$taxon, levels = tax_order)

  tax_pal <- c(
    "Amphibian" = "#1b9e77",
    "Bird"      = "#7570b3",
    "Mammal"    = "#d95f02",
    "Reptile"   = "#e7298a"
  )

  med_lines_tax <- d_eff_tax %>%
    group_by(taxon) %>%
    dplyr::summarise(med_eff = median(efficiency_s, na.rm = TRUE),
                     .groups = "drop")

  p_density_tax <- ggplot(d_eff_tax,
                          aes(x = efficiency_s, fill = taxon, color = taxon)) +
    geom_density(alpha = 0.35, linewidth = 0.5) +
    geom_vline(
      data = med_lines_tax,
      aes(xintercept = med_eff, color = taxon),
      linetype = "dashed", linewidth = 0.4, show.legend = FALSE
    ) +
    scale_fill_manual(values = tax_pal, name = NULL) +
    scale_color_manual(values = tax_pal, guide = "none") +
    coord_cartesian(xlim = c(0.80, 1.00)) +
    labs(
      x = "Species-level efficiency (1 = on frontier)",
      y = "Density"
    ) +
    fig_theme +
    theme(legend.position = "bottom")

  ggsave(paste0(paths$figures_dir, "fig.sfa.efficiency-gap.bytaxa.density.jpg"),
         plot = p_density_tax, width = 14, height = 9, units = "cm", dpi = 300)
  cat("    Saved: fig.sfa.efficiency-gap.bytaxa.density.jpg\n")
}


cat("=== Done: all SFA figures saved ===\n")
