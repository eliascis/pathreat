##################################################
### pathreat.analysis.bytaxa.fig.R ###############
### Taxon-level median species delta figure ######
##################################################
#
# Creates by-taxon figures from the centralized medians and bootstrap
# intervals produced by pathreat.analysis.bytaxa.R.
#
# Input:
#   results/pathreat.bytaxa.est.Rds
#   results/pathreat.global.fe-adjusted-levels.est.Rds
#
# Output:
#   pub/figures/fig.bytaxa.delta.est.jpg
#   pub/figures/fig.bytaxa.global-coef.est.jpg
#   pub/figures/fig.coverage-and-coef.jpg
#

library(dplyr)
library(ggplot2)
library(patchwork)

source("code/pathreat.analysis.config.R")

############################
### load data ##############
############################
{
cat("=== Loading data ===\n")

taxon_bundle_file <- file.path(paths$results_dir, "pathreat.bytaxa.est.Rds")
if (!file.exists(taxon_bundle_file)) {
  stop("Taxon estimates not found. Run code/pathreat.analysis.bytaxa.R first.")
}
taxon_bundle <- readRDS(taxon_bundle_file)
required_components <- c("settings", "counts", "effects", "levels")
missing_components <- setdiff(required_components, names(taxon_bundle))
if (length(missing_components) > 0) {
  stop("Taxon bundle is missing: ", paste(missing_components, collapse = ", "))
}

if (taxon_bundle$settings$seed != 42L ||
    taxon_bundle$settings$bootstrap_draws != 2000L) {
  stop("Taxon figures require the production seed 42 and 2,000 bootstrap draws.")
}

taxon_effects <- taxon_bundle$effects
taxon_levels <- taxon_bundle$levels
cat(sprintf(
  "  Taxon effects: %s rows; marginal levels: %s rows\n",
  format(nrow(taxon_effects), big.mark = ","),
  format(nrow(taxon_levels), big.mark = ",")
))

required_effect_columns <- c(
  "taxon", "weighting", "n_species", "pct_median", "pct_ci_low", "pct_ci_high"
)
required_level_columns <- c(
  "taxon", "n_species", "control_median", "control_approx_se",
  "protected_median", "protected_approx_se", "raw_effect_sign_p"
)
missing_effect_columns <- setdiff(required_effect_columns, names(taxon_effects))
missing_level_columns <- setdiff(required_level_columns, names(taxon_levels))
if (length(missing_effect_columns) > 0 || length(missing_level_columns) > 0) {
  stop(
    "Taxon bundle is missing required figure fields: ",
    paste(c(missing_effect_columns, missing_level_columns), collapse = ", ")
  )
}
}

##########################################
### figure: global weights + taxon #######
##########################################
{
cat("\nCreating combined global-weights + taxon figure...\n")

taxon_labels <- c(
  amphibian = "Amphibians",
  bird = "Birds",
  mammal = "Mammals",
  reptile = "Reptiles"
)

taxon_colors <- c(
  Amphibians = "#E69F00",
  Birds      = "#56B4E9",
  Mammals    = "#009E73",
  Reptiles   = "#CC79A7"
)

# --- global weighting scheme estimates (from global-threat.R) ---
global_threat_path <- file.path(paths$results_dir, "pathreat.global-threat.est.Rds")
if (!file.exists(global_threat_path)) {
  stop("Global-threat results not found. Run code/pathreat.analysis.global-threat.R first.")
}
global_threat <- readRDS(global_threat_path)$estimates

shape_map <- c(
  "Equal weights"   = "Equal weight",
  "Threatened spp." = "Spp. richness",
  "TBL"             = "TBL-weighted",
  "ED"              = "ED-weighted",
  "EDGE2"           = "EDGE-weighted",
  "EDGE"            = "EDGE-weighted"
)

global_data <- global_threat %>%
  transmute(
    label = "Global ATT",
    estimate = pct_of_control,
    ci_low = ci_low / control_mean * 100,
    ci_high = ci_high / control_mean * 100,
    estimate_type = shape_map[model],
    color_group = "Global"
  )

# --- centralized taxon median estimates (% scale, bootstrap CIs) ---
weight_labels <- c(
  equal = "Equal weight",
  tbl = "TBL-weighted",
  ed = "ED-weighted",
  edge = "EDGE-weighted"
)

taxon_data <- taxon_effects %>%
  filter(
    taxon %in% names(taxon_labels),
    weighting %in% names(weight_labels)
  ) %>%
  transmute(
    label = unname(taxon_labels[taxon]),
    group = "Taxon-level species median",
    estimate = pct_median,
    ci_low = pct_ci_low,
    ci_high = pct_ci_high,
    estimate_type = unname(weight_labels[weighting]),
    color_group = unname(taxon_labels[taxon])
  )

# combine with a blank spacer between groups
row_levels <- c(
  unname(taxon_labels),
  "",
  "Global ATT"
)

estimate_levels <- c(
  "EDGE-weighted",
  "ED-weighted",
  "TBL-weighted",
  "Spp. richness",
  "Equal weight"
)

plot_data <- bind_rows(global_data, taxon_data) %>%
  mutate(
    label = factor(label, levels = row_levels),
    estimate_type = factor(estimate_type, levels = estimate_levels)
  )

# colors
all_colors <- c("Global" = "grey30", taxon_colors)

# clip CIs for axis limits
ci_clipped <- pmin(pmax(c(plot_data$ci_low, plot_data$ci_high), -100), 100)
x_lim <- c(
  floor(min(c(plot_data$estimate, ci_clipped), na.rm = TRUE) / 10) * 10,
  ceiling(max(c(plot_data$estimate, ci_clipped), na.rm = TRUE) / 10) * 10
)

dodge <- position_dodge(width = 0.7)

fig_taxon <- ggplot(plot_data, aes(y = label,
                                   x = estimate,
                                   color = color_group,
                                   group = estimate_type)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey40") +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high),
                 position = dodge,
                 linewidth = 0.8,
                 show.legend = FALSE) +
  geom_point(aes(shape = estimate_type),
             position = dodge,
             size = 2.6,
             stroke = 0.8) +
  scale_color_manual(values = all_colors) +
  scale_shape_manual(
    values = c(
      "Equal weight"  = 16,
      "Spp. richness" = 1,
      "TBL-weighted"  = 17,
      "ED-weighted"   = 15,
      "EDGE-weighted" = 18
    ),
    na.translate = FALSE
  ) +
  scale_x_continuous(
    breaks = seq(-60, 10, by = 10),
    labels = function(x) paste0(x, "%"),
    limits = x_lim
  ) +
  scale_y_discrete(drop = FALSE) +
  labs(x = "Change in threat index (% of control mean)",
       y = NULL,
       shape = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(10, 15, 10, 10)
  ) +
  guides(color = "none")

figpath <- file.path(paths$figures_dir, "fig.bytaxa.delta.est.jpg")
ggsave(figpath, fig_taxon, width = 7.5, height = 5, dpi = 300)
cat("Saved:", figpath, "\n")
}

##########################################
### figure: global + taxon levels plot ###
##########################################
{
cat("\nCreating global + taxon levels figure (control vs. protected means)...\n")

# --- global estimate: FE-adjusted levels + country-pairs bootstrap SEs ---
global_fe_levels <- readRDS(paths$est_global_fe_levels)
if (!all(c("settings", "summary", "bootstrap_draws") %in%
         names(global_fe_levels))) {
  stop("FE-adjusted level bundle is missing required components")
}
if (global_fe_levels$settings$bootstrap_reps != 499L) {
  stop("Global level bars require the 499-draw production bootstrap")
}
required_global_level_columns <- c(
  "variable",
  "adjusted_control_mean",
  "adjusted_protected_mean",
  "bootstrap_control_se",
  "bootstrap_protected_se",
  "bootstrap_reps_valid"
)
missing_global_level_columns <- setdiff(
  required_global_level_columns,
  names(global_fe_levels$summary)
)
if (length(missing_global_level_columns) > 0) {
  stop(paste(
    "FE-adjusted level summary is missing:",
    paste(missing_global_level_columns, collapse = ", ")
  ))
}
tc_gm <- global_fe_levels$summary[
  global_fe_levels$summary$variable == "threat_composite",
]
if (nrow(tc_gm) != 1 || tc_gm$bootstrap_reps_valid != 499L) {
  stop("Expected one threat_composite row with 499 valid bootstrap draws")
}
if (any(!is.finite(as.matrix(tc_gm[, c(
  "adjusted_control_mean",
  "adjusted_protected_mean",
  "bootstrap_control_se",
  "bootstrap_protected_se"
)]))) || tc_gm$bootstrap_control_se < 0 ||
    tc_gm$bootstrap_protected_se < 0) {
  stop("Threat-composite FE levels contain invalid means or standard errors")
}
tc_draws <- global_fe_levels$bootstrap_draws[["threat_composite"]]
if (is.null(tc_draws) || nrow(tc_draws) != 499L) {
  stop("Expected 499 retained threat_composite bootstrap draws")
}
global_est <- readRDS(file.path(paths$results_dir, "pathreat.global.est.Rds"))
tc_global <- global_est[global_est$variable == "threat_composite" &
                          global_est$estimate_type == "unscaled", ]

global_row <- data.frame(
  label        = "Global PA effect",
  control_mean = tc_gm$adjusted_control_mean,
  control_se   = tc_gm$bootstrap_control_se,
  treated_mean = tc_gm$adjusted_protected_mean,
  treated_se   = tc_gm$bootstrap_protected_se,
  pval         = tc_global$pval,
  stringsAsFactors = FALSE
)

# --- centralized taxon estimates: marginal medians + bootstrap pseudo-SEs ---
taxon_rows <- taxon_levels %>%
  transmute(
    label = taxon,
    control_mean = control_median,
    control_se = control_approx_se,
    treated_mean = protected_median,
    treated_se = protected_approx_se,
    pval = raw_effect_sign_p
  )
taxon_rows$label <- dplyr::recode(
  taxon_rows$label,
  amphibian = "Amphibians",
  bird = "Birds",
  mammal = "Mammals",
  reptile = "Reptiles"
)

global_pa_label <- "Global PA effect"

level_data <- bind_rows(global_row, taxon_rows) %>%
  mutate(
    label = factor(label, levels = rev(c(
      global_pa_label, "Mammals", "Birds", "Amphibians", "Reptiles"
    ))),
    sig_label = dplyr::case_when(
      pval <= 0.01 ~ "***",
      pval <= 0.05 ~ "**",
      pval <= 0.10 ~ "*",
      TRUE         ~ "ns"
    )
  )

row_colors <- c(
  "Global PA effect" = "grey30",
  "Mammals"          = "#009E73",
  "Birds"            = "#56B4E9",
  "Amphibians"       = "#E69F00",
  "Reptiles"         = "#CC79A7"
)

# long format: one row per (label, group)
bar_long <- bind_rows(
  level_data %>%
    transmute(label, group = "Control",
              mean = control_mean, se = control_se),
  level_data %>%
    transmute(label, group = "Protected",
              mean = treated_mean, se = treated_se)
) %>%
  mutate(
    # reversed so after coord_flip Control sits on top, Protected below
    group = factor(group, levels = c("Protected", "Control")),
    ymin = mean - se,
    ymax = mean + se
  )

# shared x-axis across all 5 rows (head-room for significance bracket)
x_max <- max(bar_long$ymax, na.rm = TRUE) * 1.30

# standalone version: 5 rows stacked via facet_grid
# use top-down factor (global PA effect on top) just for the standalone facet
bar_long_standalone <- bar_long %>%
  mutate(label = factor(as.character(label),
                        levels = c(global_pa_label, "Mammals", "Birds",
                                   "Amphibians", "Reptiles")))

fig_levels <- ggplot(bar_long_standalone, aes(x = group, y = mean, fill = label)) +
  geom_bar(stat = "identity", width = 0.65) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, color = "gray40") +
  scale_fill_manual(values = row_colors) +
  scale_y_continuous(limits = c(0, x_max), expand = c(0, 0)) +
  coord_flip() +
  facet_grid(label ~ ., switch = "y") +
  labs(x = NULL, y = "Threat composite [0–1]") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 9),
    plot.margin = margin(10, 15, 10, 10)
  )

figpath2 <- file.path(paths$figures_dir, "fig.bytaxa.global-coef.est.jpg")
ggsave(figpath2, fig_levels, width = 4.5, height = 5, dpi = 300)
cat("Saved:", figpath2, "\n")
}

##########################################
### combined figure: coverage + coef #####
##########################################
{
cat("\nCreating combined coverage + coefficient figure...\n")

# load panel A plot object from species-pa.fig.R
p_cover_path <- file.path(paths$results_dir, "p_combined_cover.Rds")
if (!file.exists(p_cover_path)) {
  stop("Panel A plot object not found. Run code/pathreat.analysis.species-pa.fig.R first.")
}
p_panel_a <- readRDS(p_cover_path)

# --- split panel B into 5 row-plots matching panel A's structure ---
# row order matches panel A: Global, Mammals, Birds, Amphibians, Reptiles
coef_row_order <- c(global_pa_label, "Mammals", "Birds", "Amphibians", "Reptiles")

# shared x-axis across rows; head-room for significance bracket
x_lim <- c(0, max(bar_long$ymax, na.rm = TRUE) * 1.30)

# helper: one row = two horizontal bars (Control, Protected) + sig bracket
make_level_row <- function(row_label, subtitle_i = NULL, show_axis = FALSE) {
  d  <- bar_long[bar_long$label == row_label, , drop = FALSE]
  lv <- level_data[level_data$label == row_label, ]
  fill_col <- unname(row_colors[row_label])

  bracket_y    <- max(d$ymax, na.rm = TRUE) * 1.08
  bracket_tick <- bracket_y * 0.97
  sig_label    <- lv$sig_label

  p <- ggplot(d, aes(x = group, y = mean)) +
    geom_bar(stat = "identity", width = 0.65, fill = fill_col) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, color = "gray40") +
    annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y, linewidth = 0.3) +
    annotate("segment", x = 1, xend = 1, y = bracket_y, yend = bracket_tick, linewidth = 0.3) +
    annotate("segment", x = 2, xend = 2, y = bracket_y, yend = bracket_tick, linewidth = 0.3) +
    annotate("text", x = 1.5, y = bracket_y * 1.10, label = sig_label, size = 2.8) +
    scale_y_continuous(limits = x_lim, expand = c(0, 0)) +
    coord_flip() +
    labs(x = NULL, y = NULL, subtitle = subtitle_i)

  base_theme <- theme_minimal(base_size = 10) +
    theme(
      panel.grid          = element_blank(),
      panel.grid.major.x  = element_line(color = "grey90", linewidth = 0.3),
      axis.text.y         = element_text(size = 8),
      axis.ticks.y        = element_blank(),
      plot.subtitle       = element_text(size = 9, face = "plain")
    )

  if (show_axis) {
    p + base_theme +
      labs(y = "Threat composite [0–1]") +
      theme(plot.margin = margin(0, 10, 5, 5))
  } else {
    p + base_theme +
      theme(axis.text.x = element_blank(),
            plot.margin = margin(0, 10, 0, 5))
  }
}

# build 5 right-side panels
# row 1: global PA effect (with subtitle matching panel A's "All terrestrial land")
p_coef_1 <- make_level_row(global_pa_label, subtitle_i = global_pa_label)
# row 2: Mammals (with subtitle matching panel A's "Threatened species ranges")
p_coef_2 <- make_level_row("Mammals", subtitle_i = "Mammals")
# row 3: Birds
p_coef_3 <- make_level_row("Birds", subtitle_i = "Birds")
# row 4: Amphibians
p_coef_4 <- make_level_row("Amphibians", subtitle_i = "Amphibians")
# row 5: Reptiles (with x-axis)
p_coef_5 <- make_level_row("Reptiles", subtitle_i = "Reptiles", show_axis = TRUE)

# stack right panels
p_panel_b <- p_coef_1 / p_coef_2 / p_coef_3 / p_coef_4 / p_coef_5 +
  plot_layout(heights = rep(1, 5))

# combine side by side
fig_combined <- (p_panel_a | p_panel_b) +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(tag_levels = list(c("A", "", "", "", "", "B", "", "", "", "")))

figpath3 <- file.path(paths$figures_dir, "fig.coverage-and-coef.jpg")
ggsave(figpath3, fig_combined, width = 8, height = 4.5, dpi = 300)
cat("Saved:", figpath3, "\n")
}

cat("\n=== By-taxon figure complete ===\n")
