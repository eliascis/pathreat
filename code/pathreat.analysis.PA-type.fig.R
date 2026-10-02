##########################################
### pathreat.analysis.PA-type.fig.R ######
### PA characteristic effects figures ####
##########################################

library(dplyr)
library(ggplot2)
library(readxl)
library(patchwork)

source("code/pathreat.analysis.config.R")

############################
### load data ##############
############################
{
cat("Loading PA-type estimation results...\n")
results_file <- paste0(paths$results_dir, "pathreat.PA-type.est.Rds")
est <- readRDS(results_file)
# filter to unscaled estimates (normalized added in 2c for heatmap)
if ("estimate_type" %in% names(est)) {
  est <- est[est$estimate_type == "unscaled", ]
}
cat(sprintf("  Loaded %d rows (%d subsamples x %d outcomes)\n",
            nrow(est), length(unique(est$subsample)), length(unique(est$variable))))

# 95% CI
stopifnot(all(c("ci_low", "ci_high") %in% names(est)))
}

############################
### labels & ordering ######
############################
{
# x-axis: IUCN classes first, then size classes (matches table column order)
subsample_order <- c("Strict", "Multi-use", "Not Reported",
                     "Class 1", "Class 2", "Class 3", "Class 4")
subsample_labels <- c("Strict\n(Ia-IV)", "Multi-\nuse\n(V-VI)", "Not\nReported",
                       "Class 1\n(<10 km\u00B2)", "Class 2\n(10-500 km\u00B2)", "Class 3\n(500-2k km\u00B2)", "Class 4\n(>2k km\u00B2)")
names(subsample_labels) <- subsample_order

est$subsample <- factor(est$subsample, levels = subsample_order)

# group indicator for visual separator
est$group <- ifelse(est$subsample %in% c("Strict", "Multi-use", "Not Reported"),
                    "IUCN", "Size")

# threat labels from Data_summary.xlsx
x.labels <- read_excel(paths$data_summary) %>%
  data.frame() %>%
  filter(!is.na(threat.no)) %>%
  select(variable, threat.no, threat.category.no, threat.category, threat.label, scale.label)

# category color mapping (same as global.fig.R)
var_to_cat <- c(
  "built" = "1", "cropland" = "2", "planted" = "2", "pasture" = "2",
  "oil" = "3", "mining" = "3", "renewables" = "3",
  "roads" = "4", "powerlines" = "4",
  "fires" = "7", "swu" = "7", "dams" = "7",
  "light" = "9", "def0120_parea" = "5", "any_threat" = "99", "threat_composite" = "99"
)
cat_colors <- c(
  "1" = "chocolate1", "2" = "darkolivegreen3", "3" = "gold3",
  "4" = "dodgerblue2", "5" = "darkseagreen4", "7" = "firebrick2",
  "9" = "gray45", "99" = "gray30"
)

# significance labels
sig_label <- function(p) {
  if (is.na(p)) return("")
  if (p <= 0.01) return("***")
  if (p <= 0.05) return("**")
  if (p <= 0.1)  return("*")
  return("")
}
}

############################
### create plots ###########
############################
{
cat("\nCreating individual panels...\n")

fig_width  <- 7   # cm per panel
fig_height <- 6   # cm per panel
plot_list  <- list()

for (v in deplist) {
  d.v <- est[est$variable == v & !is.na(est$subsample), ]

  # panel title from scale.label

  label_info <- x.labels[x.labels$variable == v, ]
  if (nrow(label_info) == 0) {
    panel_title <- v
  } else {
    sl <- label_info$scale.label[1]
    if (is.na(sl)) sl <- v
    panel_title <- trimws(sub("\\s*\\[.*", "", sl))
  }

  # color
  cat_no <- var_to_cat[v]
  pt_color <- unname(cat_colors[cat_no])

  # significance labels
  d.v$sig <- sapply(d.v$pval, sig_label)

  # y-axis range (symmetric around 0 if possible, with room for labels)
  y_abs <- max(abs(c(d.v$ci_low, d.v$ci_high)), na.rm = TRUE)
  y_pad <- y_abs * 0.35
  y_lim <- c(-y_abs - y_pad, y_abs + y_pad)

  # label position: above CI for positive, below CI for negative
  d.v$label_y <- ifelse(d.v$coef >= 0,
                        d.v$ci_high + y_abs * 0.08,
                        d.v$ci_low  - y_abs * 0.08)
  d.v$label_vjust <- ifelse(d.v$coef >= 0, 0, 1)

  # separator x position between IUCN (3) and Size (4) groups
  sep_x <- 3.5

  p <- ggplot(d.v, aes(x = subsample, y = coef)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 6, angle = 0, hjust = 0.5, lineheight = 0.85),
      axis.text.y = element_text(size = 7),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.y = element_line(color = "black"),
      axis.line.x = element_blank(),
      plot.margin = margin(5, 3, 3, 3, "pt")
    ) +
    labs(title = panel_title) +
    # zero reference line
    geom_hline(yintercept = 0, linewidth = 0.4) +
    # group separator
    geom_vline(xintercept = sep_x, linetype = "dashed", color = "grey70", linewidth = 0.3) +
    # point estimates + CI
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, color = "grey40", linewidth = 0.4) +
    geom_point(size = 2, color = pt_color) +
    # significance stars
    geom_text(aes(y = label_y, label = sig, vjust = label_vjust),
              size = 2.5, color = "black") +
    scale_x_discrete(labels = subsample_labels) +
    scale_y_continuous(limits = y_lim)

  plot_list[[v]] <- p

  # save individual
  fig_ind <- paste0(paths$figures_dir, "fig.PA-type.", v, ".est.jpg")
  ggsave(plot = p, fig_ind, units = "cm", width = fig_width, height = fig_height, dpi = 300)
  cat(sprintf("  Saved: %s  (coef range: %s to %s)\n", v,
              formatC(min(d.v$coef, na.rm = TRUE), format = "f", digits = 3),
              formatC(max(d.v$coef, na.rm = TRUE), format = "f", digits = 3)))
}
}

##########################################
### combined figure ######################
##########################################
{
cat("\nCombining panels...\n")

# horizontal legend for threat categories (below grid)
legend_data <- data.frame(
  category = c(
    "Residential &\ncommercial dev.",
    "Agriculture &\naquaculture",
    "Energy prod.\n& mining",
    "Transportation &\nservice corridors",
    "Natural system\nmodifications",
    "Pollution",
    "Biological\nresource use",
    "Overall"
  ),
  color = c("chocolate1", "darkolivegreen3", "gold3", "dodgerblue2",
            "firebrick2", "gray45", "darkseagreen4", "gray30"),
  x = 1:8,
  stringsAsFactors = FALSE
)

legend_plot <- ggplot(legend_data, aes(x = x, y = 1, fill = category)) +
  geom_tile(width = 0.85, height = 0.3) +
  geom_text(aes(y = 0.3, label = category), size = 2.8, vjust = 1, lineheight = 0.85) +
  scale_fill_manual(values = setNames(legend_data$color, legend_data$category)) +
  scale_x_continuous(limits = c(0, 9), expand = c(0, 0)) +
  ylim(-0.6, 1.2) +
  theme_void() +
  theme(legend.position = "none")

# 5x3 grid with horizontal legend below
row1 <- plot_list[["built"]] + plot_list[["cropland"]] + plot_list[["planted"]] +
        plot_list[["pasture"]] + plot_list[["light"]] + plot_layout(ncol = 5)
row2 <- plot_list[["oil"]] + plot_list[["mining"]] + plot_list[["renewables"]] +
        plot_list[["roads"]] + plot_list[["powerlines"]] + plot_layout(ncol = 5)
row3 <- plot_list[["fires"]] + plot_list[["swu"]] + plot_list[["dams"]] +
        plot_list[["def0120_parea"]] + plot_list[["threat_composite"]] + plot_layout(ncol = 5)

grid <- row1 / row2 / row3
combined <- grid / legend_plot + plot_layout(heights = c(1, 1, 1, 0.22))

fig_combined <- paste0(paths$figures_dir, "fig.PA-type.combined.est.jpg")
ggsave(plot = combined, fig_combined, units = "cm", width = 36, height = 22, dpi = 300)
cat("Saved combined figure:", fig_combined, "\n")
}

##########################################
### standalone threat_composite coef + CI #
### (forest plot, matching by-biome style)
##########################################
{
cat("\nCreating standalone threat_composite forest plot...\n")

# load group means by PA subcategory
gm <- readRDS(paths$est_PA_type_group_means)

# --- subsample estimates with control_mean for normalization ---
d.at <- est[est$variable == "threat_composite" & !is.na(est$subsample), ]
gm_at <- gm[gm$variable == "threat_composite", c("subsample", "control_mean", "treated_n")]
d.at <- merge(d.at, gm_at, by = "subsample")

# normalize: coef / control_mean * 100 = % of control mean
d.at$norm_coef     <- d.at$coef    / d.at$control_mean * 100
d.at$norm_ci_lower <- d.at$ci_low  / d.at$control_mean * 100
d.at$norm_ci_upper <- d.at$ci_high / d.at$control_mean * 100

# significance tiers
d.at$sig_tier <- "Not significant"
d.at$sig_tier[d.at$pval < 0.05] <- "p < 0.05"
d.at$sig_tier[d.at$pval < 0.01] <- "p < 0.01"
d.at$sig_tier <- factor(d.at$sig_tier,
                        levels = c("p < 0.01", "p < 0.05", "Not significant"))

# order by group: size classes (bottom), then IUCN classes (top)
# gap between size and IUCN groups (skip rank 5)
group_ranks <- c("Class 4" = 1, "Class 3" = 2, "Class 2" = 3, "Class 1" = 4,
                 "Not Reported" = 6, "Multi-use" = 7, "Strict" = 8)
d.at$rank <- group_ranks[as.character(d.at$subsample)]

# clean labels
subsample_labels_forest <- c(
  "Strict" = "Strict (Ia\u2013IV)", "Multi-use" = "Multi-use (V\u2013VI)",
  "Not Reported" = "Not Reported",
  "Class 1" = "Size <10 km\u00B2", "Class 2" = "Size 10\u2013500 km\u00B2",
  "Class 3" = "Size 500\u20132k km\u00B2", "Class 4" = "Size >2k km\u00B2")
d.at$label <- subsample_labels_forest[as.character(d.at$subsample)]

# --- global estimate at top ---
est_global <- readRDS(paths$est_global)
eg <- est_global[est_global$variable == "threat_composite" &
                   est_global$estimate_type == "unscaled", ]
global_n_treat <- eg$n_obs / 2
global_rank <- max(d.at$rank) + 2
global_row <- data.frame(
  norm_coef     = eg$coef / eg$control_mean * 100,
  norm_ci_lower = eg$ci_low / eg$control_mean * 100,
  norm_ci_upper = eg$ci_high / eg$control_mean * 100,
  pval = eg$pval,
  sig_tier = factor(
    if (eg$pval < 0.01) "p < 0.01" else if (eg$pval < 0.05) "p < 0.05" else "Not significant",
    levels = c("p < 0.01", "p < 0.05", "Not significant")),
  rank = global_rank,
  n_treat = global_n_treat,
  label = "Global"
)
global_row$is_global <- TRUE

# combine
keep_cols <- c("norm_coef", "norm_ci_lower", "norm_ci_upper", "pval",
               "sig_tier", "rank", "n_treat", "label")
d_sub <- d.at[, c(keep_cols[keep_cols != "n_treat"], "treated_n")]
names(d_sub)[names(d_sub) == "treated_n"] <- "n_treat"
d_sub$is_global <- FALSE
d_plot <- rbind(d_sub[, c(keep_cols, "is_global")],
                global_row[, c(keep_cols, "is_global")])

# shared y-axis limits
y_lim <- range(d_plot$rank) + c(-0.8, 0.8)

p1 <- ggplot(d_plot, aes(x = norm_coef, y = rank)) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
  ) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = norm_ci_lower, xmax = norm_ci_upper),
                 height = 0, linewidth = 0.4, color = "grey50") +
  geom_point(aes(shape = sig_tier, fill = sig_tier, size = is_global), color = "grey30") +
  scale_shape_manual(values = c("p < 0.01" = 21, "p < 0.05" = 21, "Not significant" = 21)) +
  scale_fill_manual(values = c("p < 0.01" = "black", "p < 0.05" = "grey60", "Not significant" = "white")) +
  scale_size_manual(values = c("FALSE" = 3, "TRUE" = 4.5), guide = "none") +
  scale_y_continuous(
    breaks = d_plot$rank,
    labels = d_plot$label
  ) +
  guides(shape = guide_legend(override.aes = list(size = 3)),
         fill = guide_legend(override.aes = list(size = 3))) +
  labs(
    x = "Effect on overall threat index (% of control mean)",
    y = ""
  ) +
  coord_cartesian(xlim = c(min(d_plot$norm_ci_lower, na.rm = TRUE) - 1,
                           max(d_plot$norm_ci_upper, na.rm = TRUE) + 1),
                  ylim = y_lim)

fig_coef <- paste0(paths$figures_dir, "fig.PA-type.threat_composite.coef.est.jpg")
ggsave(plot = p1, fig_coef, units = "cm", width = 22, height = 10, dpi = 300)
cat("Saved standalone threat_composite forest plot:", fig_coef, "\n")
}

##########################################
### spider plot: IUCN class effects ######
##########################################
{
cat("\nCreating spider plot for IUCN class effects...\n")

# Load group means for normalization
gm <- readRDS(paths$est_PA_type_group_means)

# Filter to IUCN classes + 14 individual threats + threat_composite
iucn_cats <- c("Strict", "Multi-use", "Not Reported")
spider_vars <- c(threat_list_raw, "threat_composite")

d_sp <- est[est$subsample %in% iucn_cats & est$variable %in% spider_vars, ]

# Merge with control means
gm_sub <- gm[gm$subsample %in% iucn_cats & gm$variable %in% spider_vars,
              c("subsample", "variable", "control_mean")]
d_sp <- merge(d_sp, gm_sub, by = c("subsample", "variable"))

# Normalize: coef / control_mean * 100 (% change relative to control mean)
d_sp$norm_coef <- d_sp$coef / d_sp$control_mean * 100

# Drop if control mean ~ 0 (normalization undefined)
d_sp <- d_sp[abs(d_sp$control_mean) > 1e-10, ]

# Signed-log transform: compresses extremes, linear near 0, preserves sign
# sign(x) * log10(1 + |x|) maps: 0→0, ±10→±1.04, ±50→±1.71, ±100→±2.00
slog <- function(x) sign(x) * log10(1 + abs(x))
slog_inv <- function(y) sign(y) * (10^abs(y) - 1)
d_sp$y_plot <- slog(d_sp$norm_coef)
cat(sprintf("  Norm range: [%.1f, %.1f]%% → signed-log: [%.2f, %.2f]\n",
            min(d_sp$norm_coef, na.rm = TRUE), max(d_sp$norm_coef, na.rm = TRUE),
            min(d_sp$y_plot, na.rm = TRUE), max(d_sp$y_plot, na.rm = TRUE)))

# Short threat labels for radial axes
short_labels <- c(
  "built"              = "Residential dev.",
  "cropland"           = "Cropland",
  "planted"            = "Plantations",
  "pasture"            = "Pasture",
  "oil"                = "Oil & gas",
  "mining"             = "Mining",
  "renewables"         = "Renewables",
  "roads"              = "Roads",
  "powerlines"         = "Power lines",
  "def0120_parea"      = "Deforestation",
  "fires"              = "Fire",
  "swu"                = "Water abstraction",
  "dams"               = "Dams",
  "light"              = "Light pollution",
  "threat_composite"   = "Threat composite"
)

# Order by IUCN threat category (same as combined figure)
threat_order <- names(short_labels)
n_threats <- length(threat_order)

# Numeric x positions (0 to n-1) for proper polygon closure in coord_polar
var_num_map <- setNames(seq_along(threat_order) - 1, threat_order)
d_sp$var_num <- var_num_map[as.character(d_sp$variable)]

# Significance indicator
d_sp$sig <- d_sp$pval < 0.05

# Factor ordering: draw Not Reported first (background), then Multi-use, then Strict (foreground)
# Legend order controlled separately via scale_color_manual(breaks=)
draw_order <- c("Not Reported", "Multi-use", "Strict")
d_sp$subsample <- factor(d_sp$subsample, levels = draw_order)

# Sort by subsample (draw order) then variable
d_sp <- d_sp[order(d_sp$subsample, d_sp$var_num), ]

# Close polygons: duplicate first point at var_num = n_threats (wraps to angle 0)
d_closed <- do.call(rbind, lapply(split(d_sp, d_sp$subsample), function(df) {
  first_row <- df[which.min(df$var_num), ]
  first_row$var_num <- n_threats
  rbind(df, first_row)
}))
rownames(d_closed) <- NULL

# Radial axis range on signed-log scale (with padding)
y_range <- c(min(d_sp$y_plot, na.rm = TRUE) * 1.1,
             max(d_sp$y_plot, na.rm = TRUE) * 1.15)

# Radial axis breaks: nice % values mapped through slog
pct_breaks <- c(-100, -50, -25, -10, 0, 10, 25, 50, 100)
pct_breaks <- pct_breaks[slog(pct_breaks) >= y_range[1] & slog(pct_breaks) <= y_range[2]]
y_breaks_trans <- slog(pct_breaks)
y_break_labels <- paste0(ifelse(pct_breaks > 0, "+", ""), pct_breaks, "%")

# IUCN category colors
iucn_colors <- c("Strict" = "#2166ac", "Multi-use" = "#b2182b", "Not Reported" = "#999999")

p_spider <- ggplot(d_closed, aes(x = var_num, y = y_plot,
                                  group = subsample, color = subsample)) +
  # shaded zones: green = threat reduction (inside 0), red = threat increase (outside 0)
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "#4daf4a", alpha = 0.07) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "#e41a1c", alpha = 0.07) +
  # reference circle at zero effect
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  # spider lines (closed via duplicated first point)
  geom_path(linewidth = 0.8) +
  # points only at actual positions (not the closing point)
  geom_point(data = d_sp, aes(x = var_num, shape = sig), size = 2.5) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), guide = "none") +
  scale_color_manual(values = iucn_colors, breaks = iucn_cats) +
  scale_x_continuous(breaks = 0:(n_threats - 1),
                     labels = sapply(threat_order, function(v) {
                       lab <- unname(short_labels[v])
                       if (v == "threat_composite") paste0("**", lab, "**") else lab
                     }),
                     limits = c(0, n_threats)) +
  scale_y_continuous(limits = y_range, expand = expansion(0),
                     breaks = y_breaks_trans, labels = y_break_labels) +
  coord_polar(start = 0, clip = "off") +
  # zone labels
  annotate("text", x = 1.2, y = y_range[2] * 0.55, label = "Threat\nincrease",
           size = 2.8, color = "#b22222", fontface = "italic", lineheight = 0.85) +
  # radial axis labels inside the circle (along top spoke)
  annotate("text", x = 0.3, y = y_breaks_trans, label = y_break_labels,
           size = 2.5, color = "grey30", hjust = 0) +
  theme_minimal() +
  theme(
    axis.text.x = ggtext::element_markdown(size = 8),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey30"),
    plot.caption = element_text(size = 7, color = "grey40")
  ) +
  labs(
    title = "PA effects by IUCN management class",
    subtitle = "ATT as % of control mean, signed-log scale (dashed circle = zero effect)",
    caption = "Filled points: p < 0.05."
  )

fig_spider <- paste0(paths$figures_dir, "fig.PA-type.spider.iucn.est.jpg")
ggsave(plot = p_spider, fig_spider, units = "cm", width = 22, height = 22, dpi = 300)
cat("Saved spider plot:", fig_spider, "\n")
}

cat("\n=== PA-type figures complete ===\n")
cat("Individual panels: pub/figures/fig.PA-type.[variable].est.jpg\n")
cat("Combined figure:   pub/figures/fig.PA-type.combined.est.jpg\n")
cat("Coef + CI:         pub/figures/fig.PA-type.threat_composite.coef.est.jpg\n")
cat("Spider plot:       pub/figures/fig.PA-type.spider.iucn.est.jpg\n")
