##########################################
### pathreat.analysis.bycountry.fig.R ####
### By-country effects figures ###########
##########################################

library(fst)
library(dplyr)
library(ggplot2)

source("code/pathreat.analysis.config.R")

######################
### load data ########
######################
{
# load by-country results
results_file <- paste0(paths$results_dir, "pathreat.bycountry.est.Rds")
if (!file.exists(results_file)) {
  stop("By-country results not found. Run pathreat.analysis.bycountry.est.R first.")
}
d <- readRDS(results_file)

# filter to any_threat, estimated countries only
d <- d[d$variable == "any_threat" & d$status == "estimated", ]

# rename columns to match plot code
d$estimate <- d$coef
d$p_value  <- d$pval
d$n_total  <- d$n_obs

d_est <- d
cat("Countries with estimates:", nrow(d_est), "\n")

# exclude extreme outlier (N < 50 with extreme effects) for cleaner plots
d_plot <- d_est[d_est$n_total >= 50, ]
n_excluded <- nrow(d_est) - nrow(d_plot)
if (n_excluded > 0) {
  cat("Excluded", n_excluded, "countries with N < 50 for plotting\n")
}
}

############################
### prepare plot data ######
############################
{
# add significance categories
d_plot$sig_cat <- "Not significant"
d_plot$sig_cat[d_plot$p_value <= 0.10] <- "p < 0.10"
d_plot$sig_cat[d_plot$p_value <= 0.05] <- "p < 0.05"
d_plot$sig_cat[d_plot$p_value <= 0.01] <- "p < 0.01"
d_plot$sig_cat <- factor(d_plot$sig_cat,
                         levels = c("p < 0.01", "p < 0.05", "p < 0.10", "Not significant"))

# add direction + significance for coloring
d_plot$effect_type <- "Not significant"
d_plot$effect_type[d_plot$p_value <= 0.05 & d_plot$estimate < 0] <- "Negative (p < 0.05)"
d_plot$effect_type[d_plot$p_value <= 0.05 & d_plot$estimate > 0] <- "Positive (p < 0.05)"
d_plot$effect_type <- factor(d_plot$effect_type,
                             levels = c("Negative (p < 0.05)", "Not significant", "Positive (p < 0.05)"))

# order by effect size for forest plot
d_plot <- d_plot[order(d_plot$estimate), ]
d_plot$rank <- 1:nrow(d_plot)

# summary stats
median_effect <- median(d_plot$estimate, na.rm = TRUE)
mean_effect <- mean(d_plot$estimate, na.rm = TRUE)
cat("Median effect:", round(median_effect, 4), "\n")
cat("Mean effect:", round(mean_effect, 4), "\n")
}

####################################
### 1. Forest plot (significance) ##
####################################
{
cat("\nCreating forest plot...\n")

# color palette for significance
sig_colors <- c("p < 0.01" = "darkred",
                "p < 0.05" = "orangered",
                "p < 0.10" = "orange",
                "Not significant" = "gray60")

p1 <- ggplot(d_plot, aes(x = estimate, y = rank, color = sig_cat)) +

  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = median_effect, linetype = "dashed", color = "blue", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.3, alpha = 0.7) +
  geom_point(size = 1.5) +
  scale_y_continuous(
    breaks = d_plot$rank,
    labels = d_plot$country_name,
    expand = c(0.01, 0.01)
  ) +
  scale_color_manual(values = sig_colors, name = "Significance") +
  labs(
    x = "Effect on threat presence (p.p.)",
    y = "",
    title = "Protected area effects by country",
    subtitle = paste0("N = ", nrow(d_plot), " countries | Dashed line = median effect (",
                      round(median_effect * 100, 1), " p.p.)")
  ) +
  coord_cartesian(xlim = c(min(d_plot$ci_low, na.rm = TRUE) - 0.02,
                           max(d_plot$ci_high, na.rm = TRUE) + 0.02))

# save forest plot
fig1_file <- paste0(paths$figures_dir, "fig.bycountry.forest.est.jpg")
ggsave(plot = p1, fig1_file, units = "cm", width = 18, height = 40, dpi = 300)
cat("Saved:", fig1_file, "\n")
}

####################################
### 2. Histogram + density #########
####################################
{
cat("\nCreating histogram/density plot...\n")

p2 <- ggplot(d_plot, aes(x = estimate)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  ) +
  geom_histogram(aes(y = after_stat(density)), bins = 25, fill = "steelblue",
                 color = "white", alpha = 0.7) +
  geom_density(color = "darkblue", linewidth = 1) +
  geom_rug(alpha = 0.5, color = "darkblue") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = median_effect, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = 0.02, y = Inf, label = "Null", hjust = 0, vjust = 2, fontface = "bold") +
  annotate("text", x = median_effect + 0.02, y = Inf, label = paste0("Median\n(", round(median_effect * 100, 1), " p.p.)"),
           hjust = 0, vjust = 1.5, color = "red", size = 3.5) +
  labs(
    x = "Effect on threat presence (p.p.)",
    y = "Density",
    title = "Distribution of country-level PA effects",
    subtitle = paste0("N = ", nrow(d_plot), " countries | ",
                      sum(d_plot$estimate < 0), " negative, ",
                      sum(d_plot$estimate > 0), " positive | ",
                      sum(d_plot$p_value <= 0.05), " significant at 5%")
  )

# save histogram
fig2_file <- paste0(paths$figures_dir, "fig.bycountry.histogram.est.jpg")
ggsave(plot = p2, fig2_file, units = "cm", width = 20, height = 14, dpi = 300)
cat("Saved:", fig2_file, "\n")
}

################################################
### 3. Funnel plot (sample size version) #######
################################################
{
cat("\nCreating funnel plot (sample size)...\n")

# effect type colors
effect_colors <- c("Negative (p < 0.05)" = "darkred",
                   "Not significant" = "gray50",
                   "Positive (p < 0.05)" = "darkblue")

# identify countries to label (largest samples + extreme effects)
d_plot$label_flag <- FALSE
# top 5 by sample size
top_n <- d_plot[order(-d_plot$n_total), ][1:5, "country_name"]
d_plot$label_flag[d_plot$country_name %in% top_n] <- TRUE
# most extreme negative effects (top 3)
extreme_neg <- d_plot[order(d_plot$estimate), ][1:3, "country_name"]
d_plot$label_flag[d_plot$country_name %in% extreme_neg] <- TRUE
# most extreme positive effects (top 3, excluding tiny samples)
extreme_pos <- d_plot[d_plot$estimate > 0 & d_plot$n_total > 1000, ]
extreme_pos <- extreme_pos[order(-extreme_pos$estimate), ][1:3, "country_name"]
d_plot$label_flag[d_plot$country_name %in% extreme_pos] <- TRUE

p3 <- ggplot(d_plot, aes(x = estimate, y = n_total, color = effect_type)) +

  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  geom_hline(yintercept = c(100, 1000, 10000, 100000, 1000000),
             linetype = "dotted", color = "gray70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text(data = d_plot[d_plot$label_flag, ],
            aes(label = country_name), hjust = -0.1, vjust = 0.5, size = 2.5,
            show.legend = FALSE) +
  scale_y_log10(
    labels = scales::label_comma(),
    breaks = c(100, 1000, 10000, 100000, 1000000)
  ) +
  scale_color_manual(values = effect_colors, name = "Effect") +
  labs(
    x = "Effect on threat presence (p.p.)",
    y = "Sample size (log scale)",
    title = "Effect size vs. sample size by country",
    subtitle = "Larger samples (more reliable) at top | Labels: largest samples + extreme effects"
  )

# save funnel plot sample size
fig3_file <- paste0(paths$figures_dir, "fig.bycountry.funnel.n.est.jpg")
ggsave(plot = p3, fig3_file, units = "cm", width = 22, height = 16, dpi = 300)
cat("Saved:", fig3_file, "\n")
}

##################################################
### 4. Coverage vs effectiveness heatmap #########
##################################################
{
cat("\nCreating coverage vs effectiveness heatmap...\n")

# load PA coverage from unmatched data
pa_coverage <- fst::read_fst(
  paths$data_unmatched,
  columns = c("country_rast", "ever_pa")
) %>%
  group_by(country_rast) %>%
  dplyr::summarise(
    pct_pa = mean(ever_pa, na.rm = TRUE) * 100,
    pa_area_km2 = sum(ever_pa == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ungroup()

# load threat_composite estimates
d_tc <- readRDS(results_file) %>%
  filter(variable == "threat_composite")

# compute % change from control
d_tc_est <- d_tc %>%
  filter(
    status == "estimated",
    control_mean >= 0.001
  ) %>%
  mutate(pct_effect = coef / control_mean * 100)

# merge PA coverage
d_tc_est <- d_tc_est %>%
  left_join(pa_coverage, by = c("country_id" = "country_rast"))

# unassessed: countries in PA coverage but not in estimated set
assessed_ids <- d_tc_est$country_id
d_unassessed <- pa_coverage %>%
  filter(!(country_rast %in% assessed_ids))

n_assessed <- nrow(d_tc_est)
n_unassessed <- nrow(d_unassessed)
cat("  Assessed countries:", n_assessed, "\n")
cat("  Unassessed countries:", n_unassessed, "\n")

# bin parameters
bw_x <- 10
bw_y <- 10
x_cap_lo <- -40
x_cap_hi <-  40

# clamp pct_effect into extreme bins (don't drop — small sample)
y_cap_hi <- 30
d_tc_est <- d_tc_est %>%
  mutate(
    pct_effect_capped = pmax(x_cap_lo, pmin(x_cap_hi - 0.001, pct_effect)),
    pct_pa_capped = ifelse(pct_pa >= y_cap_hi, y_cap_hi, pct_pa)
  )

# bin counts for main grid
heat_counts <- d_tc_est %>%
  mutate(
    xb = floor(pct_effect_capped / bw_x) * bw_x,
    yb = floor(pct_pa_capped / bw_y) * bw_y
  ) %>%
  count(xb, yb)

# unassessed column bins
unassessed_bins <- d_unassessed %>%
  mutate(yb = floor(ifelse(pct_pa >= y_cap_hi, y_cap_hi, pct_pa) / bw_y) * bw_y) %>%
  count(yb)

# fill range across both grids
fill_range <- range(c(heat_counts$n, unassessed_bins$n))

# NA column position
na_x <- x_cap_lo - 2 * bw_x

p4 <- ggplot() +
  # main heatmap grid
  geom_tile(
    data = heat_counts,
    aes(
      x = xb + bw_x / 2,
      y = yb + bw_y / 2,
      fill = n
    ),
    width = bw_x,
    height = bw_y
  ) +
  # NA column: unassessed countries
  geom_rect(
    data = unassessed_bins,
    aes(
      xmin = na_x,
      xmax = na_x + bw_x,
      ymin = yb,
      ymax = yb + bw_y,
      fill = n
    )
  ) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey45") +
  scale_x_continuous(
    breaks = c(
      na_x + bw_x / 2,
      seq(x_cap_lo + bw_x / 2, x_cap_hi - bw_x / 2, bw_x)
    ),
    labels = function(x) {
      ifelse(
        abs(x - (na_x + bw_x / 2)) < 1, "NA",
        ifelse(x == x_cap_lo + bw_x / 2, paste0("<", x_cap_lo, "%"),
        ifelse(x == x_cap_hi - bw_x / 2, paste0(">", x_cap_hi, "%"),
        paste0(x - bw_x / 2, " \u2013 ", x + bw_x / 2, "%"))))
    },
    limits = c(na_x, x_cap_hi)
  ) +
  scale_fill_viridis_c(
    option = "inferno",
    name = "Number of countries",
    limits = fill_range
  ) +
  scale_y_continuous(
    breaks = c(seq(bw_y / 2, y_cap_hi - bw_y / 2, bw_y), y_cap_hi + bw_y / 2),
    labels = function(y) {
      ifelse(y == y_cap_hi + bw_y / 2,
        paste0(">", y_cap_hi, "%"),
        paste0(y - bw_y / 2, " \u2013 ", y + bw_y / 2, "%"))
    },
    limits = c(0, y_cap_hi + bw_y)
  ) +
  labs(
    x = "Country-level PA effectiveness (% change from control)",
    y = "Share of country area protected (%)"
  ) +
  coord_fixed(ratio = bw_x / bw_y, clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barwidth = 15
  ))

fig4_file <- paste0(paths$figures_dir, "fig.bycountry.coverage-vs-effectiveness.heatmap.jpg")
ggsave(plot = p4, fig4_file, width = 10, height = 6, dpi = 300)
cat("Saved:", fig4_file, "\n")
}

##################################################
### 5. Coverage vs threat-reduction bubble plot ##
##################################################
# Country panel of manuscript Figure 4 (companion panel lives in
# pathreat.analysis.byspecies.fig.R). Y-limits are hard-coded to match the
# species panel so the two halves of Figure 4 are directly comparable.
{
cat("\nCreating coverage vs threat-reduction bubble plot...\n")

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
if (has_ggrepel) library(ggrepel)

# reuse d_tc_est from the heatmap section (country-level threat_composite ATT
# joined to all-vintage PA coverage and area). Compute threat reduction as −coef
# on the 0–1 scale. Each equal-area raster cell is 1 km², so pa_area_km2 is
# the protected union area in the same country universe used for the x-axis.
cdat <- d_tc_est %>%
  dplyr::mutate(
    reduction  = -coef,
    effect_abs = abs(coef),
    sig = dplyr::case_when(
      pval < 0.05 & coef < 0 ~ "Effective (p<0.05)",
      pval < 0.05 & coef > 0 ~ "Counterproductive (p<0.05)",
      TRUE                    ~ "Not significant"
    )
  )

if (anyNA(cdat$pct_pa) || anyNA(cdat$pa_area_km2)) {
  stop("Missing PA coverage or total PA area for an estimated country")
}
if (any(cdat$pa_area_km2 <= 0)) {
  stop("Estimated countries must have positive total PA area")
}
cat("  Countries in bubble plot:", nrow(cdat), "\n")
cat("  Total PA area range (km²):",
    paste(range(cdat$pa_area_km2), collapse = " to "), "\n")

bio_top <- c(
  "Brazil", "Colombia", "Peru", "Indonesia", "Mexico", "Ecuador",
  "Madagascar", "Democratic Republic of the Congo", "India", "China",
  "Australia", "Venezuela", "Bolivia", "Philippines"
)
cdat$label_me <- cdat$country_name %in% bio_top |
  rank(-cdat$effect_abs, ties.method = "first") <= 4 |
  (cdat$pct_pa > 30 & cdat$reduction < 0)

sig_pal <- c(
  "Effective (p<0.05)"        = "#1a9850",
  "Not significant"           = "grey60",
  "Counterproductive (p<0.05)" = "#d73027"
)

# shared y-limits with species panel
y_lim_fig4 <- c(-0.075, 0.18)

p5 <- ggplot(cdat, aes(x = pct_pa, y = reduction)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.35) +
  geom_vline(xintercept = 30, colour = "grey20", linewidth = 0.7,
             linetype = "dotdash") +
  annotate("text", x = 30, y = y_lim_fig4[2], label = "30×30 target",
           hjust = -0.05, vjust = 1, size = 3.0, fontface = "bold",
           colour = "grey20") +
  geom_point(aes(size = pa_area_km2, fill = sig),
             shape = 21, colour = "white", stroke = 0.3, alpha = 0.85) +
  scale_size_area(max_size = 14,
                  name = "Total protected area\n(km²)",
                  breaks = pretty(cdat$pa_area_km2, n = 4)) +
  scale_fill_manual(values = sig_pal, name = NULL) +
  guides(size = "none", fill = "none") +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0.04, 0.08))) +
  scale_y_continuous(limits = y_lim_fig4,
                     breaks = scales::pretty_breaks(n = 5)) +
  labs(x = "Share of country area protected",
       y = "Threat reduction from PAs (0–1 scale)") +
  theme_minimal(base_size = 10) +
  theme(legend.position   = "right",
        legend.box        = "vertical",
        legend.key.size   = unit(0.4, "cm"),
        panel.grid.minor  = element_blank(),
        plot.margin       = margin(t = 4, r = 5, b = 4, l = 5))

if (has_ggrepel) {
  p5 <- p5 + geom_text_repel(
    data = cdat[cdat$label_me, ],
    aes(label = country_name), size = 2.6, colour = "grey15",
    min.segment.length = 0, segment.size = 0.25, segment.colour = "black",
    max.overlaps = Inf, box.padding = 0.55, point.padding = 0.25,
    force = 4, force_pull = 0.5, max.iter = 20000, max.time = 5,
    seed = 42
  )
} else {
  p5 <- p5 + geom_text(data = cdat[cdat$label_me, ],
                       aes(label = country_name), size = 2.6,
                       colour = "grey15", vjust = -0.9)
}

fig5_file <- paste0(paths$figures_dir, "fig.bycountry.coverage-vs-reduction.jpg")
save_figure_data(cdat, fig5_file)
ggsave(plot = p5, fig5_file, width = 11, height = 10, units = "cm", dpi = 300)
cat("Saved:", fig5_file, "\n")
}

#############################################################
### 6. Coverage vs percentage threat-reduction bubble plot ##
#############################################################
# Presentation version of the country bubble plot. This keeps the raw-scale
# manuscript figure above unchanged and rescales each country's estimate by
# its matched-control mean.
{
cat("\nCreating coverage vs percentage threat-reduction bubble plot...\n")

cdat_pct <- cdat %>%
  dplyr::mutate(
    pct_reduction = -100 * coef / control_mean,
    pct_effect_abs = abs(pct_reduction)
  )

if (any(!is.finite(cdat_pct$pct_reduction))) {
  stop("Percentage threat reductions must be finite")
}

# Preserve the raw figure's biodiversity anchors and high-coverage labels,
# while ranking the four effect extremes on the percentage scale displayed.
cdat_pct$label_me_pct <- cdat_pct$country_name %in% bio_top |
  rank(-cdat_pct$pct_effect_abs, ties.method = "first") <= 4 |
  (cdat_pct$pct_pa > 30 & cdat_pct$pct_reduction < 0)

pct_y_step <- 10
pct_y_lim <- c(
  -50,
  ceiling(max(cdat_pct$pct_reduction) / pct_y_step) * pct_y_step
)
pct_below_axis <- cdat_pct$pct_reduction < pct_y_lim[1]
pct_label_data <- cdat_pct[
  cdat_pct$label_me_pct & !pct_below_axis,
]

cat("  Countries in percentage bubble plot:", nrow(cdat_pct), "\n")
cat(
  "  Percentage reduction range:",
  paste(round(range(cdat_pct$pct_reduction), 1), collapse = "% to "),
  "%\n"
)
cat(
  "  Displayed y-axis:",
  paste(pct_y_lim, collapse = "% to "),
  "%\n"
)
cat("  Countries below display:", sum(pct_below_axis), "\n")

p6 <- ggplot(cdat_pct, aes(x = pct_pa, y = pct_reduction)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.35) +
  geom_vline(xintercept = 30, colour = "grey20", linewidth = 0.7,
             linetype = "dotdash") +
  annotate("text", x = 30, y = pct_y_lim[2], label = "30×30 target",
           hjust = -0.05, vjust = 1.2, size = 3.0, fontface = "bold",
           colour = "grey20") +
  geom_point(aes(size = pa_area_km2, fill = sig),
             shape = 21, colour = "white", stroke = 0.3, alpha = 0.85) +
  scale_size_area(max_size = 14,
                  name = "Total protected area\n(km²)",
                  breaks = pretty(cdat_pct$pa_area_km2, n = 4)) +
  scale_fill_manual(values = sig_pal, name = NULL) +
  guides(size = "none", fill = "none") +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0.04, 0.08))) +
  scale_y_continuous(
    breaks = seq(pct_y_lim[1], pct_y_lim[2], by = 25),
    labels = scales::label_number(accuracy = 1, suffix = "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = "Share of country area protected",
       y = "Threat reduction from PAs\n(% change from matched-control mean)") +
  coord_cartesian(ylim = pct_y_lim, clip = "off") +
  theme_minimal(base_size = 10) +
  theme(legend.position   = "right",
        legend.box        = "vertical",
        legend.key.size   = unit(0.4, "cm"),
        panel.grid.minor  = element_blank(),
        plot.margin       = margin(t = 4, r = 5, b = 4, l = 5))

if (has_ggrepel) {
  p6 <- p6 + geom_text_repel(
    data = pct_label_data,
    aes(label = country_name), size = 2.6, colour = "grey15",
    min.segment.length = 0, segment.size = 0.25, segment.colour = "black",
    max.overlaps = Inf, box.padding = 0.55, point.padding = 0.25,
    force = 4, force_pull = 0.5, max.iter = 20000, max.time = 5,
    seed = 42
  )
} else {
  p6 <- p6 + geom_text(data = pct_label_data,
                       aes(label = country_name), size = 2.6,
                       colour = "grey15", vjust = -0.9)
}

fig6_file <- paste0(
  paths$figures_dir,
  "fig.bycountry.coverage-vs-pct-reduction.jpg"
)
ggsave(plot = p6, fig6_file, width = 11, height = 10, units = "cm", dpi = 300)
cat("Saved:", fig6_file, "\n")
}

####################################
### summary ########################
####################################
{
cat("\n=== By-country figures complete ===\n")
cat("1. Forest plot:", fig1_file, "\n")
cat("2. Histogram:", fig2_file, "\n")
cat("3. Funnel (N):", fig3_file, "\n")
cat("4. Coverage heatmap:", fig4_file, "\n")
cat("5. Coverage-vs-reduction:", fig5_file, "\n")
cat("6. Coverage-vs-percentage-reduction:", fig6_file, "\n")
}
