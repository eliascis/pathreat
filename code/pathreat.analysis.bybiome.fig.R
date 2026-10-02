##########################################
### pathreat.analysis.bybiome.fig.R ######
### By-biome effects figures #############
### Forest plot (threat_composite) #######
##########################################

library(dplyr)
library(ggplot2)
library(fst)

source("code/pathreat.analysis.config.R")

######################
### load data ########
######################
{
results_file <- paths$est_bybiome
if (!file.exists(results_file)) {
  stop("By-biome results not found. Run pathreat.analysis.bybiome.est.R first.")
}
results_bybiome <- readRDS(results_file)
cat(sprintf("Loaded %d rows from %s\n", nrow(results_bybiome), results_file))
}

####################################
### Forest plot (threat_composite) #
####################################
{
cat("\nCreating forest plot (threat_composite by biome)...\n")

d_forest <- results_bybiome[results_bybiome$variable == "threat_composite" &
                              results_bybiome$status == "estimated", ]

# significance tiers
d_forest$sig_tier <- "Not significant"
d_forest$sig_tier[d_forest$pval < 0.05] <- "p < 0.05"
d_forest$sig_tier[d_forest$pval < 0.01] <- "p < 0.01"
d_forest$sig_tier <- factor(d_forest$sig_tier,
                            levels = c("p < 0.01", "p < 0.05", "Not significant"))

# normalize: coef / control_mean * 100 = % of control mean
d_forest$norm_coef <- d_forest$coef / d_forest$control_mean * 100
d_forest$norm_ci_lower <- d_forest$ci_low / d_forest$control_mean * 100
d_forest$norm_ci_upper <- d_forest$ci_high / d_forest$control_mean * 100

# order by normalized effect size
d_forest <- d_forest[order(d_forest$norm_coef), ]
d_forest$rank <- seq_len(nrow(d_forest))
# clean biome labels (PA area shown in separate panel)
d_forest$biome_label <- d_forest$biome_short

# add global estimate at top
g <- readRDS(paths$est_global)
g <- g[g$variable == "threat_composite", ]
global_n_treat <- g$n_obs / 2  # matched pairs: half are treated
global_rank <- max(d_forest$rank) + 2  # +2 for visual gap
global_row <- data.frame(
  norm_coef = g$coef / g$control_mean * 100,
  norm_ci_lower = g$ci_low / g$control_mean * 100,
  norm_ci_upper = g$ci_high / g$control_mean * 100,
  pval = g$pval,
  sig_tier = factor(
    if (g$pval < 0.01) "p < 0.01" else if (g$pval < 0.05) "p < 0.05" else "Not significant",
    levels = c("p < 0.01", "p < 0.05", "Not significant")),
  rank = global_rank,
  n_treat = global_n_treat,
  biome_label = "Global"
)
# mark global row for larger point size
global_row$is_global <- TRUE
biome_cols <- c("biome_id", "norm_coef", "norm_ci_lower", "norm_ci_upper", "pval", "sig_tier", "rank", "n_treat", "biome_label")
d_biome <- d_forest[, biome_cols]
d_biome$is_global <- FALSE
global_row$biome_id <- NA
d_forest <- rbind(d_biome, global_row[, c(biome_cols, "is_global")])

y_lim <- range(d_forest$rank) + c(-0.8, 0.8)

p1 <- ggplot(d_forest, aes(x = norm_coef, y = rank)) +
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
    breaks = d_forest$rank,
    labels = d_forest$biome_label
  ) +
  guides(shape = guide_legend(override.aes = list(size = 3)),
         fill = guide_legend(override.aes = list(size = 3))) +
  labs(
    x = "Effect on overall threat index (% of control mean)",
    y = ""
  ) +
  coord_cartesian(xlim = c(min(d_forest$norm_ci_lower, na.rm = TRUE) - 1,
                           max(d_forest$norm_ci_upper, na.rm = TRUE) + 1),
                  ylim = y_lim)

fig_file <- paste0(paths$figures_dir, "fig.bybiome.forest.est.jpg")
ggsave(plot = p1, fig_file, units = "cm", width = 22, height = 14, dpi = 300)
cat("Saved:", fig_file, "\n")
}

####################################
### Protection share by biome ######
### (LaTeX table) ##################
####################################
{
cat("\nComputing protection shares by biome...\n")

# load unmatched data
if (!"d.du" %in% ls()) {
  d.du <- read_fst(paths$data_unmatched, columns = c("biome_raster", "ever_pa"))
}

# biome lookup (same as estimation script)
biome_lookup <- data.frame(
  biome_id = 1:14,
  biome_short = c(
    "Trop. moist forests", "Trop. dry forests", "Trop. conif. forests",
    "Temp. broadleaf forests", "Temp. conif. forests", "Boreal forests/taiga",
    "Trop. grasslands", "Temp. grasslands", "Flooded grasslands",
    "Montane grasslands", "Tundra", "Mediterranean",
    "Deserts \\& xeric", "Mangroves"
  ),
  stringsAsFactors = FALSE
)

# exclude non-biome pixels (98 = rock/ice, 99 = unclassified)
d_bio <- d.du[d.du$biome_raster %in% 1:14, ]

# by-biome shares
bio_stats <- d_bio %>%
  group_by(biome_raster) %>%
  summarise(
    n_land  = n(),
    n_pa    = sum(ever_pa == 1),
    .groups = "drop"
  ) %>%
  mutate(
    share_pa    = 100 * n_pa / n_land,
    share_land  = 100 * n_land / sum(n_land),
    share_pa_global = 100 * n_pa / sum(n_pa)
  ) %>%
  arrange(biome_raster)

# global totals
global_land <- sum(bio_stats$n_land)
global_pa   <- sum(bio_stats$n_pa)
global_share <- 100 * global_pa / global_land

# merge biome names
bio_stats <- merge(bio_stats, biome_lookup, by.x = "biome_raster", by.y = "biome_id", all.x = TRUE)
bio_stats <- bio_stats[order(bio_stats$biome_raster), ]

# build LaTeX table body (for \input inside tabular in main file)
# columns: Biome | Area (10^3 km^2) | PA area (10^3 km^2) | PA share (%) | Share of land (%) | Share of all PAs (%)
tex_lines <- character()

for (i in seq_len(nrow(bio_stats))) {
  r <- bio_stats[i, ]
  tex_lines <- c(tex_lines, sprintf(
    "%s & %.1f & %.1f & %.1f & %.1f & %.1f \\\\",
    r$biome_short,
    r$n_land / 1000,
    r$n_pa / 1000,
    r$share_pa,
    r$share_land,
    r$share_pa_global
  ))
}

# add global row
tex_lines <- c(tex_lines, "\\midrule")
tex_lines <- c(tex_lines, sprintf(
  "\\textbf{Global} & %.1f & %.1f & %.1f & {100.0} & {100.0}",
  global_land / 1000,
  global_pa / 1000,
  global_share
))

# write table file (no trailing \\)
outfile <- "pub/tables/tab.biome_protection.tex"
writeLines(tex_lines, outfile)
cat("Saved:", outfile, "\n")

rm(d_bio, bio_stats); gc()
}

####################################
### summary ########################
####################################
{
cat("\n=== By-biome figures complete ===\n")
cat("Forest plot:", fig_file, "\n")
cat("Protection table:", "pub/tables/tab.biome_protection.tex\n")
}
