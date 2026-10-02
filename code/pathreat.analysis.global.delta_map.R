##################################################
### pathreat.analysis.global.delta_map.R ############
### Global raster map of PA effectiveness ########
##################################################
#
# Purpose: Global rasterized map of delta (pair-level treatment effect on
#          threat_composite) for matched PA pixels, with country boundaries.
#
# Input:   data/store/pathreat.data.merge.matched.fst
#          results/pathreat.bycountry.est.Rds
#
# Output:  pub/figures/fig.delta_map.jpg
#          pub/figures/fig.delta_map.country.jpg
#
# Notes:   delta = threat_composite[treated] - threat_composite[control]
#          Negative (blue) = PA reduces threats
#          Positive (red)  = PA increases threats
#
library(dplyr)
library(ggplot2)
library(fst)
library(sf)
library(rnaturalearth)

source("code/pathreat.analysis.config.R")

cat("=== Global Delta Map ===\n\n")

#########################################
### STEP 1: Load matched data ##########
#########################################
{
cat("STEP 1: Loading matched data...\n")

# load only needed columns (memory-efficient)
d <- read_fst(paths$data_matched,
              columns = c("x", "y", "delta", "treat", "country_rast", "threat_composite"))

# balance filter pre-applied in merge.matched.fst

# compute global control mean for normalization
ctrl_mean <- mean(d$threat_composite[d$treat == 0], na.rm = TRUE)
cat("  Control mean (threat_composite):", round(ctrl_mean, 4), "\n")

# keep treated pixels, normalize delta as % of control mean
d <- d[d$treat == 1, c("x", "y", "delta")]
d$delta_pct <- (d$delta / ctrl_mean) * 100

cat("  Treated pixels:", format(nrow(d), big.mark = ","), "\n")
cat("  Delta % range: [", round(min(d$delta_pct, na.rm = TRUE), 1), ",",
    round(max(d$delta_pct, na.rm = TRUE), 1), "]\n")
cat("  Delta % mean:", round(mean(d$delta_pct, na.rm = TRUE), 1), "\n")
}

#########################################
### STEP 2: Aggregate to grid ##########
#########################################
{
cat("\nSTEP 2: Aggregating to plotting grid...\n")

# coordinates are in World Mollweide (meters), 1km original resolution
# aggregate to 25km grid for visible tiles on global map
grid_res <- 25000
d$x_grid <- round(d$x / grid_res) * grid_res
d$y_grid <- round(d$y / grid_res) * grid_res

df_plot <- d %>%
  group_by(x = x_grid, y = y_grid) %>%
  summarise(delta_pct = mean(delta_pct, na.rm = TRUE), .groups = "drop")

cat("  Grid resolution:", grid_res / 1000, "km\n")
cat("  Grid cells:", format(nrow(df_plot), big.mark = ","),
    "(from", format(nrow(d), big.mark = ","), "pixels)\n")

rm(d)
gc()
}

#########################################
### STEP 4: Country boundaries #########
#########################################
{
cat("\nSTEP 4: Loading country boundaries...\n")

# Mollweide CRS (same as data coordinates)
moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"

world <- ne_countries(scale = "small", returnclass = "sf")
world <- st_transform(world, moll_crs)
cat("  Countries:", nrow(world), "(projected to Mollweide)\n")
}

#########################################
### STEP 5: Plot #######################
#########################################
{
cat("\nSTEP 5: Creating figure...\n")

# bin delta (%) into discrete classes, split at zero
brks <- c(-Inf, -30, -20, -10, -2, 0, 2, 10, 20, 30, Inf)
labs <- c("< -30%", "-30 to -20%", "-20 to -10%", "-10 to -2%",
          "-2 to 0%", "0 to 2%", "2 to 10%", "10 to 20%", "20 to 30%", "> 30%")
# blue (all negative) | yellow-orange-red (all positive)
bin_colors <- c("#08519C", "#3182BD", "#6BAED6", "#9ECAE1", "#C6DBEF",
                "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
df_plot$delta_bin <- cut(df_plot$delta_pct, breaks = brks, labels = labs)
cat("  Bin counts:\n")
print(table(df_plot$delta_bin))

p <- ggplot() +
  # country boundaries first
  geom_sf(data = world, fill = "#E8E8E8", color = "#AAAAAA", linewidth = 0.1) +
  # delta tiles on top
  geom_tile(data = df_plot, aes(x = x, y = y, fill = delta_bin),
            width = grid_res, height = grid_res) +
  # discrete RdYlBu scale (yellow midpoint visible against gray land)
  scale_fill_manual(
    values = bin_colors,
    name = "% change from control mean",
    drop = FALSE
  ) +
  # crop to land extent (trim empty Mollweide corners)
  coord_sf(crs = moll_crs,
           xlim = c(-14000000, 17000000), ylim = c(-6500000, 8500000),
           expand = FALSE) +
  theme_void(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    legend.key.width = unit(0.8, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 7),
    plot.margin = margin(2, 2, 2, 2, "pt")
  ) +
  guides(fill = guide_legend(nrow = 1, title.position = "top"))

out_path <- paste0(paths$figures_dir, "fig.delta_map.jpg")
save_figure_data(df_plot, out_path)
ggsave(out_path, p, width = 10, height = 5, dpi = 300, bg = "white")
cat("  Saved to:", out_path, "\n")
}

##################################################
### STEP 6: Country-level ATT choropleth ########
##################################################
{
cat("\nSTEP 6: Country-level ATT choropleth...\n")

# load by-country ATT for threat_composite
est <- readRDS(paths$est_bycountry)
est_tc <- est[est$variable == "threat_composite" & est$status == "estimated",
              c("country_id", "country_name", "coef", "pval", "control_mean")]

# percentage change relative to country-level control mean
# winsorize: skip countries with near-zero control mean (< 0.001)
est_tc$coef_pct <- ifelse(
  est_tc$control_mean >= 0.001,
  (est_tc$coef / est_tc$control_mean) * 100,
  NA_real_
)
n_skip <- sum(is.na(est_tc$coef_pct))
if (n_skip > 0) {
  cat("  Skipped", n_skip, "countries with near-zero control mean:",
      paste(est_tc$country_name[is.na(est_tc$coef_pct)], collapse = ", "), "\n")
}

# build country_rast -> gid_0 lookup from matched data
d_lu <- read_fst(paths$data_matched,
                 columns = c("country_rast", "gid_0"))
country_lu <- unique(d_lu[, c("country_rast", "gid_0")])
rm(d_lu)

# merge gid_0 into estimates
est_tc <- merge(est_tc, country_lu,
                by.x = "country_id",
                by.y = "country_rast")
est_tc$sig <- est_tc$pval < 0.05

cat("  Countries estimated:", nrow(est_tc), "\n")
cat("  Significant (p<0.05):", sum(est_tc$sig), "\n")
cat("  Coef % range: [", round(min(est_tc$coef_pct, na.rm = TRUE), 1), ",",
    round(max(est_tc$coef_pct, na.rm = TRUE), 1), "]\n")
cat("  Coef % mean:", round(mean(est_tc$coef_pct, na.rm = TRUE), 1), "\n")

# join to Natural Earth polygons (adm0_a3 covers NOR/FRA where iso_a3 = -99)
world_est <- merge(world, est_tc[, c("gid_0", "coef_pct", "sig")],
                   by.x = "adm0_a3",
                   by.y = "gid_0",
                   all.x = TRUE)

# split into layers
world_bg  <- world_est[is.na(world_est$coef_pct), ]
world_sig <- world_est[!is.na(world_est$coef_pct) & world_est$sig == TRUE, ]
world_ns  <- world_est[!is.na(world_est$coef_pct) & world_est$sig == FALSE, ]

# create diagonal hatching lines for non-significant countries
hatch_df <- NULL
if (nrow(world_ns) > 0) {
  ns_union <- st_union(world_ns)
  bbox <- st_bbox(ns_union)
  span_y <- bbox[["ymax"]] - bbox[["ymin"]]
  spacing <- 150000  # 150 km in Mollweide meters
  x_seq <- seq(bbox[["xmin"]] - span_y,
               bbox[["xmax"]],
               by = spacing)
  lines_list <- lapply(x_seq, function(x0) {
    st_linestring(matrix(c(x0, bbox[["ymin"]],
                           x0 + span_y, bbox[["ymax"]]),
                         ncol = 2, byrow = TRUE))
  })
  hatch_sfc <- st_sfc(lines_list, crs = st_crs(world_ns))
  hatch_clipped <- st_intersection(hatch_sfc, ns_union)
  hatch_df <- st_sf(geometry = hatch_clipped)
  cat("  Hatching lines for", nrow(world_ns), "non-significant countries\n")
}

# same bins as pixel-level map (% change), split at zero
brks_c <- c(-Inf, -30, -20, -10, -2, 0, 2, 10, 20, 30, Inf)
labs_c <- c("< -30%", "-30 to -20%", "-20 to -10%", "-10 to -2%",
            "-2 to 0%", "0 to 2%", "2 to 10%", "10 to 20%", "20 to 30%", "> 30%")
# same colors as pixel-level map
bin_colors_c <- c("#08519C", "#3182BD", "#6BAED6", "#9ECAE1", "#C6DBEF",
                  "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")

# apply bins — include "Not in sample" and "Not significant" as fill levels
all_levels <- c(labs_c, "Not in sample", "Not significant")
all_colors <- c(bin_colors_c, "#E8E8E8", "white")
names(all_colors) <- all_levels

world_sig$coef_bin <- factor(
  cut(world_sig$coef_pct, breaks = brks_c, labels = labs_c),
  levels = all_levels)
world_ns$coef_bin <- factor(
  cut(world_ns$coef_pct, breaks = brks_c, labels = labs_c),
  levels = all_levels)
world_bg$coef_bin <- factor("Not in sample", levels = all_levels)

cat("  Bin counts (significant):\n")
print(table(world_sig$coef_bin))
cat("  Bin counts (non-significant):\n")
print(table(world_ns$coef_bin))

# custom legend key: white box with diagonal hatching lines
draw_key_hatch <- function(data, params, size) {
  grid::grobTree(
    grid::rectGrob(gp = grid::gpar(fill = "white", col = "gray30", lwd = 0.5)),
    grid::linesGrob(x = c(0, 0.5), y = c(0, 1),
                    gp = grid::gpar(col = "gray30", lwd = 0.5, lty = "dashed")),
    grid::linesGrob(x = c(0, 1), y = c(0, 1),
                    gp = grid::gpar(col = "gray30", lwd = 0.5, lty = "dashed")),
    grid::linesGrob(x = c(0.5, 1), y = c(0, 1),
                    gp = grid::gpar(col = "gray30", lwd = 0.5, lty = "dashed"))
  )
}

# dummy sf points to force all levels into the legend (invisible on map)
# one point per level with NA coordinates — ensures empty bins still get legend keys
regular_levels <- setdiff(all_levels, "Not significant")
dummy_regular <- st_sf(
  coef_bin = factor(regular_levels, levels = all_levels),
  geometry = st_sfc(rep(list(st_point(c(NA_real_, NA_real_))), length(regular_levels)),
                    crs = st_crs(world))
)
dummy_ns <- st_sf(
  coef_bin = factor("Not significant", levels = all_levels),
  geometry = st_sfc(st_point(c(NA_real_, NA_real_)), crs = st_crs(world))
)

p_country <- ggplot() +
  # background: countries not in sample
  geom_sf(data = world_bg,
          aes(fill = coef_bin),
          color = "#AAAAAA",
          linewidth = 0.1) +
  # significant countries: solid fill
  geom_sf(data = world_sig,
          aes(fill = coef_bin),
          color = "#AAAAAA",
          linewidth = 0.15) +
  # non-significant countries: same fill
  geom_sf(data = world_ns,
          aes(fill = coef_bin),
          color = "#AAAAAA",
          linewidth = 0.15) +
  # dummy points to ensure all bin levels appear in legend
  geom_sf(data = dummy_regular,
          aes(fill = coef_bin),
          color = NA,
          na.rm = TRUE,
          key_glyph = "rect") +
  # dummy point for "Not significant" legend key with hatched glyph
  geom_sf(data = dummy_ns,
          aes(fill = coef_bin),
          key_glyph = draw_key_hatch)

# add hatching overlay (no legend — handled by dummy above)
if (!is.null(hatch_df)) {
  p_country <- p_country +
    geom_sf(data = hatch_df,
            color = "gray30",
            linewidth = 0.15,
            linetype = "dashed",
            show.legend = FALSE)
}

p_country <- p_country +
  scale_fill_manual(
    values = all_colors,
    name = "% change from control mean",
    drop = FALSE
  ) +
  coord_sf(crs = moll_crs,
           xlim = c(-14000000, 17000000),
           ylim = c(-6500000, 8500000),
           expand = FALSE) +
  theme_void(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.box.spacing = unit(0, "pt"),
    legend.key.width = unit(0.8, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 7),
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0, "pt"),
    legend.justification = "left",
    plot.margin = margin(2, 2, 2, 2, "pt")
  ) +
  guides(
    fill = guide_legend(ncol = 10, byrow = TRUE, title.position = "top")
  )

out_path2 <- paste0(paths$figures_dir, "fig.delta_map.country.jpg")
est_tc$coef_bin <- as.character(cut(est_tc$coef_pct, breaks = brks_c, labels = labs_c))
save_figure_data(est_tc, out_path2)
ggsave(out_path2, p_country, width = 10, height = 5, dpi = 300, bg = "white")
cat("  Saved to:", out_path2, "\n")
}

cat("\n=== DONE ===\n")
