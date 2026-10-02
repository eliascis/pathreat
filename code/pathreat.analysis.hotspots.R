##########################################
### pathreat.analysis.hotspots.R #########
### PA effectiveness by hotspot ##########
### threat_composite: estimation + fig ###
##########################################

library(dplyr)
library(fixest)
library(ggplot2)

source("code/pathreat.analysis.config.R")


######################
### load data ########
######################
{
cat("Loading matched data...\n")
if (!"d.mbase" %in% ls()) {
  d.mbase <- load_matched_data()
}
cat(sprintf("  Loaded %s observations\n", format(nrow(d.mbase), big.mark = ",")))
}

##############################
### hotspot lookup ###########
##############################
{
# canonical complete-pair sample for threat_composite
pair_index <- build_matched_pair_index(d.mbase)
sample_mask <- matched_pair_sample_mask(
  d.mbase,
  "threat_composite",
  pair_index
)
d <- d.mbase[
  sample_mask,
  c("matched_pair_id", "treat", "country_rast", "biome_raster",
    "hotspot_id", "hotspot_name", "threat_composite"),
  drop = FALSE
]
pair_index <- build_matched_pair_index(d)

# hotspot membership is a protected-side target characteristic; propagate it
# to the matched control so hotspot subsamples retain whole pairs
target_hotspots <- d[pair_index$treated_rows,
                     c("matched_pair_id", "hotspot_id", "hotspot_name")]
target_order <- match(d$matched_pair_id, target_hotspots$matched_pair_id)
if (anyNA(target_order)) {
  stop("Could not propagate protected-side hotspot membership to every pair")
}
d$hotspot_id <- target_hotspots$hotspot_id[target_order]
d$hotspot_name <- target_hotspots$hotspot_name[target_order]

if (any(d$hotspot_id[pair_index$treated_rows] !=
        d$hotspot_id[pair_index$control_rows], na.rm = TRUE)) {
  stop("Protected and control hotspot assignments differ after propagation")
}

# extract unique protected-side hotspot_id / hotspot_name pairs
hs_lookup <- d %>%
  filter(!is.na(hotspot_id)) %>%
  distinct(hotspot_id, hotspot_name) %>%
  arrange(hotspot_id)

# add non-hotspot group
hs_lookup <- rbind(
  data.frame(hotspot_id = 0L, hotspot_name = "Non-Hotspot", stringsAsFactors = FALSE),
  hs_lookup
)

cat(sprintf("Hotspot groups: %d (%d hotspots + non-hotspot)\n",
            nrow(hs_lookup), nrow(hs_lookup) - 1))

# recode unassigned protected targets as the non-hotspot group
d$hotspot_id[is.na(d$hotspot_id)] <- 0L

# print group counts
cat("\nTreated pixels per hotspot:\n")
for (i in seq_len(nrow(hs_lookup))) {
  hid <- hs_lookup$hotspot_id[i]
  n_t <- sum(d$hotspot_id == hid & d$treat == 1)
  cat(sprintf("  [%2d] %-50s  n_treat = %s\n",
              hid, hs_lookup$hotspot_name[i], format(n_t, big.mark = ",")))
}
}

############################
### helper functions #######
############################

# format number with magnitude-adaptive decimal places
fmt_num <- function(x) {
  if (is.na(x)) return("")
  ax <- abs(x)
  if (ax == 0) return("0.000")
  if (ax >= 100)  return(formatC(x, format = "f", digits = 1))
  if (ax >= 1)    return(formatC(x, format = "f", digits = 2))
  return(formatC(x, format = "f", digits = 3))
}

# significance stars
sig_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.01) return("***")
  if (p < 0.05) return("**")
  if (p < 0.1)  return("*")
  return("")
}

############################
### estimation loop ########
############################
{
cat("\n=== Running by-hotspot estimation for threat_composite ===\n")

dep <- "threat_composite"
results_list <- list()

for (hid in hs_lookup$hotspot_id) {
  d_hs <- d[d$hotspot_id == hid, ]
  n_total <- nrow(d_hs)
  n_treat <- sum(d_hs$treat == 1, na.rm = TRUE)
  n_control <- sum(d_hs$treat == 0, na.rm = TRUE)
  n_countries <- length(unique(d_hs$country_rast))
  n_biomes <- length(unique(d_hs$biome_raster))
  hs_name <- hs_lookup$hotspot_name[hs_lookup$hotspot_id == hid]

  if (n_treat != n_control || n_total != 2L * n_treat) {
    stop(sprintf("Hotspot %s subsample splits matched pairs", hid))
  }

  cat(sprintf("  [%2d] %s: N=%s (treat=%s, ctrl=%s, countries=%d, biomes=%d)\n",
              hid, hs_name, format(n_total, big.mark = ","),
              format(n_treat, big.mark = ","), format(n_control, big.mark = ","),
              n_countries, n_biomes))

  # skip if too few obs
  if (n_total < 10 || n_treat < 5 || n_control < 5) {
    cat("    Skipped: insufficient N\n")
    results_list[[length(results_list) + 1]] <- data.frame(
      hotspot_id = hid,
      hotspot_name = hs_name,
      coef = NA,
      se = NA,
      pval = NA,
      ci_low = NA,
      ci_high = NA,
      n_obs = n_total,
      n_treat = n_treat,
      n_control = n_control,
      n_countries = n_countries,
      n_biomes = n_biomes,
      control_mean = NA,
      status = "skipped_insufficient_n",
      row.names = NULL
    )
    next
  }

  # control mean
  ctrl_mean <- mean(d_hs[d_hs$treat == 0, dep])

  tryCatch({
    # FE fallback: both country and biome dimensions
    if (n_countries > 1 & n_biomes > 1) {
      f <- as.formula(paste(dep, "~ treat | country_rast + biome_raster"))
      e <- feols(f, data = d_hs, cluster = "country_rast")
    } else if (n_countries > 1 & n_biomes == 1) {
      f <- as.formula(paste(dep, "~ treat | country_rast"))
      e <- feols(f, data = d_hs, cluster = "country_rast")
    } else if (n_countries == 1 & n_biomes > 1) {
      f <- as.formula(paste(dep, "~ treat | biome_raster"))
      e <- feols(f, data = d_hs, vcov = "hetero")
    } else {
      f <- as.formula(paste(dep, "~ treat"))
      e <- feols(f, data = d_hs, vcov = "hetero")
    }

    ct <- e$coeftable["treat", ]
      estimate <- regression_estimate(e)
    if (e$nobs != n_total) {
      stop(sprintf("Hotspot %s model changed the canonical sample", hid))
    }
    coef_val <- ct["Estimate"]
    se_val   <- ct["Std. Error"]
    pval_val <- ct["Pr(>|t|)"]

    results_list[[length(results_list) + 1]] <- data.frame(
      hotspot_id = hid,
      hotspot_name = hs_name,
      coef = coef_val,
      se = se_val,
      pval = pval_val,
      ci_low = estimate$ci_low,
      ci_high = estimate$ci_high,
      n_obs = e$nobs,
      n_treat = n_treat,
      n_control = n_control,
      n_countries = n_countries,
      n_biomes = n_biomes,
      control_mean = ctrl_mean,
      status = "estimated",
      row.names = NULL
    )

  }, error = function(err) {
    cat(sprintf("    ERROR: %s\n", err$message))
    results_list[[length(results_list) + 1]] <<- data.frame(
      hotspot_id = hid,
      hotspot_name = hs_name,
      coef = NA,
      se = NA,
      pval = NA,
      ci_low = NA,
      ci_high = NA,
      n_obs = n_total,
      n_treat = n_treat,
      n_control = n_control,
      n_countries = n_countries,
      n_biomes = n_biomes,
      control_mean = ctrl_mean,
      status = paste0("error: ", err$message),
      row.names = NULL
    )
  })
}

# combine
results_hotspots <- do.call("rbind", results_list)
rownames(results_hotspots) <- NULL
}

############################
### summary statistics #####
############################
{
cat("\n=== By-hotspot estimation summary ===\n")

# count by status
status_counts <- table(results_hotspots$status)
print(status_counts)

# significance summary
est <- results_hotspots[results_hotspots$status == "estimated", ]
n_sig <- sum(est$pval < 0.05, na.rm = TRUE)
n_neg <- sum(est$coef < 0 & est$pval < 0.05, na.rm = TRUE)
n_pos <- sum(est$coef > 0 & est$pval < 0.05, na.rm = TRUE)
cat(sprintf("\nthreat_composite: %d estimated, %d sig at 5%% (neg=%d, pos=%d)\n",
            nrow(est), n_sig, n_neg, n_pos))
}

############################
### save Rds ###############
############################
{
output_rds <- "results/pathreat.analysis.hotspots.Rds"
saveRDS(results_hotspots, output_rds)
cat(sprintf("\nSaved Rds: %s (%d rows)\n", output_rds, nrow(results_hotspots)))
}

############################
### forest plot ############
############################
{
cat("\n=== Generating forest plot ===\n")

# filter to estimated rows
d.at <- results_hotspots[results_hotspots$status == "estimated", ]

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

# order by norm_coef ascending, assign sequential ranks
d.at <- d.at[order(d.at$norm_coef), ]
d.at$rank <- seq_len(nrow(d.at))

# labels: use hotspot_name
d.at$label <- d.at$hotspot_name

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
  label = "Global",
  is_global = TRUE,
  stringsAsFactors = FALSE
)

# combine
keep_cols <- c("norm_coef", "norm_ci_lower", "norm_ci_upper", "pval",
               "sig_tier", "rank", "n_treat", "label")
d_sub <- d.at[, c(keep_cols[keep_cols != "n_treat"], "n_treat")]
d_sub$is_global <- (d_sub$label == "Non-hotspot")
d_plot <- rbind(d_sub[, c(keep_cols, "is_global")],
                global_row[, c(keep_cols, "is_global")])

# shared y-axis limits
y_lim <- range(d_plot$rank) + c(-0.8, 0.8)

# figure height scaled to number of rows (~0.7 cm per row + margins)
fig_height <- max(22, nrow(d_plot) * 0.7 + 4)

p1 <- ggplot(d_plot, aes(x = norm_coef, y = rank)) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
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

fig_path <- paste0(paths$figures_dir, "fig.hotspots.threat_composite.coef.est.jpg")
ggsave(plot = p1, fig_path, units = "cm", width = 22, height = fig_height, dpi = 300)
cat(sprintf("Saved forest plot: %s (%.0f cm tall)\n", fig_path, fig_height))
}

############################
### cleanup ################
############################
{
rm(d, results_list, d.at, d_sub, d_plot, global_row, est_global, eg, p1)
gc()
}

cat("\n=== By-hotspot estimation complete ===\n")
cat("Output Rds:  results/pathreat.analysis.hotspots.Rds\n")
cat("Output fig:  pub/figures/fig.hotspots.threat_composite.coef.est.jpg\n")
