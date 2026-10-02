###############################################
### pathreat.analysis.global.hetero.fig.R #########
### PA-level heterogeneity figures ############
### (1) delta vs biodiversity #################
### (2) delta vs land-use pressure (est/comp/access)
### (3) delta vs control-side threat level ####
###############################################

library(dplyr)
library(ggplot2)

source("code/pathreat.analysis.config.R")

#########################
### load PA-level data ##
#########################
{
# aggregate matched pixel data to PA level (treated pixels only, except tc_threat_composite)
d_matched <- fst::read_fst(paths$data_matched, columns = c(
  "treat", "wdpaid", "tc_threat_composite", "biodiversity2024",
  "lu.pressure.est", "lu.pressure.comp", "lu.pressure.access", "size"
))
cat(sprintf("Matched data: %s rows\n", format(nrow(d_matched), big.mark = ",")))

d_treated <- d_matched[d_matched$treat == 1, ]
rm(d_matched)

d_pa <- d_treated %>%
  group_by(wdpaid) %>%
  summarise(
    tc_threat_composite = mean(tc_threat_composite, na.rm = TRUE),
    biodiversity2024 = mean(biodiversity2024, na.rm = TRUE),
    lu.pressure.est = mean(lu.pressure.est, na.rm = TRUE),
    lu.pressure.comp = mean(lu.pressure.comp, na.rm = TRUE),
    lu.pressure.access = mean(lu.pressure.access, na.rm = TRUE),
    size = first(size),
    .groups = "drop"
  )
rm(d_treated)

# Coefficients come from the canonical PA regressions. Descriptive covariates
# retain their existing aggregation and the global normalization denominator.
pa_est <- readRDS("data/store/pathreat.analysis.byPA.est.Rds")
if (!all(c("wdpaid", "coef", "status") %in% names(pa_est)) ||
    anyDuplicated(pa_est$wdpaid) || !setequal(d_pa$wdpaid, pa_est$wdpaid)) {
  stop("Canonical PA regression estimates are missing or disagree with the PA set")
}
d_pa <- d_pa %>%
  left_join(select(pa_est, wdpaid, delta = coef), by = "wdpaid") %>%
  filter(is.finite(delta))

cat(sprintf("PAs: %d\n", nrow(d_pa)))
cat(sprintf("PAs with biodiversity data: %d of %d\n",
            sum(!is.na(d_pa$biodiversity2024)), nrow(d_pa)))

# global control mean and normalized effect for reference line
g_global <- readRDS(paths$est_global)
g_tc <- g_global[g_global$variable == "threat_composite", ]
ctrl_mean <- g_tc$control_mean
global_effect_norm <- g_tc$coef / ctrl_mean * 100
}

##############################################
### bin data for both figures ################
##############################################
{
# biodiversity bins
d_b <- d_pa[!is.na(d_pa$biodiversity2024), ]
d_b$delta_norm <- d_b$delta / ctrl_mean * 100
d_b <- d_b %>% mutate(bio_bin = ntile(biodiversity2024, 400))
d_b_binned <- d_b %>%
  group_by(bio_bin) %>%
  summarise(biodiversity2024 = mean(biodiversity2024), delta_norm = mean(delta_norm),
            size = mean(size), .groups = "drop")

# pressure bins (normalize to 0-1) — three variants
pressure_vars <- list(
  est    = list(var = "lu.pressure.est",    label = "Land-use pressure index [estimate] (normalized)"),
  comp   = list(var = "lu.pressure.comp",   label = "Land-use pressure index [composite] (normalized)"),
  access = list(var = "lu.pressure.access", label = "Land-use pressure index [access] (normalized)")
)

pressure_binned <- list()
for (pname in names(pressure_vars)) {
  pvar <- pressure_vars[[pname]]$var
  d_p <- d_pa[!is.na(d_pa[[pvar]]), ]
  d_p$delta_norm <- d_p$delta / ctrl_mean * 100
  # trim access tails at P5/P95
  if (pname == "access") {
    p05 <- quantile(d_p[[pvar]], 0.05, na.rm = TRUE)
    p95 <- quantile(d_p[[pvar]], 0.95, na.rm = TRUE)
    d_p <- d_p[d_p[[pvar]] >= p05 & d_p[[pvar]] <= p95, ]
  }
  pres_min <- min(d_p[[pvar]], na.rm = TRUE)
  pres_max <- max(d_p[[pvar]], na.rm = TRUE)
  d_p$pres_norm <- (d_p[[pvar]] - pres_min) / (pres_max - pres_min)
  d_p <- d_p %>% mutate(pres_bin = ntile(pres_norm, 400))
  pressure_binned[[pname]] <- d_p %>%
    group_by(pres_bin) %>%
    summarise(pres_norm = mean(pres_norm), delta_norm = mean(delta_norm),
              size = mean(size), .groups = "drop")
}

# tc_threat_composite bins
d_tc <- d_pa[!is.na(d_pa$tc_threat_composite), ]
d_tc$delta_norm <- d_tc$delta / ctrl_mean * 100
d_tc <- d_tc %>% mutate(tc_bin = ntile(tc_threat_composite, 400))
d_tc_binned <- d_tc %>%
  group_by(tc_bin) %>%
  summarise(tc_threat_composite = mean(tc_threat_composite), delta_norm = mean(delta_norm),
            size = mean(size), .groups = "drop")

# shared y-axis limits
all_deltas <- c(d_b_binned$delta_norm, d_tc_binned$delta_norm,
                unlist(lapply(pressure_binned, `[[`, "delta_norm")))
shared_ylim <- range(all_deltas) + c(-5, 5)
}

##############################################
### delta vs biodiversity ####################
##############################################
{
cat("\nCreating figure: delta vs biodiversity...\n")

p_bio <- ggplot(d_b_binned, aes(x = biodiversity2024, y = delta_norm)) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 11, face = "bold")
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = global_effect_norm, linetype = "dashed", color = "#B2182B", linewidth = 0.5) +
  geom_point(aes(size = size), alpha = 0.5, shape = 16, color = "grey30") +
  scale_size_continuous(name = expression("PA area (km"^2*")"), range = c(0.5, 4),
                        trans = "log10", breaks = c(10, 100, 1000, 10000)) +
  geom_smooth(method = "loess", span = 0.75, linewidth = 1,
              se = TRUE, alpha = 0.15, color = "#2166AC", fill = "#2166AC") +
  labs(
    x = "Species richness",
    y = "Effect on overall threat index (% of control mean)"
  ) +
  annotate("text", x = -Inf, y = global_effect_norm, label = "Global average",
           hjust = -0.05, vjust = -0.5, size = 3, color = "#B2182B", fontface = "italic") +
  annotate("text", x = max(d_b_binned$biodiversity2024, na.rm = TRUE) * 0.95,
           y = min(d_b_binned$delta_norm) * 0.8,
           label = "PA reduces threats", hjust = 1, vjust = 1,
           size = 3, color = "grey40", fontface = "italic") +
  annotate("text", x = max(d_b_binned$biodiversity2024, na.rm = TRUE) * 0.95,
           y = max(d_b_binned$delta_norm) * 0.8,
           label = "PA increases threats", hjust = 1, vjust = 0,
           size = 3, color = "grey40", fontface = "italic") +
  coord_cartesian(ylim = shared_ylim)

fig_file_bio <- paste0(paths$figures_dir, "fig.hetero.biodiversity.est.jpg")
save_figure_data(d_b_binned, fig_file_bio)
ggsave(fig_file_bio, plot = p_bio, width = 18, height = 14, units = "cm", dpi = 300)
cat("Saved:", fig_file_bio, "\n")
}

##############################################
### delta vs land-use pressure (3 variants) ##
##############################################

fig_files_pres <- list()
for (pname in names(pressure_vars)) {
  cat("\nCreating figure: delta vs", pressure_vars[[pname]]$var, "...\n")

  d_pb <- pressure_binned[[pname]]

  p_pres <- ggplot(d_pb, aes(x = pres_norm, y = delta_norm)) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 11, face = "bold")
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = global_effect_norm, linetype = "dashed", color = "#B2182B", linewidth = 0.5) +
    geom_point(aes(size = size), alpha = 0.5, shape = 16, color = "grey30") +
    scale_size_continuous(name = expression("PA area (km"^2*")"), range = c(0.5, 4),
                          trans = "log10", breaks = c(10, 100, 1000, 10000)) +
    geom_smooth(method = "loess", span = 0.75, linewidth = 1,
                se = TRUE, alpha = 0.15, color = "#2166AC", fill = "#2166AC") +
    labs(
      x = pressure_vars[[pname]]$label,
      y = "Effect on overall threat index (% of control mean)"
    ) +
    annotate("text", x = -Inf, y = global_effect_norm, label = "Global average",
             hjust = -0.05, vjust = -0.5, size = 3, color = "#B2182B", fontface = "italic") +
    annotate("text", x = max(d_pb$pres_norm, na.rm = TRUE) * 0.99,
             y = min(d_pb$delta_norm) * 0.8,
             label = "PA reduces threats", hjust = 1, vjust = 1,
             size = 3, color = "grey40", fontface = "italic") +
    annotate("text", x = max(d_pb$pres_norm, na.rm = TRUE) * 0.99,
             y = max(d_pb$delta_norm) * 0.8,
             label = "PA increases threats", hjust = 1, vjust = 0,
             size = 3, color = "grey40", fontface = "italic") +
    coord_cartesian(ylim = shared_ylim)

  suffix <- paste0("pressure.", pname)
  fig_file <- paste0(paths$figures_dir, "fig.hetero.", suffix, ".est.jpg")
  save_figure_data(d_pb, fig_file)
  ggsave(fig_file, plot = p_pres, width = 18, height = 14, units = "cm", dpi = 300)
  cat("Saved:", fig_file, "\n")
  fig_files_pres[[pname]] <- fig_file
}

##############################################
### delta vs control-side threat level #######
##############################################
{
cat("\nCreating figure: delta vs tc_threat_composite...\n")

p_tc <- ggplot(d_tc_binned, aes(x = tc_threat_composite, y = delta_norm)) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 11, face = "bold")
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = global_effect_norm, linetype = "dashed", color = "#B2182B", linewidth = 0.5) +
  geom_point(aes(size = size), alpha = 0.5, shape = 16, color = "grey30") +
  scale_size_continuous(name = expression("PA area (km"^2*")"), range = c(0.5, 4),
                        trans = "log10", breaks = c(10, 100, 1000, 10000)) +
  geom_smooth(method = "loess", span = 0.75, linewidth = 1,
              se = TRUE, alpha = 0.15, color = "#2166AC", fill = "#2166AC") +
  labs(
    x = "Control-side threat composite index",
    y = "Effect on overall threat index (% of control mean)"
  ) +
  annotate("text", x = -Inf, y = global_effect_norm, label = "Global average",
           hjust = -0.05, vjust = -0.5, size = 3, color = "#B2182B", fontface = "italic") +
  annotate("text", x = max(d_tc_binned$tc_threat_composite, na.rm = TRUE) * 0.95,
           y = min(d_tc_binned$delta_norm) * 0.8,
           label = "PA reduces threats", hjust = 1, vjust = 1,
           size = 3, color = "grey40", fontface = "italic") +
  annotate("text", x = max(d_tc_binned$tc_threat_composite, na.rm = TRUE) * 0.95,
           y = max(d_tc_binned$delta_norm) * 0.8,
           label = "PA increases threats", hjust = 1, vjust = 0,
           size = 3, color = "grey40", fontface = "italic") +
  coord_cartesian(ylim = shared_ylim)

fig_file_tc <- paste0(paths$figures_dir, "fig.hetero.tc_threat_composite.est.jpg")
save_figure_data(d_tc_binned, fig_file_tc)
ggsave(fig_file_tc, plot = p_tc, width = 18, height = 14, units = "cm", dpi = 300)
cat("Saved:", fig_file_tc, "\n")
}

####################################
### summary ########################
####################################
{
cat("\n=== Heterogeneity figures complete ===\n")
cat("Biodiversity:", fig_file_bio, "\n")
for (pname in names(fig_files_pres)) cat("Pressure", pname, ":", fig_files_pres[[pname]], "\n")
cat("TC control:  ", fig_file_tc, "\n")
}
