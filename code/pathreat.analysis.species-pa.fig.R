##################################################
### pathreat.analysis.species-pa.fig.R ############
### Terrestrial area composition bar chart #######
### Per-taxon coverage + species-range effects ###
##################################################

library(dplyr)
library(ggplot2)
library(fst)
library(tidyr)
library(patchwork)

source("code/pathreat.analysis.config.R")


############################
### load data ##############
############################
{
cat("=== Loading unmatched data ===\n")
d_um <- read_fst(paths$data_unmatched,
                 columns = c("pixel_id", "ever_pa", "pa_year_designated"))
cat("  Total pixels:", format(nrow(d_um), big.mark = ","), "\n")
}


############################
### compute shares #########
############################
{
n_total <- nrow(d_um)

# Four categories based on ever_pa + pa_year_designated
# Unknown year (pa_year_designated == 0) lumped with pre-2001
n_pa_pre2001 <- sum(d_um$ever_pa == 1 &
                    (d_um$pa_year_designated < 2001 | d_um$pa_year_designated == 0),
                    na.rm = TRUE)
n_pa_study   <- sum(d_um$ever_pa == 1 &
                    d_um$pa_year_designated >= 2001 &
                    d_um$pa_year_designated <= 2020,
                    na.rm = TRUE)
n_pa_post2020 <- sum(d_um$ever_pa == 1 &
                     d_um$pa_year_designated > 2020,
                     na.rm = TRUE)
n_unprotected <- sum(d_um$ever_pa == 0, na.rm = TRUE)

cat("\n=== Area composition ===\n")
cat("  PAs pre-2001:     ", format(n_pa_pre2001, big.mark = ","),
    sprintf(" (%.1f%%)\n", 100 * n_pa_pre2001 / n_total))
cat("  PAs 2001-2020:    ", format(n_pa_study, big.mark = ","),
    sprintf(" (%.1f%%)\n", 100 * n_pa_study / n_total))
cat("  PAs post-2020:    ", format(n_pa_post2020, big.mark = ","),
    sprintf(" (%.1f%%)\n", 100 * n_pa_post2020 / n_total))
cat("  Unprotected:      ", format(n_unprotected, big.mark = ","),
    sprintf(" (%.1f%%)\n", 100 * n_unprotected / n_total))
cat("  Total:            ", format(n_total, big.mark = ","), "\n")

rm(d_um)
gc()
}


############################
### figure #################
############################
{
# build data frame for plotting
lbl_old   <- "Pre-2001 PAs"
lbl_study <- "2001\u20132020 PAs"
lbl_post  <- "Post-2020 PAs"
lbl_unpro <- "Unprotected"

df_bar <- data.frame(
  category = factor(
    c(lbl_old, lbl_study, lbl_post, lbl_unpro),
    levels = c(lbl_unpro, lbl_post, lbl_study, lbl_old)
  ),
  share = c(n_pa_pre2001, n_pa_study, n_pa_post2020, n_unprotected) / n_total * 100
)

# colors
cols <- setNames(c("#2c7bb6", "#abd9e9", "#e0f3f8", "#bdbdbd"),
                 c(lbl_old, lbl_study, lbl_post, lbl_unpro))

# total PA share
pa_share <- 100 * (n_pa_pre2001 + n_pa_study + n_pa_post2020) / n_total

p <- ggplot(df_bar, aes(x = 1, y = share, fill = category)) +
  geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.3) +
  # arrow marking total PA share (above the bar)
  annotate("segment", y = pa_share, yend = pa_share,
           x = 1.38, xend = 1.32,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
           linewidth = 0.4) +
  annotate("text", y = pa_share, x = 1.46,
           label = sprintf("%.1f%% protected", pa_share),
           size = 3.2, hjust = 0) +
  # white line at 30% and arrow label for 30x30 goal
  annotate("segment", y = 30, yend = 30,
           x = 0.7, xend = 1.3,
           color = "white", linewidth = 0.8) +
  annotate("segment", y = 30, yend = 30,
           x = 1.38, xend = 1.32,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
           linewidth = 0.4) +
  annotate("text", y = 30, x = 1.46,
           label = "30x30 target",
           size = 3.2, hjust = 0) +
  coord_flip() +
  scale_fill_manual(values = cols, name = NULL,
                    breaks = c(lbl_old, lbl_study, lbl_post, lbl_unpro)) +
  scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  scale_x_continuous(limits = c(0.5, 1.5)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.key.size  = unit(0.4, "cm"),
    plot.margin      = margin(10, 15, 10, 15)
  )

outfile <- paste0(paths$figures_dir, "fig.PA-species-cover-global.jpg")
ggsave(outfile, p, width = 8, height = 2, dpi = 300)
cat("\nSaved:", outfile, "\n")
}


####################################
### species-level coverage shares ##
####################################
{
cat("\n=== Species-level PA coverage shares ===\n")

# load disaggregated species x PA overlap from data script
sp_pa <- readRDS("data/store/pathreat.data.species-pa-cover.Rds")
if ("area" %in% names(sp_pa) && !"overlap_area" %in% names(sp_pa)) {
  names(sp_pa)[names(sp_pa) == "area"] <- "overlap_area"
}
cat("  Loaded species-PA overlap:", format(nrow(sp_pa), big.mark = ","), "rows,",
    length(unique(sp_pa$species_id)), "species\n")

# look up pa_year_designated per wdpaid
wdpa_yr <- read_fst("data/store/pathreat.data.wdpa.fst",
                     columns = c("wdpaid", "pa_year_designated"))
wdpa_yr <- wdpa_yr[!is.na(wdpa_yr$wdpaid), ]
wdpa_yr <- wdpa_yr[!duplicated(wdpa_yr$wdpaid), ]

# join designation year; unprotected (wdpaid == 0) gets NA
sp_pa <- merge(sp_pa, wdpa_yr,
               by = "wdpaid",
               all.x = TRUE)
rm(wdpa_yr)

# classify into period categories (unknown year lumped with pre-2001)
sp_pa <- sp_pa %>%
  mutate(
    period = case_when(
      wdpaid == 0 ~ "unprotected",
      pa_year_designated >= 2001 & pa_year_designated <= 2020 ~ "pa_2001_2020",
      pa_year_designated > 2020 ~ "pa_post2020",
      TRUE ~ "pa_pre2001"
    )
  )

# aggregate area by species x period, compute shares
# include taxon and category if present
group_vars <- c("species_id", "species_name", "range_total")
has_taxon <- "taxon" %in% names(sp_pa)
has_category <- "category" %in% names(sp_pa)
if (has_taxon) group_vars <- c(group_vars, "taxon")
if (has_category) group_vars <- c(group_vars, "category")
group_vars <- c(group_vars, "period")

species_cover <- sp_pa %>%
  group_by(across(all_of(group_vars))) %>%
  dplyr::summarize(
    n = sum(overlap_area),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = period,
    values_from = n,
    values_fill = 0,
    names_prefix = "n_"
  ) %>%
  mutate(
    n_range = range_total,
    share_pa_pre2001   = n_pa_pre2001   / n_range,
    share_pa_2001_2020 = n_pa_2001_2020 / n_range,
    share_pa_post2020  = n_pa_post2020  / n_range,
    share_unprotected  = n_unprotected  / n_range
  )

# ensure all share columns exist (in case a category is absent)
for (col in c("share_pa_pre2001", "share_pa_2001_2020",
              "share_pa_post2020", "share_unprotected")) {
  if (!col %in% names(species_cover))
    species_cover[[col]] <- 0
}

# ensure n_ columns exist too
for (col in c("n_pa_pre2001", "n_pa_2001_2020", "n_pa_post2020", "n_unprotected")) {
  if (!col %in% names(species_cover))
    species_cover[[col]] <- 0
}

cat("\n=== Species coverage summary ===\n")
cat("  Species:", nrow(species_cover), "\n")
if (has_taxon) {
  cat("  By taxon:\n")
  for (tx in sort(unique(species_cover$taxon))) {
    cat(sprintf("    %-12s %d\n", tx, sum(species_cover$taxon == tx)))
  }
}
cat("  Mean shares (across species):\n")
cat(sprintf("    PA pre-2001:   %.1f%%\n", 100 * mean(species_cover$share_pa_pre2001)))
cat(sprintf("    PA 2001-2020:  %.1f%%\n", 100 * mean(species_cover$share_pa_2001_2020)))
cat(sprintf("    PA post-2020:  %.1f%%\n", 100 * mean(species_cover$share_pa_post2020)))
cat(sprintf("    Unprotected:   %.1f%%\n", 100 * mean(species_cover$share_unprotected)))

rm(sp_pa)
gc()
}


##########################################
### per-taxon coverage bars ##############
##########################################
{
cat("\n=== Per-taxon PA vintage coverage bars ===\n")

taxon_labels_bar <- c(
  mammal    = "Mammals",
  amphibian = "Amphibians",
  reptile   = "Reptiles",
  bird      = "Birds"
)

# taxon order for the combined figure
taxon_order <- c("mammal", "bird", "amphibian", "reptile")

# store per-taxon bar plots for reuse in combined figures
p_taxon_bars <- list()

# common theme for inner panels (no legend, no x-axis labels)
theme_inner <- theme_minimal(base_size = 10) +
  theme(
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.text.x      = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "none",
    plot.margin      = margin(0, 10, 0, 5)
  )

# helper to build a single bar panel
make_coverage_bar <- function(df_bar_i, pa_share_i, subtitle_i = NULL,
                              taxon_label = NULL, inner = TRUE) {
  # midpoint of unprotected segment for label placement
  unpro_share <- df_bar_i$share[df_bar_i$category == lbl_unpro]
  label_y <- 100 - unpro_share / 2

  p_i <- ggplot(df_bar_i, aes(x = 1, y = share, fill = category)) +
    geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.3) +
    annotate("segment",
             y = pa_share_i, yend = pa_share_i,
             x = 1.38, xend = 1.32,
             arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
             linewidth = 0.4) +
    annotate("text",
             y = pa_share_i, x = 1.46,
             label = sprintf("%.1f%% protected", pa_share_i),
             size = 2.8, hjust = 0) +
    coord_flip() +
    scale_fill_manual(values = cols, name = NULL,
                      breaks = c(lbl_old, lbl_study, lbl_post, lbl_unpro)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100),
                       labels = function(x) paste0(x, "%")) +
    scale_x_continuous(limits = c(0.5, 1.5)) +
    labs(x = NULL, y = NULL, subtitle = subtitle_i)

  # taxon label centred in the grey (unprotected) area

  if (!is.null(taxon_label)) {
    p_i <- p_i +
      annotate("text",
               y = label_y, x = 1,
               label = taxon_label,
               size = 3, fontface = "italic",
               color = "grey30")
  }

  if (inner) {
    p_i <- p_i + theme_inner
  } else {
    # bottom panel: show legend, show x-axis labels
    p_i <- p_i +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        panel.grid       = element_blank(),
        legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm"),
        plot.margin      = margin(0, 10, 5, 5)
      )
  }
  p_i
}

# global bar (top panel, no legend) — with section heading
p_global_inner <- make_coverage_bar(
  df_bar, pa_share,
  subtitle_i = "All terrestrial land",
  inner = TRUE
) +
  # 30x30 target only on the global bar (below the bar)
  annotate("segment",
           y = 30, yend = 30,
           x = 0.7, xend = 1.3,
           color = "white", linewidth = 0.8) +
  annotate("segment",
           y = 30, yend = 30,
           x = 0.62, xend = 0.68,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
           linewidth = 0.4) +
  annotate("text",
           y = 30, x = 0.54,
           label = "30x30 target",
           size = 2.8, hjust = 0)

# per-taxon panels
panel_list <- list(p_global_inner)

for (i in seq_along(taxon_order)) {
  tx <- taxon_order[i]
  sc_tx <- species_cover %>% filter(taxon == tx)
  tx_label <- taxon_labels_bar[tx]

  n_pre   <- sum(sc_tx$n_pa_pre2001)
  n_study <- sum(sc_tx$n_pa_2001_2020)
  n_post  <- sum(sc_tx$n_pa_post2020)
  n_unpro <- sum(sc_tx$n_unprotected)
  n_tot   <- n_pre + n_study + n_post + n_unpro

  pa_share_tx <- 100 * (n_pre + n_study + n_post) / n_tot

  df_bar_tx <- data.frame(
    category = factor(
      c(lbl_old, lbl_study, lbl_post, lbl_unpro),
      levels = c(lbl_unpro, lbl_post, lbl_study, lbl_old)
    ),
    share = c(n_pre, n_study, n_post, n_unpro) / n_tot * 100
  )

  cat(sprintf("  %s: %.1f%% protected (n_species = %d)\n",
              tx_label, pa_share_tx, nrow(sc_tx)))

  is_last <- (i == length(taxon_order))
  is_first <- (i == 1)

  p_tx <- make_coverage_bar(
    df_bar_tx, pa_share_tx,
    subtitle_i = if (is_first) "Threatened species ranges" else NULL,
    taxon_label = tx_label,
    inner = !is_last
  )

  panel_list[[length(panel_list) + 1]] <- p_tx

  # also store compact version for per-taxon combined figures
  p_taxon_bars[[tx]] <- make_coverage_bar(
    df_bar_tx, pa_share_tx,
    taxon_label = tx_label,
    inner = TRUE
  )
}

# stack all 5 panels — negative margins to reduce vertical gaps
p_combined_cover <- Reduce("/", panel_list) +
  plot_layout(heights = rep(1, length(panel_list)))

outfile <- paste0(paths$figures_dir, "fig.PA-species-cover-all.jpg")
ggsave(outfile, p_combined_cover, width = 4, height = 4.5, dpi = 300)
cat("  Saved:", outfile, "\n")

# save plot object for reuse in combined figure (bytaxa.fig.R)
saveRDS(p_combined_cover, file.path(paths$results_dir, "p_combined_cover.Rds"))
cat("  Saved plot object: p_combined_cover.Rds\n")
}

make_species_cover_plot <- function(sc_data, taxon_label = NULL) {
  # sc_data: species_cover data for one taxon, with iucn/category column

  # IUCN category factor
  sc_data <- sc_data %>%
    mutate(
      total_pa = share_pa_pre2001 + share_pa_2001_2020 + share_pa_post2020,
      iucn = factor(category, levels = c("CR", "EN", "VU"),
                    labels = c("Critically Endangered (CR)",
                               "Endangered (EN)",
                               "Vulnerable (VU)"))
    ) %>%
    filter(!is.na(iucn)) %>%
    group_by(iucn) %>%
    arrange(desc(total_pa), .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    ungroup()

  # reshape to long
  df_long <- sc_data %>%
    select(rank,
           iucn,
           share_pa_pre2001,
           share_pa_2001_2020,
           share_pa_post2020,
           share_unprotected) %>%
    pivot_longer(cols = starts_with("share_"),
                 names_to = "cat", values_to = "share") %>%
    mutate(
      cat = factor(cat,
        levels = c("share_unprotected", "share_pa_post2020",
                   "share_pa_2001_2020", "share_pa_pre2001"),
        labels = c(lbl_unpro, lbl_post, lbl_study, lbl_old)),
      pct = share * 100
    )

  # label data for inside-block annotations
  df_labels <- sc_data %>%
    group_by(iucn) %>%
    dplyr::summarize(
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      rank = 1,
      pct = 95,
      lab = as.character(iucn)
    )

  # rank annotations: first and last per facet
  df_rank_labels <- sc_data %>%
    group_by(iucn) %>%
    dplyr::summarize(
      last = max(rank),
      .groups = "drop"
    ) %>%
    tidyr::crossing(pos = c("first", "last")) %>%
    mutate(
      rank = ifelse(pos == "first", 1, last),
      lab = as.character(rank),
      pct = -1.5,
      y_fct = reorder(rank, -rank)
    )

  p <- ggplot(df_long, aes(x = pct, y = reorder(rank, -rank), fill = cat)) +
    geom_bar(stat = "identity", width = 1, color = NA) +
    geom_text(data = df_labels,
              aes(x = pct, y = n * 0.5, label = lab),
              inherit.aes = FALSE,
              size = 2.8, color = "grey30", fontface = "bold",
              lineheight = 0.85, hjust = 1) +
    geom_text(data = df_rank_labels,
              aes(x = pct, y = y_fct, label = lab),
              inherit.aes = FALSE,
              size = 1.8, color = "grey50", hjust = 1) +
    facet_grid(rows = vars(iucn), scales = "free_y", space = "free_y") +
    scale_fill_manual(values = cols, name = NULL,
                      breaks = c(lbl_old, lbl_study, lbl_post, lbl_unpro)) +
    scale_x_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      panel.grid       = element_blank(),
      panel.spacing    = unit(8, "pt"),
      strip.text       = element_blank(),
      legend.position  = "bottom",
      legend.key.size  = unit(0.4, "cm"),
      legend.text      = element_text(size = 7),
      plot.margin      = margin(0, 10, 5, 15)
    ) +
    guides(fill = guide_legend(nrow = 1))

  return(p)
}


##########################################
### helper: effectiveness plot ###########
##########################################

# effectiveness color palette
eff_cols <- c(
  "High"             = "#1a9850",
  "Medium"           = "#91cf60",
  "Low"              = "#d9ef8b",
  "Not significant"  = "#d9d9d9",
  "Harmful"          = "#d73027",
  "No inference available" = "#969696",
  "Unassessed"       = "#f0f0f0"
)

make_species_eff_plot <- function(df_fig_data) {
  # df_fig_data: canonical species-range estimates with pct_effect_s and category
  #   Only species with estimated coefficients (assessed 2001-2020 PAs)

  df_fig_data <- df_fig_data %>%
    filter(
      status == "estimated",
      !is.na(pct_effect_s),
      is.finite(pct_effect_s)
    ) %>%
    mutate(
      iucn = factor(category, levels = c("CR", "EN", "VU"),
                    labels = c("Critically Endangered (CR)",
                               "Endangered (EN)",
                               "Vulnerable (VU)"))
    ) %>%
    filter(!is.na(iucn)) %>%
    group_by(iucn) %>%
    arrange(pct_effect_s, .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    ungroup()

  # symmetric x limits; winsorize extremes to boundaries
  x_abs <- max(abs(quantile(df_fig_data$pct_effect_s,
                             c(0.01, 0.99)))) * 1.05
  x_abs <- max(x_abs, 10)
  df_fig_data$pct_effect_s <- pmax(pmin(df_fig_data$pct_effect_s,
                                        x_abs), -x_abs)

  df_labels <- df_fig_data %>%
    group_by(iucn) %>%
    dplyr::summarize(
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      rank = 1,
      x = x_abs * 0.95,
      lab = as.character(iucn)
    )

  df_rank_labels <- df_fig_data %>%
    group_by(iucn) %>%
    dplyr::summarize(
      last = max(rank),
      .groups = "drop"
    ) %>%
    tidyr::crossing(pos = c("first", "last")) %>%
    mutate(
      rank = ifelse(pos == "first", 1, last),
      lab = as.character(rank),
      x = -x_abs - x_abs * 0.03,
      y_fct = reorder(rank, -rank)
    )

  p <- ggplot(df_fig_data,
              aes(x = pct_effect_s,
                  y = reorder(rank, -rank),
                  fill = pct_effect_s)) +
    geom_bar(stat = "identity", width = 1, color = NA) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
    geom_text(data = df_labels,
              aes(x = x, y = n * 0.5, label = lab),
              inherit.aes = FALSE,
              size = 2.8, color = "grey30", fontface = "bold",
              lineheight = 0.85, hjust = 1) +
    geom_text(data = df_rank_labels,
              aes(x = x, y = y_fct, label = lab),
              inherit.aes = FALSE,
              size = 1.8, color = "grey50", hjust = 1) +
    facet_grid(rows = vars(iucn), scales = "free_y", space = "free_y") +
    scale_fill_gradient2(
      low = "#2166ac", mid = "grey90", high = "#b2182b",
      midpoint = 0,
      name = "% change",
      limits = c(-x_abs, x_abs),
      oob = scales::squish
    ) +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(-x_abs, x_abs),
      labels = function(x) paste0(x, "%")
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y       = element_blank(),
      axis.ticks.y      = element_blank(),
      panel.grid        = element_blank(),
      panel.spacing     = unit(8, "pt"),
      strip.text        = element_blank(),
      legend.position   = "bottom",
      legend.key.width  = unit(1.5, "cm"),
      legend.key.height = unit(0.3, "cm"),
      legend.text       = element_text(size = 7),
      plot.margin       = margin(5, 10, 5, 15)
    )

  return(p)
}


##########################################
### per-taxon coverage figures ###########
##########################################
{
cat("\n=== Generating per-taxon coverage figures ===\n")

# Determine taxon list
if (has_taxon) {
  taxon_list <- sort(unique(species_cover$taxon))
} else {
  species_cover$taxon <- "mammal"
  taxon_list <- "mammal"
}

# Taxon display names (for subtitles)
taxon_labels <- c(
  mammal    = "Mammals",
  amphibian = "Amphibians",
  reptile   = "Reptiles",
  bird      = "Birds"
)

for (tx in taxon_list) {
  tx_label <- ifelse(tx %in% names(taxon_labels), taxon_labels[tx], tx)
  sc_tx <- species_cover %>% filter(taxon == tx)
  cat(sprintf("  %s: %d species\n", tx_label, nrow(sc_tx)))

  if (nrow(sc_tx) == 0) next

  p_species <- make_species_cover_plot(sc_tx)

  outfile <- sprintf("%sfig.PA-species-cover-%ss.jpg",
                     paths$figures_dir, tx)
  ggsave(outfile, p_species, width = 4, height = 6, dpi = 300)
  cat("  Saved:", outfile, "\n")
}
}


##########################################
### global PA effectiveness bar ########
##########################################
{
cat("\n=== Global PA effectiveness bar chart ===\n")

# pixel counts per wdpaid (all PA pixels in the unmatched sample)
d_wdpa <- read_fst(paths$data_unmatched, columns = c("wdpaid"))
d_wdpa <- d_wdpa[!is.na(d_wdpa$wdpaid), , drop = FALSE]
pa_pixels <- as.data.frame(table(d_wdpa$wdpaid), stringsAsFactors = FALSE)
names(pa_pixels) <- c("wdpaid", "n_pixels")
pa_pixels$wdpaid <- as.numeric(pa_pixels$wdpaid)
rm(d_wdpa)
gc()

cat("  PA pixels loaded, unique PAs:", nrow(pa_pixels), "\n")

# PA-level effectiveness from matched sample
pa_eff <- readRDS("data/store/pathreat.analysis.byPA.est.Rds")
if (!all(c("coef", "status", "p_value", "effectiveness") %in% names(pa_eff)) ||
    anyDuplicated(pa_eff$wdpaid)) stop("Canonical PA regression fields are missing")
valid_inference <- pa_eff$status == "estimated" & is.finite(pa_eff$p_value)
if (any((pa_eff$effectiveness == "no_inference") != !valid_inference)) {
  stop("PA effectiveness categories disagree with regression inference")
}

# join effectiveness; unmatched PAs → "Unassessed"
pa_pixels <- merge(pa_pixels, pa_eff[, c("wdpaid", "effectiveness")],
                   by = "wdpaid", all.x = TRUE)
pa_pixels$effectiveness[is.na(pa_pixels$effectiveness)] <- "unassessed"

# clean tier labels
pa_pixels$tier <- factor(pa_pixels$effectiveness,
  levels = c("high", "medium", "low", "not_significant", "harmful", "no_inference", "unassessed"),
  labels = c("High", "Medium", "Low", "Not significant", "Harmful", "No inference available", "Unassessed")
)

# aggregate
tier_shares <- pa_pixels %>%
  group_by(tier) %>%
  dplyr::summarize(
    n = sum(n_pixels),
    .groups = "drop"
  ) %>%
  mutate(pct = 100 * n / sum(n))

cat("  Tier shares (% of PA area):\n")
for (i in seq_len(nrow(tier_shares))) {
  cat(sprintf("    %-18s %6.1f%%\n", tier_shares$tier[i], tier_shares$pct[i]))
}

# effective share = high + medium + low
eff_share <- sum(tier_shares$pct[tier_shares$tier %in% c("High", "Medium", "Low")])
cat(sprintf("  Effective share: %.1f%%\n", eff_share))

# standalone global bar (same style as fig.PA-species-cover-global.jpg)
p_global_eff <- ggplot(tier_shares, aes(x = 1, y = pct, fill = tier)) +
  geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.3,
           position = position_stack(reverse = TRUE)) +
  annotate("segment", y = eff_share, yend = eff_share,
           x = 1.38, xend = 1.32,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
           linewidth = 0.4) +
  annotate("text", y = eff_share, x = 1.46,
           label = sprintf("%.1f%% effective", eff_share),
           size = 3.2, hjust = 0) +
  coord_flip() +
  scale_fill_manual(
    values = eff_cols,
    name = NULL,
    breaks = c("High", "Medium", "Low", "Not significant", "Harmful", "No inference available", "Unassessed")
  ) +
  scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  scale_x_continuous(limits = c(0.5, 1.5)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.key.size  = unit(0.4, "cm"),
    plot.margin      = margin(10, 15, 10, 15)
  ) +
  guides(fill = guide_legend(nrow = 2))

outfile_global_eff <- paste0(paths$figures_dir, "fig.PA-species-effectiveness-global.jpg")
ggsave(outfile_global_eff, p_global_eff, width = 8, height = 2, dpi = 300)
cat("Saved:", outfile_global_eff, "\n")

# global panel (top) — no legend, no x-axis labels
p_global_eff_top <- ggplot(tier_shares, aes(x = 1, y = pct, fill = tier)) +
  geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.3,
           position = position_stack(reverse = TRUE)) +
  annotate("segment", y = eff_share, yend = eff_share,
           x = 0.62, xend = 0.68,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
           linewidth = 0.4) +
  annotate("text", y = eff_share + 1, x = 0.54,
           label = sprintf("%.1f%% effective", eff_share),
           size = 2.5, hjust = 0) +
  coord_flip() +
  scale_fill_manual(
    values = eff_cols,
    name = NULL,
    breaks = c("High", "Medium", "Low", "Not significant", "Harmful", "No inference available", "Unassessed")
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100),
                     labels = function(x) paste0(x, "%")) +
  scale_x_continuous(limits = c(0.5, 1.5)) +
  labs(x = NULL, y = NULL, subtitle = "All PA area") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.text.x      = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "none",
    plot.margin      = margin(5, 10, 0, 5)
  )

rm(pa_pixels, pa_eff, tier_shares)
gc()
}


##########################################
### per-taxon effectiveness figures ######
##########################################
{
cat("\n=== Per-taxon canonical species-range effect figures ===\n")

# load by-species coefficient estimates for Panel B
sp_est <- readRDS(file.path(paths$results_dir, "pathreat.byspecies.est.Rds"))

for (tx in taxon_list) {
  tx_label <- ifelse(tx %in% names(taxon_labels), taxon_labels[tx], tx)
  df_fig <- sp_est %>%
    filter(
      taxon == tx,
      category %in% c("CR", "EN", "VU")
    )

  cat(sprintf("  %s: %d species with estimated coefficient\n",
              tx_label,
              sum(df_fig$status == "estimated" &
                    !is.na(df_fig$pct_effect_s) &
                    is.finite(df_fig$pct_effect_s))))

  if (nrow(df_fig) == 0) next

  p_eff <- make_species_eff_plot(df_fig)

  outfile <- sprintf("%sfig.PA-species-effectiveness-%ss.jpg",
                     paths$figures_dir, tx)
  ggsave(outfile, p_eff, width = 4, height = 6, dpi = 300)
  cat("  Saved:", outfile, "\n")
}
}
