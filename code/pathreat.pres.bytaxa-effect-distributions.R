##########################################################
### pathreat.pres.bytaxa-effect-distributions.R ##########
### Taxon effect distributions ###########################
##########################################################
#
# Purpose: Plot fixed-scale taxon density panels for the canonical
#          species-range ATTs. Density estimation uses every eligible species;
#          only the displayed x-window is restricted.
#
# Inputs:  results/pathreat.byspecies.est.Rds
#          results/pathreat.bytaxa.est.Rds (centralized count validation)
#
# Output:  pub/figures/
#          fig.bytaxa.effect-distributions.jpg
#          (consumed by pathreat_pres_03.tex)
#          fig.bytaxa.effect-distributions.manuscript.jpg
#          (same estimates and density data, sized for pathreat_10.tex)
#

library(ggplot2)

source("code/pathreat.analysis.config.R")

cat("=== Taxon Effect Distributions ===\n\n")

species_file <- file.path(paths$results_dir, "pathreat.byspecies.est.Rds")
taxon_file <- file.path(paths$results_dir, "pathreat.bytaxa.est.Rds")
output_file <- paste0(
  paths$figures_dir,
  "fig.bytaxa.effect-distributions.jpg"
)

required_files <- c(species_file, taxon_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required inputs: ", paste(missing_files, collapse = ", "))
}

if (!dir.exists(paths$figures_dir)) {
  dir.create(paths$figures_dir, recursive = TRUE)
}

taxon_labels <- c(
  mammal = "Mammals",
  bird = "Birds",
  amphibian = "Amphibians",
  reptile = "Reptiles"
)

taxon_colors <- c(
  Mammals = "#009E73",
  Birds = "#56B4E9",
  Amphibians = "#E69F00",
  Reptiles = "#CC79A7"
)

############################################
### STEP 1: Validate the eligible sample ###
############################################
{
cat("STEP 1: Loading canonical species estimates...\n")

species_estimates <- readRDS(species_file)
required_columns <- c(
  "species_id",
  "taxon",
  "delta_s",
  "se_delta_s",
  "p_delta_s",
  "pct_effect_s",
  "control_mean_s",
  "status"
)

missing_columns <- setdiff(required_columns, names(species_estimates))
if (length(missing_columns) > 0) {
  stop(
    "Canonical species estimates are missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}
if (anyDuplicated(species_estimates$species_id)) {
  stop("Canonical species estimates contain duplicated species_id values.")
}

eligible <- with(
  species_estimates,
  status == "estimated" &
    is.finite(delta_s) &
    is.finite(se_delta_s) &
    is.finite(p_delta_s) &
    is.finite(pct_effect_s) &
    is.finite(control_mean_s)
)

plot_data <- species_estimates[
  eligible,
  c("species_id", "taxon", "delta_s"),
  drop = FALSE
]

if (nrow(plot_data) == 0) {
  stop("No species satisfy the centralized taxon-summary eligibility rule.")
}
if (any(!is.finite(plot_data$delta_s))) {
  stop("Eligible species include non-finite delta_s values.")
}

unexpected_taxa <- setdiff(unique(plot_data$taxon), names(taxon_labels))
if (length(unexpected_taxa) > 0) {
  stop(
    "Eligible species include unexpected taxa: ",
    paste(unexpected_taxa, collapse = ", ")
  )
}

observed_counts <- table(
  factor(
    plot_data$taxon,
    levels = names(taxon_labels)
  )
)
if (any(observed_counts == 0)) {
  stop(
    "No eligible species estimates for: ",
    paste(names(observed_counts)[observed_counts == 0], collapse = ", ")
  )
}
if (any(observed_counts < 2)) {
  stop("At least two eligible species are required in every taxon.")
}

taxon_bundle <- readRDS(taxon_file)
if (!all(c("settings", "counts") %in% names(taxon_bundle))) {
  stop("Centralized taxon bundle is missing settings or counts.")
}

central_counts <- taxon_bundle$counts
required_count_columns <- c(
  "taxon",
  "n_species_estimated",
  "n_species_equal"
)
missing_count_columns <- setdiff(
  required_count_columns,
  names(central_counts)
)
if (length(missing_count_columns) > 0) {
  stop(
    "Centralized taxon counts are missing required columns: ",
    paste(missing_count_columns, collapse = ", ")
  )
}
if (anyDuplicated(central_counts$taxon)) {
  stop("Centralized taxon counts contain duplicate taxa.")
}

central_counts <- central_counts[
  match(names(taxon_labels), central_counts$taxon),
  required_count_columns,
  drop = FALSE
]
if (any(is.na(central_counts$taxon))) {
  stop("Centralized taxon counts do not contain all four taxa.")
}
if (any(central_counts$n_species_estimated != central_counts$n_species_equal)) {
  stop("Centralized estimated and equal-weight taxon counts disagree.")
}

expected_counts <- as.integer(central_counts$n_species_equal)
if (!identical(as.integer(observed_counts), expected_counts)) {
  count_comparison <- paste(
    sprintf(
      "%s: observed %s, centralized %s",
      unname(taxon_labels),
      format(as.integer(observed_counts), big.mark = ","),
      format(expected_counts, big.mark = ",")
    ),
    collapse = "; "
  )
  stop("Eligible species counts disagree: ", count_comparison)
}

cat(sprintf(
  "  Validated %s unique eligible species across four taxa.\n",
  format(nrow(plot_data), big.mark = ",")
))
for (i in seq_along(taxon_labels)) {
  cat(sprintf(
    "  %-11s n = %s\n",
    unname(taxon_labels[i]),
    format(observed_counts[i], big.mark = ",")
  ))
}
}

####################################################
### STEP 2: Estimate full-sample distributions ######
####################################################
{
cat("\nSTEP 2: Estimating full-sample distributions...\n")

taxon_annotations_list <- lapply(names(taxon_labels), function(tx) {
  effect_values <- plot_data$delta_s[plot_data$taxon == tx]

  data.frame(
    taxon = tx,
    negative_share_pct = 100 * mean(effect_values < 0),
    stringsAsFactors = FALSE
  )
})
taxon_annotations <- do.call(rbind, taxon_annotations_list)
rownames(taxon_annotations) <- NULL

x_limits <- unname(
  quantile(
    plot_data$delta_s,
    probs = c(0.01, 0.99),
    type = 7
  )
)
if (any(!is.finite(x_limits)) || x_limits[1] >= x_limits[2]) {
  stop("Could not calculate a valid pooled P1-P99 display window.")
}

common_bandwidth <- stats::bw.nrd0(plot_data$delta_s)
if (!is.finite(common_bandwidth) || common_bandwidth <= 0) {
  stop("Could not calculate a valid common density bandwidth.")
}

# Each density uses the complete taxon vector. The from/to arguments only set
# the shared evaluation grid; they do not trim observations before estimation.
density_list <- lapply(names(taxon_labels), function(tx) {
  effect_values <- plot_data$delta_s[plot_data$taxon == tx]
  density_estimate <- stats::density(
    effect_values,
    bw = common_bandwidth,
    from = x_limits[1],
    to = x_limits[2],
    n = 2048
  )

  data.frame(
    taxon = tx,
    effect = density_estimate$x,
    density = density_estimate$y,
    stringsAsFactors = FALSE
  )
})
density_data <- do.call(rbind, density_list)
rownames(density_data) <- NULL

y_limit <- ceiling(max(density_data$density) * 1.05 / 5) * 5
if (!is.finite(y_limit) || y_limit <= 0) {
  stop("Could not calculate a valid shared density limit.")
}

panel_labels <- sprintf(
  "%s (n = %s)",
  unname(taxon_labels),
  format(
    as.integer(observed_counts),
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
)
names(panel_labels) <- names(taxon_labels)

density_data$taxon_label <- unname(taxon_labels[density_data$taxon])
density_data$panel_label <- unname(panel_labels[density_data$taxon])
taxon_annotations$taxon_label <- unname(
  taxon_labels[taxon_annotations$taxon]
)
taxon_annotations$panel_label <- unname(
  panel_labels[taxon_annotations$taxon]
)
taxon_annotations$annotation_x <- x_limits[1] + 0.03 * diff(x_limits)
taxon_annotations$annotation_y <- 0.94 * y_limit
taxon_annotations$annotation_label <- sprintf(
  "%.0f%%\nwith ATT < 0",
  taxon_annotations$negative_share_pct
)

density_data$panel_label <- factor(
  density_data$panel_label,
  levels = panel_labels
)
taxon_annotations$panel_label <- factor(
  taxon_annotations$panel_label,
  levels = panel_labels
)

cat(sprintf(
  "  Shared bandwidth: %.6f\n  Display window: [%.6f, %.6f]\n",
  common_bandwidth,
  x_limits[1],
  x_limits[2]
))
cat(
  paste(
    sprintf(
      "  %-11s %.1f%% with ATT < 0",
      taxon_annotations$taxon_label,
      taxon_annotations$negative_share_pct
    ),
    collapse = "\n"
  ),
  "\n"
)
}

#########################################
### STEP 3: Build taxon distribution figure ###
#########################################
{
cat("\nSTEP 3: Building taxon distribution figure...\n")

x_breaks <- pretty(x_limits, n = 5)
x_breaks <- x_breaks[x_breaks >= x_limits[1] & x_breaks <= x_limits[2]]
y_breaks <- pretty(c(0, y_limit), n = 4)
y_breaks <- y_breaks[y_breaks >= 0 & y_breaks <= y_limit]

figure <- ggplot(
  density_data,
  aes(
    x = effect,
    y = density,
    group = taxon_label
  )
) +
  geom_area(
    aes(fill = taxon_label),
    alpha = 0.24,
    color = NA
  ) +
  geom_line(
    aes(color = taxon_label),
    linewidth = 1.05
  ) +
  geom_vline(
    xintercept = 0,
    color = "grey70",
    linetype = "solid",
    linewidth = 0.90
  ) +
  geom_text(
    data = taxon_annotations,
    aes(
      x = annotation_x,
      y = annotation_y,
      label = annotation_label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    color = "grey25",
    fontface = "bold",
    size = 6.6,
    lineheight = 0.95
  ) +
  facet_wrap(
    vars(panel_label),
    ncol = 2,
    scales = "fixed"
  ) +
  scale_color_manual(
    values = taxon_colors,
    guide = "none"
  ) +
  scale_fill_manual(
    values = taxon_colors,
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = function(x) sprintf("%.2f", x),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = y_breaks,
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(
    xlim = x_limits,
    ylim = c(0, y_limit),
    expand = FALSE
  ) +
  labs(
    x = "Species-level ATT on threat composite [0\u20131]",
    y = "Density"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      color = "grey89",
      linewidth = 0.35
    ),
    panel.spacing = grid::unit(0.9, "lines"),
    strip.text = element_text(
      size = 15,
      face = "bold",
      color = "grey20",
      margin = margin(5, 0, 5, 0)
    ),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "grey25"),
    plot.margin = margin(8, 14, 8, 10)
  )

ggsave(
  output_file,
  figure,
  width = 11.5,
  height = 6.4,
  dpi = 320,
  bg = "white"
)

cat("  Saved: ", output_file, "\n", sep = "")

paper_figure <- figure
paper_figure$layers[[which(vapply(paper_figure$layers,
  function(layer) inherits(layer$geom, "GeomText"), logical(1)))]]$aes_params$size <- 3.4
paper_figure <- paper_figure + theme(
  strip.text = element_text(size = 10, face = "bold"),
  axis.title = element_text(size = 9),
  axis.text = element_text(size = 8),
  panel.spacing = grid::unit(0.65, "lines"),
  plot.margin = margin(4, 8, 4, 5)
)
ggsave(file.path(paths$figures_dir, "fig.bytaxa.effect-distributions.manuscript.jpg"),
       paper_figure, width = 6.2, height = 3.4, dpi = 320, bg = "white")
}

cat("\n=== Taxon Effect Distributions Complete ===\n")
