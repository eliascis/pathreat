##########################################
### pathreat.analysis.threat-weights.R ###
### Country threat-composite weights #####
##########################################

library(dplyr)
library(fst)
library(ggplot2)

source("code/pathreat.analysis.config.R")

out_rds <- "results/pathreat.threat_composite_weights.Rds"
out_csv <- "results/pathreat.threat_composite_weights.csv"
out_tex <- "pub/tables/tab.threat_composite_weights.tex"
out_fig <- paste0(paths$figures_dir, "fig.threat_composite_weights.heatmap.jpg")

short_labels <- c(
  built = "Built",
  cropland = "Crops",
  planted = "Plant.",
  pasture = "Pasture",
  oil = "Oil",
  mining = "Mining",
  renewables = "Renew.",
  roads = "Roads",
  powerlines = "Power",
  def0120_parea = "Deforest.",
  fires = "Fire",
  swu = "Water",
  dams = "Dams",
  light = "Light"
)

category_short_labels <- c(
  tc_residential = "Residential",
  tc_agriculture = "Agriculture",
  tc_energy = "Energy",
  tc_transport = "Transport",
  tc_biological = "Biological",
  tc_modification = "Modification",
  tc_pollution = "Pollution"
)

threat_category_lookup <- stack(threat_categories) %>%
  transmute(
    threat = as.character(values),
    category = as.character(ind)
  )

tex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%&_#$])", "\\\\\\1", x)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x <- gsub("~", "\\\\textasciitilde{}", x)
  x
}

cat("=== Loading merged threat data ===\n")
d <- read_fst(
  "data/store/pathreat.data.merge.fst",
  columns = c("country_rast", "country", threat_list_raw)
)

cat("=== Resolving presence thresholds ===\n")
thresholds <- dep_thresholds
for (v in thresholds$variable) {
  qp <- thresholds[v, "quantile_prob"]
  if (!is.na(qp)) {
    thresholds[v, "value"] <- quantile(d[[v]], qp, na.rm = TRUE)
  }
}

cat("=== Computing country-level prevalence weights ===\n")
countries <- sort(unique(d$country_rast))
prevalence_mat <- matrix(
  0,
  nrow = length(countries),
  ncol = length(threat_list_raw),
  dimnames = list(as.character(countries), threat_list_raw)
)

for (v in threat_list_raw) {
  thresh <- thresholds[v, "value"]
  dir <- thresholds[v, "dir"]
  x <- d[[v]]
  pres <- if (dir == "<") {
    !is.na(x) & x < thresh
  } else {
    !is.na(x) & x > thresh
  }
  cm <- tapply(pres, d$country_rast, mean, na.rm = TRUE)
  prevalence_mat[as.character(names(cm)), v] <- cm
}

cat("=== Normalizing weights within threat categories ===\n")
cat_prevalence_mat <- sapply(threat_categories, function(vars) {
  rowSums(prevalence_mat[, vars, drop = FALSE], na.rm = TRUE)
})
weight_mat <- prevalence_mat
for (cat_name in names(threat_categories)) {
  cat_vars <- threat_categories[[cat_name]]
  cat_prevalence <- cat_prevalence_mat[, cat_name]
  for (v in cat_vars) {
    weight_mat[, v] <- ifelse(cat_prevalence > 0, prevalence_mat[, v] / cat_prevalence, 0)
  }
}

within_cat_sums <- sapply(threat_categories, function(vars) {
  rowSums(weight_mat[, vars, drop = FALSE], na.rm = TRUE)
})
active_cat <- cat_prevalence_mat > 0
max_sum_error <- max(abs(within_cat_sums[active_cat] - 1), na.rm = TRUE)
if (is.finite(max_sum_error) && max_sum_error > 1e-8) {
  stop("Within-category threat weights do not sum to one for all active country-categories")
}
cat("Maximum within-category weight-sum error:", format(max_sum_error, scientific = TRUE), "\n")

country_lookup <- d %>%
  filter(!is.na(country_rast), !is.na(country)) %>%
  count(country_rast, country, name = "n") %>%
  arrange(country_rast, desc(n), country) %>%
  group_by(country_rast) %>%
  slice(1) %>%
  ungroup() %>%
  select(country_rast, country)

prevalence_wide <- as.data.frame(prevalence_mat)
prevalence_wide$country_rast <- as.integer(rownames(prevalence_mat))
prevalence_wide <- prevalence_wide %>%
  left_join(country_lookup, by = "country_rast") %>%
  mutate(
    country = ifelse(is.na(country), paste0("Country ", country_rast), country),
    total_prevalence = rowSums(across(all_of(threat_list_raw)), na.rm = TRUE),
    active_categories = rowSums(cat_prevalence_mat[match(country_rast, as.integer(rownames(prevalence_mat))), , drop = FALSE] > 0)
  ) %>%
  arrange(desc(total_prevalence), country) %>%
  select(country_rast, country, active_categories, total_prevalence, all_of(threat_list_raw))

weights_wide <- as.data.frame(weight_mat)
weights_wide$country_rast <- as.integer(rownames(weight_mat))
weights_wide <- weights_wide %>%
  left_join(country_lookup, by = "country_rast") %>%
  mutate(
    country = ifelse(is.na(country), paste0("Country ", country_rast), country),
    total_prevalence = rowSums(prevalence_mat[match(country_rast, as.integer(rownames(prevalence_mat))), , drop = FALSE], na.rm = TRUE),
    active_categories = rowSums(cat_prevalence_mat[match(country_rast, as.integer(rownames(prevalence_mat))), , drop = FALSE] > 0)
  ) %>%
  arrange(desc(total_prevalence), country) %>%
  select(country_rast, country, active_categories, total_prevalence, all_of(threat_list_raw))

# Final contribution weight of threat k in country c, as implemented in
# pathreat.data.threat-indices.R step 4:
#   (within-category share) x (category share of total prevalence)
#   = (wbar_kc / sum_{l in C} wbar_lc) x (sum_{l in C} wbar_lc / sum_l wbar_lc)
#   = wbar_kc / sum_l wbar_lc
# i.e. each threat's share of the country's total presence prevalence. The
# earlier division by the count of active categories corresponded to an
# equal-weight average over categories, which is not the published index.
final_weight_mat <- weight_mat
active_category_count <- rowSums(cat_prevalence_mat > 0)
total_prevalence_country <- rowSums(prevalence_mat, na.rm = TRUE)
for (v in threat_list_raw) {
  final_weight_mat[, v] <- ifelse(total_prevalence_country > 0,
                                  prevalence_mat[, v] / total_prevalence_country,
                                  0)
}

final_weight_sums <- rowSums(final_weight_mat[, threat_list_raw, drop = FALSE], na.rm = TRUE)
active_country <- total_prevalence_country > 0
max_final_sum_error <- max(abs(final_weight_sums[active_country] - 1), na.rm = TRUE)
if (is.finite(max_final_sum_error) && max_final_sum_error > 1e-8) {
  stop("Final threat weights do not sum to one for all active countries")
}
cat("Maximum final country weight-sum error:", format(max_final_sum_error, scientific = TRUE), "\n")

final_weights_wide <- as.data.frame(final_weight_mat)
final_weights_wide$country_rast <- as.integer(rownames(final_weight_mat))
final_weights_wide <- final_weights_wide %>%
  left_join(country_lookup, by = "country_rast") %>%
  mutate(
    country = ifelse(is.na(country), paste0("Country ", country_rast), country),
    active_categories = active_category_count[match(country_rast, as.integer(rownames(final_weight_mat)))],
    total_prevalence = rowSums(prevalence_mat[match(country_rast, as.integer(rownames(prevalence_mat))), , drop = FALSE], na.rm = TRUE)
  ) %>%
  arrange(desc(total_prevalence), country) %>%
  select(country_rast, country, active_categories, total_prevalence, all_of(threat_list_raw))

weights_long <- data.frame(
  country = rep(final_weights_wide$country, times = length(threat_list_raw)),
  total_prevalence = rep(final_weights_wide$total_prevalence, times = length(threat_list_raw)),
  threat = rep(threat_list_raw, each = nrow(weights_wide)),
  weight = as.vector(as.matrix(final_weights_wide[, threat_list_raw])),
  stringsAsFactors = FALSE
) %>%
  left_join(threat_category_lookup, by = "threat") %>%
  mutate(
    plot_weight = ifelse(weight > 0, weight, NA_real_),
    category = factor(category, levels = names(threat_categories),
                      labels = category_short_labels[names(threat_categories)]),
    threat = factor(threat, levels = threat_list_raw, labels = short_labels[threat_list_raw]),
    country = factor(country, levels = rev(final_weights_wide$country))
  )
max_plot_weight <- max(weights_long$weight, na.rm = TRUE)

saveRDS(
  list(
    weights = final_weights_wide,
    within_category_weights = weights_wide,
    prevalence = prevalence_wide,
    thresholds = thresholds,
    within_category_sums = within_cat_sums,
    final_weight_sums = final_weight_sums
  ),
  out_rds
)
write.csv(final_weights_wide, out_csv, row.names = FALSE)

cat("=== Writing LaTeX table fragment ===\n")
tex <- character()
for (i in seq_len(nrow(final_weights_wide))) {
  vals <- sprintf("%.3f", as.numeric(final_weights_wide[i, threat_list_raw]))
  tex <- c(
    tex,
    paste(
      c(tex_escape(final_weights_wide$country[i]), vals),
      collapse = " & "
    )
  )
}
tex <- paste0(tex, " \\\\")
writeLines(tex, out_tex)

cat("=== Writing heatmap ===\n")
p <- ggplot(weights_long, aes(x = threat, y = country, fill = plot_weight)) +
  geom_tile(color = "grey90", linewidth = 0.06) +
  facet_grid(. ~ category, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(
    colors = c("#fff7bc", "#fec44f", "#f03b20", "#bd0026"),
    limits = c(0, max_plot_weight),
    na.value = "white",
    name = "Final composite\nweight"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 3.6),
    strip.text.x = element_text(size = 7, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    panel.spacing.x = grid::unit(1.2, "mm"),
    legend.position = "bottom",
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave(out_fig, p, width = 18, height = 30, units = "cm", dpi = 300)

cat("Output RDS:   ", out_rds, "\n")
cat("Output CSV:   ", out_csv, "\n")
cat("Output table: ", out_tex, "\n")
cat("Output figure:", out_fig, "\n")
