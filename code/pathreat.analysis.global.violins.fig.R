###################################################
### pathreat.analysis.global.violins.fig.R ########
### Country-level violin alternative to Figure 1 ##
###################################################
###
### Per outcome, plots the distribution of country-level means
### (one observation per country, weighted by matched-pair count)
### for control and protected sides side-by-side, overlaid with
### regression-standardized global levels averaged over protected
### pixels under protection and no protection. Whiskers show 95%
### percentile confidence intervals from 499 country-cluster
### bootstrap draws. In the current complete-pair sample, these
### levels equal the observed pixel-level grand means.
###
### Inputs:
###   merge.matched.fst   (via load_matched_data; cached on first run)
###   pathreat.global.est.Rds  (significance brackets sourced from
###                             country-clustered ATTs)
###   pathreat.global.fe-adjusted-levels.est.Rds
###     (regression-standardized levels and bootstrap intervals)
###   Data_summary.xlsx   (scale labels)
###
### Outputs:
###   data/store/pathreat.data.global.violins.country-means.pair-complete.Rds
###   pub/figures/fig.global.violins.country.combined.jpg
###   pub/figures/fig.global.violins.country.<v>.jpg  (16 panels)
###   results/pathreat.global.violins.country.diagnostics.csv

library(dplyr)
library(ggplot2)
library(patchwork)
library(readxl)

source("code/pathreat.analysis.config.R")

############################
### cache helpers ##########
############################

country_violin_cache_schema_version <- 1L
country_violin_contract_version <- "pathreat-country-violin-pair-complete-v1"

fingerprint_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required cache input does not exist: ", path)
  }

  info <- file.info(path)
  checksum <- unname(tools::md5sum(path))
  if (is.na(info$size) || is.na(info$mtime) || is.na(checksum)) {
    stop("Could not fingerprint required cache input: ", path)
  }

  list(
    path = path,
    size_bytes = as.numeric(info$size),
    mtime_utc = format(
      info$mtime,
      "%Y-%m-%dT%H:%M:%OS6Z",
      tz = "UTC"
    ),
    md5 = checksum
  )
}

save_rds_atomic <- function(object, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  temporary_path <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = dirname(path)
  )
  on.exit(unlink(temporary_path), add = TRUE)

  saveRDS(object, temporary_path)
  if (!file.rename(temporary_path, path)) {
    stop("Could not atomically move RDS into place: ", path)
  }

  invisible(path)
}

validate_country_violin_data <- function(data, expected_variables) {
  required_columns <- c(
    "country_rast",
    "group",
    "country_mean",
    "n_pairs",
    "variable"
  )

  if (!is.data.frame(data)) {
    return("cached country means are not a data frame")
  }
  if (!identical(names(data), required_columns)) {
    return("cached country means do not have the required field order")
  }
  if (nrow(data) == 0L) {
    return("cached country means are empty")
  }
  if (!is.numeric(data$country_rast) ||
      !is.character(data$group) ||
      !is.numeric(data$country_mean) ||
      !is.numeric(data$n_pairs) ||
      !is.character(data$variable)) {
    return("cached country means have invalid field types")
  }
  if (anyNA(data$country_rast) ||
      anyNA(data$group) ||
      anyNA(data$country_mean) ||
      anyNA(data$n_pairs) ||
      anyNA(data$variable)) {
    return("cached country means contain missing values")
  }
  if (any(!is.finite(data$country_mean)) ||
      any(!is.finite(data$n_pairs)) ||
      any(data$n_pairs <= 0) ||
      any(data$n_pairs != floor(data$n_pairs))) {
    return("cached country means contain invalid means or pair counts")
  }
  if (!setequal(unique(data$group), c("Control", "Protected"))) {
    return("cached country means have invalid treatment-side labels")
  }
  if (!setequal(unique(data$variable), expected_variables)) {
    return("cached country means do not cover the configured outcomes")
  }

  country_variable_key <- paste(
    data$variable,
    data$country_rast,
    sep = "\r"
  )
  row_key <- paste(country_variable_key, data$group, sep = "\r")
  if (anyDuplicated(row_key)) {
    return("cached country means contain duplicated country/outcome/side rows")
  }
  if (any(table(country_variable_key) != 2L)) {
    return("cached country means do not contain both sides of every pair group")
  }

  control <- data[data$group == "Control", , drop = FALSE]
  protected <- data[data$group == "Protected", , drop = FALSE]
  control_key <- paste(control$variable, control$country_rast, sep = "\r")
  protected_key <- paste(protected$variable, protected$country_rast, sep = "\r")
  protected_order <- match(control_key, protected_key)
  if (anyNA(protected_order) ||
      any(control$n_pairs != protected$n_pairs[protected_order])) {
    return("cached country means have unequal control/protected pair counts")
  }

  NULL
}

validate_country_violin_cache <- function(cache,
                                          matched_fst_fingerprint,
                                          expected_variables) {
  required_components <- c(
    "cache_schema_version",
    "aggregation_contract_version",
    "matched_fst",
    "variables",
    "country_means"
  )

  if (!is.list(cache) || !identical(names(cache), required_components)) {
    return("old or malformed cache bundle")
  }
  if (!identical(
        cache$cache_schema_version,
        country_violin_cache_schema_version
      )) {
    return("cache schema version mismatch")
  }
  if (!identical(
        cache$aggregation_contract_version,
        country_violin_contract_version
      )) {
    return("aggregation contract version mismatch")
  }
  if (!identical(cache$matched_fst, matched_fst_fingerprint)) {
    return("matched-FST fingerprint mismatch")
  }
  if (!identical(cache$variables, expected_variables)) {
    return("configured outcome list mismatch")
  }

  validate_country_violin_data(cache$country_means, expected_variables)
}

############################
### country-level means ####
############################
{
expected_variables <- unname(as.character(deplist))
matched_fst_fingerprint <- fingerprint_file(paths$data_matched)
rebuild_cache <- TRUE

if (file.exists(paths$country_violin_means)) {
  cache <- tryCatch(
    readRDS(paths$country_violin_means),
    error = function(err) err
  )

  if (inherits(cache, "error")) {
    cat(
      "Rebuilding unreadable country-level means cache:",
      cache$message,
      "\n"
    )
  } else {
    cache_problem <- validate_country_violin_cache(
      cache,
      matched_fst_fingerprint,
      expected_variables
    )
    if (is.null(cache_problem)) {
      cat(
        "Loading fingerprinted country-level means from:",
        paths$country_violin_means,
        "\n"
      )
      ctry_means <- cache$country_means
      rebuild_cache <- FALSE
    } else {
      cat(
        "Rebuilding incompatible country-level means cache:",
        cache_problem,
        "\n"
      )
    }
  }
}

if (rebuild_cache) {
  cat("Aggregating matched data to country level...\n")
  d.mbase <- load_matched_data()
  pair_index <- build_matched_pair_index(d.mbase)

  ctry_mean_list <- lapply(deplist, function(v) {
    cat(sprintf("  %-20s ", v))

    sample_mask <- matched_pair_sample_mask(d.mbase, v, pair_index)
    d_v <- d.mbase[
      sample_mask,
      c("country_rast", "treat", v),
      drop = FALSE
    ]
    n_treated <- sum(d_v$treat == 1)
    n_control <- sum(d_v$treat == 0)
    if (n_treated != n_control || nrow(d_v) != 2L * n_treated) {
      stop(sprintf("%s: violin sample is not pair-complete", v))
    }

    agg <- d_v %>%
      group_by(country_rast, treat) %>%
      summarize(
        mean_value = mean(.data[[v]]),
        n_pixels = n(),
        .groups = "drop"
      )

    ctrl <- agg %>%
      filter(treat == 0) %>%
      select(country_rast,
             mean_control = mean_value,
             n_control = n_pixels)
    prot <- agg %>%
      filter(treat == 1) %>%
      select(country_rast,
             mean_treated = mean_value,
             n_treated = n_pixels)

    if (!setequal(ctrl$country_rast, prot$country_rast)) {
      stop(sprintf("%s: control/protected country sets differ", v))
    }
    ctry <- inner_join(ctrl, prot, by = "country_rast")
    if (any(ctry$n_control != ctry$n_treated)) {
      stop(sprintf("%s: at least one country has unequal pair-side counts", v))
    }
    ctry <- ctry %>%
      mutate(n_pairs = n_control)

    long <- bind_rows(
      ctry %>%
        transmute(
          country_rast,
          group = "Control",
          country_mean = mean_control,
          n_pairs
        ),
      ctry %>%
        transmute(
          country_rast,
          group = "Protected",
          country_mean = mean_treated,
          n_pairs
        )
    )
    long$variable <- v

    cat(sprintf("%4d countries, %s matched pairs total\n",
                nrow(ctry), format(sum(ctry$n_pairs), big.mark = ",")))
    long
  })
  ctry_means <- bind_rows(ctry_mean_list)

  cache_problem <- validate_country_violin_data(
    ctry_means,
    expected_variables
  )
  if (!is.null(cache_problem)) {
    stop("Refusing to cache invalid country-level means: ", cache_problem)
  }

  cache <- list(
    cache_schema_version = country_violin_cache_schema_version,
    aggregation_contract_version = country_violin_contract_version,
    matched_fst = matched_fst_fingerprint,
    variables = expected_variables,
    country_means = ctry_means
  )
  save_rds_atomic(cache, paths$country_violin_means)
  cat(sprintf("Saved cache: %s (%s rows)\n",
              paths$country_violin_means,
              format(nrow(ctry_means), big.mark = ",")))
  rm(d.mbase, pair_index, ctry_mean_list)
}

ctry_means$group <- factor(ctry_means$group, levels = c("Control", "Protected"))
}

############################
### load supporting data ###
############################
{
est_results  <- readRDS(paths$est_global)
est_unscaled <- est_results[est_results$estimate_type == "unscaled", ]

x.labels <- read_excel(paths$data_summary) %>%
  data.frame() %>%
  filter(!is.na(threat.no)) %>%
  select(variable, threat.no, threat.category.no, threat.category, threat.label, scale.label)
}

############################
### colors (match bars) ####
############################
{
var_to_cat <- c(
  "built" = "1",
  "cropland" = "2", "planted" = "2", "pasture" = "2",
  "oil" = "3", "mining" = "3", "renewables" = "3",
  "roads" = "4", "powerlines" = "4",
  "fires" = "7", "swu" = "7", "dams" = "7",
  "light" = "9",
  "def0120_parea" = "5",
  "any_threat" = "99", "threat_composite" = "99"
)
cat_colors <- c(
  "1"  = "chocolate1",
  "2"  = "darkolivegreen1",
  "3"  = "gold",
  "4"  = "dodgerblue1",
  "5"  = "darkseagreen",
  "7"  = "firebrick1",
  "9"  = "gray50",
  "99" = "gray40"
)
}

############################
### per-panel builder ######
############################
build_panel <- function(v) {
  d_v <- ctry_means %>% filter(variable == v)

  diag <- d_v %>%
    group_by(group) %>%
    summarize(
      n_country = n(),
      zero_share = mean(country_mean == 0),
      iqr = IQR(country_mean),
      .groups = "drop"
    ) %>%
    mutate(variable = v)

  levels <- readRDS(paths$est_global_fe_levels)$summary
  levels <- levels[levels$variable == v, , drop = FALSE]
  if (nrow(levels) != 1L) stop("Missing regression-standardized levels: ", v)
  stats <- data.frame(
    group = factor(c("Control", "Protected"), levels = c("Control", "Protected")),
    mean = c(levels$adjusted_control_mean, levels$adjusted_protected_mean),
    se = c(levels$bootstrap_control_se, levels$bootstrap_protected_se),
    ci_lo = c(levels$bootstrap_control_ci_low, levels$bootstrap_protected_ci_low),
    ci_hi = c(levels$bootstrap_control_ci_high, levels$bootstrap_protected_ci_high)
  )

  label_info <- x.labels[x.labels$variable == v, ]
  if (nrow(label_info) == 0) {
    panel_title <- v
    unit_label  <- ""
  } else {
    sl <- label_info$scale.label[1]
    if (is.na(sl)) sl <- v
    if (grepl("\\[", sl)) {
      panel_title <- trimws(sub("\\s*\\[.*", "", sl))
      unit_label  <- unit_label_plot(sub(".*\\[(.*)\\].*", "[\\1]", sl))
    } else {
      panel_title <- sl
      unit_label  <- ""
    }
  }

  cat_no <- var_to_cat[v]
  if (is.na(cat_no)) cat_no <- "99"
  fill_color <- unname(cat_colors[cat_no])

  est_row <- est_unscaled[est_unscaled$variable == v, ]
  pval <- if (nrow(est_row)) est_row$pval[1] else NA
  sig_label <- "ns"
  if (!is.na(pval)) {
    if (pval <= 0.1)  sig_label <- "*"
    if (pval <= 0.05) sig_label <- "**"
    if (pval <= 0.01) sig_label <- "***"
  }

  # Adaptive y-cap: q95 default, drop to q90 for highly right-skewed.
  y_q50 <- quantile(d_v$country_mean, 0.50, na.rm = TRUE)
  y_q90 <- quantile(d_v$country_mean, 0.90, na.rm = TRUE)
  y_q95 <- quantile(d_v$country_mean, 0.95, na.rm = TRUE)
  skew_ratio <- if (y_q50 > 0) y_q95 / y_q50 else Inf
  y_top_data <- if (is.finite(skew_ratio) && skew_ratio > 5) y_q90 else y_q95
  y_min_data <- min(d_v$country_mean, na.rm = TRUE)
  y_lo <- min(0, y_min_data, min(stats$ci_lo, na.rm = TRUE))
  ci_hi_max <- max(stats$ci_hi, na.rm = TRUE)
  y_hi <- max(y_top_data, ci_hi_max * 1.5, 1e-6)

  # Per-outcome upper-bound overrides (visual; data not dropped).
  y_hi_override <- c(
    "renewables" = 0.05,
    "dams"       = 0.03
  )
  if (v %in% names(y_hi_override)) {
    # Keep model intervals below the significance bracket.
    y_hi <- max(y_hi_override[[v]], ci_hi_max / 0.84)
  }

  if (y_lo == y_hi) y_hi <- y_lo + 1
  bracket_y    <- y_lo + (y_hi - y_lo) * 0.92
  bracket_tick <- y_lo + (y_hi - y_lo) * 0.88

  ctrl_mean <- stats$mean[stats$group == "Control"]
  prot_mean <- stats$mean[stats$group == "Protected"]

  iqr_protected <- diag$iqr[diag$group == "Protected"]
  use_strip <- !is.na(iqr_protected) && iqr_protected == 0

  base <- ggplot(d_v, aes(x = group, y = country_mean)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 8),
      axis.text.x  = element_text(size = 8, angle = 0, hjust = 0.5),
      axis.text.y  = element_text(size = 7),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.line.y = element_line(color = "black"),
      axis.line.x = element_blank(),
      plot.margin = margin(5, 5, 5, 5, "pt")
    ) +
    labs(title = panel_title, y = unit_label)

  if (use_strip) {
    d_zero    <- d_v[d_v$country_mean == 0, ]
    d_nonzero <- d_v[d_v$country_mean > 0, ]
    p <- base +
      geom_jitter(data = d_zero,
                  aes(x = group, y = country_mean),
                  inherit.aes = FALSE,
                  width = 0.18, height = 0,
                  alpha = 0.15, size = 0.5, color = fill_color) +
      geom_jitter(data = d_nonzero,
                  aes(x = group, y = country_mean),
                  inherit.aes = FALSE,
                  width = 0.18, height = 0,
                  alpha = 0.7, size = 1.1, color = fill_color)
  } else {
    p <- base +
      geom_violin(aes(weight = n_pairs),
                  fill = fill_color, alpha = 0.55,
                  trim = TRUE, scale = "area", adjust = 1, color = NA)
  }

  p <- p +
    annotate("segment",
             x = 1, xend = 2,
             y = ctrl_mean, yend = prot_mean,
             linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_errorbar(data = stats,
                  aes(x = group, ymin = ci_lo, ymax = ci_hi),
                  inherit.aes = FALSE,
                  width = 0.12, color = "black", linewidth = 0.5) +
    geom_point(data = stats,
               aes(x = group, y = mean),
               inherit.aes = FALSE,
               shape = 23, size = 2.2,
               fill = fill_color, color = "black", stroke = 0.5) +
    coord_cartesian(ylim = c(y_lo, y_hi)) +
    scale_y_continuous(expand = c(0.01, 0)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "gray60") +
    annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y) +
    annotate("segment", x = 1, xend = 1, y = bracket_y, yend = bracket_tick) +
    annotate("segment", x = 2, xend = 2, y = bracket_y, yend = bracket_tick) +
    annotate("text", x = 1.5,
             y = y_lo + (y_hi - y_lo) * 0.97,
             label = sig_label, size = 3)

  list(plot = p, diag = diag, used_strip = use_strip)
}

############################
### build all panels #######
############################
{
fig_width  <- 4
fig_height <- 6

panels <- lapply(names(deplist), build_panel)
names(panels) <- names(deplist)

plot_list <- lapply(panels, `[[`, "plot")
diag_df   <- bind_rows(lapply(panels, `[[`, "diag"))
diag_df$used_strip <- vapply(panels[diag_df$variable], `[[`, logical(1), "used_strip")

diag_csv <- file.path(paths$results_dir,
                      "pathreat.global.violins.country.diagnostics.csv")
write.csv(diag_df, diag_csv, row.names = FALSE)
cat("Saved diagnostics:", diag_csv, "\n")

for (v in names(plot_list)) {
  ggsave(file.path(paths$figures_dir,
                   paste0("fig.global.violins.country.", v, ".jpg")),
         plot_list[[v]],
         units = "cm", width = fig_width, height = fig_height, dpi = 300)
}
}

############################
### combined layout ########
############################
{
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
  color = c(
    "chocolate1", "darkolivegreen1", "gold", "dodgerblue1",
    "firebrick1", "gray50", "darkseagreen", "gray40"
  ),
  x = 1:8,
  stringsAsFactors = FALSE
)
legend_plot <- ggplot(legend_data, aes(x = x, y = 1, fill = category)) +
  geom_tile(width = 0.85, height = 0.3) +
  geom_text(aes(y = 0.3, label = category),
            size = 2.8, vjust = 1, lineheight = 0.85) +
  scale_fill_manual(values = setNames(legend_data$color, legend_data$category)) +
  scale_x_continuous(limits = c(0, 9), expand = c(0, 0)) +
  ylim(-0.6, 1.2) +
  theme_void() +
  theme(legend.position = "none")

row1 <- plot_list[["built"]]    + plot_list[["cropland"]]   + plot_list[["planted"]] +
        plot_list[["pasture"]]  + plot_list[["oil"]]        + plot_layout(ncol = 5)
row2 <- plot_list[["mining"]]   + plot_list[["renewables"]] + plot_list[["roads"]] +
        plot_list[["powerlines"]] + plot_list[["def0120_parea"]] + plot_layout(ncol = 5)
row3 <- plot_list[["fires"]]    + plot_list[["swu"]]        + plot_list[["dams"]] +
        plot_list[["light"]]    + plot_list[["threat_composite"]] + plot_layout(ncol = 5)

grid <- row1 / row2 / row3
combined_plot <- grid / legend_plot + plot_layout(heights = c(1, 1, 1, 0.22))

combined_path <- file.path(paths$figures_dir,
                           "fig.global.violins.country.combined.jpg")
save_figure_data(ctry_means, combined_path)
ggsave(combined_path, combined_plot,
       units = "cm", width = 34, height = 22, dpi = 300)
cat("Saved combined:", combined_path, "\n")
}
