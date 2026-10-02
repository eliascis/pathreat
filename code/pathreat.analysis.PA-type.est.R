##########################################
### pathreat.analysis.PA-type.est.R ######
### Subsample estimation by PA type #####
### Produces Table A.6 (PA_char_table) ###
##########################################

library(dplyr)
library(fixest)
library(fst)

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

####################################
### prepare subsample indicators ###
####################################
{
cat("\nPreparing subsample indicators...\n")

table(d.mbase$pa_iucn_cat)
table(d.mbase$iucn_class)
table(d.mbase$iucn_class, d.mbase$pa_iucn_cat)

table(d.mbase$size_class)
# summary(d.mbase[d.mbase$size_class=="class 1","size"])
# summary(d.mbase[d.mbase$size_class=="class 2","size"])
# summary(d.mbase[d.mbase$size_class=="class 3","size"])
# summary(d.mbase[d.mbase$size_class=="class 4","size"])

# Assign PA characteristics from treated pixel to both members of each pair
# (controls inherit their matched treated pixel's PA characteristic)
treated_chars <- d.mbase[d.mbase$treat == 1, c("matched_pair_id", "iucn_class", "size_class")]
names(treated_chars) <- c("matched_pair_id", "pair_iucn", "pair_size")
stopifnot(!anyDuplicated(treated_chars$matched_pair_id))

d.mbase <- merge(d.mbase, treated_chars, by = "matched_pair_id", all.x = TRUE)

cat("\n  IUCN class distribution (treated pixels):\n")
print(table(d.mbase$pair_iucn[d.mbase$treat == 1], useNA = "always"))
cat("\n  Size class distribution (treated pixels):\n")
print(table(d.mbase$pair_size[d.mbase$treat == 1], useNA = "always"))

rm(treated_chars); gc()
}

############################
### helper functions #######
############################
{
# format number with fixed 3 decimal places
fmt_num <- function(x) {
  if (is.na(x)) return("")
  formatC(x, format = "f", digits = 3, big.mark = "")
}

# significance stars
sig_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.01) return("***")
  if (p < 0.05) return("**")
  if (p < 0.1)  return("*")
  return("")
}

# format coefficient cell for LaTeX (value + stars, plain minus for S columns)
fmt_coef <- function(coef, pval) {
  if (is.na(coef)) return("")
  val <- fmt_num(coef)
  stars <- sig_stars(pval)
  # plain minus sign (compatible with siunitx S column type)
  if (coef < 0) val <- sub("^-", "-", val)
  if (nchar(stars) > 0) {
    paste0(val, "\\sym{", stars, "}")
  } else {
    val
  }
}

# format SE cell for LaTeX
fmt_se <- function(se) {
  if (is.na(se)) return("")
  paste0("(", fmt_num(se), ")")
}

}

############################
### define subsamples ######
############################

# detect iucn_class level names
iucn_levels <- sort(unique(d.mbase$pair_iucn[!is.na(d.mbase$pair_iucn)]))
cat("\nDetected IUCN levels:", paste(iucn_levels, collapse = ", "), "\n")

# map levels: "Strict protection" -> Strict, "Less strict" -> Multi-use, "Not Reported" -> Not Reported
iucn_map <- c(
  "Strict"       = "Strict protection",
  "Multi-use"    = "Less strict",
  "Not Reported" = "Not Reported"
)
# verify all expected levels exist
missing <- setdiff(iucn_map, iucn_levels)
if (length(missing) > 0) {
  warning("IUCN levels not found in data: ", paste(missing, collapse = ", "),
          "\n  Available: ", paste(iucn_levels, collapse = ", "))
}

subsample_defs <- list(
  list(name = "Strict",       col = "pair_iucn", val = iucn_map["Strict"]),
  list(name = "Multi-use",    col = "pair_iucn", val = iucn_map["Multi-use"]),
  list(name = "Not Reported", col = "pair_iucn", val = iucn_map["Not Reported"]),
  list(name = "Class 1",      col = "pair_size", val = "class 1"),
  list(name = "Class 2",      col = "pair_size", val = "class 2"),
  list(name = "Class 3",      col = "pair_size", val = "class 3"),
  list(name = "Class 4",      col = "pair_size", val = "class 4")
)

cat("Subsample mapping:\n")
for (ss in subsample_defs) cat(sprintf("  %s -> %s == '%s'\n", ss$name, ss$col, ss$val))

# merge() above can reorder rows, so build the canonical pair index afterwards
pair_index <- build_matched_pair_index(d.mbase)

############################
### run estimations ########
############################
{
cat("\n=== Running subsample estimations ===\n")

subsample_names <- vapply(subsample_defs, `[[`, character(1), "name")
results <- setNames(lapply(subsample_names, function(x) list()), subsample_names)
group_mean_rows <- list()

for (dep in deplist) {
  outcome_mask <- matched_pair_sample_mask(d.mbase, dep, pair_index)
  n_outcome_treated <- sum(outcome_mask & d.mbase$treat == 1)
  n_outcome_control <- sum(outcome_mask & d.mbase$treat == 0)
  if (n_outcome_treated != n_outcome_control ||
      sum(outcome_mask) != 2L * n_outcome_treated) {
    stop(sprintf("%s: canonical sample is not pair-complete", dep))
  }
  cat(sprintf("\n--- Outcome: %s (%s eligible pairs) ---\n",
              dep, format(n_outcome_treated, big.mark = ",")))

  for (ss in subsample_defs) {
    ss_filter <- d.mbase[[ss$col]] == ss$val & !is.na(d.mbase[[ss$col]])
    combined_filter <- outcome_mask & ss_filter
    n_obs <- sum(combined_filter, na.rm = TRUE)
    n_treated <- sum(combined_filter & d.mbase$treat == 1)
    n_control <- sum(combined_filter & d.mbase$treat == 0)
    n_countries <- length(unique(d.mbase$country_rast[combined_filter]))

    if (n_treated != n_control || n_obs != 2L * n_treated) {
      stop(sprintf("%s / %s: PA-characteristic subsample splits pairs",
                   dep, ss$name))
    }

    control_mean <- if (n_control > 0) {
      mean(d.mbase[[dep]][combined_filter & d.mbase$treat == 0])
    } else {
      NA_real_
    }
    treated_mean <- if (n_treated > 0) {
      mean(d.mbase[[dep]][combined_filter & d.mbase$treat == 1])
    } else {
      NA_real_
    }
    group_mean_rows[[length(group_mean_rows) + 1L]] <- data.frame(
      subsample = ss$name,
      variable = dep,
      control_mean = control_mean,
      control_n = n_control,
      treated_mean = treated_mean,
      treated_n = n_treated,
      row.names = NULL
    )

    if (n_obs < 20 || n_countries < 2) {
      cat(sprintf("  [%s] skipped (n=%d, countries=%d)\n",
                  ss$name, n_obs, n_countries))
      results[[ss$name]][[dep]] <- list(
        coef = NA,
        se = NA,
        pval = NA,
        ci_low = NA_real_,
        ci_high = NA_real_,
        n = n_obs,
        n_countries = n_countries,
        r2 = NA
      )
      next
    }

    f <- as.formula(paste(dep, "~ treat | country_rast + biome_raster"))

    fit_result <- tryCatch({
      e <- feols(f, data = d.mbase, cluster = "country_rast", subset = combined_filter)
      if (e$nobs != n_obs) {
        stop(sprintf("Model used %d rows but canonical sample contains %d", e$nobs, n_obs))
      }
      ct <- e$coeftable["treat", ]
      estimate <- regression_estimate(e)
      result <- list(
        coef = ct["Estimate"],
        se   = ct["Std. Error"],
        pval = estimate$pval,
        ci_low = estimate$ci_low,
        ci_high = estimate$ci_high,
        n    = e$nobs,
        n_countries = n_countries,
        r2   = fitstat(e, "ar2")$ar2
      )
      cat(sprintf("  [%s] coef=%s se=%s n=%s\n",
                  ss$name, fmt_num(ct["Estimate"]), fmt_num(ct["Std. Error"]),
                  format(e$nobs, big.mark = ",")))
      result
    }, error = function(err) {
      cat(sprintf("  [%s] ERROR: %s\n", ss$name, err$message))
      list(
        coef = NA,
        se = NA,
        pval = NA,
        ci_low = NA_real_,
        ci_high = NA_real_,
        n = n_obs,
        n_countries = n_countries,
        r2 = NA
      )
    })
    results[[ss$name]][[dep]] <- fit_result
  }

  rm(outcome_mask)
  invisible(gc())
}
}

##############################################
### group means ##############################
##############################################
{
cat("\n=== Computing group means by PA subcategory ===\n")

group_means_PA <- do.call(rbind, group_mean_rows)
group_means_PA <- group_means_PA[
  order(
    match(group_means_PA$subsample, subsample_names),
    match(group_means_PA$variable, deplist)
  ),
]
rownames(group_means_PA) <- NULL
rm(group_mean_rows)
gc()

cat(sprintf("  %d rows (%d subsamples x %d outcomes)\n",
            nrow(group_means_PA), length(unique(group_means_PA$subsample)),
            length(unique(group_means_PA$variable))))
}

############################
### format LaTeX table #####
############################
{
cat("\n=== Generating LaTeX table ===\n")

# variable labels with units from config (Data_summary.xlsx via deplist_scale_labels)
tex_label <- function(x) { x <- gsub("%", "\\\\%", x); gsub("&", "\\\\&", x) }
threat_labels <- setNames(tex_label(unit_label_tex(deplist_scale_labels[deplist])), deplist)

subsample_names <- sapply(subsample_defs, function(x) x$name)

# build table lines
lines <- character()
## Note: \begin{tabular}, \toprule, \bottomrule, \end{tabular} are in pathreat_03.tex
lines <- c(lines, " & \\multicolumn{3}{c}{IUCN class} & \\multicolumn{4}{c}{Size class} \\\\")
lines <- c(lines, "\\cmidrule(lr){2-4} \\cmidrule(lr){5-8}")
lines <- c(lines, " & {Strict} & {Multi-use} & {Not} & {Class 1} & {Class 2} & {Class 3} & {Class 4} \\\\")
lines <- c(lines, " & {(Ia-IV)} & {(V-VI)} & {Reported} & {(<10 km\\textsuperscript{2})} & {(10-500 km\\textsuperscript{2})} & {(500-2000 km\\textsuperscript{2})} & {(>2000 km\\textsuperscript{2})} \\\\")
lines <- c(lines, "\\midrule")

for (dep in deplist) {
  label <- threat_labels[dep]

  # coefficient row
  coef_cells <- sapply(subsample_names, function(ss) {
    r <- results[[ss]][[dep]]
    fmt_coef(r$coef, r$pval)
  })
  lines <- c(lines, paste0(label, " & ", paste(coef_cells, collapse = " & "), " \\\\"))

  # SE row
  se_cells <- sapply(subsample_names, function(ss) {
    r <- results[[ss]][[dep]]
    fmt_se(r$se)
  })
  lines <- c(lines, paste0(" & ", paste(se_cells, collapse = " & "), " \\\\"))

  # control mean row
  cmean_cells <- sapply(subsample_names, function(ss) {
    gm <- group_means_PA[group_means_PA$subsample == ss & group_means_PA$variable == dep, ]
    if (nrow(gm) == 0 || is.na(gm$control_mean[1])) return("")
    paste0("{[", fmt_num(gm$control_mean[1]), "]}")
  })
  lines <- c(lines, paste0(" & ", paste(cmean_cells, collapse = " & "), " \\\\"))
}

## \bottomrule, \\, and \end{tabular} are in pathreat_03.tex
## Strip trailing \\ from last line (main file provides \input{table}\\)
lines[length(lines)] <- sub(" *\\\\\\\\$", "", lines[length(lines)])

# write to file
outfile <- "pub/tables/tab.PA_char_table.tex"
writeLines(lines, outfile)
cat(sprintf("Saved LaTeX table to: %s\n", outfile))
}

##############################################
### PA category summary table (unmatched) ####
##############################################
{
cat("\n=== Generating PA category summary table ===\n")

# Load unmatched data — full PA universe (no country exclusions)
d.all <- fst::read_fst(paths$data_unmatched, columns = c("treat", "iucn_class", "size_class"))
d.pa <- d.all[d.all$treat == 1, ]
rm(d.all); gc()

cat(sprintf("  Total treated pixels (full unmatched): %s\n", format(nrow(d.pa), big.mark = ",")))

# IUCN class tabulation
iucn_tab <- as.data.frame(table(d.pa$iucn_class, useNA = "no"), stringsAsFactors = FALSE)
names(iucn_tab) <- c("raw_level", "n")
iucn_label_map <- c("Strict protection" = "Strict (Ia--IV)",
                     "Less strict"       = "Multi-use (V--VI)",
                     "Not Reported"      = "Not Reported")
iucn_tab$label <- iucn_label_map[iucn_tab$raw_level]
iucn_tab <- iucn_tab[!is.na(iucn_tab$label), ]
iucn_tab$area_1000 <- iucn_tab$n / 1000
iucn_tab$share <- iucn_tab$n / sum(iucn_tab$n) * 100
# order: Strict, Multi-use, Not Reported
iucn_tab <- iucn_tab[match(names(iucn_label_map), iucn_tab$raw_level), ]

# Size class tabulation
size_tab <- as.data.frame(table(d.pa$size_class, useNA = "no"), stringsAsFactors = FALSE)
names(size_tab) <- c("raw_level", "n")
size_label_map <- c("class 1" = "Class 1 ($<$10 km\\textsuperscript{2})",
                     "class 2" = "Class 2 (10--500 km\\textsuperscript{2})",
                     "class 3" = "Class 3 (500--2000 km\\textsuperscript{2})",
                     "class 4" = "Class 4 ($>$2000 km\\textsuperscript{2})")
size_tab$label <- size_label_map[size_tab$raw_level]
size_tab <- size_tab[!is.na(size_tab$label), ]
size_tab$area_1000 <- size_tab$n / 1000
size_tab$share <- size_tab$n / sum(size_tab$n) * 100
size_tab <- size_tab[match(names(size_label_map), size_tab$raw_level), ]

rm(d.pa); gc()

# Build LaTeX table body
fmt1 <- function(x) formatC(x, format = "f", digits = 1)
fmt0 <- function(x) formatC(x, format = "f", digits = 1)

tex <- character()
tex <- c(tex, " & {N pixels (10\\textsuperscript{3} km\\textsuperscript{2})} & {Share (\\%)} \\\\")
tex <- c(tex, "\\midrule")

# Panel A: IUCN class
tex <- c(tex, "\\textit{IUCN management class} & & \\\\")
for (i in seq_len(nrow(iucn_tab))) {
  tex <- c(tex, sprintf("%s & %s & %s \\\\",
                        iucn_tab$label[i], fmt1(iucn_tab$area_1000[i]), fmt0(iucn_tab$share[i])))
}

tex <- c(tex, "\\midrule")

# Panel B: Size class
tex <- c(tex, "\\textit{Size class} & & \\\\")
for (i in seq_len(nrow(size_tab))) {
  tex <- c(tex, sprintf("%s & %s & %s \\\\",
                        size_tab$label[i], fmt1(size_tab$area_1000[i]), fmt0(size_tab$share[i])))
}

# Strip trailing \\ from last line
tex[length(tex)] <- sub(" *\\\\\\\\$", "", tex[length(tex)])

outfile_summary <- "pub/tables/tab.PA_category_summary.tex"
writeLines(tex, outfile_summary)
cat(sprintf("Saved PA category summary table to: %s\n", outfile_summary))
}

############################
### save results ###########
############################
{
results_df <- do.call(rbind, lapply(names(results), function(ss) {
  do.call(rbind, lapply(names(results[[ss]]), function(dep) {
    r <- results[[ss]][[dep]]
    data.frame(
      subsample = ss,
      variable = dep,
      coef = r$coef,
      se = r$se,
      pval = r$pval,
      ci_low = r$ci_low,
      ci_high = r$ci_high,
      n_obs = r$n,
      n_countries = r$n_countries,
      r2 = r$r2,
      row.names = NULL
    )
  }))
}))

results_file <- paste0(paths$results_dir, "pathreat.PA-type.est.Rds")
saveRDS(results_df, results_file)
cat(sprintf("Saved results to: %s (%d rows)\n", results_file, nrow(results_df)))

rm(results, results_df); gc()
}

############################
### save group means #######
############################
{
gm_file <- paths$est_PA_type_group_means
saveRDS(group_means_PA, gm_file)
cat(sprintf("Saved group means to: %s\n", gm_file))

rm(d.mbase, group_means_PA); gc()
}

cat("\n=== Subsample estimation complete ===\n")
cat("Output table: pub/tables/tab.PA_char_table.tex\n")
cat("Output data:  results/pathreat.PA-type.est.Rds\n")
cat("Output data:  ", paths$est_PA_type_group_means, "\n")
