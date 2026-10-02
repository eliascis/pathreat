##########################################
### pathreat.analysis.bycountry.est.R ####
### By-country effects estimation ########
### All outcomes, biome FE, robust SEs ###
##########################################

library(dplyr)
library(fixest)

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

# format coefficient cell for LaTeX
fmt_coef_tex <- function(coef, pval) {
  if (is.na(coef)) return("NA")
  val <- fmt_num(coef)
  stars <- sig_stars(pval)
  if (nchar(stars) > 0) {
    paste0("\\mbox{", val, "\\textsuperscript{", stars, "}}")
  } else {
    paste0("\\mbox{", val, "}")
  }
}

############################
### prepare data ###########
############################
{
d <- d.mbase

# country lookup (country_rast -> country name)
country_lookup <- unique(d[, c("country_rast", "country")])
names(country_lookup) <- c("country_id", "country_name")

countries <- sort(unique(d$country_rast))
cat(sprintf("Number of countries: %d\n", length(countries)))
cat(sprintf("Outcomes: %s\n", paste(deplist, collapse = ", ")))

}

############################
### by-country estimation ##
############################
{
cat("\n=== Running by-country estimation for all outcomes ===\n")

results_list <- list()
pair_index <- build_matched_pair_index(d)

for (dep in deplist) {
  outcome_mask <- matched_pair_sample_mask(d, dep, pair_index)
  dep_columns <- unique(c(
    "country_rast",
    "country",
    "biome_raster",
    "treat",
    dep
  ))
  d_dep <- d[outcome_mask, dep_columns, drop = FALSE]
  n_dep_treat <- sum(d_dep$treat == 1)
  n_dep_control <- sum(d_dep$treat == 0)
  if (n_dep_treat != n_dep_control || nrow(d_dep) != 2L * n_dep_treat) {
    stop(sprintf("%s: canonical sample is not pair-complete", dep))
  }
  cat(sprintf("  [%s] %s eligible pairs\n",
              dep, format(n_dep_treat, big.mark = ",")))

  for (ctry in countries) {
    d_ctry <- d_dep[d_dep$country_rast == ctry, , drop = FALSE]
    n_total <- nrow(d_ctry)
    n_treat <- sum(d_ctry$treat == 1)
    n_control <- sum(d_ctry$treat == 0)
    ctry_name <- country_lookup$country_name[country_lookup$country_id == ctry][1]

    if (n_treat != n_control || n_total != 2L * n_treat) {
      stop(sprintf("%s / %s: country subsample splits matched pairs",
                   dep, ctry_name))
    }

    # skip if too few eligible pairs for this outcome
    if (n_total < 10 || n_treat < 5 || n_control < 5) {
      results_list[[length(results_list) + 1]] <- data.frame(
        country_id = ctry, variable = dep,
        coef = NA, se = NA, pval = NA,
        ci_low = NA, ci_high = NA,
        n_obs = n_total, n_treat = n_treat, n_control = n_control,
        control_mean = NA, status = "skipped_insufficient_n",
        row.names = NULL
      )
      next
    }

    # control mean
    ctrl_mean <- mean(d_ctry[d_ctry$treat == 0, dep])

    # biome FE if >1 biome, otherwise no FE
    n_biomes_sub <- length(unique(d_ctry$biome_raster))

    tryCatch({
      if (n_biomes_sub > 1) {
        f <- as.formula(paste(dep, "~ treat | biome_raster"))
      } else {
        f <- as.formula(paste(dep, "~ treat"))
      }
      e <- feols(f, data = d_ctry, vcov = "hetero")
      if (e$nobs != n_total) {
        stop(sprintf("Model used %d rows but canonical sample contains %d",
                     e$nobs, n_total))
      }

      ct <- e$coeftable["treat", ]
      estimate <- regression_estimate(e)
      coef_val <- ct["Estimate"]
      se_val   <- ct["Std. Error"]
      pval_val <- ct["Pr(>|t|)"]

      results_list[[length(results_list) + 1]] <- data.frame(
        country_id = ctry, variable = dep,
        coef = coef_val, se = se_val, pval = pval_val,
        ci_low = estimate$ci_low,
        ci_high = estimate$ci_high,
        n_obs = e$nobs, n_treat = n_treat, n_control = n_control,
        control_mean = ctrl_mean, status = "estimated",
        row.names = NULL
      )

    }, error = function(err) {
      cat(sprintf("  [%s / %s] ERROR: %s\n", ctry_name, dep, err$message))
      results_list[[length(results_list) + 1]] <<- data.frame(
        country_id = ctry, variable = dep,
        coef = NA, se = NA, pval = NA,
        ci_low = NA, ci_high = NA,
        n_obs = n_total, n_treat = n_treat, n_control = n_control,
        control_mean = ctrl_mean, status = paste0("error: ", err$message),
        row.names = NULL
      )
    })
    # progress (every 20 countries)
    idx <- which(countries == ctry)
    if (idx %% 20 == 0) {
      cat(sprintf("    ... processed %d / %d countries\n",
                  idx, length(countries)))
    }
  }

  rm(d_dep, outcome_mask)
  invisible(gc())
}

# combine
results_bycountry <- do.call("rbind", results_list)
rownames(results_bycountry) <- NULL

# merge country names
results_bycountry <- merge(results_bycountry, country_lookup, by = "country_id", all.x = TRUE)

# reorder columns
col_order <- c("country_id", "country_name", "variable",
               "coef", "se", "pval", "ci_low", "ci_high",
               "n_obs", "n_treat", "n_control", "control_mean", "status")
results_bycountry <- results_bycountry[, col_order]

# sort by country name, then variable
results_bycountry <- results_bycountry[order(results_bycountry$country_name,
                                             match(results_bycountry$variable, deplist)), ]
rownames(results_bycountry) <- NULL
}

############################
### summary statistics #####
############################
{
cat("\n=== By-country estimation summary ===\n")

# count by status
status_counts <- table(results_bycountry$status)
print(status_counts)

# per-outcome summary
cat("\nSignificant results by outcome (p < 0.05):\n")
for (dep in deplist) {
  sub <- results_bycountry[results_bycountry$variable == dep &
                            results_bycountry$status == "estimated", ]
  n_sig <- sum(sub$pval < 0.05, na.rm = TRUE)
  n_neg <- sum(sub$coef < 0 & sub$pval < 0.05, na.rm = TRUE)
  n_pos <- sum(sub$coef > 0 & sub$pval < 0.05, na.rm = TRUE)
  cat(sprintf("  %-15s: %3d estimated, %3d sig (neg=%d, pos=%d)\n",
              dep, nrow(sub), n_sig, n_neg, n_pos))
}
}

############################
### save Rds ###############
############################
{
output_rds <- paste0(paths$results_dir, "pathreat.bycountry.est.Rds")
saveRDS(results_bycountry, output_rds)
cat(sprintf("\nSaved Rds: %s (%d rows)\n", output_rds, nrow(results_bycountry)))
}


############################
### LaTeX table ############
############################
{
cat("\n=== Generating LaTeX table (Table A.7) ===\n")

table_vars <- deplist

# column labels from deplist_scale_labels (from Data_summary.xlsx via config)
tex_escape <- function(x) { x <- gsub("%", "\\\\%", x); gsub("&", "\\\\&", x) }
col_labels <- tex_escape(unit_label_tex(deplist_scale_labels[table_vars]))

# countries in alphabetical order (only those with at least one estimated result)
estimated_countries <- unique(results_bycountry$country_name[results_bycountry$status == "estimated"])
table_countries <- sort(estimated_countries)

# build column spec (longtable for multi-page)
ncols <- length(table_vars)
col_spec <- paste0(">{\\raggedright} p{2cm} ", paste(rep("p{1.2cm}", ncols), collapse = " "))

# build longtable
lines <- character()
lines <- c(lines, sprintf("\\begin{longtable}{%s}", col_spec))
lines <- c(lines, "\\caption{Protected areas effect by country} \\label{tab.ATT_by_country} \\\\")
lines <- c(lines, "\\toprule")

# header rows
header_cells <- paste(col_labels, collapse = " & ")
col_numbers <- paste(sprintf("(%d)", seq_len(ncols)), collapse = " & ")
lines <- c(lines, paste0("         & ", header_cells, " \\\\"))
lines <- c(lines, paste0("         & ", col_numbers, " \\\\"))
lines <- c(lines, "\\midrule")
lines <- c(lines, "\\endfirsthead")

# continuation header
lines <- c(lines, paste0("\\multicolumn{", ncols + 1, "}{l}{\\textit{Table \\ref{tab.ATT_by_country} continued}} \\\\"))
lines <- c(lines, "\\toprule")
lines <- c(lines, paste0("         & ", header_cells, " \\\\"))
lines <- c(lines, paste0("         & ", col_numbers, " \\\\"))
lines <- c(lines, "\\midrule")
lines <- c(lines, "\\endhead")

# continuation footer
lines <- c(lines, "\\midrule")
lines <- c(lines, paste0("\\multicolumn{", ncols + 1, "}{r}{\\textit{Continued on next page}} \\\\"))
lines <- c(lines, "\\endfoot")

# final footer
lines <- c(lines, "\\bottomrule")
lines <- c(lines, "\\endlastfoot")

# data rows
for (ctry_name in table_countries) {
  ctry_data <- results_bycountry[results_bycountry$country_name == ctry_name, ]

  cells <- sapply(table_vars, function(v) {
    row <- ctry_data[ctry_data$variable == v, ]
    if (nrow(row) == 0 || is.na(row$coef[1])) return("NA")
    fmt_coef_tex(row$coef[1], row$pval[1])
  })

  # escape ampersand in country names
  safe_name <- gsub("&", "\\\\&", ctry_name)
  lines <- c(lines, paste0(safe_name, " & ", paste(cells, collapse = " & "), " \\\\"))
}

lines <- c(lines, "\\end{longtable}")

# write
outfile <- "pub/tables/tab.ATT_by_country.tex"
writeLines(lines, outfile)
cat(sprintf("Saved LaTeX table: %s (%d countries x %d outcomes)\n",
            outfile, length(table_countries), ncols))
}

############################
### cleanup ################
############################
{
rm(d, results_list); gc()
}

cat("\n=== By-country estimation complete ===\n")
cat("Output Rds:  ", paste0(paths$results_dir, "pathreat.bycountry.est.Rds"), "\n")
cat("Output tex:  pub/tables/tab.ATT_by_country.tex\n")
