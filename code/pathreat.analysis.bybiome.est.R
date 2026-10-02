##########################################
### pathreat.analysis.bybiome.est.R ######
### By-biome effects estimation ##########
### All outcomes, country FE, clustered ##
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

############################
### biome lookup ###########
############################
{
# WWF terrestrial biome codes (1-14); 98 = rock/ice, 99 = unclassified
biome_lookup <- data.frame(
  biome_id = 1:14,
  biome_name = c(
    "Tropical & subtropical moist broadleaf forests",
    "Tropical & subtropical dry broadleaf forests",
    "Tropical & subtropical coniferous forests",
    "Temperate broadleaf & mixed forests",
    "Temperate coniferous forests",
    "Boreal forests/taiga",
    "Tropical & subtropical grasslands, savannas & shrublands",
    "Temperate grasslands, savannas & shrublands",
    "Flooded grasslands & savannas",
    "Montane grasslands & shrublands",
    "Tundra",
    "Mediterranean forests, woodlands & scrub",
    "Deserts & xeric shrublands",
    "Mangroves"
  ),
  biome_short = c(
    "Trop. moist forests",
    "Trop. dry forests",
    "Trop. conif. forests",
    "Temp. broadleaf forests",
    "Temp. conif. forests",
    "Boreal forests/taiga",
    "Trop. grasslands",
    "Temp. grasslands",
    "Flooded grasslands",
    "Montane grasslands",
    "Tundra",
    "Mediterranean",
    "Deserts & xeric",
    "Mangroves"
  ),
  stringsAsFactors = FALSE
)
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

# biomes present in data (exclude 98/99)
biomes <- sort(unique(d$biome_raster))
biomes <- biomes[biomes %in% biome_lookup$biome_id]
cat(sprintf("Number of biomes: %d\n", length(biomes)))
cat(sprintf("Outcomes: %s\n", paste(deplist, collapse = ", ")))
}

############################
### by-biome estimation ####
############################
{
cat("\n=== Running by-biome estimation for all outcomes ===\n")

results_list <- list()
pair_index <- build_matched_pair_index(d)

for (dep in deplist) {
  outcome_mask <- matched_pair_sample_mask(d, dep, pair_index)
  dep_columns <- unique(c("country_rast", "biome_raster", "treat", dep))
  d_dep <- d[outcome_mask, dep_columns, drop = FALSE]
  n_dep_treat <- sum(d_dep$treat == 1)
  n_dep_control <- sum(d_dep$treat == 0)
  if (n_dep_treat != n_dep_control || nrow(d_dep) != 2L * n_dep_treat) {
    stop(sprintf("%s: canonical sample is not pair-complete", dep))
  }
  cat(sprintf("  [%s] %s eligible pairs\n",
              dep, format(n_dep_treat, big.mark = ",")))

  for (bio in biomes) {
    d_bio <- d_dep[d_dep$biome_raster == bio, , drop = FALSE]
    n_total <- nrow(d_bio)
    n_treat <- sum(d_bio$treat == 1)
    n_control <- sum(d_bio$treat == 0)
    n_countries <- length(unique(d_bio$country_rast))
    bio_name <- biome_lookup$biome_name[biome_lookup$biome_id == bio]
    bio_short <- biome_lookup$biome_short[biome_lookup$biome_id == bio]

    if (n_treat != n_control || n_total != 2L * n_treat) {
      stop(sprintf("%s / biome %s: subsample splits matched pairs", dep, bio))
    }

    cat(sprintf("    [%d] %s: N=%s (pairs=%s, countries=%d)\n",
                bio, bio_short, format(n_total, big.mark = ","),
                format(n_treat, big.mark = ","), n_countries))

    # skip if too few eligible pairs for this outcome
    if (n_total < 10 || n_treat < 5 || n_control < 5) {
      results_list[[length(results_list) + 1]] <- data.frame(
        biome_id = bio, variable = dep,
        coef = NA, se = NA, pval = NA,
        ci_low = NA, ci_high = NA,
        n_obs = n_total, n_treat = n_treat, n_control = n_control,
        n_countries = n_countries, control_mean = NA,
        status = "skipped_insufficient_n",
        row.names = NULL
      )
      next
    }

    # control mean
    ctrl_mean <- mean(d_bio[d_bio$treat == 0, dep])

    tryCatch({
      # country FE if >1 country, otherwise no FE
      if (n_countries > 1) {
        f <- as.formula(paste(dep, "~ treat | country_rast"))
        e <- feols(f, data = d_bio, cluster = "country_rast")
      } else {
        f <- as.formula(paste(dep, "~ treat"))
        e <- feols(f, data = d_bio, vcov = "hetero")
      }
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
        biome_id = bio, variable = dep,
        coef = coef_val, se = se_val, pval = pval_val,
        ci_low = estimate$ci_low,
        ci_high = estimate$ci_high,
        n_obs = e$nobs, n_treat = n_treat, n_control = n_control,
        n_countries = n_countries, control_mean = ctrl_mean,
        status = "estimated",
        row.names = NULL
      )

    }, error = function(err) {
      cat(sprintf("    [%s] ERROR: %s\n", dep, err$message))
      results_list[[length(results_list) + 1]] <<- data.frame(
        biome_id = bio, variable = dep,
        coef = NA, se = NA, pval = NA,
        ci_low = NA, ci_high = NA,
        n_obs = n_total, n_treat = n_treat, n_control = n_control,
        n_countries = n_countries, control_mean = ctrl_mean,
        status = paste0("error: ", err$message),
        row.names = NULL
      )
    })
  }

  rm(d_dep, outcome_mask)
  invisible(gc())
}

# combine
results_bybiome <- do.call("rbind", results_list)
rownames(results_bybiome) <- NULL

# merge biome names
results_bybiome <- merge(results_bybiome, biome_lookup, by = "biome_id", all.x = TRUE)

# reorder columns
col_order <- c("biome_id", "biome_name", "biome_short", "variable",
               "coef", "se", "pval", "ci_low", "ci_high",
               "n_obs", "n_treat", "n_control", "n_countries",
               "control_mean", "status")
results_bybiome <- results_bybiome[, col_order]

# sort by biome id, then variable
results_bybiome <- results_bybiome[order(results_bybiome$biome_id,
                                         match(results_bybiome$variable, deplist)), ]
rownames(results_bybiome) <- NULL
}

############################
### summary statistics #####
############################
{
cat("\n=== By-biome estimation summary ===\n")

# count by status
status_counts <- table(results_bybiome$status)
print(status_counts)

# per-outcome summary
cat("\nSignificant results by outcome (p < 0.05):\n")
for (dep in deplist) {
  sub <- results_bybiome[results_bybiome$variable == dep &
                          results_bybiome$status == "estimated", ]
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
output_rds <- paths$est_bybiome
saveRDS(results_bybiome, output_rds)
cat(sprintf("\nSaved Rds: %s (%d rows)\n", output_rds, nrow(results_bybiome)))
}


############################
### LaTeX table ############
############################
{
cat("\n=== Generating LaTeX table ===\n")

table_vars <- deplist

# column labels from deplist_scale_labels (from Data_summary.xlsx via config)
tex_escape <- function(x) { x <- gsub("%", "\\\\%", x); gsub("&", "\\\\&", x) }
col_labels <- tex_escape(unit_label_tex(deplist_scale_labels[table_vars]))

# biomes with at least one estimated result
estimated_biomes <- unique(results_bybiome$biome_id[results_bybiome$status == "estimated"])
table_biomes <- biome_lookup[biome_lookup$biome_id %in% estimated_biomes, ]
table_biomes <- table_biomes[order(table_biomes$biome_id), ]

# build column spec. The by-biome table is short enough to fit on one
# landscape page, so emit a boxed tabular rather than longtable.
ncols <- length(table_vars)
col_spec <- paste0(">{\\raggedright\\arraybackslash}p{3.5cm} ",
                   paste(rep("p{1.2cm}", ncols), collapse = " "))

# build tabular
lines <- character()
lines <- c(lines, sprintf("\\begin{tabular}{%s}", col_spec))
lines <- c(lines, "\\toprule")

# header rows
header_cells <- paste(col_labels, collapse = " & ")
col_numbers <- paste(sprintf("(%d)", seq_len(ncols)), collapse = " & ")
lines <- c(lines, paste0("         & ", header_cells, " \\\\"))
lines <- c(lines, paste0("         & ", col_numbers, " \\\\"))
lines <- c(lines, "\\midrule")

# data rows
for (i in seq_len(nrow(table_biomes))) {
  bio_id <- table_biomes$biome_id[i]
  bio_short <- table_biomes$biome_short[i]
  bio_data <- results_bybiome[results_bybiome$biome_id == bio_id, ]

  cells <- sapply(table_vars, function(v) {
    row <- bio_data[bio_data$variable == v, ]
    if (nrow(row) == 0 || is.na(row$coef[1])) return("NA")
    fmt_coef_tex(row$coef[1], row$pval[1])
  })

  # escape ampersand in biome names
  safe_name <- gsub("&", "\\\\&", bio_short)
  lines <- c(lines, paste0(safe_name, " & ", paste(cells, collapse = " & "), " \\\\"))
}

lines <- c(lines, "\\bottomrule")
lines <- c(lines, "\\end{tabular}")

# write
outfile <- "pub/tables/tab.ATT_by_biome.tex"
writeLines(lines, outfile)
cat(sprintf("Saved LaTeX table: %s (%d biomes x %d outcomes)\n",
            outfile, nrow(table_biomes), ncols))
}

############################
### cleanup ################
############################
{
rm(d, results_list); gc()
}

cat("\n=== By-biome estimation complete ===\n")
cat("Output Rds:  ", paths$est_bybiome, "\n")
cat("Output tex:  pub/tables/tab.ATT_by_biome.tex\n")
