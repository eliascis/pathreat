##########################################
### pathreat.analysis.covbalance.R #######
### Covariate balance + VIF diagnostics ##
##########################################

library(dplyr)
library(ggplot2)
library(fst)

source("code/pathreat.analysis.config.R")

vif_sample_n <- 200000L

w_mean <- function(x, w) {
  stats::weighted.mean(x, w, na.rm = TRUE)
}

w_var <- function(x, w) {
  ok <- !is.na(x) & !is.na(w)
  x <- x[ok]
  w <- w[ok]
  if (length(x) < 2 || sum(w) <= 0) return(NA_real_)
  mu <- stats::weighted.mean(x, w)
  sum(w * (x - mu)^2) / sum(w)
}

######################
### load data ########
######################
{
if (!"d.du" %in% ls()) {
  d.du <- read_fst(paths$data_unmatched)
}
if (!"d.mbase" %in% ls()) {
  d.mbase <- load_matched_data()
}
}

##################################
### matching covariate balance ###
##################################
{
du <- d.du
du$weight <- 1
dm <- d.mbase

# statistics across samples
a <- lapply(
  mlist,
  function(cv) {
    print(cv)

    # unmatched differences
    mean.u.t <- w_mean(du[du$treat == 1, cv], du[du$treat == 1, "weight"])
    mean.u.c <- w_mean(du[du$treat == 0, cv], du[du$treat == 0, "weight"])
    sd.u.t   <- w_var(du[du$treat == 1, cv], du[du$treat == 1, "weight"])
    sd.u.c   <- w_var(du[du$treat == 0, cv], du[du$treat == 0, "weight"])
    sdiff.u  <- (mean.u.t - mean.u.c) / sqrt((sd.u.t + sd.u.c) / 2)
    ttest.u  <- t.test(du[du$treat == 1, cv], du[du$treat == 0, cv], paired = F)

    # matched differences
    mean.m.t <- w_mean(dm[dm$treat == 1, cv], dm[dm$treat == 1, "weight"])
    mean.m.c <- w_mean(dm[dm$treat == 0, cv], dm[dm$treat == 0, "weight"])
    sd.m.t   <- w_var(dm[dm$treat == 1, cv], dm[dm$treat == 1, "weight"])
    sd.m.c   <- w_var(dm[dm$treat == 0, cv], dm[dm$treat == 0, "weight"])
    sdiff.m  <- (mean.m.t - mean.m.c) / sqrt((sd.u.t + sd.u.c) / 2)
    ttest.m  <- t.test(dm[dm$treat == 1, cv], dm[dm$treat == 0, cv], paired = F)

    # smd % improvement
    smd.improve <- 100 * (abs(sdiff.u) - abs(sdiff.m)) / abs(sdiff.u)
    # var.ratio
    var.r.u <- sd.u.t / sd.u.c
    var.r.m <- sd.m.t / sd.m.c
    var.r.improve <- 100 * (abs(log(var.r.u)) - abs(log(var.r.m))) / abs(log(var.r.u))

    s <- data.frame(
      cov = cv,
      mean.u.t = round(mean.u.t, 3),
      mean.m.t = round(mean.m.t, 3),
      mean.u.c = round(mean.u.c, 3),
      mean.m.c = round(mean.m.c, 3),
      diff.u = round(mean.u.t - mean.u.c, 3),
      diff.u.sd = round(ttest.u$stderr, 3),
      diff.u.pval = round(ttest.u$p.value),
      diff.m = round(mean.m.t - mean.m.c),
      diff.m.sd = round(ttest.m$stderr, 3),
      diff.m.pval = round(ttest.m$p.value),
      sdiff.u = round(sdiff.u, 3),
      sdiff.m = round(sdiff.m, 3),
      sdiff.improve = round(smd.improve, 3),
      var.ratio.u = round(var.r.u, 3),
      var.ratio.m = round(var.r.m, 3),
      var.ratio.improve = round(var.r.improve, 3),
      stringsAsFactors = F
    )
    return(s)
  }
)
a <- do.call("rbind", a)
a$varnames <- mlist_labels[a$cov]
a <- a %>% relocate(varnames, .after = cov)

# add averages row
x <- data.frame(t(rep(NA, length(a))))
names(x) <- names(a)
x$cov <- "all"
x$sdiff.u <- mean(a$sdiff.u)
x$sdiff.m <- mean(a$sdiff.m)
a <- rbind(a, x)

# save
options(xtable.include.rownames = FALSE)
write.csv(a, paste0(paths$results_dir, "pathreat.covbalance.csv"), row.names = FALSE)
d.covbal <- a
}

##################################
### variance inflation factors ###
##################################
{
cat("\n=== Computing VIF diagnostics ===\n")

compute_vif <- function(d, sample_n = vif_sample_n) {
  set.seed(42)
  if (nrow(d) > sample_n) {
    d <- d[sample.int(nrow(d), sample_n), ]
  }

  x <- as.data.frame(d[, mlist])
  x <- x[complete.cases(x), , drop = FALSE]
  cat("  VIF sample size:", format(nrow(x), big.mark = ","), "\n")

  r <- stats::cor(x)
  vif <- diag(solve(r))
  names(vif) <- colnames(r)
  round(vif, 2)
}

cat("  Unmatched sample...\n")
vif_unmatched <- compute_vif(d.du)

cat("  Matched sample...\n")
vif_matched <- compute_vif(d.mbase)

d.vif <- data.frame(
  cov = mlist,
  vif.u = as.numeric(vif_unmatched[mlist]),
  vif.m = as.numeric(vif_matched[mlist]),
  stringsAsFactors = FALSE
)
}

##########################
### visualizing balance ##
##########################
{
s <- d.covbal
i <- which(s$cov == "all")
s <- s[-i, ]
s$varnames <- with(s, reorder(varnames, sdiff.u))

p <- ggplot() +
  geom_vline(xintercept = 0.1, linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = -0.1, linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 1) +
  geom_hline(yintercept = 0.5, linetype = "solid", linewidth = 1) +
  geom_point(data = s, aes(x = sdiff.u, y = varnames, colour = "Unmatched"), size = 5) +
  geom_point(data = s, aes(x = sdiff.m, y = varnames, colour = "Matched"), size = 5) +
  xlim(c(-0.2, 0.5)) +
  labs(x = "Std. Mean Differences", y = "") +
  theme(text = element_text(size = 20), panel.background = element_blank(),
        legend.position = "bottom") +
  scale_colour_manual("", values = c("Unmatched" = "#5e3c99", "Matched" = "#e66101"))

ggsave(plot = p, paste0(paths$figures_dir, "fig.match.balance.jpg"),
       units = "cm", width = 24, height = 10, dpi = 300)
}

##################################
### LaTeX diagnostics table ######
##################################
{
cat("\n=== Writing matching diagnostics table ===\n")

bal <- d.covbal[d.covbal$cov != "all", ]
diag_table <- bal %>%
  left_join(d.vif, by = "cov")

tex_lines <- c(
  "& \\multicolumn{3}{c}{Std.\\ mean difference} & \\multicolumn{3}{c}{Variance ratio} & \\multicolumn{2}{c}{VIF} \\\\",
  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7} \\cmidrule(lr){8-9}",
  "Predictor variable & {Unmatched} & {Matched} & {\\% impr.} & {Unmatched} & {Matched} & {\\% impr.} & {Unmatched} & {Matched} \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(diag_table))) {
  tex_lines <- c(tex_lines,
    sprintf("%s & %s & %s & %d & %.3f & %.3f & %d & %.2f & %.2f \\\\",
            diag_table$varnames[i],
            format(diag_table$sdiff.u[i], nsmall = 3),
            format(diag_table$sdiff.m[i], nsmall = 3),
            round(diag_table$sdiff.improve[i]),
            diag_table$var.ratio.u[i],
            diag_table$var.ratio.m[i],
            round(diag_table$var.ratio.improve[i]),
            diag_table$vif.u[i],
            diag_table$vif.m[i]))
}

tex_lines[length(tex_lines)] <- sub(" *\\\\\\\\$", "", tex_lines[length(tex_lines)])

writeLines(tex_lines, "pub/tables/tab.matching_diagnostics.tex")
write.csv(diag_table, paste0(paths$results_dir, "pathreat.matching_diagnostics.csv"),
          row.names = FALSE)
cat("Output: pub/tables/tab.matching_diagnostics.tex\n")
}

cat("\n=== Covariate diagnostics complete ===\n")
cat("Output: ", paths$results_dir, "pathreat.covbalance.csv\n")
cat("Output: ", paths$results_dir, "pathreat.matching_diagnostics.csv\n")
cat("Figure: ", paths$figures_dir, "fig.match.balance.jpg\n")
cat("Table:  pub/tables/tab.matching_diagnostics.tex\n")
