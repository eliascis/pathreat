##########################################
### pathreat.analysis.stats.fig.pa-by-year.R ###
### Bar charts: PAs and pixels by year ###
##########################################

library(dplyr)
library(ggplot2)

source("code/pathreat.analysis.config.R")

############################
### load data ##############
############################

d <- load_matched_data()

# keep treated pixels only
d_treat <- d %>%
  filter(treat == 1)

############################
### panel a: PAs by year ###
############################

pa_by_year <- d_treat %>%
  distinct(wdpaid, .keep_all = TRUE) %>%
  count(pa_year_designated) %>%
  rename(year = pa_year_designated)

p_pa <- ggplot(pa_by_year, aes(x = year, y = n)) +
  geom_col(fill = "steelblue", width = 0.7) +
  scale_x_continuous(breaks = seq(2001, 2020, by = 2)) +
  labs(
    x = "Designation year",
    y = "Number of PAs"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

################################
### panel b: pixels by year ####
################################

pix_by_year <- d_treat %>%
  count(pa_year_designated) %>%
  rename(year = pa_year_designated) %>%
  mutate(n_thousands = n / 1000)

p_pix <- ggplot(pix_by_year, aes(x = year, y = n_thousands)) +
  geom_col(fill = "steelblue", width = 0.7) +
  scale_x_continuous(breaks = seq(2001, 2020, by = 2)) +
  labs(
    x = "Designation year",
    y = "Number of pixels (thousands)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

############################
### save figures ###########
############################

fig_pa_file <- file.path(paths$figures_dir, "fig.sample.pa_by_year.jpg")
fig_pix_file <- file.path(paths$figures_dir, "fig.sample.pixels_by_year.jpg")

save_figure_data(pa_by_year, fig_pa_file)
ggsave(
  plot = p_pa,
  filename = fig_pa_file,
  units = "cm",
  width = 12,
  height = 9,
  dpi = 300
)

save_figure_data(pix_by_year, fig_pix_file)
ggsave(
  plot = p_pix,
  filename = fig_pix_file,
  units = "cm",
  width = 12,
  height = 9,
  dpi = 300
)

cat("Saved figures to", paths$figures_dir, "\n")
