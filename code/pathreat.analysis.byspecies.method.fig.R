##################################################
### pathreat.analysis.byspecies.method.fig.R #####
### Schematic: species-range FE estimation #######
##################################################
#
# Shows the canonical species estimator: range membership is evaluated at the
# protected pixel, the matched control is retained regardless of its own range
# membership, and the selected pairs enter a species-specific adaptive-FE
# regression.
#
# Output:
#   pub/figures/fig.byspecies.method.jpg
#

library(ggplot2)

############################
### geometry ###############
############################

range_theta <- seq(0, 2 * pi, length.out = 500)
range_radius <- 3.4 +
  0.45 * sin(2 * range_theta + 0.2) +
  0.24 * cos(3 * range_theta + 0.8) +
  0.10 * sin(6 * range_theta)
species_range <- data.frame(
  x = 4.7 + range_radius * cos(range_theta),
  y = 5.6 + range_radius * sin(range_theta)
)

pa_inside <- data.frame(
  x = c(3.55, 5.55, 5.40, 3.45),
  y = c(4.45, 4.55, 7.35, 7.20)
)
pa_outside <- data.frame(
  x = c(8.25, 9.40, 9.50, 8.35),
  y = c(7.55, 7.65, 9.15, 9.05)
)

included_pairs <- data.frame(
  pair = c("Pair 1", "Pair 2"),
  treat_x = c(4.20, 5.05),
  treat_y = c(6.35, 5.15),
  control_x = c(1.15, 8.65),
  control_y = c(3.65, 4.00)
)
excluded_pair <- data.frame(
  treat_x = 8.85,
  treat_y = 8.35,
  control_x = 7.65,
  control_y = 9.25
)

############################
### figure #################
############################

p <- ggplot() +
  geom_polygon(
    data = species_range,
    aes(x, y),
    fill = "#D6EAF8",
    color = "#2980B9",
    linewidth = 1.0,
    alpha = 0.75
  ) +
  geom_polygon(
    data = pa_inside,
    aes(x, y),
    fill = "#A9DFBF",
    color = "#1E8449",
    linewidth = 0.8,
    alpha = 0.80
  ) +
  geom_polygon(
    data = pa_outside,
    aes(x, y),
    fill = "#E5E8E8",
    color = "#7F8C8D",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  geom_segment(
    data = included_pairs,
    aes(x = treat_x, y = treat_y, xend = control_x, yend = control_y),
    color = "#95A5A6",
    linewidth = 0.55,
    linetype = "dotted"
  ) +
  geom_segment(
    data = excluded_pair,
    aes(x = treat_x, y = treat_y, xend = control_x, yend = control_y),
    color = "#BFC9CA",
    linewidth = 0.45,
    linetype = "dashed"
  ) +
  geom_point(
    data = included_pairs,
    aes(treat_x, treat_y),
    color = "#CB4335",
    shape = 15,
    size = 4.0
  ) +
  geom_point(
    data = included_pairs,
    aes(control_x, control_y),
    color = "#2E86C1",
    shape = 15,
    size = 4.0
  ) +
  geom_point(
    data = excluded_pair,
    aes(treat_x, treat_y),
    color = "#AAB7B8",
    shape = 15,
    size = 3.6
  ) +
  geom_point(
    data = excluded_pair,
    aes(control_x, control_y),
    color = "#D5DBDB",
    shape = 15,
    size = 3.6
  ) +
  geom_vline(xintercept = 10.1, color = "#D5D8DC", linewidth = 0.7) +
  annotate(
    "text",
    x = 0.25,
    y = 10.65,
    label = "A  Select complete matched pairs",
    hjust = 0,
    fontface = "bold",
    size = 4.4,
    color = "#2C3E50"
  ) +
  annotate(
    "text",
    x = 1.0,
    y = 9.65,
    label = expression(bold("Species ") * italic(s) * bold(" range")),
    hjust = 0,
    size = 4.0,
    color = "#2980B9"
  ) +
  annotate(
    "text",
    x = 4.45,
    y = 7.75,
    label = "Protected area",
    fontface = "bold",
    size = 3.5,
    color = "#1E8449"
  ) +
  annotate(
    "label",
    x = 3.1,
    y = 1.35,
    label = "Include when the protected pixel is in range",
    size = 3.15,
    color = "#2C3E50",
    fill = "white",
    linewidth = 0.25
  ) +
  annotate(
    "segment",
    x = 3.65,
    xend = 4.1,
    y = 1.75,
    yend = 5.85,
    color = "#CB4335",
    linewidth = 0.45,
    arrow = arrow(length = grid::unit(0.12, "cm"), type = "closed")
  ) +
  annotate(
    "label",
    x = 7.25,
    y = 2.55,
    label = "Keep its matched control\neven outside the range",
    size = 3.05,
    color = "#2E86C1",
    fill = "white",
    linewidth = 0.25,
    lineheight = 0.9
  ) +
  annotate(
    "segment",
    x = 7.85,
    xend = 8.45,
    y = 2.95,
    yend = 3.82,
    color = "#2E86C1",
    linewidth = 0.45,
    arrow = arrow(length = grid::unit(0.12, "cm"), type = "closed")
  ) +
  annotate(
    "text",
    x = 9.75,
    y = 6.95,
    label = "Protected pixel outside range:\nexclude the whole pair",
    hjust = 1,
    color = "#7F8C8D",
    fontface = "italic",
    size = 2.85,
    lineheight = 0.9
  ) +
  annotate(
    "text",
    x = 10.55,
    y = 10.65,
    label = "B  Estimate separately for species s",
    hjust = 0,
    fontface = "bold",
    size = 4.4,
    color = "#2C3E50"
  ) +
  annotate(
    "text",
    x = 13.1,
    y = 8.95,
    label = expression(
      Y[i]^s == alpha[s] + delta[s] * T[i] + lambda[c] + mu[b] + epsilon[i]
    ),
    size = 4.8,
    color = "#2C3E50"
  ) +
  annotate(
    "label",
    x = 13.1,
    y = 7.45,
    label = expression(atop(
      "Species-specific sample " * italic(I)[s],
      "One protected and one control row per pair"
    )),
    size = 3.25,
    color = "#2C3E50",
    fill = "#F8F9F9",
    linewidth = 0.25,
    lineheight = 0.95
  ) +
  annotate(
    "text",
    x = 10.65,
    y = 6.10,
    label = "Adaptive fixed effects",
    hjust = 0,
    fontface = "bold",
    size = 3.65,
    color = "#2C3E50"
  ) +
  annotate(
    "text",
    x = 10.85,
    y = 4.75,
    label = paste(
      "Country + biome if both vary",
      "Country only if only country varies",
      "Biome only if only biome varies",
      "No FE if neither varies",
      sep = "\n"
    ),
    hjust = 0,
    vjust = 0.5,
    size = 3.15,
    color = "#34495E",
    lineheight = 1.15
  ) +
  annotate(
    "label",
    x = 13.1,
    y = 2.25,
    label = "SE: country-clustered when countries vary;\nheteroskedasticity-robust otherwise",
    size = 3.05,
    color = "#2C3E50",
    fill = "white",
    linewidth = 0.25,
    lineheight = 0.95
  ) +
  annotate("point", x = 0.65, y = 0.35, color = "#CB4335", shape = 15, size = 3.3) +
  annotate("text", x = 0.95, y = 0.35, label = "Protected pixel", hjust = 0, size = 2.8) +
  annotate("point", x = 3.15, y = 0.35, color = "#2E86C1", shape = 15, size = 3.3) +
  annotate("text", x = 3.45, y = 0.35, label = "Matched control", hjust = 0, size = 2.8) +
  annotate("segment", x = 5.80, xend = 6.45, y = 0.35, yend = 0.35,
           color = "#95A5A6", linetype = "dotted", linewidth = 0.55) +
  annotate("text", x = 6.65, y = 0.35, label = "Retained pair", hjust = 0, size = 2.8) +
  coord_fixed(xlim = c(0, 16.3), ylim = c(0, 11.2), clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(8, 8, 8, 8),
    plot.background = element_rect(fill = "white", color = NA)
  )

############################
### save ###################
############################

ggsave(
  "pub/figures/fig.byspecies.method.jpg",
  p,
  width = 8.5,
  height = 7,
  dpi = 300,
  bg = "white"
)
cat("Saved: pub/figures/fig.byspecies.method.jpg\n")
