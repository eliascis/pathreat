##################################################
### pathreat.analysis.byspecies.fig.R ###########
### Top and bottom species effects ###############
##################################################
#
# Plots the 10 most negative and 10 most positive range-based
# species effects from pathreat.analysis.byspecies.R as ranked
# pointrange panels with 95% confidence intervals on a shared x-axis.
# The script also attempts to download Wikimedia
# species photos for the plotted rows, creates circular center-cropped
# thumbnails, and adds them in a right-hand strip.
#
# Input:
#   results/pathreat.byspecies.est.Rds
#
# Output:
#   pub/figures/fig.byspecies.delta.top-bottom.est.jpg
#   pub/figures/fig.byspecies.coverage-vs-pct.heatmap.jpg
#   pub/figures/fig.byspecies.coverage-vs-pct.heatmap-noNA.jpg
#   pub/figures/fig.byspecies.coverage-vs-reduction.jpg
#

library(dplyr)
library(ggplot2)
library(httr)
library(jsonlite)
library(cowplot)
library(png)
library(jpeg)
library(grid)
library(patchwork)

source("code/pathreat.analysis.config.R")


############################
### helpers ################
############################

sanitize_file_stub <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- tolower(gsub("[^a-zA-Z0-9]+", "_", x))
  x <- gsub("^_+|_+$", "", x)

  ifelse(nchar(x) == 0, "species", x)
}

species_initials <- function(x) {
  bits <- strsplit(x, " +")[[1]]
  bits <- bits[nzchar(bits)]

  if (length(bits) == 0) {
    return("?")
  }
  if (length(bits) == 1) {
    return(toupper(substr(bits[1], 1, 2)))
  }

  paste0(
    toupper(substr(bits[1], 1, 1)),
    toupper(substr(bits[2], 1, 1))
  )
}

photo_overrides <- c(
  "Presbytis chrysomelas" = "https://upload.wikimedia.org/wikipedia/commons/c/c6/Presbytis_chrysomelas.jpg",
  "Pelophylax caralitanus" = "https://upload.wikimedia.org/wikipedia/commons/c/c2/Pelophylax_caralitanus_Bey%C5%9Fehir_%2801%29.jpg"
)

thumb_overrides <- list(
  "Presbytis chrysomelas" = list(
    crop_frac = 0.84,
    focus_x = 0.42,
    focus_y = 0.44
  )
)

bad_photo_patterns <- c(
  "(^|[^a-z])(dis|distribution|map|range|locator|location|area)([^a-z]|$)",
  "blank[_ -]?map",
  "loc[_ -]?map",
  "figure[_ -]?[0-9]",
  "plate",
  "diagram",
  "chart",
  "catalogue",
  "catalog",
  "bulletin",
  "\\.pdf$"
)

is_bad_photo_candidate <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(x)) {
    return(TRUE)
  }

  lower <- tolower(x)

  any(vapply(
    bad_photo_patterns,
    function(pattern) grepl(pattern, lower, perl = TRUE),
    logical(1)
  ))
}

score_photo_candidate <- function(x, species_name) {
  if (is_bad_photo_candidate(x)) {
    return(-Inf)
  }

  lower <- tolower(x)
  ext <- tolower(tools::file_ext(sub("\\?.*$", "", lower)))
  species_bits <- strsplit(sanitize_file_stub(species_name), "_")[[1]]
  species_bits <- species_bits[nzchar(species_bits)]

  score <- 0

  if (ext %in% c("jpg", "jpeg")) {
    score <- score + 20
  } else if (ext == "png") {
    score <- score + 8
  }

  score <- score + 5 * sum(vapply(
    species_bits,
    function(bit) grepl(bit, lower, fixed = TRUE),
    logical(1)
  ))

  if (grepl("commons|wikipedia", lower)) {
    score <- score + 2
  }

  score
}

wiki_get_json <- function(url) {
  resp <- try(
    httr::GET(
      url,
      httr::user_agent("pathreat-byspecies-fig/1.0"),
      httr::timeout(20)
    ),
    silent = TRUE
  )

  if (inherits(resp, "try-error") || httr::status_code(resp) >= 300) {
    return(NULL)
  }

  txt <- httr::content(resp, as = "text", encoding = "UTF-8")
  out <- try(jsonlite::fromJSON(txt, simplifyVector = TRUE), silent = TRUE)

  if (inherits(out, "try-error")) {
    return(NULL)
  }

  out
}

wiki_get_summary <- function(title) {
  url <- paste0(
    "https://en.wikipedia.org/api/rest_v1/page/summary/",
    URLencode(gsub(" ", "_", title), reserved = TRUE)
  )

  wiki_get_json(url)
}

wiki_search_titles <- function(query, limit = 5L) {
  url <- paste0(
    "https://en.wikipedia.org/w/api.php?action=query&format=json&list=search",
    "&srlimit=", limit,
    "&srsearch=", URLencode(paste0("\"", query, "\""), reserved = TRUE)
  )

  out <- wiki_get_json(url)
  if (is.null(out) || is.null(out$query$search$title)) {
    return(character(0))
  }

  hits <- out$query$search$title
  hits <- as.character(hits)
  hits[nzchar(hits)]
}

commons_search_files <- function(query, limit = 8L) {
  url <- paste0(
    "https://commons.wikimedia.org/w/api.php?action=query&format=json&generator=search",
    "&gsrnamespace=6",
    "&gsrlimit=", limit,
    "&gsrsearch=", URLencode(query, reserved = TRUE),
    "&prop=imageinfo",
    "&iiprop=url"
  )

  out <- wiki_get_json(url)
  if (is.null(out) || is.null(out$query$pages)) {
    return(data.frame())
  }

  pages <- out$query$pages
  rows <- lapply(pages, function(page) {
    image_url <- NA_character_
    if (!is.null(page$imageinfo[[1]]$url) && nzchar(page$imageinfo[[1]]$url)) {
      image_url <- page$imageinfo[[1]]$url
    }

    data.frame(
      title = if (!is.null(page$title)) page$title else NA_character_,
      image_url = image_url,
      index = if (!is.null(page$index)) page$index else Inf,
      stringsAsFactors = FALSE
    )
  })

  bind_rows(rows) %>%
    filter(is.finite(index), !is.na(image_url), nzchar(image_url))
}

resolve_species_photo <- function(species_name) {
  if (species_name %in% names(photo_overrides)) {
    return(
      list(
        page_title = species_name,
        image_url = unname(photo_overrides[[species_name]])
      )
    )
  }

  candidate_titles <- c(species_name, wiki_search_titles(species_name, limit = 5L))
  candidate_titles <- unique(candidate_titles[nzchar(candidate_titles)])

  for (title in candidate_titles) {
    summary <- wiki_get_summary(title)
    if (is.null(summary)) {
      next
    }

    image_url <- NULL
    if (!is.null(summary$thumbnail$source) && nzchar(summary$thumbnail$source)) {
      image_url <- summary$thumbnail$source
    } else if (!is.null(summary$originalimage$source) &&
               nzchar(summary$originalimage$source)) {
      image_url <- summary$originalimage$source
    }

    if (!is.null(image_url) &&
        !is_bad_photo_candidate(paste(title, image_url))) {
      return(
        list(
          page_title = if (!is.null(summary$title)) summary$title else title,
          image_url = image_url
        )
      )
    }
  }

  commons_hits <- commons_search_files(species_name, limit = 8L)
  if (nrow(commons_hits) > 0) {
    commons_hits <- commons_hits %>%
      mutate(
        score = vapply(
          paste(title, image_url),
          score_photo_candidate,
          numeric(1),
          species_name = species_name
        )
      ) %>%
      filter(is.finite(score)) %>%
      arrange(desc(score), index, title)

    if (nrow(commons_hits) > 0) {
      return(
        list(
          page_title = commons_hits$title[1],
          image_url = commons_hits$image_url[1]
        )
      )
    }
  }

  NULL
}

pick_best_cached_photo <- function(paths, species_name) {
  if (length(paths) == 0) {
    return(NA_character_)
  }

  scores <- vapply(
    basename(paths),
    score_photo_candidate,
    numeric(1),
    species_name = species_name
  )

  if (!any(is.finite(scores))) {
    return(paths[1])
  }

  paths[order(scores, decreasing = TRUE)][1]
}

download_species_photo <- function(species_id, species_name, cache_dir) {
  stub <- sprintf("%s_%s", species_id, sanitize_file_stub(species_name))
  photo_meta <- NULL
  if (species_name %in% names(photo_overrides)) {
    photo_meta <- list(
      page_title = species_name,
      image_url = unname(photo_overrides[[species_name]])
    )
  }

  existing <- list.files(
    cache_dir,
    pattern = paste0("^", stub, "\\.(jpg|jpeg|png)$"),
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (is.null(photo_meta) && length(existing) > 0) {
    return(pick_best_cached_photo(existing, species_name))
  }

  if (is.null(photo_meta)) {
    photo_meta <- resolve_species_photo(species_name)
    if (is.null(photo_meta)) {
      return(NA_character_)
    }
  }

  ext <- tolower(tools::file_ext(sub("\\?.*$", "", photo_meta$image_url)))
  if (!ext %in% c("jpg", "jpeg", "png")) {
    ext <- "jpg"
  }

  out_path <- file.path(cache_dir, paste0(stub, ".", ext))
  if (file.exists(out_path)) {
    return(out_path)
  }

  tmp_path <- paste0(out_path, ".tmp")

  resp <- try(
    httr::GET(
      photo_meta$image_url,
      httr::user_agent("pathreat-byspecies-fig/1.0"),
      httr::timeout(30),
      httr::write_disk(tmp_path, overwrite = TRUE)
    ),
    silent = TRUE
  )

  if (inherits(resp, "try-error") || httr::status_code(resp) >= 300 || !file.exists(tmp_path)) {
    if (file.exists(tmp_path)) {
      unlink(tmp_path)
    }
    return(NA_character_)
  }

  file.rename(tmp_path, out_path)
  out_path
}

read_bitmap <- function(path) {
  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("jpg", "jpeg")) {
    img <- jpeg::readJPEG(path, native = FALSE)
  } else if (ext == "png") {
    img <- png::readPNG(path, native = FALSE)
  } else {
    stop(sprintf("Unsupported image format: %s", path))
  }

  if (length(dim(img)) == 2) {
    h <- nrow(img)
    w <- ncol(img)
    tmp <- array(0, dim = c(h, w, 3))
    tmp[, , 1] <- img
    tmp[, , 2] <- img
    tmp[, , 3] <- img
    img <- tmp
  }

  h <- dim(img)[1]
  w <- dim(img)[2]
  channels <- dim(img)[3]

  rgba <- array(1, dim = c(h, w, 4))
  rgba[, , 1:3] <- img[, , 1:3]
  if (channels >= 4) {
    rgba[, , 4] <- img[, , 4]
  }

  rgba
}

trim_bitmap <- function(img,
                        alpha_threshold = 0.02,
                        white_threshold = 0.97,
                        pad_frac = 0.04) {
  alpha <- img[, , 4]
  near_white <- img[, , 1] > white_threshold &
    img[, , 2] > white_threshold &
    img[, , 3] > white_threshold

  content_mask <- alpha > alpha_threshold & !near_white
  if (!any(content_mask)) {
    content_mask <- alpha > alpha_threshold
  }
  if (!any(content_mask)) {
    return(img)
  }

  rows <- which(rowSums(content_mask) > 0)
  cols <- which(colSums(content_mask) > 0)

  row_pad <- max(1L, floor(length(rows) * pad_frac))
  col_pad <- max(1L, floor(length(cols) * pad_frac))

  row_min <- max(1L, min(rows) - row_pad)
  row_max <- min(dim(img)[1], max(rows) + row_pad)
  col_min <- max(1L, min(cols) - col_pad)
  col_max <- min(dim(img)[2], max(cols) + col_pad)

  img[row_min:row_max, col_min:col_max, , drop = FALSE]
}

create_round_thumbnail <- function(src_path,
                                   thumb_path,
                                   crop_frac = 0.76,
                                   focus_x = 0.5,
                                   focus_y = 0.5,
                                   overwrite = FALSE) {
  if (!overwrite && file.exists(thumb_path)) {
    return(thumb_path)
  }

  img <- trim_bitmap(read_bitmap(src_path))
  h <- dim(img)[1]
  w <- dim(img)[2]

  crop_side <- max(10L, floor(min(h, w) * crop_frac))
  crop_side <- min(crop_side, min(h, w))

  row_center <- 1 + (h - 1) * focus_y
  col_center <- 1 + (w - 1) * focus_x
  row_start <- round(row_center - crop_side / 2)
  col_start <- round(col_center - crop_side / 2)
  row_start <- max(1L, min(row_start, h - crop_side + 1L))
  col_start <- max(1L, min(col_start, w - crop_side + 1L))

  img_crop <- img[
    row_start:(row_start + crop_side - 1L),
    col_start:(col_start + crop_side - 1L),
    ,
    drop = FALSE
  ]

  yy <- seq(-1, 1, length.out = dim(img_crop)[1])
  xx <- seq(-1, 1, length.out = dim(img_crop)[2])
  dist <- outer(yy, xx, function(y, x) sqrt(x^2 + y^2))
  alpha_mask <- pmin(1, pmax(0, (1.0 - dist) / 0.015))
  border_mask <- dist > 0.97 & dist <= 1.0

  img_crop[, , 4] <- img_crop[, , 4] * alpha_mask

  # Draw the circular border directly into the PNG so the visible edge
  # matches the image mask exactly.
  img_crop[, , 1][border_mask] <- 0.78
  img_crop[, , 2][border_mask] <- 0.78
  img_crop[, , 3][border_mask] <- 0.78
  img_crop[, , 4][border_mask] <- 1

  png::writePNG(img_crop, target = thumb_path)
  thumb_path
}

make_effect_panel <- function(panel_data, taxon_colors, x_limits) {
  y_breaks <- panel_data$y_pos
  y_labels <- as.character(panel_data$species_name)

  ggplot(panel_data, aes(x = delta_s, y = y_pos)) +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey45") +
    geom_segment(
      aes(
        x = ci_low,
        xend = ci_high,
        y = y_pos,
        yend = y_pos,
        color = taxon_label
      ),
      linewidth = 0.9
    ) +
    geom_point(
      aes(color = taxon_label),
      size = 2.8
    ) +
    scale_y_continuous(
      limits = c(0.5, nrow(panel_data) + 0.5),
      breaks = y_breaks,
      labels = function(x) y_labels[match(x, y_breaks)],
      expand = c(0, 0)
    ) +
    scale_x_continuous(limits = x_limits, expand = expansion(mult = c(0, 0))) +
    scale_color_manual(values = taxon_colors, name = NULL, drop = FALSE) +
    labs(
      title = unique(panel_data$panel),
      x = expression(hat(delta)[s] ~ "(95% CI)"),
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 9, face = "italic"),
      axis.title.x = element_text(margin = margin(t = 8)),
      plot.title = element_text(face = "bold", margin = margin(b = 8)),
      legend.position = "none",
      plot.margin = margin(10, 5, 10, 15)
    )
}

make_placeholder_grob <- function(label) {
  grid::grobTree(
    grid::circleGrob(
      r = unit(0.485, "npc"),
      gp = grid::gpar(fill = "white", col = "grey75", lwd = 1.1)
    ),
    grid::textGrob(
      label,
      gp = grid::gpar(col = "grey45", fontsize = 11, fontface = "bold")
    )
  )
}

make_photo_panel <- function(panel_data) {
  thumb_half_width <- 0.355
  thumb_half_height <- 0.655

  p <- ggplot(panel_data, aes(x = 0.5, y = y_pos)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(
      limits = c(0.5, nrow(panel_data) + 0.5),
      breaks = panel_data$y_pos,
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
      plot.margin = margin(28, 15, 18, 0)
    )

  for (i in seq_len(nrow(panel_data))) {
    row <- panel_data[i, ]

    if (!is.na(row$thumb_path) && file.exists(row$thumb_path)) {
      grob <- grid::rasterGrob(
        png::readPNG(row$thumb_path),
        interpolate = TRUE
      )

      p <- p + annotation_custom(
        grob = grob,
        xmin = row$thumb_x - thumb_half_width,
        xmax = row$thumb_x + thumb_half_width,
        ymin = row$y_pos - thumb_half_height,
        ymax = row$y_pos + thumb_half_height
      )
    } else {
      p <- p + annotation_custom(
        grob = make_placeholder_grob(row$photo_label),
        xmin = row$thumb_x - thumb_half_width,
        xmax = row$thumb_x + thumb_half_width,
        ymin = row$y_pos - thumb_half_height,
        ymax = row$y_pos + thumb_half_height
      )
    }
  }

  p
}

make_legend_panel <- function(taxon_colors) {
  legend_data <- data.frame(
    x = c(0.06, 0.32, 0.58, 0.82),
    y = 0.5,
    label = names(taxon_colors),
    stringsAsFactors = FALSE
  )

  ggplot(legend_data, aes(x = x, y = y)) +
    geom_point(aes(color = label), size = 3.2) +
    geom_text(
      aes(x = x + 0.045, label = label),
      hjust = 0,
      size = 3.5,
      color = "grey20"
    ) +
    scale_color_manual(values = taxon_colors, drop = FALSE) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 20, 0, 20)
    )
}


############################
### load data ##############
############################
{
cat("=== Loading species estimates ===\n")

sp_path <- file.path(paths$results_dir, "pathreat.byspecies.est.Rds")
if (!file.exists(sp_path)) {
  stop(sprintf("Missing %s. Run code/pathreat.analysis.byspecies.R first.", sp_path))
}

sp_est <- readRDS(sp_path)
cat(sprintf("  Species estimates: %s rows\n", format(nrow(sp_est), big.mark = ",")))

required_cols <- c(
  "species_id",
  "species_name",
  "taxon",
  "delta_s",
  "se_delta_s",
  "ci_low",
  "ci_high",
  "pct_effect_s",
  "share_protected",
  "status"
)
missing_cols <- setdiff(required_cols, names(sp_est))
if (length(missing_cols) > 0) {
  stop(sprintf(
    "Missing required columns in %s: %s",
    sp_path,
    paste(missing_cols, collapse = ", ")
  ))
}
}


############################
### prepare data ###########
############################
{
taxon_labels <- c(
  amphibian = "Amphibians",
  bird = "Birds",
  mammal = "Mammals",
  reptile = "Reptiles"
)

taxon_colors <- c(
  Amphibians = "#E69F00",
  Birds      = "#56B4E9",
  Mammals    = "#009E73",
  Reptiles   = "#CC79A7"
)

sp_plot <- sp_est %>%
  filter(status == "estimated", is.finite(delta_s), is.finite(se_delta_s)) %>%
  mutate(
    taxon_label = factor(taxon_labels[taxon], levels = taxon_labels)
  )

if (nrow(sp_plot) == 0) {
  stop("No regression-estimated species have finite delta_s and se_delta_s.")
}

bottom_species <- sp_plot %>%
  arrange(delta_s, species_name) %>%
  slice_head(n = min(10L, nrow(sp_plot))) %>%
  mutate(panel = "10 largest threat reductions")

top_species <- sp_plot %>%
  arrange(desc(delta_s), species_name) %>%
  slice_head(n = min(10L, nrow(sp_plot))) %>%
  mutate(panel = "10 largest threat increases")

plot_data <- bind_rows(bottom_species, top_species) %>%
  distinct(species_id, .keep_all = TRUE) %>%
  mutate(
    photo_label = vapply(species_name, species_initials, character(1))
  )

cat("\nSelected species:\n")
plot_data %>%
  arrange(panel, delta_s) %>%
  select(panel, species_name, taxon, delta_s, se_delta_s) %>%
  mutate(
    delta_s = round(delta_s, 4),
    se_delta_s = round(se_delta_s, 4)
  ) %>%
  as.data.frame() %>%
  print(row.names = FALSE)
}


############################
### photos #################
############################
{
  photo_dir <- "data/store/species-pics/photos"
  thumb_dir <- "data/store/species-pics/thumbs"
  dir.create(photo_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(thumb_dir, recursive = TRUE, showWarnings = FALSE)

  photo_paths <- rep(NA_character_, nrow(plot_data))
  thumb_paths <- rep(NA_character_, nrow(plot_data))

cat("\n=== Fetching species photos ===\n")
for (i in seq_len(nrow(plot_data))) {
  sid <- plot_data$species_id[i]
  sname <- plot_data$species_name[i]
  cat(sprintf("  [%02d/%02d] %s\n", i, nrow(plot_data), sname))

  photo_paths[i] <- download_species_photo(sid, sname, photo_dir)

  if (!is.na(photo_paths[i]) && file.exists(photo_paths[i])) {
    thumb_name <- paste0(
      tools::file_path_sans_ext(basename(photo_paths[i])),
      ".png"
    )
    thumb_args <- thumb_overrides[[sname]]
    if (is.null(thumb_args)) {
      thumb_args <- list()
    }
    thumb_paths[i] <- create_round_thumbnail(
      src_path = photo_paths[i],
      thumb_path = file.path(thumb_dir, thumb_name),
      overwrite = TRUE,
      crop_frac = if (!is.null(thumb_args$crop_frac)) thumb_args$crop_frac else 0.76,
      focus_x = if (!is.null(thumb_args$focus_x)) thumb_args$focus_x else 0.5,
      focus_y = if (!is.null(thumb_args$focus_y)) thumb_args$focus_y else 0.5
    )
  }
}

plot_data$photo_path <- photo_paths
plot_data$thumb_path <- thumb_paths

cat(sprintf(
  "  Photo thumbnails available for %d of %d species\n",
  sum(!is.na(plot_data$thumb_path) & file.exists(plot_data$thumb_path)),
  nrow(plot_data)
))
}


############################
### figure #################
############################
{
legend_panel <- make_legend_panel(taxon_colors)

reduction_data <- plot_data %>%
  filter(panel == "10 largest threat reductions") %>%
  arrange(delta_s, species_name) %>%
  mutate(
    y_pos = rev(seq_len(n())),
    thumb_x = ifelse(seq_len(n()) %% 2 == 1, 0.34, 0.66)
  )

increase_data <- plot_data %>%
  filter(panel == "10 largest threat increases") %>%
  arrange(desc(delta_s), species_name) %>%
  mutate(
    y_pos = rev(seq_len(n())),
    thumb_x = ifelse(seq_len(n()) %% 2 == 1, 0.34, 0.66)
  )

x_min <- min(plot_data$ci_low, na.rm = TRUE)
x_max <- max(plot_data$ci_high, na.rm = TRUE)
x_pad <- max(0.01, 0.05 * (x_max - x_min))
x_limits <- c(min(x_min - x_pad, 0), max(x_max + x_pad, 0))

reduction_panel <- cowplot::plot_grid(
  make_effect_panel(reduction_data, taxon_colors, x_limits),
  make_photo_panel(reduction_data),
  nrow = 1,
  rel_widths = c(4.9, 1.95),
  align = "h",
  axis = "tb"
)

increase_panel <- cowplot::plot_grid(
  make_effect_panel(increase_data, taxon_colors, x_limits),
  make_photo_panel(increase_data),
  nrow = 1,
  rel_widths = c(4.9, 1.95),
  align = "h",
  axis = "tb"
)

fig <- cowplot::plot_grid(
  reduction_panel,
  increase_panel,
  legend_panel,
  ncol = 1,
  rel_heights = c(1, 1, 0.14)
)

figpath <- file.path(paths$figures_dir, "fig.byspecies.delta.top-bottom.est.jpg")
ggsave(figpath, fig, width = 11.2, height = 8.8, dpi = 300)
cat(sprintf("\nSaved figure: %s\n", figpath))
}

######################################
### coverage vs pct effectiveness ####
######################################
{
cat("\n=== Coverage vs effectiveness heatmap (% change) ===\n")

# Keep all finite effects. Values outside [-40, 40] enter explicit tail bins.
x_tail_lo <- -40
x_tail_hi <- 40
x_plot_lo <- -50
x_plot_hi <- 50
x_tail_center_lo <- -45
x_tail_center_hi <- 45
x_edge_eps <- 1e-6
cat(sprintf(
  "  Grouping finite pct_effect_s below %.0f%% and above %.0f%% into tail bins\n",
  x_tail_lo,
  x_tail_hi
))

sp_scatter <- sp_est %>%
  mutate(
    pct_protected = pmin(share_protected * 100, 99.999),
    pct_effect_plot = case_when(
      !is.finite(pct_effect_s) ~ NA_real_,
      pct_effect_s < x_tail_lo ~ x_tail_center_lo,
      pct_effect_s > x_tail_hi ~ x_tail_center_hi,
      TRUE ~ pmin(
        pmax(pct_effect_s, x_tail_lo + x_edge_eps),
        x_tail_hi - x_edge_eps
      )
    )
  )

n_unprotected <- sum(sp_scatter$share_protected == 0, na.rm = TRUE) +
  sum(is.na(sp_scatter$share_protected))
n_unassessed <- sum(!is.finite(sp_scatter$pct_effect_s) & sp_scatter$share_protected > 0,
                    na.rm = TRUE)
n_assessed <- sum(is.finite(sp_scatter$pct_effect_s))
cat(sprintf("  Assessed: %d, Unprotected: %d, Protected unassessed: %d\n",
            n_assessed, n_unprotected, n_unassessed))

# compute bin counts to set shared fill limits
bw_x <- 10
bw_y <- 20
heat_counts <- sp_scatter %>%
  filter(is.finite(pct_effect_plot), is.finite(pct_protected)) %>%
  mutate(
    xb = floor(pct_effect_plot / bw_x) * bw_x,
    yb = floor(pct_protected / bw_y) * bw_y
  ) %>%
  count(xb, yb)

# unassessed column: protected species with NA pct_effect_s, binned by pct_protected
unassessed_bins <- sp_scatter %>%
  filter(!is.finite(pct_effect_s), share_protected > 0, is.finite(pct_protected)) %>%
  mutate(yb = floor(pct_protected / bw_y) * bw_y) %>%
  count(yb)

fill_range <- range(c(heat_counts$n, unassessed_bins$n))

# place the NA column one bin-width left of the data, with a gap
na_x <- x_plot_lo - 2 * bw_x

pct_bin_labels <- function(x) {
  ifelse(
    x == x_tail_center_lo,
    paste0("<", x_tail_lo, "%"),
    ifelse(
      x == x_tail_center_hi,
      paste0(">", x_tail_hi, "%"),
      paste0(x - bw_x / 2, " \u2013 ", x + bw_x / 2, "%")
    )
  )
}

p_heat <- ggplot(sp_scatter, aes(x = pct_effect_plot, y = pct_protected)) +
  geom_bin2d(binwidth = c(bw_x, bw_y)) +
  # NA column: unassessed species by pct_protected bin
  geom_rect(
    data = unassessed_bins,
    aes(
      xmin = na_x,
      xmax = na_x + bw_x,
      ymin = yb,
      ymax = yb + bw_y,
      fill = n
    ),
    inherit.aes = FALSE
  ) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey45") +
  scale_x_continuous(
    breaks = c(
      na_x + bw_x / 2,
      seq(x_plot_lo + bw_x / 2, x_plot_hi - bw_x / 2, bw_x)
    ),
    labels = function(x) {
      ifelse(
        abs(x - (na_x + bw_x / 2)) < 1, "NA",
        pct_bin_labels(x)
      )
    },
    limits = c(na_x, x_plot_hi)
  ) +
  scale_fill_viridis_c(
    option = "inferno",
    name = "Number of species",
    trans = "log10",
    limits = fill_range
  ) +
  scale_y_continuous(
    breaks = seq(bw_y / 2, 100 - bw_y / 2, bw_y),
    labels = function(y) paste0(y - bw_y / 2, " \u2013 ", y + bw_y / 2, "%"),
    limits = c(0, 100)
  ) +
  labs(
    x = "Species-level PA effectiveness (% change from control)",
    y = "Share of species range protected (%)"
  ) +
  coord_fixed(ratio = bw_x / bw_y, clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barwidth = 15
  ))

heatpath <- file.path(paths$figures_dir, "fig.byspecies.coverage-vs-pct.heatmap.jpg")
ggsave(heatpath, p_heat, width = 10, height = 6, dpi = 300)
cat(sprintf("Saved figure: %s\n", heatpath))

sp_heat_estimated <- sp_est %>%
  filter(
    status == "estimated",
    is.finite(pct_effect_s),
    is.finite(share_protected)
  ) %>%
  transmute(
    coverage_pct = pmin(pmax(share_protected * 100, 0), 99.999),
    reduction_pct = -pct_effect_s,
    reduction_plot = case_when(
      reduction_pct < x_tail_lo ~ x_tail_center_lo,
      reduction_pct > x_tail_hi ~ x_tail_center_hi,
      TRUE ~ pmin(
        pmax(reduction_pct, x_tail_lo + x_edge_eps),
        x_tail_hi - x_edge_eps
      )
    )
  )

n_estimated <- nrow(sp_heat_estimated)
n_positive_reductions <- sum(sp_heat_estimated$reduction_pct > 0)
n_negative_reductions <- sum(sp_heat_estimated$reduction_pct < 0)
n_reduction_below <- sum(sp_heat_estimated$reduction_pct < x_tail_lo)
n_reduction_interior <- sum(
  sp_heat_estimated$reduction_pct >= x_tail_lo &
    sp_heat_estimated$reduction_pct <= x_tail_hi
)
n_reduction_above <- sum(sp_heat_estimated$reduction_pct > x_tail_hi)
cat(sprintf(
  paste0(
    "  Regression-estimated heatmap sample: %s species ",
    "(%s positive reductions; %s negative reductions)\n",
    "(<%.0f%%: %s; [%.0f%%, %.0f%%]: %s; >%.0f%%: %s)\n"
  ),
  format(n_estimated, big.mark = ","),
  format(n_positive_reductions, big.mark = ","),
  format(n_negative_reductions, big.mark = ","),
  x_tail_lo,
  format(n_reduction_below, big.mark = ","),
  x_tail_lo,
  x_tail_hi,
  format(n_reduction_interior, big.mark = ","),
  x_tail_hi,
  format(n_reduction_above, big.mark = ",")
))

heat_counts_estimated <- sp_heat_estimated %>%
  mutate(
    coverage_bin = floor(coverage_pct / bw_y) * bw_y,
    reduction_bin = floor(reduction_plot / bw_x) * bw_x
  ) %>%
  count(coverage_bin, reduction_bin, name = "n") %>%
  mutate(
    coverage_mid = coverage_bin + bw_y / 2,
    reduction_mid = reduction_bin + bw_x / 2,
    direction = if_else(
      reduction_mid > 0,
      "Positive reduction",
      "Negative reduction"
    ),
    direction = factor(
      direction,
      levels = c("Positive reduction", "Negative reduction")
    )
  )

if (sum(heat_counts_estimated$n) != n_estimated) {
  stop("Binned species counts do not sum to the regression-estimated sample.")
}

count_limits <- range(heat_counts_estimated$n)
count_breaks <- c(3, 10, 30, 100, 300)
count_breaks <- count_breaks[
  count_breaks >= count_limits[1] & count_breaks <= count_limits[2]
]

direction_colors <- c(
  "Positive reduction" = "#238B45",
  "Negative reduction" = "#CB181D"
)

p_heat_no_na <- ggplot(
  heat_counts_estimated,
  aes(x = coverage_mid, y = reduction_mid)
) +
  geom_tile(
    aes(fill = direction, alpha = n),
    width = bw_y,
    height = bw_x
  ) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey35") +
  scale_x_continuous(
    breaks = seq(bw_y / 2, 100 - bw_y / 2, bw_y),
    labels = function(x) paste0(x - bw_y / 2, " \u2013 ", x + bw_y / 2, "%"),
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(x_plot_lo + bw_x / 2, x_plot_hi - bw_x / 2, bw_x),
    labels = pct_bin_labels,
    limits = c(x_plot_lo, x_plot_hi),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = direction_colors,
    guide = "none"
  ) +
  scale_alpha_continuous(
    name = "Number of\nspecies\n(log scale)",
    trans = "log10",
    range = c(0.25, 1),
    limits = count_limits,
    breaks = count_breaks
  ) +
  labs(
    x = "Share of species range protected (%)",
    y = "Species-level threat reduction\n(% change from control)"
  ) +
  coord_fixed(ratio = bw_y / bw_x, clip = "off") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  guides(
    alpha = guide_legend(
      title.position = "top",
      position = "right",
      override.aes = list(fill = "grey25")
    )
  )

heatpath_no_na <- file.path(paths$figures_dir, "fig.byspecies.coverage-vs-pct.heatmap-noNA.jpg")
ggsave(
  heatpath_no_na,
  p_heat_no_na,
  width = 11,
  height = 10,
  units = "cm",
  dpi = 300,
  bg = "white"
)
# The manuscript stacks both panels at full width. Only the cell aspect ratio
# and output dimensions differ from the existing square heatmap.
p_heat_paper <- p_heat_no_na + coord_cartesian(clip = "off")
ggsave(file.path(paths$figures_dir,
                "fig.byspecies.coverage-vs-pct.heatmap-manuscript.jpg"),
       p_heat_paper, width = 15.8, height = 7.6, units = "cm", dpi = 320,
       bg = "white")
cat(sprintf("Saved figure: %s\n", heatpath_no_na))
}

######################################
### coverage vs threat reduction #####
######################################
# Species panel of manuscript Figure 4 (country panel lives in
# pathreat.analysis.bycountry.fig.R). Y-limits match the country panel so
# the two halves of Figure 4 are directly comparable.
{
cat("\n=== Coverage vs threat-reduction (species panel of Fig 4) ===\n")

y_lim_fig4 <- c(-0.075, 0.18)

tax_pal <- c(
  "Amphibian" = "#1b9e77",
  "Bird"      = "#7570b3",
  "Mammal"    = "#d95f02",
  "Reptile"   = "#e7298a"
)

sp_reduction <- sp_est %>%
  dplyr::filter(
    status == "estimated",
    !is.na(delta_s), !is.na(share_protected)
  ) %>%
  dplyr::mutate(
    share_pct = share_protected * 100,
    reduction = -delta_s,
    taxon     = tools::toTitleCase(taxon)
  )

tax_means <- sp_reduction %>%
  dplyr::group_by(taxon) %>%
  dplyr::summarise(
    mean_x = mean(share_pct, na.rm = TRUE),
    mean_y = mean(reduction, na.rm = TRUE),
    .groups = "drop"
  )

p_red <- ggplot(sp_reduction, aes(x = share_pct, y = reduction)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_point(aes(fill = taxon),
             shape = 21, size = 1.4, colour = "white",
             stroke = 0.1, alpha = 0.55) +
  geom_vline(data = tax_means,
             aes(xintercept = mean_x),
             colour = "grey25", linewidth = 0.3, linetype = "dashed") +
  geom_hline(data = tax_means,
             aes(yintercept = mean_y),
             colour = "grey25", linewidth = 0.3, linetype = "dashed") +
  scale_fill_manual(values = tax_pal, guide = "none") +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     limits = c(0, 100),
                     breaks = seq(0, 100, 25)) +
  scale_y_continuous(limits = y_lim_fig4,
                     breaks = scales::pretty_breaks(n = 5)) +
  facet_wrap(~ taxon, ncol = 2) +
  labs(x = "Share of species range protected",
       y = "Threat reduction from PAs (0–1 scale)") +
  theme_minimal(base_size = 10) +
  theme(strip.text        = element_text(face = "bold"),
        panel.grid.minor  = element_blank(),
        legend.position   = "none",
        plot.margin       = margin(t = 4, r = 5, b = 4, l = 5))

redpath <- file.path(paths$figures_dir, "fig.byspecies.coverage-vs-reduction.jpg")
ggsave(redpath, p_red, width = 11, height = 10, units = "cm", dpi = 300)
cat(sprintf("Saved figure: %s\n", redpath))
}

cat("\n=== Done ===\n")
