##################################################
### pathreat.data.species.R ######################
### Species metadata reference table #############
### All four vertebrate taxa (threatened only) ###
##################################################

library(sf)
library(readxl)
library(foreign)


############################
### helper functions #######
############################
taxon_lookup <- c(
  Amphibia = "amphibian",
  Aves = "bird",
  Mammalia = "mammal",
  Reptilia = "reptile"
)

build_key <- function(taxon, value) {
  value_chr <- as.character(value)
  out <- paste(taxon, value_chr, sep = "::")
  out[is.na(value) | is.na(value_chr) | value_chr == ""] <- NA_character_
  out
}

normalize_priority_sheet <- function(sheet, taxon, species_col = "Species",
                                     id_col = "RL.ID", genus_col = NULL) {
  d <- read_excel(file.path("data/external/EDGE", "2024_EDGE_species_external.xlsx"),
                  sheet = sheet)

  species_name <- d[[species_col]]
  if (!is.null(genus_col)) {
    species_name <- paste(d[[genus_col]], d[[species_col]])
  }

  id_no <- rep(NA_integer_, nrow(d))
  if (id_col %in% names(d)) {
    id_no <- suppressWarnings(as.integer(d[[id_col]]))
  }

  rank_col <- grep("^EDGE\\.rank$|^EDGE\\.Rank$", names(d), value = TRUE)
  if (length(rank_col) == 0) {
    edge_rank <- rep(NA_real_, nrow(d))
  } else {
    edge_rank <- d[[rank_col[1]]]
  }

  out <- data.frame(
    taxon = taxon,
    id_no = id_no,
    species_name = species_name,
    tbl_median = d[["TBL.median"]],
    ed_median = d[["ED.median"]],
    edge_median = d[["EDGE.median"]],
    edge_species = TRUE,
    edge_rank = edge_rank,
    edge_tier = d[["Tier"]],
    stringsAsFactors = FALSE
  )

  out <- out[!duplicated(out[, c("taxon", "id_no", "species_name")]), ]
  rownames(out) <- NULL
  out
}

fill_core_from_source <- function(target, source, key_col, source_name, match_name) {
  source <- source[!is.na(source[[key_col]]) & !duplicated(source[[key_col]]), ]
  idx <- match(target[[key_col]], source[[key_col]])
  hit <- !is.na(idx) & is.na(target$edge_median)

  if (!any(hit)) {
    return(target)
  }

  target$tbl_median[hit] <- source$tbl_median[idx[hit]]
  target$ed_median[hit] <- source$ed_median[idx[hit]]
  target$edge_median[hit] <- source$edge_median[idx[hit]]
  target$edge_species[hit] <- source$edge_species[idx[hit]]
  target$edge_source[hit] <- source_name
  target$edge_match_key[hit] <- match_name

  target
}

fill_optional_from_source <- function(target, source, key_col, cols) {
  source <- source[!is.na(source[[key_col]]) & !duplicated(source[[key_col]]), ]
  idx <- match(target[[key_col]], source[[key_col]])

  for (col in cols) {
    if (!col %in% names(source)) {
      next
    }
    replace <- !is.na(idx) & is.na(target[[col]]) & !is.na(source[[col]][idx])
    target[[col]][replace] <- source[[col]][idx[replace]]
  }

  target
}


############################
### species metadata #######
############################
{
  cat("=== Building species metadata for all taxa ===\n")

  ext_dir <- "data/external"

  # Mammals, amphibians, and reptiles from IUCN shapefile DBFs
  shp_taxa <- list(
    mammal = list(
      dir = "MAMMALS",
      parts = c("MAMMALS_PART1.dbf", "MAMMALS_PART2.dbf")
    ),
    amphibian = list(
      dir = "AMPHIBIANS",
      parts = c("AMPHIBIANS_PART1.dbf", "AMPHIBIANS_PART2.dbf")
    ),
    reptile = list(
      dir = "REPTILES",
      parts = c("REPTILES_PART1.dbf", "REPTILES_PART2.dbf")
    )
  )

  sp_meta_list <- list()
  for (tx in names(shp_taxa)) {
    cfg <- shp_taxa[[tx]]
    parts <- lapply(cfg$parts, function(p) {
      foreign::read.dbf(file.path(ext_dir, cfg$dir, p), as.is = TRUE)[,
        c("id_no", "sci_name", "category")]
    })
    meta <- do.call(rbind, parts)
    meta <- meta[!duplicated(meta$id_no), ]
    meta$taxon <- tx
    meta$species_id <- meta$id_no
    meta$sisid <- NA_integer_
    names(meta)[names(meta) == "sci_name"] <- "species_name"
    sp_meta_list[[tx]] <- meta
    cat(sprintf("  %-12s %d species\n", tx, nrow(meta)))
  }

  # Birds from GPKG attributes plus Red List checklist
  bird_attrs <- st_drop_geometry(
    st_read(file.path(ext_dir, "BIRDS", "BOTW_2025.gpkg"),
            layer = "all_species", quiet = TRUE)
  )
  bird_attrs <- bird_attrs[!duplicated(bird_attrs$sisid),
                           c("sisid", "sci_name")]

  checklist <- read_excel(
    file.path(
      ext_dir,
      "BIRDS",
      "Handbook of the Birds of the World and BirdLife International Digital Checklist of the Birds of the World_Version_10.xlsx"
    ),
    sheet = 1,
    skip = 3
  )
  checklist <- checklist[!is.na(checklist$SISRecID), ]
  checklist$category <- checklist[["2025 IUCN Red List category"]]
  checklist$SISRecID <- as.integer(checklist$SISRecID)
  checklist$category <- sub("^CR \\(PE\\)$", "CR", checklist$category)
  checklist$category <- sub("^CR \\(PEW\\)$", "CR", checklist$category)

  bird_meta <- merge(
    bird_attrs,
    checklist[, c("SISRecID", "category")],
    by.x = "sisid",
    by.y = "SISRecID",
    all.x = TRUE
  )
  bird_meta$taxon <- "bird"
  bird_meta$id_no <- NA_integer_
  bird_meta$species_id <- bird_meta$sisid
  names(bird_meta)[names(bird_meta) == "sci_name"] <- "species_name"
  sp_meta_list[["bird"]] <- bird_meta
  cat(sprintf("  %-12s %d species\n", "bird", nrow(bird_meta)))

  rm(bird_attrs, checklist, bird_meta)

  sp_meta <- do.call(rbind, sp_meta_list)
  rownames(sp_meta) <- NULL
  rm(sp_meta_list)

  cat(sprintf("  Total metadata: %d species\n", nrow(sp_meta)))
}


############################
### filter species #########
############################
{
  sp_meta <- sp_meta[sp_meta$category %in% c("CR", "EN", "VU"), ]
  sp_meta$raster_file <- paste0("sp_", sp_meta$taxon, "_", sp_meta$species_id, ".tif")
  rownames(sp_meta) <- NULL

  cat(sprintf("  Threatened only: %d species\n", nrow(sp_meta)))
  cat("  By taxon:\n")
  for (tx in c("mammal", "amphibian", "reptile", "bird")) {
    cat(sprintf("    %-12s %d\n", tx, sum(sp_meta$taxon == tx)))
  }
}


############################
### merge EDGE data ########
############################
{
  cat("\n=== Merging EDGE metadata ===\n")

  sp_meta$id_no_chr <- as.character(sp_meta$id_no)
  sp_meta$edge_rl_key <- build_key(sp_meta$taxon, sp_meta$id_no_chr)
  sp_meta$edge_name_key <- build_key(sp_meta$taxon, sp_meta$species_name)

  # Primary source: full 2024.1 tetrapod file
  edge_full <- read_excel(
    file.path("data/external/EDGE", "Full_EDGE_2024.1_data_tetrapods.xlsx"),
    sheet = 1
  )
  edge_full$taxon <- unname(taxon_lookup[edge_full$class])
  edge_full$tbl_median <- edge_full$tbl
  edge_full$ed_median <- edge_full$ed
  edge_full$edge_median <- edge_full$edge
  edge_full$edge_species <- edge_full$edgeSpecies == "YES"
  edge_full$edge_rl_key <- build_key(edge_full$taxon, edge_full$rlTaxonId)
  edge_full$edge_name_key <- build_key(edge_full$taxon, edge_full$species)
  edge_full <- edge_full[, c(
    "edge_rl_key", "edge_name_key",
    "tbl_median", "ed_median", "edge_median", "edge_species"
  )]

  # 2024 priority file for rank/tier and fallback matches not found in full 2024.1
  edge_priority <- do.call(
    rbind,
    list(
      normalize_priority_sheet("Amphibans", "amphibian"),
      normalize_priority_sheet("Birds", "bird", id_col = "__none__"),
      normalize_priority_sheet("Mammals", "mammal"),
      normalize_priority_sheet("Lepidosaurs", "reptile"),
      normalize_priority_sheet("Crocodylians", "reptile", species_col = "Species", genus_col = "Genus"),
      normalize_priority_sheet("Testudines", "reptile", species_col = "Species", genus_col = "Genus")
    )
  )
  edge_priority$edge_rl_key <- build_key(edge_priority$taxon, edge_priority$id_no)
  edge_priority$edge_name_key <- build_key(edge_priority$taxon, edge_priority$species_name)
  edge_priority <- edge_priority[, c(
    "edge_rl_key", "edge_name_key",
    "tbl_median", "ed_median", "edge_median", "edge_species",
    "edge_rank", "edge_tier"
  )]

  # 2023 mammals-only file for uncertainty ranges and mammal fallback matches
  edge_mammal_2023 <- read_excel(
    file.path("data/external/EDGE", "journal.pbio.3001991.s003.xlsx"),
    sheet = "All EDGE2 Scores"
  )
  edge_mammal_2023 <- data.frame(
    taxon = "mammal",
    id_no = suppressWarnings(as.integer(edge_mammal_2023[["RL.ID"]])),
    species_name = edge_mammal_2023[["Species"]],
    tbl_median = edge_mammal_2023[["TBL.median"]],
    ed_median = edge_mammal_2023[["ED.median"]],
    edge_median = edge_mammal_2023[["EDGE.median"]],
    edge_species = edge_mammal_2023[["EDGE.species"]] == "YES",
    tbl_iqr_low = edge_mammal_2023[["TBL.IQR.low"]],
    tbl_iqr_high = edge_mammal_2023[["TBL.IQR.high"]],
    ed_iqr_low = edge_mammal_2023[["ED.IQR.low"]],
    ed_iqr_high = edge_mammal_2023[["ED.IQR.high"]],
    edge_iqr_low = edge_mammal_2023[["EDGE.IQR.low"]],
    edge_iqr_high = edge_mammal_2023[["EDGE.IQR.high"]],
    edge_no_above_median = edge_mammal_2023[["no.above.median"]],
    stringsAsFactors = FALSE
  )
  edge_mammal_2023$edge_rl_key <- build_key(edge_mammal_2023$taxon, edge_mammal_2023$id_no)
  edge_mammal_2023$edge_name_key <- build_key(edge_mammal_2023$taxon, edge_mammal_2023$species_name)

  sp_meta$tbl_median <- NA_real_
  sp_meta$ed_median <- NA_real_
  sp_meta$edge_median <- NA_real_
  sp_meta$edge_species <- NA
  sp_meta$edge_rank <- NA_real_
  sp_meta$edge_tier <- NA_real_
  sp_meta$tbl_iqr_low <- NA_real_
  sp_meta$tbl_iqr_high <- NA_real_
  sp_meta$ed_iqr_low <- NA_real_
  sp_meta$ed_iqr_high <- NA_real_
  sp_meta$edge_iqr_low <- NA_real_
  sp_meta$edge_iqr_high <- NA_real_
  sp_meta$edge_no_above_median <- NA_real_
  sp_meta$edge_source <- NA_character_
  sp_meta$edge_match_key <- NA_character_

  sp_meta <- fill_core_from_source(sp_meta, edge_full, "edge_rl_key", "full_2024_1", "rl_id")
  sp_meta <- fill_core_from_source(sp_meta, edge_full, "edge_name_key", "full_2024_1", "species_name")
  sp_meta <- fill_core_from_source(sp_meta, edge_priority, "edge_rl_key", "priority_2024", "rl_id")
  sp_meta <- fill_core_from_source(sp_meta, edge_priority, "edge_name_key", "priority_2024", "species_name")
  sp_meta <- fill_core_from_source(sp_meta, edge_mammal_2023, "edge_rl_key", "mammals_2023", "rl_id")
  sp_meta <- fill_core_from_source(sp_meta, edge_mammal_2023, "edge_name_key", "mammals_2023", "species_name")

  sp_meta <- fill_optional_from_source(
    sp_meta,
    edge_priority,
    "edge_rl_key",
    c("edge_rank", "edge_tier")
  )
  sp_meta <- fill_optional_from_source(
    sp_meta,
    edge_priority,
    "edge_name_key",
    c("edge_rank", "edge_tier")
  )
  sp_meta <- fill_optional_from_source(
    sp_meta,
    edge_mammal_2023,
    "edge_rl_key",
    c("tbl_iqr_low", "tbl_iqr_high", "ed_iqr_low", "ed_iqr_high",
      "edge_iqr_low", "edge_iqr_high", "edge_no_above_median")
  )
  sp_meta <- fill_optional_from_source(
    sp_meta,
    edge_mammal_2023,
    "edge_name_key",
    c("tbl_iqr_low", "tbl_iqr_high", "ed_iqr_low", "ed_iqr_high",
      "edge_iqr_low", "edge_iqr_high", "edge_no_above_median")
  )

  matched_edge <- !is.na(sp_meta$edge_median)
  cat(sprintf("  EDGE matched: %d / %d species (%.1f%%)\n",
              sum(matched_edge), nrow(sp_meta), 100 * mean(matched_edge)))
  cat("  By taxon:\n")
  for (tx in c("mammal", "amphibian", "reptile", "bird")) {
    keep <- sp_meta$taxon == tx
    cat(sprintf("    %-12s %4d / %4d (%.1f%%)\n",
                tx,
                sum(matched_edge[keep]),
                sum(keep),
                100 * mean(matched_edge[keep])))
  }

  cat("  Core EDGE source:\n")
  source_counts <- sort(table(sp_meta$edge_source), decreasing = TRUE)
  for (nm in names(source_counts)) {
    cat(sprintf("    %-14s %d\n", nm, source_counts[[nm]]))
  }

  mammal_uncertainty <- sp_meta$taxon == "mammal" & !is.na(sp_meta$edge_iqr_low)
  cat(sprintf("  Mammals with 2023 uncertainty ranges: %d / %d\n",
              sum(mammal_uncertainty),
              sum(sp_meta$taxon == "mammal")))

  if (any(!matched_edge)) {
    missing_sp <- sp_meta[!matched_edge, c("taxon", "species_name")]
    cat(sprintf("  WARNING: %d species still unmatched after EDGE merge\n",
                nrow(missing_sp)))
    print(utils::head(missing_sp, 20), row.names = FALSE)
  }

  sp_meta$edge_species <- as.logical(sp_meta$edge_species)

  sp_meta$id_no_chr <- NULL
  sp_meta$edge_rl_key <- NULL
  sp_meta$edge_name_key <- NULL
}


############################
### reorder + save #########
############################
{
  sp_meta <- sp_meta[, c(
    "species_id", "id_no", "sisid", "species_name", "taxon", "category",
    "raster_file",
    "tbl_median", "ed_median", "edge_median", "edge_species",
    "edge_rank", "edge_tier",
    "tbl_iqr_low", "tbl_iqr_high",
    "ed_iqr_low", "ed_iqr_high",
    "edge_iqr_low", "edge_iqr_high",
    "edge_no_above_median",
    "edge_source", "edge_match_key"
  )]
  rownames(sp_meta) <- NULL

  outfile <- "data/store/pathreat.data.species.Rds"
  saveRDS(sp_meta, outfile)
  cat(sprintf("\nSaved: %s (%.1f KB)\n",
              outfile, file.size(outfile) / 1e3))
}


######have a look #######
d<-readRDS("data/store/pathreat.data.species.Rds")
nrow(d)
head(d)
table(d$category,useNA="always")
table(is.na(d$species_id))
table(is.na(d$edge_median))
