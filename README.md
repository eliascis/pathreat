# PA-Threat

Analysis code and derived estimates for

> Vincent, C., Cisneros, E., and Rondinini, C. (2026). *Improving existing protected areas could double their threat reduction impact.* Working paper under review.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23108161.svg)](https://doi.org/10.5281/zenodo.23108161)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Claire Vincent (Sapienza University of Rome), Elías Cisneros (The University of Texas at Dallas), Carlo Rondinini (Sapienza University of Rome).

## What the study does

The paper evaluates terrestrial protected areas (PAs) designated between 2001 and 2020 against 14 mapped anthropogenic threats on a global 1-km grid. Grid cells inside PAs are matched 1:1 (Mahalanobis distance on eight biophysical and accessibility covariates) to unprotected cells within country-biome strata, with large strata subdivided by subnational boundaries. Treatment effects are estimated with country and biome fixed effects and country-clustered standard errors on a composite threat index, an any-threat indicator and each of the 14 indicators, and are broken down by country, biome, biodiversity hotspot and PA type. Protected-area-level regressions feed a stochastic frontier model that separates achieved threat reduction from estimated underperformance. Species-range fixed-effects regressions estimate effects for threatened terrestrial vertebrates (mammals, birds, amphibians, reptiles), summarised per taxon with a bootstrap.

The balanced matched sample covers 5,711,528 matched pairs (11,423,056 observations) in 107 countries and 22,090 PAs.

## What this repository contains

| Path | Content |
|---|---|
| `code/` | R scripts (and two Python raster helpers) for data preparation, matching, estimation, and the manuscript figures and tables. `code/tests/` holds a data-free unit test of the PA regression helper. |
| `results/` | Canonical estimate tables as `.Rds` with CSV twins, plus `MANIFEST.csv` (file, rows, columns, md5). |
| `data/` | Directory skeleton only. `data/README.md` lists every input, its source, version and licence. A few small derived files are tracked under `data/store/`: variable metadata, per-PA estimates, PA-by-country membership, and per-figure source data. |
| `renv.lock`, `SESSIONINFO.txt` | Exact R package versions and the geospatial library versions used. |
| `requirements.txt` | Python dependency for the raster helpers. |

Not included: the raw input datasets (about 100 GB; most are not redistributable), the intermediate pixel-level datasets (about 100 GB), the matched sample (1.6 GB, available on request), and the manuscript.

## Pipeline

All scripts are run from the repository root with `Rscript code/<script>` and source `code/pathreat.analysis.config.R` for paths, variable lists, labels and shared estimation helpers. Scripts write to `results/`, `pub/figures/`, `pub/tables/` and `data/store/`; the output directories are created on demand.

### 1. Data preparation (needs the raw inputs in `data/external/`)

| Script | Purpose |
|---|---|
| `pathreat.data.prepdata.R` | Reads the compiled 1-km pixel table (`data/store/pre-process/global_prematch_data.Rds`), basic cleaning |
| `pathreat.data.rasterbase.R` | 1-km Mollweide grid template defining the project geometry |
| `pathreat.data.deforestation.R` | Yearly forest-loss counts to cumulative 2001-2020 area |
| `pathreat.data.hotspots.R` | Rasterises the 36 biodiversity hotspots onto the grid |
| `pathreat.data.land-vulnerability-download.py`, `-resample.py`, `pathreat.data.land-vulnerability.R` | Download, resample and extract the land-cover vulnerability 2050 layer |
| `pathreat.data.wdpa.union.R`, `wdpa.R`, `wdpa.points.R` | Unify public and restricted PA layers, rasterise polygons (oldest PA wins), flag point-derived PAs |
| `pathreat.data.countries.R`, `admin1.R` | GADM level-0 and Natural Earth admin-1 boundaries rasterised onto the grid |
| `pathreat.data.wdpa-buffer.R` | 5-km buffer indicator around study-period PAs |
| `pathreat.data.merge.R` | Merges all pixel-level inputs; derives `treat` and `ever_pa`. Needs `R_MAX_VSIZE=128Gb` |
| `pathreat.data.threat-indices.R` | Presence thresholds, composite threat index and any-threat indicator |
| `pathreat.data.lu.pressure.calc.R` | Three land-use pressure indices |
| `pathreat.data.species.R`, `species.raster.R`, `species-count.R`, `species-pa-cover.R` | Threatened-species reference table, per-species range rasters, pixel-level species counts and score sums, species-PA overlap |
| `pathreat.data.merge.unmatched.R` | Adds species counts, indices and pressure; applies country and PA exclusions; writes the analysis-ready unmatched sample |

### 2. Matching and assembly

| Script | Purpose |
|---|---|
| `pathreat.data.outcome-sample-filters.R` | Outcome-specific eligibility rules |
| `pathreat.data.matching.R` | 1:1 Mahalanobis matching within country-biome-region strata (668 effective strata), with a manifest of inputs, covariates and outputs |
| `pathreat.data.merge.matched.R` | Assembles the validated matched pairs, applies the country balance filter (average absolute SMD at most 0.10), writes `data/store/pathreat.data.merge.matched.fst` |
| `pathreat.randomrob.data.R` | Matching-stage robustness sweep on 1% subsamples (nine matching variants, five seeds) |

### 3. Estimation

| Script | Output |
|---|---|
| `pathreat.analysis.covbalance.R` | Covariate balance and VIF diagnostics |
| `pathreat.analysis.global.est.R` | Global ATT for all 16 outcomes with country and biome FE; standardised levels with a country-pairs bootstrap |
| `pathreat.analysis.global-threat.R` | Biodiversity-weighted variants of the global composite ATT |
| `pathreat.analysis.threat-weights.R` | Country-specific prevalence weights behind the composite index |
| `pathreat.analysis.PA-type.est.R`, `bycountry.est.R`, `bybiome.est.R`, `hotspots.R` | Subgroup ATTs by PA type, country, biome and biodiversity hotspot |
| `pathreat.analysis.robustness.est.R`, `sensitivity.R`, `randomrob.analysis.R` | Specification robustness, E-values and Rosenbaum bounds, matching-variant robustness |
| `pathreat.analysis.byPA.est.R` | One composite-outcome regression per PA (22,090 PAs) |
| `pathreat.analysis.sfa.est.R`, `sfa.aggregates.R` | PA-level stochastic cost frontier and area-weighted aggregates |
| `pathreat.analysis.byspecies.R`, `bytaxa.R` | Species-range FE regressions and taxon medians with 2,000-draw bootstrap |
| `pathreat.analysis.statistics.R` | Audited headline statistics quoted in the paper |
| `pathreat.analysis.export-csv.R` | CSV twins of the `.Rds` tables and `results/MANIFEST.csv` |

### 4. Figures and tables

Every figure script writes the data frame it plots to `data/store/<figure root>.csv` when that data is otherwise computed only in memory (`save_figure_data()` in the config). The table below maps each asset in the manuscript to its producing script and the data behind it.

| Manuscript asset | Script | Data |
|---|---|---|
| Fig. 1 `fig.global.violins.country.combined.jpg` | `pathreat.analysis.global.violins.fig.R` | `results/pathreat.global.est.Rds`, `results/pathreat.global.fe-adjusted-levels.est.Rds`, `data/store/fig.global.violins.country.combined.csv` |
| Fig. 2 `fig.delta_map.jpg` | `pathreat.analysis.global.delta_map.R` | `data/store/fig.delta_map.csv` (25-km grid means of pair differences) |
| Fig. 3 `fig.bycountry.coverage-vs-reduction.jpg` | `pathreat.analysis.bycountry.fig.R` | `results/pathreat.bycountry.est.Rds`, `data/store/fig.bycountry.coverage-vs-reduction.csv` |
| Fig. 4 `fig.sfa.closable-gap.jpg` | `pathreat.analysis.sfa.fig.R` | `results/pathreat.sfa.est.Rds`, `results/pathreat.sfa.country.est.Rds`, `results/pathreat.sfa.models.est.Rds`, `results/pathreat.global.est.Rds` |
| Fig. 5A `fig.bytaxa.effect-distributions.manuscript.jpg` | `pathreat.pres.bytaxa-effect-distributions.R` | `results/pathreat.byspecies.est.Rds`, `results/pathreat.bytaxa.est.Rds` |
| Fig. 5B `fig.byspecies.coverage-vs-pct.heatmap-manuscript.jpg` | `pathreat.analysis.byspecies.fig.R` | `results/pathreat.byspecies.est.Rds` |
| Extended Data: `fig.bybiome.forest.est.jpg` | `pathreat.analysis.bybiome.fig.R` | `results/pathreat.bybiome.est.Rds` |
| Extended Data: `fig.hotspots.threat_composite.coef.est.jpg` | `pathreat.analysis.hotspots.R` | `results/pathreat.analysis.hotspots.Rds` |
| Extended Data: `fig.PA-type.combined.est.jpg`, `fig.PA-type.threat_composite.coef.est.jpg` | `pathreat.analysis.PA-type.fig.R` | `results/pathreat.PA-type.est.Rds`, `results/pathreat.PA-type.group_means.est.Rds` |
| Extended Data: `fig.hetero.{biodiversity,pressure.est,pressure.access,tc_threat_composite}.est.jpg` | `pathreat.analysis.global.hetero.fig.R` | `data/store/pathreat.analysis.byPA.est.Rds`, `data/store/fig.hetero.*.csv` |
| Extended Data: `fig.sfa.levels.frontier-curve.uhet.{country,pa}.jpg`, `fig.sfa.efficiency-gap.density.jpg`, `fig.sfa.country_lollipop-level.jpg` | `pathreat.analysis.sfa.fig.R` | SFA results as for Fig. 4 |
| Supplementary: `fig.sample.pa_by_year.jpg`, `fig.sample.pixels_by_year.jpg` | `pathreat.analysis.stats.fig.pa-by-year.R` | `data/store/fig.sample.*.csv` |
| Supplementary: `fig.threat_composite_weights.heatmap.jpg` | `pathreat.analysis.threat-weights.R` | `results/pathreat.threat_composite_weights.Rds` |
| Supplementary: `fig.match.balance.jpg` | `pathreat.analysis.covbalance.R` | `results/pathreat.covbalance.csv` |
| Supplementary: `fig.randomrob.coef.est.jpg` | `pathreat.randomrob.analysis.R` | `results/pathreat.randomrob.est.Rds` |
| Supplementary: `fig.bytaxa.delta.est.jpg` | `pathreat.analysis.bytaxa.fig.R` | `results/pathreat.bytaxa.est.Rds`, `results/pathreat.global-threat.est.Rds`, `results/p_combined_cover.Rds` |
| Supplementary: `fig.PA-species-{cover,effectiveness}-{mammals,amphibians,birds,reptiles}.jpg` | `pathreat.analysis.species-pa.fig.R` | species-PA overlap and matched data (not distributed), `results/pathreat.byspecies.est.Rds` |
| Supplementary: `fig.pres.byspecies.method.jpg` | `pathreat.analysis.byspecies.method.fig.R` | schematic, no data |
| Tables `tab.ATT_by_threat_table` | `pathreat.analysis.global.est.R` | `results/pathreat.global.est.Rds` |
| `tab.sfa.frontier` | `pathreat.analysis.sfa.est.R` | `results/pathreat.sfa.models.est.Rds` |
| `tab.PA_category_summary`, `tab.PA_char_table` | `pathreat.analysis.PA-type.est.R` | `results/pathreat.PA-type.est.Rds`; category summary needs the unmatched sample |
| `tab.matching_diagnostics` | `pathreat.analysis.covbalance.R` | `results/pathreat.matching_diagnostics.csv` |
| `tab.rob.controls`, `tab.sensitivity`, `tab.randomrob` | `pathreat.analysis.robustness.est.R`, `sensitivity.R`, `randomrob.analysis.R` | corresponding `results/*.Rds` |
| `tab.biome_protection`, `tab.ATT_by_biome`, `tab.ATT_by_country` | `pathreat.analysis.bybiome.fig.R`, `bybiome.est.R`, `bycountry.est.R` | `results/pathreat.bybiome.est.Rds`, `results/pathreat.bycountry.est.Rds`; protection shares need the unmatched sample |
| `tab.Dataset_Table` | hand-written in the manuscript | - |

## Reproducing the results

Three levels of reproduction are possible, depending on which inputs you have.

**From this repository alone** (estimate tables and the shipped `data/store/` files): `pathreat.analysis.sfa.fig.R`, `sfa.aggregates.R`, `PA-type.fig.R`, `byspecies.fig.R` (the ranked-species panel downloads species photos from Wikimedia), `byspecies.method.fig.R`, `bytaxa.fig.R`, `pres.bytaxa-effect-distributions.R`, `statistics.R`, `export-csv.R`, and `code/tests/pathreat.analysis.byPA.R`.

**With the matched sample** (`data/store/pathreat.data.merge.matched.fst`, 1.6 GB, available from the corresponding author subject to the input providers' terms): all estimation scripts except `covbalance.R`, `threat-weights.R`, `PA-type.est.R`, `robustness.est.R`, `byspecies.R` and `randomrob.*`, plus `global.violins.fig.R`, `global.hetero.fig.R`, `global.delta_map.R`, `stats.fig.pa-by-year.R` and `bytaxa.R`. Loading the matched sample needs several GB of RAM.

**With the unmatched analysis sample** (`data/store/pathreat.data.merge.unmatched.fst`, about 127 million pixels, not distributed; rebuilt by the data-preparation stage): `pathreat.analysis.bycountry.fig.R` (Fig. 3 country coverage and PA area), `bybiome.fig.R` (biome protection table), `species-pa.fig.R`, `PA-type.est.R` (category summary), `covbalance.R`, `robustness.est.R` and `randomrob.data.R`.

**With the raw inputs** listed in `data/README.md`: the full pipeline. Data preparation needs about 200 GB of disk and 128 GB of RAM for `pathreat.data.merge.R`; species rasterisation and matching take hours on 8 cores.

Set the number of worker processes with `PATHREAT_CORES` (default: all but one core, at most 8). The matching and species scripts also read `PATHREAT_MATCHING_CORES` and `PATHREAT_BYSPECIES_CORES`. Three scripts access the network: `pathreat.data.countries.R` downloads GADM boundaries, `pathreat.data.land-vulnerability-download.py` downloads the land-cover vulnerability raster, and `pathreat.analysis.byspecies.fig.R` fetches species photos. `pathreat.data.admin1.R` reads Natural Earth boundaries bundled in the `rnaturalearthhires` data package, which is installed once.

## Software environment

- R 4.6.1. Package versions are pinned in `renv.lock`; restore them with

  ```r
  install.packages("renv")
  options(pkgType = "binary")   # macOS / Windows: avoid compiling sf and terra
  renv::restore()
  ```

  CRAN keeps binaries only for the current version of each package, so locked versions that have since been superseded are built from source. On macOS this needs the Xcode command-line tools and the GNU Fortran compiler from [mac.r-project.org/tools](https://mac.r-project.org/tools/) (RcppArmadillo, and through it `sfaR` and `MatchIt`, need Fortran). If compilers are not available, restore everything else at the locked versions and install current binaries of the superseded packages. This route was tested on macOS (arm64) on 2026-10-02: the unit test and the results-only scripts ran, and the global estimates from the matched sample reproduced the archived values to 1e-8.

  ```r
  superseded <- c("MatchIt", "RcppArmadillo", "curl", "ggtext", "gridtext", "hpa",
                  "httr", "mnorm", "rnaturalearth", "sfaR", "texreg")
  renv::restore(exclude = superseded)
  install.packages(superseded, type = "binary")
  ```

- `rnaturalearthhires` is distributed through rOpenSci's r-universe rather than CRAN; the lockfile records that repository, and `renv::restore()` resolves it. Manual install: `install.packages("rnaturalearthhires", repos = c("https://ropensci.r-universe.dev", "https://cloud.r-project.org"))`.
- System libraries for `sf` and `terra`: GDAL 3.8.5, GEOS 3.14.1, PROJ 9.5.1 were used (see `SESSIONINFO.txt`). `rsvg` needs librsvg.
- Python 3 with `requests` (see `requirements.txt`) and the GDAL command-line tools (`gdalbuildvrt`, `gdal_translate`, `gdalwarp`) for the two land-vulnerability helpers.
- Main R packages: `dplyr`, `tidyr`, `fixest`, `Matching`, `MatchIt`, `sfaR`, `sf`, `terra`, `fst`, `ggplot2`, `patchwork`, `ggrepel`, `rnaturalearth`, `readxl`, `texreg`.

## Results directory

| File (`.Rds` and `.csv`) | Unit of observation |
|---|---|
| `pathreat.global.est`, `pathreat.global.group_means.est`, `pathreat.global.fe-adjusted-levels.est` | Outcome (16): ATT, inference, group means, standardised levels |
| `pathreat.global-threat.est` | Weighting scheme (5): biodiversity-weighted composite ATT |
| `pathreat.bycountry.est` | Country by outcome |
| `pathreat.bybiome.est` | Biome by outcome |
| `pathreat.hotspots.est`, `pathreat.analysis.hotspots` | Biodiversity hotspot (31 + non-hotspot), composite ATT |
| `pathreat.PA-type.est`, `pathreat.PA-type.group_means.est` | PA type (IUCN class, size class) by outcome |
| `pathreat.analysis.byPA.est` (CSV here; `.Rds` under `data/store/`) | Protected area (22,090): composite coefficient, inference, counts, effectiveness class |
| `pathreat.sfa.est`, `pathreat.sfa.country.est`, `pathreat.sfa.aggregates`, `pathreat.sfa.models.est` | Protected area and country: frontier, inefficiency, efficiency; fitted `sfaR` models |
| `pathreat.byspecies.est` | Threatened species (7,279): range FE coefficient, inference, coverage |
| `pathreat.bytaxa.est` (`counts`, `effects`, `levels`) | Taxon (4): medians under equal/TBL/ED/EDGE weighting with bootstrap intervals |
| `pathreat.randomrob.est`, `pathreat.robustness.est`, `pathreat.sensitivity.*` | Specification or matching variant |
| `pathreat.covbalance.*`, `pathreat.matching_diagnostics.csv`, `pathreat.matched-sample.summary.csv` | Matching diagnostics and sample counts |
| `pathreat.threat_composite_weights.*` | Country by threat: composite-index weights, prevalence, thresholds |
| `pathreat.analysis.statistics.csv` / `.txt` | Headline statistics quoted in the paper |
| `MANIFEST.csv` | Rows, columns, size and md5 of the canonical `.Rds` tables and their CSV twins |

The CSV of `pathreat.sfa.est` is a subset that omits the PA size and IUCN category columns copied from the WDPA; the `.Rds` versions of the SFA tables retain them because the aggregation scripts need PA area, and the fitted models in `pathreat.sfa.models.est.Rds` carry log PA area in their model frames.

## Citation

Please cite the paper and the software record (see `CITATION.cff`). Zenodo concept DOI (always the latest release): `10.5281/zenodo.23108161`. The snapshot submitted with the manuscript is release `v1.0.1`; each release has its own version DOI, listed on the Zenodo record.

## Licence

Code and derived estimate tables: MIT (see `LICENSE`). Input datasets remain under their providers' licences and are not redistributed.

## Contact

Claire Vincent (corresponding author) or Elías Cisneros; see the author affiliations in the paper.
