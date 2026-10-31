# Stable isotopes and diet of *Chelonia mydas* – reef food web

R code and data for the analysis of the trophic ecology of the green turtle
(*Chelonia mydas*, juveniles and sub-adults) and its potential food sources
on reefs of northeastern Brazil (Maragogi, Japaratinga and São Miguel dos
Milagres). The analyses combine stable isotopes (δ¹³C and δ¹⁵N) and diet data.

All code and data are public in this repository.

## Repository structure

```
.
├── R/
│   ├── 01_biplot_sia.R            # δ13C–δ15N biplot (means ± SD)
│   ├── 02_trophic_position.R      # Bayesian trophic position
│   ├── 03_siber_niche.R           # isotopic niche (SIBER)
│   └── 04_diet_analysis.R         # diet: PERMANOVA, nMDS, niche indices, Manly, Figure S2
├── data/
│   ├── dataset_sia.xlsx           # isotope data (sheets: dataset, trophicposition)
│   ├── diet_data.xlsx             # diet data
│   └── benthic_cover.xlsx         # benthic cover of the sampled reefs
├── output/                        # created by the scripts (tables, figures)
│   └── figures/
├── .gitignore
├── LICENSE
└── README.md
```

## Scripts

Run the scripts from the project root, each one in a fresh R session. They
are independent of each other.

| Script | What it does | Data used |
|---|---|---|
| `01_biplot_sia.R` | Biplot of δ¹³C vs δ¹⁵N with individual samples and group means ± SD | `dataset_sia.xlsx`, sheet `dataset` |
| `02_trophic_position.R` | Trophic position of each consumer with `tRophicPosition` (one-baseline Bayesian model), credibility intervals and pairwise comparisons | `dataset_sia.xlsx`, sheet `trophicposition` |
| `03_siber_niche.R` | Standard ellipse areas (SEAc, SEA-B), Layman metrics and ellipse overlap between juvenile and sub-adult turtles with `SIBER` | `dataset_sia.xlsx`, sheet `dataset` |
| `04_diet_analysis.R` | PERMANOVA and nMDS of diet composition, Levins niche breadth, Morisita-Horn overlap, Manly selection ratios, and Figure S2 (benthic composition) | `diet_data.xlsx`, `benthic_cover.xlsx` |

Each script has a header listing its inputs and outputs. Tables are saved in
`output/` and figures in `output/figures/`.

## Requirements

- R (>= 4.0)
- [JAGS](https://mcmc-jags.sourceforge.io/) (needed by `02` and `03`)
- R packages:

```r
install.packages(c("here", "readxl", "dplyr", "tidyr", "forcats", "ggplot2",
                   "vegan", "writexl", "adehabitatHS", "SIBER", "hdrcde",
                   "remotes"))

# tRophicPosition is not on CRAN
remotes::install_github("clquezada/tRophicPosition")
```

The R and package versions used are saved by each script in
`output/sessionInfo*.txt`.

## How to run

1. Clone this repository and open the `.Rproj` file in RStudio (this makes
   `here()` find the project root).
2. Check that the files in `data/` are present.
3. Open a script from `R/` and run it from top to bottom.

Bayesian models (JAGS) are not fully controlled by `set.seed()`, so small
differences between runs are expected in `02` and `03`.

## Data

| File | Sheet / range | Used by | Contents |
|---|---|---|---|
| `data/dataset_sia.xlsx` | `dataset` | `01`, `03` | One row per isotope sample: `location`, `species`, `type`, `ccc`, `size_class`, `sample_id`, `d13C`, `d15N`, `C`, `N` |
| `data/dataset_sia.xlsx` | `trophicposition` | `02` | `study`, `location`, `spp`, `FG` (functional group; `Benthic_BL` = benthic baseline), `d13C`, `d15N`, `NS`, `grouping` |
| `data/diet_data.xlsx` | `Sheet1`, range `A1:AK37` | `04` | One row per turtle: `Code`, `Stage`, `Class`, `CCC`, `EST_g`, `Triado_g`, then the dietary categories (%) |
| `data/benthic_cover.xlsx` | `benthic_cover` | `04` | `Division`, `family`, `item`, then benthic cover (%) at each site: `Maragogi`, `Japaratinga`, `SM_Milagres` |

**Units and conventions**

- δ¹³C and δ¹⁵N are in ‰ (δ¹³C relative to VPDB, δ¹⁵N relative to atmospheric N₂).
- Diet composition and benthic cover are in percent.
- `type` / `grouping` identify the sample group (for example `primary_producer`,
  `fish`, `sea-urchin`, `amphipod`, `turtle_chemyd_J` for juvenile and
  `turtle_chemyd_SA` for sub-adult *C. mydas*).
- TODO: describe `ccc` / `CCC` (unit), `size_class`, `Class`, `EST_g`,
  `Triado_g`, `NS`, `C` and `N`, and how missing values are coded.

Sample collection and laboratory methods are described in the associated
manuscript. The resource availability and use values for the Manly selection
ratios are given in `R/04_diet_analysis.R`; availability corresponds to the
benthic composition shown in Figure S2.

## Citation

The manuscript describing this work has been submitted and is not yet
published. Citation information will be added here after publication.

## License

- Code: MIT License (see `LICENSE`).
- Data and figures: CC BY 4.0.

## Contact

Gabriel C. Cardozo-Ferreira – gabrielccf@gmail.com
