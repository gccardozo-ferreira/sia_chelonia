#################################################################
# Title:    Bayesian trophic position (tRophicPosition) - Chelonia mydas
#           and other consumers of the reef food web
# Purpose:  Estimate trophic position (TP) of each consumer using the
#           benthic baseline, compare TPs between consumers, and plot
#           credibility intervals and a food web biplot.
# Input:    data/dataset_sia.xlsx, sheet "trophicposition"
#           Columns used: study, location, spp, FG, d13C, d15N, grouping
# Outputs:  output/TP_models_df.csv
#           output/TP_quantiles.csv
#           output/pairwiseTP.csv
#           output/figures/TP_credibility_intervals.png
#           output/figures/foodweb_biplot.png
#           output/sessionInfo_trophic_position.txt
# Author:   Gabriel C. Cardozo-Ferreira
# Date:     2026-10-26
# Usage:    Open the .Rproj file (so here() finds the project root),
#           start a fresh R session, then run the script top to bottom.
#
# Requirements:
#   - JAGS must be installed on the system (needed by rjags).
#   - tRophicPosition is not on CRAN. Install it once with:
#       install.packages("remotes")
#       remotes::install_github("clquezada/tRophicPosition")
#################################################################

# ---- Packages ---------------------------------------------------
library(here)
library(readxl)
library(dplyr)
library(forcats)
library(ggplot2)
library(tRophicPosition)

set.seed(3)

# ---- Paths ------------------------------------------------------
input_file <- here("data", "dataset_sia.xlsx")
out_dir    <- here("output")
fig_dir    <- here("output", "figures")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Parameters -------------------------------------------------
baseline_label     <- "Benthic_BL"   # value of column FG used as baseline
show_isotope_plots <- TRUE           # print/plot each isotopeData object

# MCMC settings
n_adapt  <- 10000
n_iter   <- 10000
n_burnin <- 10000
n_chains <- 2

# Colours for the credibility interval plot
my_cols <- c("#EC663B", "#66CC66", "#727171", "#E4E255",
             "#00B0F6", "#727171", "#E76BF3")

# ---- Consumers and trophic discrimination factors (TDFs) -------
# One row per consumer. Samples are selected where
# <filter_col> == <filter_val>; baseline samples (FG == baseline_label)
# are always included. TDFs are simulated with simulateTDF().
consumer_config <- tribble(
  ~label,                      ~filter_col, ~filter_val,         ~nN, ~meanN, ~sdN,  ~nC, ~meanC, ~sdC,  ~tdf_source,
  "Chelonia mydas (juvenile)",  "grouping", "turtle_chemyd_J",    20,  4.15,  0.47,   20,  2.18,  0.59,  "<add reference>",
  "Chelonia mydas (subadult)",  "grouping", "turtle_chemyd_SA",   20,  4.93,  0.59,   20,  2.58,  1.19,  "<add reference>",
  "Amphipods",                  "grouping", "amphipod",            5,  2.60,  0.15,    5, -0.11,  0.14,  "Caut et al. 2009 (check reference)",
  "Echinometra lucunter",       "spp",      "Echinometra lucunter", 11, 0.62,  0.15,   11,  0.84,  0.51,  "L. variegatus, Prado et al. 2012",
  "Sparisoma axillare",         "spp",      "Sparisoma axillare", 20,  3.58,  0.09,   20,  0.27,  0.25,  "Caut et al. 2009",
  "Ophioblennius trinitatis",   "spp",      "Ophioblennius trinitatis", 10, 3.75, 0.13, 10, -0.24, 0.11, "Caut et al. 2009",
  "Stegastes fuscus",           "spp",      "Stegastes fuscus",   25,  3.02,  0.06,   25,  0.06,  0.13,  "Caut et al. 2009"
)

# ---- Import data ------------------------------------------------
# No fixed cell range: read_excel() detects the used range, so rows
# added later are not silently dropped.
dataset <- read_excel(input_file, sheet = "trophicposition") %>%
  as.data.frame()

str(dataset)

required_cols <- c("study", "location", "spp", "FG", "d13C", "d15N", "grouping")
missing_cols  <- setdiff(required_cols, names(dataset))
if (length(missing_cols) > 0) {
  stop("Missing column(s) in dataset: ", paste(missing_cols, collapse = ", "))
}

# Number of samples per species and group (use it to check the filters below)
print(table(dataset$spp, dataset$grouping, useNA = "ifany"))

# ---- Trophic position: one model per consumer ------------------
TP_dfs  <- list()   # summary data frame of each model
TP_post <- list()   # posterior TP distributions of each model

for (i in seq_len(nrow(consumer_config))) {

  cfg <- consumer_config[i, ]

  # Baseline samples + the samples of this consumer only
  is_consumer <- dataset[[cfg$filter_col]] == cfg$filter_val
  n_consumer  <- sum(is_consumer, na.rm = TRUE)
  if (n_consumer == 0) {
    stop("No samples found for '", cfg$label, "' (", cfg$filter_col,
         " == '", cfg$filter_val, "'). Check the table printed above.")
  }

  mydata.TP <- dataset[which(dataset$FG == baseline_label | is_consumer), ]
  message(cfg$label, ": ", n_consumer, " consumer samples, ",
          sum(mydata.TP$FG == baseline_label, na.rm = TRUE), " baseline samples")

  # Consumer TDFs
  TDFn <- simulateTDF(nN = cfg$nN, meanN = cfg$meanN, sdN = cfg$sdN)
  TDFc <- simulateTDF(nC = cfg$nC, meanC = cfg$meanC, sdC = cfg$sdC)

  TPList <- extractIsotopeData(mydata.TP,
                               b1 = baseline_label, b2 = NULL,
                               baselineColumn = "FG",
                               consumersColumn = "spp",
                               groupsColumn = "location",
                               d13C = "d13C", d15N = "d15N",
                               deltaC = TDFc, deltaN = TDFn,
                               seed = 3)

  stopifnot(length(TPList) > 0)

  if (show_isotope_plots) {
    for (community in TPList) {
      print(summary(community))
      plot(community)
    }
  }

  TP_run <- multiSpeciesTP(TPList, model = "oneBaseline",
                           n.adapt = n_adapt, n.iter = n_iter,
                           burnin = n_burnin, n.chains = n_chains,
                           lambda = 1, print = FALSE)

  # Summary data frame, labelled with the consumer name
  run_df          <- TP_run$df
  run_df$consumer <- cfg$label
  TP_dfs[[i]]     <- run_df

  # Posterior distributions, with unique names across consumers
  post <- as.list(TP_run$TPs)
  names(post) <- if (length(post) == 1) {
    cfg$label
  } else {
    paste(cfg$label, names(post), sep = " - ")
  }
  TP_post <- c(TP_post, post)
}

# ---- Combine results --------------------------------------------
TP_df <- bind_rows(TP_dfs)
print(TP_df)

write.csv(TP_df, file.path(out_dir, "TP_models_df.csv"), row.names = FALSE)

TP_quantiles <- sapply(TP_post, quantile, probs = c(0.025, 0.5, 0.975)) %>%
  round(3)
print(TP_quantiles)
write.csv(t(TP_quantiles), file.path(out_dir, "TP_quantiles.csv"))

TP_E <- getPosteriorMode(TP_post)
print(TP_E)

# ---- Credibility intervals plot ---------------------------------
TP_df %>%
  mutate(consumer = fct_reorder(consumer, mode)) %>%
  credibilityIntervals(
    x = "consumer",
    xlab = "Consumers",
    plotAlpha = FALSE,
    legend = NULL,
    y1lim = c(1, 3.1),
    group_by = "consumer",
    scale_colour_manual = my_cols
  )

# ggsave() without 'plot' saves the last ggplot created
ggsave(file.path(fig_dir, "TP_credibility_intervals.png"),
       width = 8, height = 5, dpi = 300)

# ---- Food web biplot (d13C vs d15N) -----------------------------
dplyr::glimpse(dataset)

foodweb_biplot <- screenFoodWeb(dataset,
                                grouping = c("spp", "study"),
                                order = TRUE, printSummary = TRUE)

ggsave(file.path(fig_dir, "foodweb_biplot.png"),
       width = 8, height = 6, dpi = 300)

# ---- Pairwise comparison of posterior trophic positions ---------
pairwiseTP <- pairwiseComparisons(TP_post, print = TRUE)
print(pairwiseTP)

write.csv(pairwiseTP, file.path(out_dir, "pairwiseTP.csv"))

# ---- Session info (reproducibility) -----------------------------
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo_trophic_position.txt"))
