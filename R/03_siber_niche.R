#################################################################
# Title:    Isotopic niche analysis with SIBER - Chelonia mydas and
#           other consumers of the reef food web
# Purpose:  Standard ellipse areas (SEA, SEAc, SEA-B), Layman metrics
#           and ellipse overlap between juvenile and sub-adult
#           Chelonia mydas
# Input:    data/dataset_sia.xlsx, sheet "dataset"
#           Columns used: species, type, d13C, d15N
#           (d13C and d15N in per mil, vs. VPDB and air N2)
# Outputs:  output/TA_SEAc.csv
#           output/overlap_bayes95.csv
#           output/figures/siber_ellipses.png
#           output/figures/siber_SEA_density.png
#           output/figures/siber_layman_community<k>.png
#           output/figures/siber_overlap_area_hist.png
#           output/figures/siber_overlap_prop_hist.png
#           output/sessionInfo_siber.txt
# Author:   Gabriel C. Cardozo-Ferreira
# Date:     2026-10-26
# Usage:    Open the .Rproj file (so here() finds the project root),
#           start a fresh R session, then run the script top to bottom.
#
# Requirements:
#   - JAGS must be installed on the system (needed by rjags/SIBER).
#   - R packages: SIBER, hdrcde, readxl, dplyr, here
#   - Save this file with UTF-8 encoding (the data contains "Anfípodas").
#################################################################

# ---- Packages ---------------------------------------------------
library(here)
library(readxl)
library(dplyr)
library(SIBER)

# MCMC runs in JAGS are not fully controlled by set.seed(), so small
# differences between runs of the Bayesian estimates are expected.
set.seed(3)

# ---- Paths ------------------------------------------------------
input_file  <- if (file.exists(here("data", "dataset_sia.xlsx"))) {
  here("data", "dataset_sia.xlsx")
} else {
  here("dataset_sia.xlsx")
}
input_sheet <- "dataset"
out_dir    <- here("output")
fig_dir    <- here("output", "figures")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Copies the plot currently on screen to a PNG file
save_current_plot <- function(file, width = 7, height = 5) {
  dev.copy(png, filename = file, width = width, height = height,
           units = "in", res = 300)
  invisible(dev.off())
}

# ---- Parameters -------------------------------------------------
# Sample types left out of the niche analysis
excluded_types <- "primary_producer"

# Ellipses compared in the overlap analysis ("<community>.<group>")
ellipse_juv <- "turtle_chemyd_J.Chelonia mydas"
ellipse_sa  <- "turtle_chemyd_SA.Chelonia mydas"

# Community used for the Layman metrics (1 = first community, amphipods;
# note that Layman metrics need several groups within a community)
layman_community <- 1

# Number of posterior draws used for the Bayesian overlap
overlap_draws <- 100

# Legend position in the ellipse plot (data coordinates)
legend_x <- -20.2
legend_y <- 12

# Display names for the SEA density plot, keyed by "<community>.<group>"
display_names <- c(
  "amphipod.Anfípodas"              = "Amphipod",
  "fish.Ophioblennius trinitatis"   = "Ophioblennius trinitatis",
  "fish.Sparisoma axillare"         = "Sparisoma axillare",
  "fish.Stegastes fuscus"           = "Stegastes fuscus",
  "sea-urchin.Echinometra lucunter" = "Echinometra lucunter",
  "turtle_chemyd_J.Chelonia mydas"  = "Chelonia mydas - Juvenile",
  "turtle_chemyd_SA.Chelonia mydas" = "Chelonia mydas - Sub-adult"
)

# Colours (applied through palette()) and plot symbols
my_clrs <- c("#66CC66", "#727171", "#FF6633",
             "#00B0F6", "#A3A550", "#E76BF3",
             "#B26BF3", "#00A1A1", "#16FF66")

# ---- Import data ------------------------------------------------
# No fixed cell range: read_excel() detects the used range, so rows
# added later are not silently dropped.
dataset <- read_excel(input_file, sheet = input_sheet)

head(dataset)
str(dataset)

required_cols <- c("species", "type", "d13C", "d15N")
missing_cols  <- setdiff(required_cols, names(dataset))
if (length(missing_cols) > 0) {
  stop("Missing column(s) in dataset: ", paste(missing_cols, collapse = ", "))
}

# ---- Prepare data -----------------------------------------------
data <- dataset %>%
  filter(!is.na(d13C), !is.na(d15N)) %>%
  transmute(
    iso1      = d13C,
    iso2      = d15N,
    group     = as.factor(species),
    community = as.factor(type)
  ) %>%
  as.data.frame()

head(data)

turtles <- data %>%
  filter(!community %in% excluded_types) %>%
  droplevels()                             # remove empty factor levels

turtles <- turtles[order(turtles$group), ]

unique(turtles$group)

# Lookup table (community / group codes), used in the legend
lookup <- turtles %>%
  distinct(community, group) %>%
  mutate(community_id = as.integer(community),
         group_id     = as.integer(group))

# Sample size and isotopic range per group
# (SIBER needs at least 3 samples per group to fit an ellipse)
group_check <- turtles %>%
  group_by(community, group) %>%
  summarise(
    n        = n(),
    dX_range = diff(range(iso1)),
    dY_range = diff(range(iso2)),
    .groups  = "drop"
  )
print(group_check)

if (any(group_check$n < 3)) {
  warning("Some groups have fewer than 3 samples; SIBER ellipses may fail.")
}

# ---- Create SIBER object ----------------------------------------
siber.turtles <- createSiberObject(turtles)

community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.4, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "black")

palette(my_clrs)

par(mfrow = c(1, 1))

plotSiberObject(siber.turtles,
                ax.pad = .5,
                hulls = FALSE,
                community.hulls.args = community.hulls.args,
                ellipses = TRUE,
                group.ellipses.args = group.ellipses.args,
                group.hulls = FALSE,
                group.hull.args = group.hull.args,
                bty = "L",
                pch = rep(16, 9),
                iso.order = c(1, 2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)

combo <- lookup %>% arrange(community_id, group_id)

legend(legend_x, legend_y,
       legend = paste0(combo$group, " (", combo$community, ")"),
       pch    = combo$community,           # symbol follows community
       col    = my_clrs[combo$group_id],   # colour follows group
       cex = .8, bty = "n", ncol = 1,
       pt.cex = 1.3, pt.lwd = 2)

save_current_plot(file.path(fig_dir, "siber_ellipses.png"), width = 8, height = 6)

palette("default")   # restore the default palette

# ---- Maximum likelihood metrics: TA, SEA and SEAc ----------------
group.ML <- groupMetricsML(siber.turtles)
print(group.ML)

write.csv(group.ML, file.path(out_dir, "TA_SEAc.csv"))

# ---- Bayesian ellipses (JAGS) -----------------------------------
# options for running JAGS
parms <- list()
parms$n.iter   <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3   # discard the first set of values
parms$n.thin   <- 10         # thin the posterior by this many
parms$n.chains <- 2          # run this many chains

# define the priors
priors <- list()
priors$R       <- 1 * diag(2)
priors$k       <- 2
priors$tau.mu  <- 1.0E-3

ellipses.posterior <- siberMVN(siber.turtles, parms, priors)

# ---- Compare groups: Standard Ellipse Area (SEA-B) ---------------
SEA.B <- siberEllipses(ellipses.posterior)

# x tick labels, taken from the group names in the same order as SEA.B
group_keys <- colnames(group.ML)
sea_labels <- unname(ifelse(group_keys %in% names(display_names),
                            display_names[group_keys],
                            group_keys))

# Some groups can have too little variation to support a density estimate.
# Skip those columns before calling SIBER's density plot and hdr().
has_density_data <- function(x) {
  x <- x[is.finite(x)]
  length(x) >= 3 && length(unique(x)) > 1 && stats::sd(x, na.rm = TRUE) > 0
}

SEA.B.valid <- as.data.frame(SEA.B)
SEA.B.valid_cols <- vapply(SEA.B.valid, has_density_data, logical(1))

if (sum(SEA.B.valid_cols) > 0) {
  SEA.B.plot <- SEA.B[, SEA.B.valid_cols, drop = FALSE]
  sea_labels_plot <- sea_labels[SEA.B.valid_cols]

  siberDensityPlot(SEA.B.plot,
                   xticklabels = sea_labels_plot,
                   xlab = "",
                   ylab = expression("Standard Ellipse Area " ('\u2030' ^2)),
                   bty = "L",
                   ylims = c(0, 11),
                   las = 1,
                   main = "Isotopic niche")

  # Add red x's for the ML estimated SEAc only for columns that were plotted
  points(seq_along(sea_labels_plot), group.ML[3, SEA.B.valid_cols],
         col = "red", pch = "x", lwd = 2)

  save_current_plot(file.path(fig_dir, "siber_SEA_density.png"), width = 9, height = 6)
} else {
  warning("No SEA-B columns had enough variation for density plotting; skipping siber_SEA_density.png.")
}

# Credible intervals (highest density regions, in %) and modes of SEA-B.
# 'prob' and 'all.modes' are passed explicitly to hdr().
cr.p <- c(50, 95, 99)

SEA.B.credibles <- lapply(
  SEA.B.valid[SEA.B.valid_cols],
  function(x) {
    if (!has_density_data(x)) return(NA_real_)
    tryCatch(hdrcde::hdr(x, prob = cr.p)$hdr,
             error = function(e) NA_real_)
  }
)

# all.modes = TRUE picks up multimodal posterior distributions, if present
SEA.B.modes <- lapply(
  SEA.B.valid[SEA.B.valid_cols],
  function(x) {
    if (!has_density_data(x)) return(NA_real_)
    tryCatch(hdrcde::hdr(x, all.modes = TRUE)$mode,
             error = function(e) NA_real_)
  }
)

print(SEA.B.credibles)
print(SEA.B.modes)

# ---- Compare communities: Layman metrics -------------------------
# extract the posterior means
mu.post <- extractPosteriorMeans(siber.turtles, ellipses.posterior)

# calculate the corresponding distribution of Layman metrics
layman.B <- bayesianLayman(mu.post)

# Visualise one community; skip any metrics with insufficient posterior spread.
layman_valid <- as.data.frame(layman.B[[layman_community]])
layman_valid_cols <- vapply(layman_valid, has_density_data, logical(1))

if (sum(layman_valid_cols) > 0) {
  layman_plot <- layman.B[[layman_community]][, layman_valid_cols, drop = FALSE]

  siberDensityPlot(layman_plot,
                   xticklabels = colnames(layman_plot),
                   bty = "L", ylim = c(0, 20))

  # add the ML estimates, extracted from the same community
  comm.layman.ml <- laymanMetrics(
    siber.turtles$ML.mu[[layman_community]][1, 1, ],
    siber.turtles$ML.mu[[layman_community]][1, 2, ]
  )
  points(seq_along(colnames(layman_plot)),
         comm.layman.ml$metrics[layman_valid_cols],
         col = "red", pch = "x", lwd = 2)

  save_current_plot(
    file.path(fig_dir, paste0("siber_layman_community", layman_community, ".png")))
} else {
  warning("No Layman metrics had enough variation for density plotting; skipping the community density plot.")
}

# ---- Ellipse overlap --------------------------------------------
stopifnot(all(c(ellipse_juv, ellipse_sa) %in% colnames(group.ML)))

# Maximum likelihood overlap between the SEAs (p.interval = NULL)
sea.overlap <- maxLikOverlap(ellipse_juv, ellipse_sa, siber.turtles,
                             p.interval = NULL, n = 100)
print(sea.overlap)

# Overlap between the corresponding 95% prediction ellipses
ellipse95.overlap <- maxLikOverlap(ellipse_juv, ellipse_sa, siber.turtles,
                                   p.interval = 0.95, n = 100)

# Overlap as a proportion of the total area covered by the two ellipses
prop.95.over <- ellipse95.overlap[3] / (ellipse95.overlap[2] +
                                          ellipse95.overlap[1] -
                                          ellipse95.overlap[3])
print(prop.95.over)

# Bayesian overlap between the 95% prediction ellipses
bayes95.overlap <- bayesianOverlap(ellipse_juv, ellipse_sa, ellipses.posterior,
                                   draws = overlap_draws, p.interval = 0.95,
                                   n = 100)

hist(bayes95.overlap[, 3], 10,
     main = "Bayesian overlap (95% ellipses)",
     xlab = expression("Overlap area " ('\u2030' ^2)))
save_current_plot(file.path(fig_dir, "siber_overlap_area_hist.png"))

# Same as above, as a proportion of the total area of the two ellipses
bayes.prop.95.over <- bayes95.overlap[, 3] / (bayes95.overlap[, 2] +
                                                bayes95.overlap[, 1] -
                                                bayes95.overlap[, 3])
print(bayes.prop.95.over)

hist(bayes.prop.95.over, 10,
     main = "Bayesian overlap (proportion)",
     xlab = "Overlap / total area")
save_current_plot(file.path(fig_dir, "siber_overlap_prop_hist.png"))

write.csv(data.frame(bayes95.overlap, prop.overlap = bayes.prop.95.over),
          file.path(out_dir, "overlap_bayes95.csv"), row.names = FALSE)

# ---- Session info (reproducibility) -----------------------------
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo_siber.txt"))
