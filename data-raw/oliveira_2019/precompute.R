## Precompute the worked-example results shown in
## vignettes/rcisignal.Rmd §12.
##
## This script reads the Oliveira et al. (2019) Study 1 data (OSF
## DOI: 10.17605/osf.io/hr5pd) and produces:
##
##   - Tables (RDS) under vignettes/precomputed/oliveira_2019/
##   - Figures (PNG) under vignettes/figures/oliveira_2019/
##
## The vignette loads these on render via readRDS() and
## knitr::include_graphics() so the user sees real numbers and
## figures rather than placeholder text. Re-run this script after
## any package change that affects the metric outputs.
##
## Outputs:
##   - vignettes/precomputed/oliveira_2019/*.rds   (small tables)
##   - vignettes/figures/oliveira_2019/*.png       (figures)
##
## Data paths are configurable via the rcisignal_OLV_DATA env var.
## The default points at the development location used during the
## v0.3 build; adjust to your local copy of the OSF data.

suppressPackageStartupMessages({
  library(devtools)
  library(dplyr)
  library(png)
})

devtools::load_all(".")

# ---- paths ----------------------------------------------------------

data_dir <- Sys.getenv(
  "rcisignal_OLV_DATA",
  unset = "temp/rcicrdiagnostics/temp/oliveira_study"
)
stopifnot(dir.exists(data_dir))

fig_dir <- "vignettes/figures/oliveira_2019"
cache_dir <- "vignettes/precomputed/oliveira_2019"
dir.create(fig_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

img_dims <- c(256L, 256L)

# ---- read data ------------------------------------------------------

message("Reading responses, noise matrix, base face ...")

raw <- read.csv2(file.path(data_dir, "data/study1data.csv"),
                 stringsAsFactors = FALSE)
raw <- raw |>
  rename(participant_id = subject) |>
  mutate(
    participant_id = as.character(participant_id),
    trait          = tolower(trait),
    response       = as.integer(response),
    stimulus       = as.integer(stimulus),
    trial          = as.integer(trial)
  )

noise_matrix <- readRDS(
  file.path(data_dir, "diagnostics/noise_matrix.rds")
)
storage.mode(noise_matrix) <- "double"
dimnames(noise_matrix) <- NULL

env <- new.env(parent = emptyenv())
load(file.path(data_dir, "diagnostics/study1_stimuli_modernised.RData"),
     envir = env)
base_face <- env$base_faces[["male"]]

# Save base face as PNG (used by plot_ci_overlay).
base_face_path <- file.path(fig_dir, "base_face.png")
png::writePNG(base_face, base_face_path)

traits <- sort(unique(raw$trait))
message("Traits: ", paste(traits, collapse = ", "))

# ---- 1. diagnostic report ------------------------------------------

message("Running run_diagnostics() ...")
diag_report <- run_diagnostics(
  raw |> select(participant_id, stimulus, response),
  method     = "2ifc",
  rdata      = file.path(data_dir, "diagnostics/study1_stimuli_modernised.RData"),
  expected_n = 300L
)
saveRDS(diag_report, file.path(cache_dir, "diagnostic_report.rds"))

# Per-trait response-bias check.
trait_bias <- lapply(traits, function(tr) {
  sub <- raw |> filter(trait == tr) |>
    select(participant_id, stimulus, response)
  check_response_bias(sub, method = "2ifc")
})
names(trait_bias) <- traits
saveRDS(trait_bias, file.path(cache_dir, "trait_bias.rds"))

# ---- 2. signal matrices per trait ----------------------------------

message("Building per-trait signal matrices ...")

build_signal_matrix <- function(label) {
  trials <- raw |> filter(trait == label) |>
    arrange(participant_id, trial)
  ids <- unique(trials$participant_id)
  m <- matrix(NA_real_, nrow = nrow(noise_matrix),
              ncol = length(ids),
              dimnames = list(NULL, ids))
  for (i in seq_along(ids)) {
    pr <- trials |> filter(participant_id == ids[i])
    m[, i] <- (noise_matrix[, pr$stimulus] %*% pr$response) /
                nrow(pr)
  }
  attr(m, "img_dims") <- img_dims
  attr(m, "source")   <- "raw"
  m
}

sm <- stats::setNames(lapply(traits, build_signal_matrix), traits)
saveRDS(sm, file.path(cache_dir, "signal_matrices.rds"))

# ---- 3. within-condition reliability per trait ---------------------

message("Computing reliability per trait ...")

rel_table <- data.frame(
  trait    = traits,
  r_sb     = NA_real_,
  icc_3_1  = NA_real_,
  icc_3_k  = NA_real_,
  stringsAsFactors = FALSE
)
for (i in seq_along(traits)) {
  rep <- run_reliability(sm[[traits[i]]],
                         n_permutations = 2000L,
                         seed = 1L, progress = FALSE)
  rel_table$r_sb[i]    <- rep$results$split_half$r_sb
  rel_table$icc_3_1[i] <- rep$results$icc$icc_3_1
  rel_table$icc_3_k[i] <- rep$results$icc$icc_3_k
}
saveRDS(rel_table, file.path(cache_dir, "reliability_table.rds"))

# ---- 4. infoVal per trait, full-face oval --------------------------

message("Computing infoVal per trait (full-face oval) ...")

trial_counts_for <- function(label) {
  trials <- raw |> filter(trait == label)
  ids <- unique(trials$participant_id)
  out <- as.integer(table(trials$participant_id)[ids])
  stats::setNames(out, ids)
}

oval_mask <- make_face_mask(img_dims, region = "full")

iv_full <- data.frame(
  trait                 = traits,
  median_z              = NA_real_,
  n_above_1.96          = NA_integer_,
  group_z               = NA_real_,
  stringsAsFactors      = FALSE
)
for (i in seq_along(traits)) {
  label <- traits[i]
  tc <- trial_counts_for(label)
  iv <- infoval(sm[[label]], noise_matrix, tc,
                iter = 1000L, mask = oval_mask,
                seed = 1L, progress = FALSE)
  iv_full$median_z[i]      <- median(iv$infoval)
  iv_full$n_above_1.96[i]  <- sum(iv$infoval >= 1.96)
  # Group-mean CI z, with a reference distribution matched to the
  # group construction (N random producers averaged together, each
  # at its real trial count). Uses the same internal helper that
  # diagnose_infoval() relies on, so the two stay consistent.
  # iter=500 here (not 1000) because each iteration is N=20x more
  # expensive than per-producer infoval; 500 gives single-z stability
  # of about +/- 0.07, which is fine for a publishable headline.
  iv_full$group_z[i] <- rcisignal:::group_mean_z(
    sm[[label]], noise_matrix, tc,
    iter = 500L, mask = oval_mask,
    seed = 1L, progress = FALSE
  )
}
saveRDS(iv_full, file.path(cache_dir, "infoval_full.rds"))

# ---- 5. region-by-region cluster tests -----------------------------

message("Computing region-by-region cluster tests for three contrasts ...")

contrasts <- list(
  "Trust vs Friendly"     = list(a = "trust",     b = "friendly"),
  "Competent vs Dominant" = list(a = "competent", b = "dominant"),
  "Trust vs Dominant"     = list(a = "trust",     b = "dominant")
)
regions <- c("full", "eyes", "mouth", "upper_face")

cluster_grid <- expand.grid(
  contrast = names(contrasts),
  region   = regions,
  stringsAsFactors = FALSE
)
cluster_grid$n_clusters    <- NA_integer_
cluster_grid$n_significant <- NA_integer_
cluster_grid$min_p         <- NA_real_

for (i in seq_len(nrow(cluster_grid))) {
  cname  <- cluster_grid$contrast[i]
  region <- cluster_grid$region[i]
  m  <- make_face_mask(img_dims, region = region)
  ct <- rel_cluster_test(
    sm[[contrasts[[cname]]$a]], sm[[contrasts[[cname]]$b]],
    img_dims          = img_dims,
    mask              = m,
    cluster_threshold = 2.0,
    n_permutations    = 2000L,
    seed              = 1L,
    progress          = FALSE
  )
  cl <- ct$clusters
  cluster_grid$n_clusters[i]    <- if (is.null(cl)) 0L else nrow(cl)
  cluster_grid$n_significant[i] <- sum(cl$significant, na.rm = TRUE)
  cluster_grid$min_p[i] <-
    if (is.null(cl) || nrow(cl) == 0L) NA_real_ else
      min(cl$p_value, na.rm = TRUE)
}
saveRDS(cluster_grid, file.path(cache_dir, "cluster_grid.rds"))

# ---- 6. dissimilarity per contrast ---------------------------------

message("Computing dissimilarity per contrast (full face) ...")

dissim_full <- list()
for (cname in names(contrasts)) {
  dr <- rel_dissimilarity(
    sm[[contrasts[[cname]]$a]], sm[[contrasts[[cname]]$b]],
    n_boot = 2000L, seed = 1L, progress = FALSE
  )
  dissim_full[[cname]] <- dr
}
saveRDS(dissim_full, file.path(cache_dir, "dissim_full.rds"))

# ---- 7. infoVal per trait per region -------------------------------

message("Computing infoVal per trait per region ...")

iv_focal_traits <- c("trust", "friendly", "competent", "dominant")

iv_grid <- expand.grid(
  trait  = iv_focal_traits,
  region = regions,
  stringsAsFactors = FALSE
)
iv_grid$median_z     <- NA_real_
iv_grid$n_above_1.96 <- NA_integer_

for (i in seq_len(nrow(iv_grid))) {
  label  <- iv_grid$trait[i]
  region <- iv_grid$region[i]
  tc     <- trial_counts_for(label)
  m      <- make_face_mask(img_dims, region = region)
  iv     <- infoval(sm[[label]], noise_matrix, tc,
                    iter = 1000L, mask = m,
                    seed = 1L, progress = FALSE)
  iv_grid$median_z[i]     <- median(iv$infoval)
  iv_grid$n_above_1.96[i] <- sum(iv$infoval >= 1.96)
}
saveRDS(iv_grid, file.path(cache_dir, "infoval_grid.rds"))

# ---- 8. pairwise cluster-agreement maps ---------------------------

message("Generating pairwise cluster-agreement maps ...")

build_sig_mask <- function(ct) {
  sig_pos <- ct$clusters$cluster_id[
    ct$clusters$direction == "pos" & ct$clusters$significant
  ]
  sig_neg <- ct$clusters$cluster_id[
    ct$clusters$direction == "neg" & ct$clusters$significant
  ]
  as.vector((ct$pos_labels %in% sig_pos) |
              (ct$neg_labels %in% sig_neg))
}

pairwise <- list(
  list(name = "trust_vs_friendly",
       a = "trust", b = "friendly",
       title = "Trust minus Friendly (FWER-controlled clusters)"),
  list(name = "dominant_vs_competent",
       a = "dominant", b = "competent",
       title = "Dominant minus Competent (FWER-controlled clusters)")
)

cluster_pairwise <- list()
for (p in pairwise) {
  message("  cluster test: ", p$name)
  ct <- rel_cluster_test(
    sm[[p$a]], sm[[p$b]],
    img_dims          = img_dims,
    cluster_threshold = 2.0,
    n_permutations    = 2000L,
    seed              = 1L,
    progress          = FALSE
  )
  cluster_pairwise[[p$name]] <- ct

  diff_signal <- rowMeans(sm[[p$a]]) - rowMeans(sm[[p$b]])

  desc_path <- file.path(fig_dir,
                         sprintf("pairwise_descriptive_%s.png", p$name))
  desc_title <- sub(" \\(FWER-controlled clusters\\)$",
                    " (descriptive)", p$title)
  grDevices::png(desc_path, width = 600, height = 600, res = 100)
  plot_ci_overlay(
    diff_signal,
    base_image = base_face_path,
    img_dims   = img_dims,
    mask       = oval_mask,
    main       = desc_title
  )
  grDevices::dev.off()

  sig_mask <- build_sig_mask(ct)
  out_path <- file.path(fig_dir, sprintf("pairwise_%s.png", p$name))
  grDevices::png(out_path, width = 600, height = 600, res = 100)
  plot_ci_overlay(
    diff_signal,
    base_image = base_face_path,
    img_dims   = img_dims,
    mask       = sig_mask,
    main       = p$title
  )
  grDevices::dev.off()
}

saveRDS(cluster_pairwise,
        file.path(cache_dir, "cluster_pairwise.rds"))

# ---- 9. dissimilarity grid figure ---------------------------------

message("Generating plot_dissimilarity_grid figure ...")
grDevices::png(file.path(fig_dir, "dissim_grid.png"),
               width = 900, height = 360, res = 100)
plot_dissimilarity_grid(
  "Trust vs Friendly"     = dissim_full[["Trust vs Friendly"]],
  "Competent vs Dominant" = dissim_full[["Competent vs Dominant"]],
  "Trust vs Dominant"     = dissim_full[["Trust vs Dominant"]]
)
grDevices::dev.off()

# ---- 10. face-mask reference figures ------------------------------

message("Generating face-mask reference figures ...")
for (region in regions) {
  m <- make_face_mask(img_dims, region = region)
  grDevices::png(
    file.path(fig_dir, sprintf("mask_%s.png", region)),
    width = 360, height = 360, res = 100
  )
  plot_face_mask(m, img_dims = img_dims,
                 base_image = base_face_path,
                 main = paste("region =", region))
  grDevices::dev.off()
}

message("Precompute complete. Outputs:")
message("  ", cache_dir)
message("  ", fig_dir)
