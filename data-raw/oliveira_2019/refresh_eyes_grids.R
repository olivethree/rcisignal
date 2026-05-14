## Re-run only the worked-example region grids that depend on the
## "eyes" face-mask geometry: cluster_grid and iv_grid. The new
## rectangle "eyes" geometry (v0.1.4) yields different numbers from
## the previous two-ellipses definition, so the cached RDS files
## must be refreshed. Everything else in
## data-raw/oliveira_2019/precompute.R is independent of the eye
## mask and is left intact.

suppressPackageStartupMessages({
  library(devtools)
  library(dplyr)
})
devtools::load_all(".")

data_dir <- Sys.getenv(
  "rcisignal_OLV_DATA",
  unset = "temp/rcicrdiagnostics/temp/oliveira_study"
)
stopifnot(dir.exists(data_dir))

cache_dir <- "vignettes/precomputed/oliveira_2019"
img_dims  <- c(256L, 256L)

message("Reading responses + noise matrix ...")
raw <- read.csv2(file.path(data_dir, "data/study1data.csv"),
                 stringsAsFactors = FALSE) |>
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

sm <- readRDS(file.path(cache_dir, "signal_matrices.rds"))

trial_counts_for <- function(trait_label) {
  sub <- raw[raw$trait == trait_label, , drop = FALSE]
  pid <- factor(sub$participant_id, levels = colnames(sm[[trait_label]]))
  as.integer(table(pid))
}

contrasts <- list(
  "Trust vs Friendly"     = list(a = "trust",     b = "friendly"),
  "Competent vs Dominant" = list(a = "competent", b = "dominant"),
  "Trust vs Dominant"     = list(a = "trust",     b = "dominant")
)
regions <- c("full", "eyes", "mouth", "upper_face")

# ---- cluster_grid ---------------------------------------------------

message("Refreshing cluster_grid (12 cells) ...")
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
  message(sprintf("  %s @ %s ...", cname, region))
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

# ---- iv_grid --------------------------------------------------------

message("Refreshing iv_grid (16 cells) ...")
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
  message(sprintf("  %s @ %s ...", label, region))
  tc     <- trial_counts_for(label)
  m      <- make_face_mask(img_dims, region = region)
  iv     <- infoval(sm[[label]], noise_matrix, tc,
                    iter = 1000L, mask = m,
                    seed = 1L, progress = FALSE)
  iv_grid$median_z[i]     <- median(iv$infoval)
  iv_grid$n_above_1.96[i] <- sum(iv$infoval >= 1.96)
}
saveRDS(iv_grid, file.path(cache_dir, "infoval_grid.rds"))

message("Refresh complete.")
