## Resume the iv_grid portion after the cluster_grid refresh
## (refresh_eyes_grids.R) already completed.

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

trial_counts_for <- function(label) {
  trials <- raw |> filter(trait == label)
  ids <- unique(trials$participant_id)
  out <- as.integer(table(trials$participant_id)[ids])
  stats::setNames(out, ids)
}

regions <- c("full", "eyes", "mouth", "upper_face")
iv_focal_traits <- c("trust", "friendly", "competent", "dominant")

message("Refreshing iv_grid (16 cells) ...")
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
message("iv_grid refresh complete.")
