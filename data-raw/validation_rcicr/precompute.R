## Precompute per-trait per-participant infoVal side-by-side
## (rcisignal::infoval() vs rcicr::computeInfoVal2IFC()) for the
## Oliveira et al. (2019) Study 1 data, as shown in vignette
## "validation_rcicr".
##
## Reads the same data the oliveira_2019 precompute uses (OSF DOI:
## 10.17605/osf.io/hr5pd) and writes:
##
##   vignettes/precomputed/validation_rcicr/rcicr_per_trait.rds
##
## Output shape: named list (one element per trait) of data.frames
## with columns:
##   - participant_id (character)
##   - n_trials (integer)
##   - rcisignal_z (double)  per-trial-count reference
##   - rcicr_z    (double)   pool-size reference (rcicr default)
##
## Runtime (sequential): about 30-45 minutes (10 traits x ~3-4 min
## at iter=10000). With trait-level parallelism on a modern Mac
## (8 workers), expect ~6-10 minutes wall time.
##
## Data paths are configurable via the rcisignal_OLV_DATA env var,
## matching data-raw/oliveira_2019/precompute.R. Worker count via
## the rcisignal_VAL_CORES env var (defaults to detectCores() - 1).
## If the rdata references base-face image files that don't live
## under `data_dir`, point the script at the right directory via
## the rcisignal_OLV_BASE_DIR env var.

## Cap BLAS threads BEFORE loading any package that initialises a
## BLAS handle. Trait-level parallelism via mclapply forks N workers
## each running heavy matrix multiplies; without this cap, Accelerate
## (macOS) / OpenBLAS / MKL spawn additional threads per worker and
## oversubscribe the CPU, making the parallel run slower than serial.
Sys.setenv(
  VECLIB_MAXIMUM_THREADS = "1",  # Apple Accelerate
  OPENBLAS_NUM_THREADS   = "1",
  MKL_NUM_THREADS        = "1",
  OMP_NUM_THREADS        = "1"
)

suppressPackageStartupMessages({
  library(devtools)
  library(dplyr)
  library(parallel)
})

# Belt-and-braces: also set runtime BLAS threads to 1 if RhpcBLASctl
# is available. Env vars only bite if read at BLAS init; this catches
# the case where R has already initialised BLAS in this session.
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}

devtools::load_all(".")

stopifnot(requireNamespace("rcicr", quietly = TRUE))

# rcicr's batchGenerateCI2IFC() uses %dopar% and, if a foreach
# parallel backend is registered (or it tries to spawn its own
# PSOCK cluster), each mclapply-forked worker ends up nesting
# clusters. PSOCK sub-workers cannot reliably open sockets back to
# a forked-child master and die with "cannot open the connection".
# Force foreach to run %dopar% sequentially inside each fork, so
# rcicr stays single-threaded per worker and parallelism happens
# only at the trait level.
if (requireNamespace("foreach", quietly = TRUE)) {
  foreach::registerDoSEQ()
}

# ---- paths ----------------------------------------------------------

data_dir <- Sys.getenv(
  "rcisignal_OLV_DATA",
  unset = "temp/rcicrdiagnostics/temp/oliveira_study"
)
stopifnot(dir.exists(data_dir))

cache_dir <- "vignettes/precomputed/validation_rcicr"
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

rdata_path <- file.path(data_dir,
                        "diagnostics/study1_stimuli_modernised.RData")
stopifnot(file.exists(rdata_path))

img_dims <- c(256L, 256L)

# The Oliveira rdata stores base faces under labels "male" / "female".
# ci_from_responses_2ifc()'s default base_image branch expects "base"
# (the canonical rcicr label); pick "male" explicitly here.
base_label <- "male"

# Patch the Oliveira rdata once before forking workers. The file
# pre-dates the rcicr version that started persisting noise_type /
# nscales / sigma into the rdata environment, but
# rcicr::generateReferenceDistribution2IFC() does `load(rdata)` and
# then references these names directly. Without them the call aborts
# with `object 'noise_type' not found`. Defaults match rcicr's own
# generateStimuli2IFC() defaults. Saved as a new tempfile; the
# canonical data file stays untouched.
patch_rdata <- function(src, dst, data_dir) {
  env <- new.env(parent = emptyenv())
  load(src, envir = env)
  added <- character()
  if (!exists("noise_type", envir = env, inherits = FALSE)) {
    env$noise_type <- "sinusoid"
    added <- c(added, "noise_type")
  }
  if (!exists("nscales", envir = env, inherits = FALSE)) {
    env$nscales <- 5L
    added <- c(added, "nscales")
  }
  if (!exists("sigma", envir = env, inherits = FALSE)) {
    env$sigma <- 25
    added <- c(added, "sigma")
  }

  # Resolve relative base-face image paths against the data dir.
  # rcicr reads these via base R file I/O, which uses the worker's
  # CWD; relative paths fail when CWD is not the data dir.
  if (exists("base_face_files", envir = env, inherits = FALSE)) {
    bff <- env$base_face_files
    # User-provided override takes precedence over the search.
    extra <- Sys.getenv("rcisignal_OLV_BASE_DIR", unset = "")
    candidates <- unique(c(
      if (nzchar(extra)) extra else character(0L),
      data_dir,
      dirname(src),
      file.path(data_dir, "base_stimuli"),
      file.path(dirname(src), "base_stimuli"),
      file.path(data_dir, "baseimages"),
      file.path(dirname(src), "baseimages"),
      file.path(data_dir, "diagnostics"),
      file.path(data_dir, "diagnostics", "baseimages")
    ))
    # Preserve original type (list vs named character vector); rcicr
    # may rely on either via `[[name]]` indexing.
    resolved <- bff
    for (i in seq_along(bff)) {
      p <- bff[[i]]
      if (file.exists(p)) {
        resolved[[i]] <- normalizePath(p, mustWork = TRUE)
      } else {
        hit <- NA_character_
        for (d in candidates) {
          cand <- file.path(d, p)
          if (file.exists(cand)) {
            hit <- normalizePath(cand, mustWork = TRUE)
            break
          }
        }
        if (is.na(hit)) {
          stop("Could not locate base-face image '", p,
               "' under ", data_dir,
               ". Candidate dirs tried:\n  ",
               paste(candidates, collapse = "\n  "))
        }
        resolved[[i]] <- hit
      }
    }
    env$base_face_files <- resolved
    added <- c(added, "base_face_files (absolutised)")
  }

  save(list = ls(env), envir = env, file = dst)
  added
}

rdata_fixed <- tempfile(fileext = ".RData")
added <- patch_rdata(rdata_path, rdata_fixed, data_dir)
if (length(added) > 0L) {
  message("Patched rdata with: ", paste(added, collapse = "; "))
}

# ---- read data ------------------------------------------------------

message("Reading responses + noise matrix ...")

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

traits <- sort(unique(raw$trait))
message("Traits: ", paste(traits, collapse = ", "))

# ---- helpers --------------------------------------------------------

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

process_trait <- function(tr) {
  trial_data <- raw |> filter(trait == tr)

  sm_tr <- build_signal_matrix(tr)
  ids   <- colnames(sm_tr)
  tc    <- stats::setNames(
    as.integer(table(trial_data$participant_id)[ids]),
    ids
  )

  # rcicr writes its reference-distribution cache back into the
  # rdata file on first call. Give each trait its own tempfile
  # copy so parallel workers never write to the same path and the
  # canonical data file stays untouched.
  rdata_local <- tempfile(fileext = ".RData")
  on.exit(unlink(rdata_local), add = TRUE)
  file.copy(rdata_fixed, rdata_local, overwrite = TRUE)

  iv_rs <- infoval(sm_tr, noise_matrix, tc,
                   iter = 10000L, seed = 1L, progress = FALSE)

  # rcicr::computeInfoVal2IFC() calls bare write(...) which writes
  # a file literally named "data" (and similar) to CWD. With N
  # workers all hitting the same CWD simultaneously, this races on
  # the same path and leaks debris into the package root. Run rcicr
  # from a per-worker tempdir so any leaked files land somewhere
  # disposable.
  worker_wd <- tempfile("rcicr_wd_")
  dir.create(worker_wd, recursive = TRUE)
  old_wd <- getwd()
  setwd(worker_wd)
  on.exit({ setwd(old_wd); unlink(worker_wd, recursive = TRUE) },
          add = TRUE, after = FALSE)

  # rcicr's native 2IFC infoVal: build per-producer CIs via its
  # batch helper, then call computeInfoVal2IFC() one producer at a
  # time. This is what compute_infoval_summary() used to wrap; with
  # that wrapper gone, the validation precompute hits rcicr directly
  # so the engine-equivalence vignette still compares like with like.
  for (pkg in c("foreach", "tibble", "dplyr")) {
    if (!paste0("package:", pkg) %in% search()) {
      suppressPackageStartupMessages(attachNamespace(pkg))
    }
  }
  rc_cis <- rcicr::batchGenerateCI2IFC(
    data        = trial_data |>
      select(participant_id, stimulus, response) |>
      as.data.frame(),
    by          = "participant_id",
    stimuli     = "stimulus",
    responses   = "response",
    baseimage   = base_label,
    rdata       = rdata_local,
    save_as_png = FALSE,
    targetpath  = tempfile("rc_targets_")
  )
  rc_z <- vapply(
    names(rc_cis),
    function(nm) rcicr::computeInfoVal2IFC(rc_cis[[nm]],
                                           rdata_local,
                                           iter = 10000L),
    numeric(1L)
  )
  rc_ids <- sub(paste0("^", base_label, "_participant_id_"), "",
                names(rc_cis))
  rc_df <- data.frame(participant_id = rc_ids,
                      rcicr_z        = unname(rc_z),
                      stringsAsFactors = FALSE)

  rs_df <- data.frame(
    participant_id = ids,
    n_trials       = unname(tc),
    rcisignal_z    = unname(iv_rs$infoval),
    stringsAsFactors = FALSE
  )

  merge(rs_df, rc_df, by = "participant_id",
        all = FALSE, sort = FALSE)
}

# ---- run (parallel by trait) ---------------------------------------

n_cores <- as.integer(
  Sys.getenv("rcisignal_VAL_CORES",
             unset = max(1L, parallel::detectCores() - 1L))
)
n_cores <- min(n_cores, length(traits))

message(sprintf("Running %d traits across %d worker(s) ...",
                length(traits), n_cores))
t0 <- Sys.time()

if (n_cores > 1L && .Platform$OS.type == "unix") {
  per_trait <- parallel::mclapply(
    traits, process_trait,
    mc.cores       = n_cores,
    mc.preschedule = FALSE
  )
} else {
  if (n_cores > 1L) {
    message("Note: trait-level parallelism uses fork (mclapply); ",
            "falling back to sequential on this platform.")
  }
  per_trait <- lapply(traits, function(tr) {
    message("Trait: ", tr)
    process_trait(tr)
  })
}
names(per_trait) <- traits

# Surface worker errors (mclapply returns try-error objects in
# place of results when a worker fails).
errs <- vapply(per_trait, inherits, logical(1L), what = "try-error")
if (any(errs)) {
  bad <- names(per_trait)[errs]
  stop("Worker error(s) for trait(s): ", paste(bad, collapse = ", "),
       "\nFirst error: ",
       attr(per_trait[[bad[1L]]], "condition")$message)
}

dt <- difftime(Sys.time(), t0, units = "mins")
message(sprintf("Done in %.1f min.", as.numeric(dt)))

saveRDS(per_trait, file.path(cache_dir, "rcicr_per_trait.rds"))

message("Output:")
message("  ", file.path(cache_dir, "rcicr_per_trait.rds"))
