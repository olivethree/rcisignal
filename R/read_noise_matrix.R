#' Read a noise matrix from any supported source, with caching
#'
#' @description
#' Returns the `pixels x pool_size` numeric matrix of basis noise
#' patterns that [ci_from_responses_briefrc()], [infoval()], and the
#' diagnostic checks need. The function detects the input format
#' from the file extension (with a magic-byte fallback) and
#' transparently caches slow-to-parse formats to a sibling `.rds`
#' so the second read is near-instant.
#'
#' Supported sources, with no `format` flag — detection is
#' automatic:
#'
#' * **`.rds`** (R serialized object). Loaded directly. This is the
#'   fast path; everything else gets cached to one of these.
#' * **`.txt` / `.csv` / `.tsv`** (delimited text matrix, one column
#'   per stimulus). Read via [data.table::fread()]; the parsed
#'   matrix is materialised to `<basename>.rds` next to the source
#'   for next time.
#' * **`.RData` / `.Rdata`** (rcicr stimulus generation output).
#'   The rdata does not store the noise matrix directly — it
#'   stores `stimuli_params` (per-trial parameter rows) and `p`
#'   (the sinusoid / gabor basis). This loader reconstructs each
#'   trial's noise pattern via [rcicr::generateNoiseImage()] and
#'   caches the result to `<basename>.rds`. Requires `rcicr` to be
#'   installed for the first read; subsequent reads from the cache
#'   do not.
#'
#' @details
#' The cache invalidates automatically. The `.rds` records the
#' source file's size and modification time alongside the matrix;
#' on the next call, if either differs, the source is reparsed and
#' the cache overwritten. A once-per-session `cli` line announces
#' "cache built" or "cache reused" so caching is never silent.
#' Silence with `options(rcisignal.silence_cache_messages = TRUE)`.
#'
#' @section Noise matrix vs signal matrix:
#' This function returns the **noise matrix** (`pixels x pool_size`,
#' the basis patterns from stimulus generation, *input* to CI
#' computation). The **signal matrix** ([read_signal_matrix()]) is a
#' different object: `pixels x participants`, one column per
#' producer's CI after base subtraction, *output* of CI computation
#' and *input* to reliability metrics. Both readers are exported
#' under `read_*` names and operate independently.
#'
#' @section What is in an rcicr `.RData`:
#' `load("rcic_stimuli.Rdata")` adds the following objects:
#' \describe{
#'   \item{`base_face_files`}{Named list of file paths to the
#'     original face images. List names are the labels passed
#'     downstream as `baseimage = "..."`.}
#'   \item{`base_faces`}{Base images themselves, loaded as numeric
#'     matrices of grayscale pixels in `[0, 1]`.}
#'   \item{`img_size`}{Side length of the (square) images in pixels.}
#'   \item{`n_trials`}{Number of stimulus pairs generated.}
#'   \item{`p`}{The noise basis. List with `$patches` (sinusoidal
#'     dictionary) and `$patchIdx` (index map).}
#'   \item{`stimuli_params`}{Per-trial recipe: list of matrices, one
#'     per base label. Each row is one trial's contrast weights.}
#'   \item{`seed`, `label`, `stimulus_path`, `trial`,
#'     `generator_version`, `use_same_parameters`}{Bookkeeping.}
#'   \item{`reference_norms`}{Random-responder Frobenius norms,
#'     inserted in place by `rcicr::computeInfoVal2IFC()`. Not
#'     present at first.}
#' }
#' Trial-level noise images are not stored in the rdata; this
#' loader reconstructs them on demand.
#'
#' @param path Path to the noise-matrix source (any of the formats
#'   listed in `Description`).
#' @param baseimage For rcicr `.RData` inputs: which base label to
#'   reconstruct noise for. Defaults to the only label if the rdata
#'   contains exactly one; aborts with a list of options otherwise.
#' @param stimuli_object For `.RData` files that contain a pre-saved
#'   `stimuli` matrix object (rare; rcicr does not save one by
#'   default), the object name to look up. Defaults to `"stimuli"`.
#' @param cache Logical. When `TRUE` (default), slow sources
#'   (text, rdata) are materialised to a sibling `.rds`. Set to
#'   `FALSE` for read-only filesystems or when you want to force a
#'   fresh parse without writing.
#' @param cache_path Optional explicit cache path. Default `NULL`
#'   places the cache next to the source (`<source>.rds`).
#' @param header Logical. When the source is a delimited text file,
#'   forwarded to [data.table::fread()]. Default `FALSE` (matches the
#'   Schmitz et al. (2024) Brief-RC noise-matrix text convention).
#'   Ignored for `.rds` and `.RData` sources.
#' @return A numeric `pixels x pool_size` matrix.
#' @seealso [read_signal_matrix()], [validate_noise_matrix()],
#'   [ci_from_responses_briefrc()]
#' @export
#' @examples
#' \dontrun{
#' # Plain text (Schmitz et al. 2024 OSF format).
#' # First call parses + writes data/noise_matrix.rds.
#' nm <- read_noise_matrix("data/noise_matrix.txt")
#' # Second call loads from the cache.
#' nm <- read_noise_matrix("data/noise_matrix.txt")
#'
#' # rcicr .RData (reconstructs each trial; slow first time, cached after).
#' nm <- read_noise_matrix("data/rcicr_stimuli.Rdata")
#' }
read_noise_matrix <- function(path,
                              baseimage      = NULL,
                              stimuli_object = "stimuli",
                              cache          = TRUE,
                              cache_path     = NULL,
                              header         = FALSE) {
  validate_path(path, "path")

  ext <- tolower(tools::file_ext(path))
  if (identical(ext, "rds")) {
    return(load_noise_rds_direct(path))
  }

  src_info <- file.info(path)
  cache_target <- resolve_noise_cache_path(path, cache_path)

  if (isTRUE(cache) && !is.null(cache_target) &&
        file.exists(cache_target)) {
    cached <- try_load_noise_cache(cache_target, src_info, path)
    if (!is.null(cached)) {
      inform_cache_reused(cache_target)
      return(cached)
    }
  }

  parsed <- parse_noise_source(path, ext, baseimage,
                               stimuli_object, header)
  matrix_out <- parsed$matrix

  if (isTRUE(cache) && !is.null(cache_target) &&
        !identical(parsed$source_kind, "fresh_rds")) {
    write_noise_cache(cache_target, matrix_out, path, src_info,
                      parsed$source_kind, baseimage)
    inform_cache_built(cache_target)
  }

  matrix_out
}

# ----- internal helpers ----------------------------------------------

#' @keywords internal
#' @noRd
load_noise_rds_direct <- function(path) {
  obj <- readRDS(path)
  if (is.list(obj) && identical(obj$.rcisignal_noise_cache, TRUE)) {
    return(coerce_to_noise_matrix(obj$matrix, path))
  }
  coerce_to_noise_matrix(obj, path)
}

#' @keywords internal
#' @noRd
coerce_to_noise_matrix <- function(x, path) {
  if (is.matrix(x) && is.numeric(x)) {
    out <- x
    dimnames(out) <- NULL
    storage.mode(out) <- "double"
    return(out)
  }
  if (is.data.frame(x)) {
    out <- as.matrix(x)
    if (!is.numeric(out)) {
      cli::cli_abort(c(
        "Object in {.path {path}} is not numeric after coercion.",
        "i" = "Got column classes: \\
               {.val {vapply(x, function(z) class(z)[1L], character(1L))}}"
      ))
    }
    dimnames(out) <- NULL
    storage.mode(out) <- "double"
    return(out)
  }
  cli::cli_abort(c(
    "Cannot interpret object in {.path {path}} as a noise matrix.",
    "i" = "Expected a numeric matrix or data.frame; \\
           got {.cls {class(x)}}."
  ))
}

#' @keywords internal
#' @noRd
resolve_noise_cache_path <- function(source_path, cache_path = NULL) {
  if (!is.null(cache_path)) return(cache_path)
  base <- tools::file_path_sans_ext(source_path)
  paste0(base, ".rds")
}

#' @keywords internal
#' @noRd
try_load_noise_cache <- function(cache_target, src_info, source_path) {
  obj <- tryCatch(readRDS(cache_target), error = function(e) NULL)
  if (is.null(obj)) return(NULL)
  if (!is.list(obj) || !identical(obj$.rcisignal_noise_cache, TRUE)) {
    return(NULL)
  }
  if (!identical(normalizePath(obj$source_path, mustWork = FALSE),
                 normalizePath(source_path, mustWork = FALSE))) {
    return(NULL)
  }
  if (!identical(as.numeric(obj$source_size),
                 as.numeric(src_info$size))) {
    return(NULL)
  }
  if (!identical(as.numeric(obj$source_mtime),
                 as.numeric(src_info$mtime))) {
    return(NULL)
  }
  coerce_to_noise_matrix(obj$matrix, cache_target)
}

#' @keywords internal
#' @noRd
write_noise_cache <- function(cache_target, mat, source_path,
                              src_info, source_kind,
                              baseimage = NULL) {
  payload <- list(
    .rcisignal_noise_cache = TRUE,
    rcisignal_version      = utils::packageVersion("rcisignal"),
    matrix                = mat,
    source_path           = normalizePath(source_path,
                                          mustWork = FALSE),
    source_size           = src_info$size,
    source_mtime          = src_info$mtime,
    source_kind           = source_kind,
    baseimage             = baseimage,
    cached_at             = Sys.time()
  )
  saveRDS(payload, cache_target)
  invisible(cache_target)
}

#' @keywords internal
#' @noRd
parse_noise_source <- function(path, ext, baseimage,
                               stimuli_object, header = FALSE) {
  if (ext %in% c("rdata", "RData")) {
    return(parse_noise_rdata(path, baseimage, stimuli_object))
  }
  if (ext %in% c("txt", "csv", "tsv", "dat")) {
    return(parse_noise_text(path, header))
  }
  # No extension or unknown: sniff the magic bytes.
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, "raw", n = 4L)
  is_gzip <- length(magic) >= 2L &&
    magic[1] == as.raw(0x1f) && magic[2] == as.raw(0x8b)
  is_rdx  <- length(magic) >= 4L &&
    rawToChar(magic[1:4]) %in%
      c("RDX2", "RDX3", "RDA2", "RDA3")
  if (is_rdx || is_gzip) {
    rds_obj <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!is.null(rds_obj)) {
      return(list(
        matrix      = coerce_to_noise_matrix(rds_obj, path),
        source_kind = "fresh_rds"
      ))
    }
    return(parse_noise_rdata(path, baseimage, stimuli_object))
  }
  parse_noise_text(path, header)
}

#' @keywords internal
#' @noRd
parse_noise_text <- function(path, header = FALSE) {
  tbl <- data.table::fread(path, header = header)
  if (ncol(tbl) == 0L) {
    cli::cli_abort(
      "Text source {.path {path}} parsed to zero columns."
    )
  }
  if (!all(vapply(tbl, is.numeric, logical(1L)))) {
    cli::cli_abort(c(
      "Text source {.path {path}} has non-numeric columns.",
      "i" = "Column classes: \\
             {.val {vapply(tbl, function(x) class(x)[1L], character(1L))}}"
    ))
  }
  mat <- as.matrix(tbl)
  dimnames(mat) <- NULL
  storage.mode(mat) <- "double"
  list(matrix = mat, source_kind = "text")
}

#' @keywords internal
#' @noRd
parse_noise_rdata <- function(path, baseimage, stimuli_object) {
  env <- new.env(parent = emptyenv())
  load_ok <- tryCatch({ load(path, envir = env); TRUE },
                     error = function(e) FALSE)
  if (!load_ok) {
    cli::cli_abort(
      "Could not {.fn load} {.path {path}} as an .RData file."
    )
  }
  objs <- ls(env)

  if (stimuli_object %in% objs) {
    mat <- as.matrix(env[[stimuli_object]])
    dimnames(mat) <- NULL
    storage.mode(mat) <- "double"
    return(list(matrix = mat, source_kind = "rdata_prebuilt"))
  }

  if (all(c("stimuli_params", "p", "img_size") %in% objs)) {
    if (!requireNamespace("rcicr", quietly = TRUE)) {
      cli::cli_abort(c(
        "Reconstructing the noise matrix from an rcicr rdata \\
         requires the {.pkg rcicr} package.",
        "i" = "Install: \\
               {.run remotes::install_github(\"rdotsch/rcicr\")}"
      ))
    }
    labels <- names(env$stimuli_params)
    if (is.null(baseimage)) {
      if (length(labels) != 1L) {
        cli::cli_abort(c(
          "Multiple base labels in rdata; pick one via \\
           {.arg baseimage}.",
          "i" = "Available: {.val {labels}}"
        ))
      }
      baseimage <- labels[1L]
    }
    if (!baseimage %in% labels) {
      cli::cli_abort(c(
        "{.arg baseimage} = {.val {baseimage}} not in rdata.",
        "i" = "Available: {.val {labels}}"
      ))
    }
    params  <- env$stimuli_params[[baseimage]]
    p_basis <- env$p
    img_sz  <- env$img_size
    n_trials <- nrow(params)
    n_pix    <- as.integer(img_sz) * as.integer(img_sz)
    mat <- matrix(NA_real_, nrow = n_pix, ncol = n_trials)
    for (i in seq_len(n_trials)) {
      mat[, i] <- as.vector(rcicr::generateNoiseImage(
        params[i, ], p_basis
      ))
    }
    return(list(matrix = mat, source_kind = "rdata_reconstructed"))
  }

  cli::cli_abort(c(
    "Unrecognised objects in {.path {path}}.",
    "i" = "Expected either {.var {stimuli_object}} (pre-saved \\
           matrix) or {.var stimuli_params} + {.var p} + \\
           {.var img_size} (rcicr rdata).",
    "*" = "Found: {.val {objs}}"
  ))
}
