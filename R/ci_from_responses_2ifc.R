#' Compute individual 2IFC CIs from trial-level responses
#'
#' @description
#' Thin wrapper around [rcicr::batchGenerateCI2IFC()] that handles
#' the known rcicr-integration gotchas and exposes a uniform return
#' shape across the 2IFC and Brief-RC paths.
#'
#' Use this when you have 2IFC trial-level responses and want the
#' package to compute individual CIs ready for the reliability
#' metrics.
#'
#' @details
#' What the wrapper does for you:
#' - Attaches `foreach`, `tibble`, `dplyr` at runtime
#'   (`%dopar%` / `tribble` / `%>%` are not namespace-prefixed
#'   inside rcicr).
#' - Matches the `.Rdata` extension case-insensitively (rcicr
#'   writes lowercase on some filesystems).
#' - Extracts the per-participant CI noise component (`$ci`) from
#'   the rcicr result list and stacks it into a pixels x
#'   participants signal matrix, already base-subtracted, ready
#'   for `rel_*()`.
#' - When `keep_rendered = TRUE`, also extracts the rendered
#'   `$combined` image and stacks it as `$rendered_ci` for
#'   visualization (display only, not for downstream stats).
#'
#' `rcicr` must be installed (it is a Suggests dep; install with
#' `remotes::install_github("rdotsch/rcicr")` if needed).
#'
#' @section Reading the result:
#' * `$signal_matrix` is the **raw mask** (rcicr's `$ci` per
#'   producer). This is what every `rel_*` function expects.
#' * `$rendered_ci`, present when `keep_rendered = TRUE`, is the
#'   `base + scaling(mask)` image rcicr would have written to PNG.
#'   **Visualization only.**
#' * `$rcicr_result` is the raw return value of
#'   [rcicr::batchGenerateCI2IFC()].
#'
#' @section Common mistakes:
#' * Passing `$rendered_ci` to `rel_*` functions; it carries the
#'   scaling step, which distorts variance-based metrics.
#' * Forgetting that the warning about Mode 1 also applies if you
#'   later read `$rendered_ci` back from PNG instead of using
#'   `$signal_matrix` directly.
#'
#' @param responses A data frame with one row per trial. Must
#'   contain `col_participant` (producer id), the column named by
#'   `col_stimulus` (stimulus id, integer, index into the rcicr
#'   noise pool), and `col_response` with values in `{-1, +1}`.
#' @param rdata_path Path to the `.Rdata` file produced by
#'   [rcicr::generateStimuli2IFC()]. Either `rdata_path` or
#'   `stimuli` must be supplied; if both are given `stimuli` wins.
#' @param stimuli In-memory stimuli list as returned in
#'   `$stimuli` by [simulate_2ifc_data()]. Use this in place of
#'   `rdata_path` when the simulation has been saved with
#'   [saveRDS()] and reloaded in a different R session: the path
#'   stored on `$rdata_path` no longer resolves, but `$stimuli`
#'   is self-contained. Internally the list is written to a
#'   fresh tempdir-backed `.Rdata` before the call into rcicr.
#' @param base_image Base face image. Three forms are accepted,
#'   matching [ci_from_responses_briefrc()]:
#'   * a numeric matrix in `[0, 1]` (used directly),
#'   * a single string path to a PNG / JPEG, or
#'   * a single string label naming an entry in the rdata's
#'     `base_faces` list (the historical 2IFC form).
#'   When `NULL`, the rdata's single base is used; if the rdata
#'   contains more than one base, the call aborts and lists the
#'   available labels. Matrix and path forms are injected into a
#'   temporary copy of the rdata under a synthetic label so the
#'   rcicr call sees the same structure it always has.
#' @param col_participant,col_stimulus,col_response Column names in
#'   `responses`.
#' @param scaling rcicr scaling option; one of `"autoscale"`,
#'   `"independent"`, `"constant"`, `"none"`. Passed through to
#'   `rcicr::batchGenerateCI2IFC()`. Scaling only affects the
#'   `$combined` (rendered) image, not `$ci`, so the returned
#'   `$signal_matrix` is the raw mask regardless of this argument.
#' @param keep_rendered If `TRUE`, also extract rcicr's `$combined`
#'   image (base + scaled noise) per producer and stack into
#'   `$rendered_ci`. Default `FALSE`. Visualization only.
#' @param targetpath Where rcicr writes PNGs. Defaults to an
#'   auto-cleaned tempdir so the working directory is not polluted.
#' @param save_as_png Whether rcicr writes individual CI PNGs.
#'   Defaults to `FALSE` for speed in a pure reliability pipeline.
#' @param group_by Optional character vector of column names in
#'   `responses` to group producers by. When supplied, the return
#'   list additionally contains `$group_ci`, a pixels x n_groups
#'   matrix built by calling [group_ci()] internally with the same
#'   `col_participant`. Single column (e.g. `"condition"`) gives
#'   one CI per level; length 2+ gives a factorial grouping with
#'   cell labels joined by `"_"`. Default `NULL` (no group CI
#'   built; only the per-producer `$signal_matrix` is returned).
#' @return A list with `signal_matrix`, optionally `rendered_ci`,
#'   `participants`, `img_dims`, `scaling`, `rcicr_result`, and
#'   optionally `group_ci` (when `group_by` is supplied).
#' @seealso [ci_from_responses_briefrc()],
#'   [rcicr::batchGenerateCI2IFC()],
#'   [group_ci()],
#'   [save_ci_images()]
#' @export
#' @examples
#' \dontrun{
#' # In your study, you replace the simulated inputs with your own:
#' #   - `sim$data` -> your responses data frame (CSV via read_responses())
#' #   - `sim$rdata_path` -> path to the rcicr stimuli .Rdata file used to
#' #     generate your stimuli.
#' # The `simulate_2ifc_data()` helper just produces the same shapes so the
#' # example below is runnable end-to-end.
#' sim <- simulate_2ifc_data(n_per_condition = 10, n_trials = 60, seed = 1)
#' res <- ci_from_responses_2ifc(sim$data, rdata_path = sim$rdata_path)
#' dim(res$signal_matrix)   # n_pixels x n_producers
#' rel_split_half(res$signal_matrix, n_permutations = 200L, seed = 1)
#'
#' # Same call, but using the portable in-memory stimuli list -- the
#' # form that survives saveRDS()/readRDS() across R sessions.
#' res2 <- ci_from_responses_2ifc(sim$data, stimuli = sim$stimuli)
#'
#' # Supply a base image as a path or matrix (same shapes as Brief-RC).
#' res3 <- ci_from_responses_2ifc(sim$data, rdata_path = sim$rdata_path,
#'                                base_image = sim$base_image_path)
#' res4 <- ci_from_responses_2ifc(sim$data, rdata_path = sim$rdata_path,
#'                                base_image = sim$base_face)
#' }
ci_from_responses_2ifc <- function(responses,
                                   rdata_path      = NULL,
                                   stimuli         = NULL,
                                   base_image      = NULL,
                                   col_participant = "participant_id",
                                   col_stimulus    = "stimulus",
                                   col_response    = "response",
                                   scaling         = "autoscale",
                                   keep_rendered   = FALSE,
                                   targetpath      = tempfile("rcisignal_2ifc_"),
                                   save_as_png     = FALSE,
                                   group_by        = NULL) {
  rdata_path <- resolve_rdata_input(rdata_path, stimuli,
                                    method = "2ifc")
  if (!requireNamespace("rcicr", quietly = TRUE)) {
    cli::cli_abort(c(
      "The 2IFC CI path requires the {.pkg rcicr} package.",
      "i" = "Install via: \\
             {.run remotes::install_github(\"rdotsch/rcicr\")}"
    ))
  }

  responses <- as.data.frame(responses)
  required <- c(col_participant, col_stimulus, col_response)
  missing_cols <- setdiff(required, colnames(responses))
  if (length(missing_cols) > 0L) {
    n_missing <- length(missing_cols)
    cli::cli_abort(c(
      "Missing {n_missing} column{?s} in {.arg responses}:",
      "*" = "{.val {missing_cols}}",
      "i" = "Have: {.val {colnames(responses)}}"
    ))
  }

  unique_resp <- sort(unique(as.numeric(responses[[col_response]])))
  if (!identical(unique_resp, c(-1, 1))) {
    msg <- c(
      "Column {.var {col_response}} must contain only values in \\
       {.val {c(-1, 1)}}.",
      "*" = "Got: {.val {unique_resp}}"
    )
    if (identical(unique_resp, c(0, 1)) ||
          identical(unique_resp, c(0, 1L))) {
      msg <- c(
        msg,
        "i" = "Did you mean {.code -1 / +1}? The {.code {{0, 1}}} \\
               coding is the most common silent failure in 2IFC \\
               pipelines (often produced by experiment software \\
               that records 'left' / 'right' as 0 / 1).",
        "i" = "Recode in one line: \\
               {.code responses${col_response} <- \\
                       2 * responses${col_response} - 1}"
      )
    }
    cli::cli_abort(msg)
  }

  resolved   <- resolve_2ifc_base_image(base_image, rdata_path)
  baseimage  <- resolved$label
  rdata_path <- resolved$rdata_path
  img_dims   <- resolved$img_dims

  ensure_attached(c("foreach", "tibble", "dplyr"))

  dir.create(targetpath, recursive = TRUE, showWarnings = FALSE)

  cis <- rcicr::batchGenerateCI2IFC(
    data        = responses,
    by          = col_participant,
    stimuli     = col_stimulus,
    responses   = col_response,
    baseimage   = baseimage,
    rdata       = rdata_path,
    save_as_png = save_as_png,
    targetpath  = targetpath,
    scaling     = scaling
  )

  participants <- unique(as.character(responses[[col_participant]]))
  n_pix <- prod(img_dims)
  signal_matrix <- matrix(
    NA_real_,
    nrow = n_pix,
    ncol = length(participants),
    dimnames = list(NULL, participants)
  )
  rendered <- if (isTRUE(keep_rendered)) {
    matrix(NA_real_, nrow = n_pix, ncol = length(participants),
           dimnames = list(NULL, participants))
  } else NULL

  for (pid in participants) {
    key <- grep(paste0(col_participant, "_", pid, "$"),
                names(cis), value = TRUE)
    if (length(key) == 0L) {
      cli::cli_abort(c(
        "rcicr produced no CI for producer {.val {pid}}.",
        "i" = "batchGenerateCI2IFC keys: {.val {names(cis)}}"
      ))
    }
    ci_mat <- cis[[key[1L]]]$ci
    if (is.null(ci_mat)) {
      cli::cli_abort(
        "rcicr CI element missing {.var $ci} component."
      )
    }
    signal_matrix[, pid] <- as.vector(ci_mat)
    if (!is.null(rendered)) {
      combined <- cis[[key[1L]]]$combined
      if (is.null(combined)) {
        cli::cli_abort(
          "rcicr CI element missing {.var $combined} component, \\
           cannot fill {.var $rendered_ci}."
        )
      }
      rendered[, pid] <- as.vector(combined)
    }
  }
  attr(signal_matrix, "img_dims") <- img_dims
  attr(signal_matrix, "source")   <- "raw"
  attr(signal_matrix, "ci_level") <- "individual"
  if (!is.null(rendered)) {
    attr(rendered, "img_dims") <- img_dims
    attr(rendered, "source")   <- "rendered"
  }

  out <- list(
    signal_matrix = signal_matrix,
    participants  = participants,
    img_dims      = img_dims,
    scaling       = scaling,
    rcicr_result  = cis
  )
  if (!is.null(rendered)) out$rendered_ci <- rendered

  if (!is.null(group_by)) {
    out$group_ci <- group_ci(
      signal_matrix   = signal_matrix,
      responses       = responses,
      by              = group_by,
      col_participant = col_participant
    )
  }

  out
}

#' Resolve the `base_image` argument for 2IFC functions
#'
#' Accepts a numeric matrix, an image path, a label string, or NULL.
#' For matrix / path forms, injects the resolved base into a fresh
#' temp copy of the rdata under a synthetic label so the rcicr call
#' downstream sees the same `base_faces` / `stimuli_params` layout it
#' always has. Returns the label and (possibly new) rdata path
#' downstream code should use, plus the resolved image dims.
#'
#' Disambiguation rule for a character input: if it ends in an image
#' extension (`.png`, `.jpe?g`, `.tiff?`, `.bmp`) it is treated as a
#' path; otherwise as a label lookup in `env$base_faces`.
#'
#' @keywords internal
#' @noRd
resolve_2ifc_base_image <- function(base_image, rdata_path) {
  env <- new.env(parent = emptyenv())
  load(rdata_path, envir = env)
  if (!"base_faces" %in% ls(env)) {
    cli::cli_abort(c(
      "No {.var base_faces} in {.path {rdata_path}}.",
      "i" = "Expected the output of {.fn rcicr::generateStimuli2IFC}."
    ))
  }
  labels <- names(env$base_faces)

  if (is.null(base_image)) {
    if (length(labels) != 1L) {
      cli::cli_abort(c(
        "Multiple base images in rdata; pick one via \\
         {.arg base_image}.",
        "i" = "Available: {.val {labels}}"
      ))
    }
    label <- labels[1L]
    img_dims <- as.integer(dim(env$base_faces[[label]]))
    return(list(label = label, rdata_path = rdata_path,
                img_dims = img_dims))
  }

  if (is.character(base_image) && length(base_image) == 1L) {
    is_path <- grepl("\\.(png|jpe?g|tiff?|bmp)$", base_image,
                     ignore.case = TRUE)
    if (!is_path) {
      if (!base_image %in% labels) {
        cli::cli_abort(c(
          "{.arg base_image} = {.val {base_image}} not found in rdata.",
          "i" = "Available labels: {.val {labels}}",
          "i" = "To supply an image file, give a path ending in \\
                 {.val .png} / {.val .jpg} / {.val .jpeg}."
        ))
      }
      img_dims <- as.integer(dim(env$base_faces[[base_image]]))
      return(list(label = base_image, rdata_path = rdata_path,
                  img_dims = img_dims))
    }
    base_mat <- read_image_as_gray(base_image)
  } else if (is.matrix(base_image) && is.numeric(base_image)) {
    if (any(!is.finite(base_image) & !is.na(base_image))) {
      cli::cli_abort("{.arg base_image} matrix contains non-finite values.")
    }
    if (!all(base_image >= 0 & base_image <= 1, na.rm = TRUE)) {
      cli::cli_abort(c(
        "{.arg base_image} matrix values must lie in {.code [0, 1]}.",
        "*" = "Range observed: \\
               [{min(base_image, na.rm = TRUE)}, \\
               {max(base_image, na.rm = TRUE)}]"
      ))
    }
    base_mat <- base_image
  } else {
    cli::cli_abort(c(
      "{.arg base_image} must be one of:",
      "*" = "{.code NULL} (use the rdata's single base);",
      "*" = "a numeric matrix in {.code [0, 1]};",
      "*" = "a single string path to a PNG / JPEG; or",
      "*" = "a single string label naming an entry in \\
             {.code base_faces}."
    ))
  }

  primary <- labels[1L]
  primary_dims <- as.integer(dim(env$base_faces[[primary]]))
  supplied_dims <- as.integer(dim(base_mat))
  if (!identical(primary_dims, supplied_dims)) {
    cli::cli_abort(c(
      "Supplied {.arg base_image} dimensions \\
       {.val {supplied_dims}} do not match the rdata's base \\
       {.val {primary_dims}}.",
      "i" = "The noise pool was generated for the rdata's image \\
             size and is not interchangeable."
    ))
  }

  synth_label <- "_rcisignal_supplied_base"
  env$base_faces[[synth_label]] <- base_mat
  if ("stimuli_params" %in% ls(env) &&
        is.list(env$stimuli_params) &&
        primary %in% names(env$stimuli_params)) {
    env$stimuli_params[[synth_label]] <- env$stimuli_params[[primary]]
  }
  if ("base_face_files" %in% ls(env) &&
        is.list(env$base_face_files) &&
        primary %in% names(env$base_face_files)) {
    env$base_face_files[[synth_label]] <- env$base_face_files[[primary]]
  }

  tmp_rdata <- tempfile("rcisignal_2ifc_basesub_", fileext = ".Rdata")
  save(list = ls(env), envir = env, file = tmp_rdata)
  list(label = synth_label, rdata_path = tmp_rdata,
       img_dims = supplied_dims)
}
