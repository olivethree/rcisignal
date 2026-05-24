# Read a noise matrix from any supported source, with caching

Returns the `pixels x pool_size` numeric matrix of basis noise patterns
that
[`ci_from_responses_briefrc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_briefrc.md),
[`infoval()`](https://olivethree.github.io/rcisignal/reference/infoval.md),
and the diagnostic checks need. The function detects the input format
from the file extension (with a magic-byte fallback) and transparently
caches slow-to-parse formats to a sibling `.rds` so the second read is
near-instant.

Supported sources, with no `format` flag — detection is automatic:

- **`.rds`** (R serialized object). Loaded directly. This is the fast
  path; everything else gets cached to one of these.

- **`.txt` / `.csv` / `.tsv`** (delimited text matrix, one column per
  stimulus). Read via
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html);
  the parsed matrix is materialised to `<basename>.rds` next to the
  source for next time.

- **`.RData` / `.Rdata`** (rcicr stimulus generation output). The rdata
  does not store the noise matrix directly — it stores `stimuli_params`
  (per-trial parameter rows) and `p` (the sinusoid / gabor basis). This
  loader reconstructs each trial's noise pattern via
  [`rcicr::generateNoiseImage()`](https://rdrr.io/pkg/rcicr/man/generateNoiseImage.html)
  and caches the result to `<basename>.rds`. Requires `rcicr` to be
  installed for the first read; subsequent reads from the cache do not.

## Usage

``` r
read_noise_matrix(
  path,
  base_image = NULL,
  stimuli_object = "stimuli",
  cache = TRUE,
  cache_path = NULL,
  header = FALSE
)
```

## Arguments

- path:

  Path to the noise-matrix source (any of the formats listed in
  `Description`).

- base_image:

  For rcicr `.RData` inputs: which base label to reconstruct noise for.
  Defaults to the only label if the rdata contains exactly one; aborts
  with a list of options otherwise.

- stimuli_object:

  For `.RData` files that contain a pre-saved `stimuli` matrix object
  (rare; rcicr does not save one by default), the object name to look
  up. Defaults to `"stimuli"`.

- cache:

  Logical. When `TRUE` (default), slow sources (text, rdata) are
  materialised to a sibling `.rds`. Set to `FALSE` for read-only
  filesystems or when you want to force a fresh parse without writing.

- cache_path:

  Optional explicit cache path. Default `NULL` places the cache next to
  the source (`<source>.rds`).

- header:

  Logical. When the source is a delimited text file, forwarded to
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html).
  Default `FALSE` (matches the Schmitz et al. (2024) Brief-RC
  noise-matrix text convention). Ignored for `.rds` and `.RData`
  sources.

## Value

A numeric `pixels x pool_size` matrix.

## Details

The cache invalidates automatically. The `.rds` records the source
file's size and modification time alongside the matrix; on the next
call, if either differs, the source is reparsed and the cache
overwritten. A once-per-session `cli` line announces "cache built" or
"cache reused" so caching is never silent. Silence with
`options(rcisignal.silence_cache_messages = TRUE)`.

## Noise matrix vs signal matrix

This function returns the **noise matrix** (`pixels x pool_size`, the
basis patterns from stimulus generation, *input* to CI computation). The
**signal matrix**
([`read_signal_matrix()`](https://olivethree.github.io/rcisignal/reference/read_signal_matrix.md))
is a different object: `pixels x participants`, one column per
producer's CI after base subtraction, *output* of CI computation and
*input* to reliability metrics. Both readers are exported under `read_*`
names and operate independently.

## What is in an rcicr `.RData`

`load("rcic_stimuli.Rdata")` adds the following objects:

- `base_face_files`:

  Named list of file paths to the original face images. List names are
  the labels passed downstream as `base_image = "..."`.

- `base_faces`:

  Base images themselves, loaded as numeric matrices of grayscale pixels
  in `[0, 1]`.

- `img_size`:

  Side length of the (square) images in pixels.

- `n_trials`:

  Number of stimulus pairs generated.

- `p`:

  The noise basis. List with `$patches` (sinusoidal dictionary) and
  `$patchIdx` (index map).

- `stimuli_params`:

  Per-trial recipe: list of matrices, one per base label. Each row is
  one trial's contrast weights.

- `seed`, `label`, `stimulus_path`, `trial`, `generator_version`,
  `use_same_parameters`:

  Bookkeeping.

- `reference_norms`:

  Random-responder Frobenius norms, inserted in place by
  [`rcicr::computeInfoVal2IFC()`](https://rdrr.io/pkg/rcicr/man/computeInfoVal2IFC.html).
  Not present at first.

Trial-level noise images are not stored in the rdata; this loader
reconstructs them on demand.

## See also

[`read_signal_matrix()`](https://olivethree.github.io/rcisignal/reference/read_signal_matrix.md),
[`validate_noise_matrix()`](https://olivethree.github.io/rcisignal/reference/validate_noise_matrix.md),
[`ci_from_responses_briefrc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_briefrc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Plain text (Schmitz et al. 2024 OSF format).
# First call parses + writes data/noise_matrix.rds.
nm <- read_noise_matrix("data/noise_matrix.txt")
# Second call loads from the cache.
nm <- read_noise_matrix("data/noise_matrix.txt")

# rcicr .RData (reconstructs each trial; slow first time, cached after).
nm <- read_noise_matrix("data/rcicr_stimuli.Rdata")
} # }
```
