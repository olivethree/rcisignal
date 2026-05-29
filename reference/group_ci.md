# Collapse a per-producer signal matrix into per-group means

Convenience helper that takes a per-producer `signal_matrix` (pixels x
n_producers, the object returned by
[`ci_from_responses_briefrc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_briefrc.md)
/
[`ci_from_responses_2ifc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_2ifc.md))
and a trial-level `responses` data frame, and collapses producers into
per-group means via [`rowMeans()`](https://rdrr.io/r/base/colSums.html).

Uses the same "data frame plus column name" idiom as every other
responses-consuming function in the package: pass `responses` plus the
name of the column(s) you want to group by. Producer-to-group alignment
happens internally via `colnames(signal_matrix)`.

`ci_from_responses_*()` will call this for you if you pass them a
`group_by =` argument, so most users do not need to call `group_ci()`
directly. Reach for it when you want to build group CIs from a signal
matrix you already have in hand (e.g. one read back from disk).

## Usage

``` r
group_ci(
  signal_matrix,
  responses,
  by,
  col_participant = "participant_id",
  drop = TRUE
)
```

## Arguments

- signal_matrix:

  Numeric matrix of pixels x n_producers, as returned by
  [`ci_from_responses_briefrc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_briefrc.md)
  or
  [`ci_from_responses_2ifc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_2ifc.md)
  in their `$signal_matrix` field. Must have non-empty column names (the
  producer ids that map into `responses[[col_participant]]`).

- responses:

  Trial-level data frame containing one row per trial. Must contain the
  column named by `col_participant` (mapping to
  `colnames(signal_matrix)`) and every column named in `by`. Each
  producer's `by` value(s) must be consistent across all of their rows;
  an inconsistency aborts with a teaching message naming the offending
  producer.

- by:

  Character vector of column names in `responses`. Length 1 selects a
  single grouping column (e.g. `by = "condition"`). Length 2+ produces a
  factorial grouping where cell labels are the levels joined with `"_"`
  in the given order (e.g. `by = c("condition", "sex")` yields cells
  like `"happy_F"`).

- col_participant:

  Name of the participant-id column in `responses`. Default
  `"participant_id"`, matching the rest of the package.

- drop:

  Logical. When `TRUE` (default), empty cells are dropped from the
  output. When `FALSE`, empty cells are present as `NA` columns with
  `n = 0L`.

## Value

A numeric matrix of pixels x n_groups. Carries:

- column names = group labels;

- `attr(., "n")` = named integer vector of per-group producer counts;

- `attr(., "img_dims")` = inherited from
  `attr(signal_matrix, "img_dims")` if present;

- `attr(., "by_name")` = text description of the grouping;

- `attr(., "ci_level") = "group"` = a one-string hint used by
  [`save_ci_images()`](https://olivethree.github.io/rcisignal/reference/save_ci_images.md)
  to pick the filename prefix. Not a class; no S3 dispatch attached.

## Details

For each group, the corresponding output column is
`rowMeans(signal_matrix[, producers_in_group, drop = FALSE])`. For a
factorial `by` (length 2+), the cell label is the levels joined by `"_"`
in the column order given.

## See also

[`ci_from_responses_briefrc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_briefrc.md),
[`ci_from_responses_2ifc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_2ifc.md)
(their `group_by =` argument calls this internally);
[`plot_ci_distance_matrix()`](https://olivethree.github.io/rcisignal/reference/plot_ci_distance_matrix.md),
[`plot_ci_mds()`](https://olivethree.github.io/rcisignal/reference/plot_ci_mds.md),
[`plot_ci_correlogram()`](https://olivethree.github.io/rcisignal/reference/plot_ci_correlogram.md)
(compare multiple CIs, accept either per-producer or group-level
matrices);
[`save_ci_images()`](https://olivethree.github.io/rcisignal/reference/save_ci_images.md).

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_briefrc_data(
  n_per_condition = 10, n_trials = 60,
  conditions = c("A", "B"), seed = 1
)
cis  <- ci_from_responses_briefrc(sim$data,
                                  noise_matrix = sim$noise_matrix)
gcis <- group_ci(cis$signal_matrix, sim$data, by = "condition")
gcis                              # n_pixels x 2 (A and B)
attr(gcis, "n")                   # per-group producer counts

# Or skip the standalone call: pass group_by to the generator and
# get both per-producer and group matrices in one shot.
both <- ci_from_responses_briefrc(sim$data,
                                  noise_matrix = sim$noise_matrix,
                                  group_by     = "condition")
identical(both$group_ci, gcis)
} # }
```
