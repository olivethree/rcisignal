# Compute only the empirical-null distribution for split-half

Standalone access to the null-distribution simulation used by
`rel_split_half(null = ...)`. Useful when the user wants to precompute
and reuse a null across multiple observed analyses (e.g., when comparing
across conditions with the same producer count).

## Usage

``` r
rel_split_half_null(
  n_participants,
  n_pixels,
  null = c("permutation", "random_responders"),
  noise_matrix = NULL,
  n_permutations = 2000L,
  seed = NULL,
  pid = NULL
)
```

## Arguments

- n_participants:

  Integer.

- n_pixels:

  Integer.

- null:

  Either `"permutation"` or `"random_responders"`.

- noise_matrix:

  Required for `"random_responders"`.

- n_permutations:

  Integer. Default 2000.

- seed:

  Optional integer.

- pid:

  Optional internal progress-bar id; end users should leave at `NULL`.

## Value

Numeric vector of length `n_permutations`.
