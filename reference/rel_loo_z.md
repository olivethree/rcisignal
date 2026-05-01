# Z-scored leave-one-out influence (accessor)

Convenience accessor that returns a data frame of producer ids and their
z-scored LOO influence, ordered from most-influential (lowest, most
negative `z_score`) to least. Accepts either a signal matrix (runs
[`rel_loo()`](https://olivethree.github.io/rcisignal/reference/rel_loo.md)
under the hood) or an existing `rcisignal_rel_loo` result object (cheap,
no recomputation).

## Usage

``` r
rel_loo_z(x, ...)
```

## Arguments

- x:

  Either a `pixels x participants` signal matrix or an object of class
  `rcisignal_rel_loo` (as returned by
  [`rel_loo()`](https://olivethree.github.io/rcisignal/reference/rel_loo.md)).

- ...:

  Passed to
  [`rel_loo()`](https://olivethree.github.io/rcisignal/reference/rel_loo.md)
  when `x` is a signal matrix (e.g. `flag_threshold`, `flag_method`).

## Value

A data frame with columns `participant_id`, `correlation`, `z_score`,
`flag`, sorted by `z_score` ascending.

## See also

[`rel_loo()`](https://olivethree.github.io/rcisignal/reference/rel_loo.md)
