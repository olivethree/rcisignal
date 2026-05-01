# Construct an `rcisignal_diag_report` object

Collects the outputs of multiple diagnostic checks into a single
printable summary, with the method that was run and the names of checks
that were skipped.

## Usage

``` r
rcisignal_diag_report(
  results,
  skipped_checks = character(),
  method = c("2ifc", "briefrc")
)
```

## Arguments

- results:

  A named list of
  [`rcisignal_diag_result()`](https://olivethree.github.io/rcisignal/reference/rcisignal_diag_result.md)
  objects.

- skipped_checks:

  Character vector of check names that were not executed. Defaults to
  empty.

- method:

  `"2ifc"` or `"briefrc"`.

## Value

An object of class `"rcisignal_diag_report"`.

## Examples

``` r
r <- rcisignal_diag_report(
  results = list(
    a = rcisignal_diag_result("pass", "Check A", "Looks fine.")
  ),
  method = "2ifc"
)
print(r)
#> == Data-quality report (2ifc) ==
#> 
#> [PASS] Check A
#>   Looks fine.
#> 
#> Summary: pass=1, warn=0, fail=0, skip=0
```
