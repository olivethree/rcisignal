# Construct an `rcisignal_diag_result` object

Every diagnostic check function returns an object of class
`"rcisignal_diag_result"`. This constructor is the single source of
truth for the result shape.

## Usage

``` r
rcisignal_diag_result(status, label, detail, data = list())
```

## Arguments

- status:

  One of `"pass"`, `"warn"`, `"fail"`, or `"skip"`.

- label:

  Short (one-line) description of the check.

- detail:

  Human-readable explanation. Multiple lines allowed.

- data:

  Optional list of supporting data (data frames, flagged participant
  ids, summary statistics). Defaults to empty.

## Value

An object of class `"rcisignal_diag_result"`: a list with elements
`status`, `label`, `detail`, and `data`.

## Examples

``` r
rcisignal_diag_result("pass", "Response coding",
                     "All responses coded {-1, 1}.")
#> [PASS] Response coding
#>   All responses coded {-1, 1}.
```
