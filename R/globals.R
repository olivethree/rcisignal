## Package-level declarations.
##
## Both lines below are required for a package that uses
## `[.data.table` syntax but namespaces `data.table::` rather than
## `import(data.table)`-ing it:
##
## - `utils::globalVariables()` silences R CMD check NOTEs about
##   unbound sentinels used inside data.table calls.
## - `.datatable.aware <- TRUE` flips the runtime flag that makes
##   `cedta` recognize this package as data.table-aware, so
##   `[.data.table` dispatches correctly instead of falling through
##   to `[.data.frame`.
##
## Do NOT add short common names like "n" here — that silences real
## typos. Aggregations producing such columns must be renamed at the
## site (e.g. `n_obs`, `n_trials`).

utils::globalVariables(c(".N", ".SD", "n_obs", "response"))

.datatable.aware <- TRUE
