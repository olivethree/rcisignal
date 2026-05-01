# rcisignal: Quality Checks for Reverse-Correlation Data and Classification Images

A consolidated toolkit for assessing the quality of reverse correlation
(RC) experiments end-to-end. Input-side diagnostics
([`run_diagnostics()`](https://olivethree.github.io/rcisignal/reference/run_diagnostics.md)
and the `check_*` family) catch silent data-processing errors before
classification image (CI) computation. Output-side analyses
([`run_reliability()`](https://olivethree.github.io/rcisignal/reference/run_reliability.md),
[`run_discriminability()`](https://olivethree.github.io/rcisignal/reference/run_discriminability.md),
[`infoval()`](https://olivethree.github.io/rcisignal/reference/infoval.md),
[`agreement_map_test()`](https://olivethree.github.io/rcisignal/reference/agreement_map_test.md))
quantify CI quality after computation.

Works with both the standard 2IFC pipeline (via the upstream `rcicr`
package) and the Brief-RC pipeline (Schmitz, Rougier & Yzerbyt, 2024 —
implemented natively here). Operates directly on the pixel-level signal
produced by the original producers, so it does not depend on a
second-phase trait-rating study.

## See also

Useful links:

- <https://github.com/olivethree/rcisignal>

- <https://olivethree.github.io/rcisignal/>

- Report bugs at <https://github.com/olivethree/rcisignal/issues>

## Author

**Maintainer**: Manuel Oliveira <m.j.barbosa.de.oliveira@tue.nl>
([ORCID](https://orcid.org/0000-0002-6220-0695))
