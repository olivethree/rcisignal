<!-- README.md (rendered at https://olivethree.github.io/rcisignal/) -->

# rcisignal <img src="man/figures/logo.png" align="left" width="180" hspace="20"/>

An R toolkit for quality assessment of reverse correlation (RC)
classification images. It bundles per-face-region `infoval()` analyses
and cross-producer agreement maps across the face (both descriptive
and FWER-controlled inferential tests) into one workflow, for 2IFC
(Dotsch, 2016) and Brief-RC (Schmitz, Rougier & Yzerbyt, 2024)
designs.

<br clear="left"/>

<!-- badges: start -->
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.19961180-1682d4.svg)](https://doi.org/10.5281/zenodo.19961180)
[![R-CMD-check](https://github.com/olivethree/rcisignal/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/olivethree/rcisignal/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/olivethree/rcisignal/actions/workflows/pkgdown.yaml/badge.svg)](https://olivethree.github.io/rcisignal/)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

> Package is not on CRAN, distribution is GitHub-only.

**[User guide](https://olivethree.github.io/rcisignal/articles/rcisignal.html) · [Installation](#installation) · [Using rcisignal on a real data set](#using-rcisignal-on-a-real-data-set) · [Citation](#citation)**

## Installation

From GitHub:

``` r
# install.packages("remotes")
remotes::install_github("olivethree/rcisignal")
```

If you ran a 2IFC study, you will also need `rcicr` (used to compute
the individual CIs):

``` r
remotes::install_github("rdotsch/rcicr")
```

Brief-RC users can skip `rcicr`; the Brief-RC code is fully native to
`rcisignal`.

> **Re-install before each fresh analysis.** `rcisignal` is in an
> experimental stage and exported functions are still being refined.
> Re-running `remotes::install_github("olivethree/rcisignal")` at the
> start of an analysis session pulls the latest version; the
> [user guide](https://olivethree.github.io/rcisignal/articles/rcisignal.html)
> is kept in sync with new and updated functions.

## Using rcisignal on a real data set

A re-analysis of the open data from [Oliveira, Garcia-Marques,
Dotsch & Garcia-Marques (2019)](https://doi.org/10.1002/ejsp.2569),
2IFC reverse correlation on a 256 x 256 grayscale male base face,
20 producers per trait, 300 trials each.

### Cross-producer agreement maps (fire palette)

![Cross-producer agreement maps for four traits (Trustworthy,
Friendly, Competent, Dominant), fire palette, restricted to the
full-face oval](vignettes/figures/oliveira_2019/agreement_maps_fire.png)

For each trait, brighter pixels mark image regions where the 20
producers in the RC task agreed more strongly on the location of
the target representation in their perception. The colorbar shows
the *degree of agreement (|t|)* from a per-pixel one-sample test of
producer CIs against zero, computed within the full-face oval. Hot
regions indicate consistent across-producer structure regardless of
sign; pale regions indicate weak or inconsistent agreement. For the
directional view (where the signal sits above or below the base),
see the user guide.

### `infoval()` per face region

Per (trait, region) cell. Values are the **median producer z-score**
and (in parentheses) the **number of producers (out of 20) clearing
`z >= 1.96`**, using face-region masks and a
trial-count-matched reference distribution.

| Trait       | Full face    | Upper face   | Eyes         | Mouth        |
|-------------|--------------|--------------|--------------|--------------|
| Trustworthy | +0.50 (3/20) | +0.30 (2/20) | +0.50 (1/20) | +0.53 (3/20) |
| Friendly    | +0.97 (5/20) | +0.23 (2/20) | +0.34 (2/20) | +0.75 (4/20) |
| Competent   | +0.70 (3/20) | +0.21 (3/20) | +0.36 (2/20) | +0.25 (5/20) |
| Dominant    | +0.89 (6/20) | +0.37 (5/20) | +0.91 (3/20) | +0.38 (2/20) |

The [**user guide**](https://olivethree.github.io/rcisignal/articles/rcisignal.html)
walks through every exported function, the interpretation of these
maps and the per-region grid, and the full *Worked example:
Oliveira et al. (2019), Study 1* section.

## Validation status

Several metrics in this package are package-level extensions whose
behavior on social-face RC data has not been independently
validated. See [*§1.2 Validation
status*](https://olivethree.github.io/rcisignal/articles/rcisignal.html#validation-status)
in the user guide for the breakdown of validated versus unvalidated
metrics, and treat the unvalidated ones as exploratory in published
work. The
[*rcisignal vs rcicr*](https://olivethree.github.io/rcisignal/articles/validation_rcicr.html)
vignette has the engine-equivalence receipts for the per-producer
2IFC infoVal claim and a Brief-RC signal-recovery sanity check.

## Citation

If `rcisignal` helps your research, please cite it:

```
Oliveira, M. (2026). rcisignal: Quality checks for reverse-correlation
data and classification images (Version 0.3.0) [R package]. Zenodo.
https://doi.org/10.5281/zenodo.19961180
```

Run `citation("rcisignal")` in R for a BibTeX entry.

Please also cite the methodological sources appropriate to your
pipeline:

-   **2IFC**: Dotsch (2016, 2023) for the `rcicr` package; Brinkman et
    al. (2019) for infoVal.
-   **Brief-RC**: Schmitz, Rougier, and Yzerbyt (2024).

## References

- Brinkman, L., Goffin, S., van de Schoot, R., van Haren, N. E. M., Dotsch, R., & Aarts, H. (2019). Quantifying the informational value of classification images. *Behavior Research Methods*, *51*(5), 2059-2073. [https://doi.org/10.3758/s13428-019-01232-2](https://doi.org/10.3758/s13428-019-01232-2)
- Dotsch, R. (2016, 2023). *rcicr: Reverse-correlation image-classification toolbox* [R package]. [https://github.com/rdotsch/rcicr](https://github.com/rdotsch/rcicr)
- Oliveira, M., Garcia-Marques, T., Dotsch, R., & Garcia-Marques, L. (2019). Dominance and competence face to face: Dissociations obtained with a reverse correlation approach. *European Journal of Social Psychology*. [https://doi.org/10.1002/ejsp.2569](https://doi.org/10.1002/ejsp.2569)
- Schmitz, M., Rougier, M., & Yzerbyt, V. (2024). Introducing the brief reverse correlation: An improved tool to assess visual representations. *European Journal of Social Psychology*. [https://doi.org/10.1002/ejsp.3100](https://doi.org/10.1002/ejsp.3100)

## License

Released under the [MIT License](LICENSE.md).

## Credits

Designed by [Manuel Oliveira](https://manueloliveira.nl/) <a href="https://orcid.org/0000-0002-6220-0695"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD"/></a>

Code and documentation were co-built with Claude (Opus 4.6, Anthropic;
April-May 2026).

### Acknowledgements to the community

This package builds on the excellent foundational work by
Ron Dotsch, Loek Brinkman, Alex Todorov, Mathias Schmitz, Marine
Rougier, Vincent Yzerbyt, and their many collaborators across the
years. Reverse correlation is no longer a central part of my research,
but I still find a lot of enjoyment in working on these side projects.
The inspiration to build these tools and tutorials comes mostly from
occasional collaborations with my PhD supervisors (Teresa
Garcia-Marques, Leonel Garcia-Marques) and all the warm and competent
colleagues from the research groups in Lisbon (Goncalo Oliveira, Rui
Costa-Lopes and their teams) with whom I have been greatly enjoying
working together on this stuff. Hopefully this toolkit will come in
handy to all the RC research enthusiasts out there! :)
