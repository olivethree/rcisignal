# rcisignal 0.1.1

## New features

* `ci_from_responses_briefrc()` now accepts `method = "briefrc20"`
  in addition to `"briefrc12"`. Both Brief-RC variants are
  validated in Schmitz, Rougier, & Yzerbyt (2024). The CI
  computation is identical for both (the `genMask()` formula is
  symmetric in the per-trial split); the `method` argument is
  recorded as `$method` in the result list as provenance metadata.

## Documentation

* Vignette §13.1 rewritten to document Brief-RC 12 and Brief-RC 20
  as the two validated variants. Removed the inaccurate reference
  to a Brief-RC 6 variant (Schmitz et al. mention 4 / 6 / 8 / 10
  only as future research directions, not as published variants).
* Vignette §15 paragraph 5 (group-mean infoVal interpretation)
  rewritten to remove an unsourced "5-10x per-producer median"
  claim and to flag that the sqrt(N)-style inflation of group-mean
  z is conditional on producer-level signal alignment, not
  structural. Brinkman et al. (2019) numerical claims and
  recommendations are retained with explicit page citations.
* Vignette §15 paragraph 1 reframed to drop unsourced
  "10-30% / 70-90%" pixel-fraction figures while keeping the
  Frobenius-norm dilution mechanism.
* README validation section trimmed to a brief pointer; full
  validated-vs-unvalidated breakdown now lives in vignette §1.2.
* New educational background paragraphs added to §7 (split-half /
  Spearman-Brown, ICC), §8 (multiple-comparisons problem, cluster
  permutation, k-connectivity, TFCE, Frobenius/Euclidean
  dissimilarity), and §9 (per-producer and group-mean infoVal).

# rcisignal 0.1.0

First release of the package. Provides a consolidated toolkit for
quality assessment of reverse-correlation data and classification
images: input-side diagnostics (response coding, RT, alignment,
balance) and output-side reliability (pixel-level reliability,
cluster inference, per-producer infoVal).
