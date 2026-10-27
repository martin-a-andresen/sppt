# sppt 1.1.0

## Changes to defaults

* `conf_level` now defaults to 0.85 in `sppt()` and `sppt_longitudinal()`
  (was 0.683). Two marginal 85% intervals fail to overlap when the
  difference exceeds about 2.04 standard errors, which corresponds to a
  two-sample test at about `alpha = 0.05`, as in Andresen (2026). Results
  with `calculate_ci = TRUE`, and all envelope results of
  `sppt_longitudinal()`, differ from version 1.0.0. P-value results are
  unchanged.

## New features

* `sppt()` assigns each flagged unit a confidence tier (`confidence_tier`):
  core (inside a local Moran cluster), edge (bordering a core unit with the
  same sign) or outside. Tier counts are printed and stored in the summary,
  and a new tier map is available as `sppt_map(res, "tiers")`.
* `sppt()` notes a widespread shift (most flagged units share one direction
  and lie outside any cluster) and adds it to the two-way verdict.
* New `sppt_clusters()` tests clusters of flagged units against the
  clusters that false positives alone produce on the same map, under the
  one-gap or contiguous rule.
* `plot()` methods for `sppt()` and `sppt_longitudinal()` results.
* `sppt_longitudinal()` results now have class `sppt_long_result`, with
  `summary()` and `print()` methods.
* `sppt()` gains `correction` for multiple testing corrections, kept for
  reproducing earlier analyses and not recommended. Raw p-values stay in
  `pvalue`, adjusted ones go in `pvalue_adj`, and the sameness test is not
  reported under a correction.

## Bug fixes

* Local Moran outlines now require the matching sign: a unit is outlined
  as "test > base" only if its difference is positive, and as
  "base > test" only if it is negative. Local Moran classifies around the
  mean of the signed surface, so in one-sided comparisons version 1.0.0
  could outline units under the wrong direction.

## Other changes

* Cluster outlines on the maps are drawn as one dissolved boundary per
  direction, with lighter unit borders, so clusters read at fine scales.
* The decision tier text follows the current guidance: `alpha = 0.01` for
  false discovery priority only at coarse scales (`lambda_min` above
  roughly 10 to 15).
* Help pages note where the implementation differs from the published
  descriptions of the method.

# sppt 1.0.0

* Initial CRAN release.
* `sppt()`: two-pattern test with bounded Monte Carlo bootstrap p-values,
  optional marginal confidence interval overlap, the robust S index with a
  sameness test, and an arrangement test on the signed evidence surface.
* `sppt_longitudinal()`: longitudinal extension with envelope and p-value
  stability, trend classification and Jenks trend tables.
* `sppt_sameness()`: sameness test for a robust S index.
* `sppt_map()`: redraw any map with a chosen legend, north arrow and scale
  bar layout.
* Data: `vancouver_ct` and `vancouver_da`, residential and commercial
  burglary in Vancouver, 2003-2025.
