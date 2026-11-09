# sppt 1.1.0

## Changes to defaults

* `conf_level` now defaults to 0.85 in `sppt()` and `sppt_longitudinal()`
  (was 0.683). Two marginal 85% intervals fail to overlap when the
  difference exceeds about 2.04 standard errors, which corresponds to a
  two-sample test at about `alpha = 0.05`, as in Andresen (2026). Results
  with `calculate_ci = TRUE`, and all envelope results of
  `sppt_longitudinal()`, differ from version 1.0.0. P-value results are
  unchanged.
* `B` now defaults to `"auto"` in `sppt()` and `sppt_longitudinal()`. The
  number of replicates is the largest of a minimum of 4999, enough for a
  Monte Carlo standard error at `alpha` of at most 10% of `alpha`, and,
  under a correction, enough for the most significant unit to survive it.
  At `alpha = 0.05` without a correction this is 4999, as before, so default
  results are unchanged. An integer can still be given.

## New features

* `sppt()` assigns each flagged unit a confidence tier (`confidence_tier`):
  core (inside a local Moran cluster), edge (bordering a core unit with the
  same sign) or outside. Tier counts are printed and stored in the summary,
  and a new tier map is available as `sppt_map(res, "tiers")`.
* `sppt()` notes a widespread shift (most flagged units share one direction
  and lie outside any cluster) and adds it to the two-way verdict.
* `sppt()` reports a degree of similarity, S_D: one minus half the sum of
  absolute differences in shares over the flagged units, alongside the index
  of dissimilarity D (Duncan and Duncan 1955), the share of D held by
  flagged units and the number of events of each pattern that would have to
  move. A note is printed when flagged units hold less than 60% of D. S tests
  sameness; S_D measures how similar two patterns are. New
  `sppt_similarity()` returns these values from a result.
* `sppt()` tests a single pattern against a random distribution:
  `count_col = "BNER2021"` or `c("random", "BNER2021")` compares the column
  with random allocations of its own total to the units, in proportion to
  each unit's exposure: area for polygons, length for lines such as street
  segments (`random_by = "exposure"`, the default). `random_by = "units"`
  gives equal probability per unit; naming a numeric column uses it as the
  exposure.
  The test is a Monte Carlo test against the random envelope, and all the
  usual output follows (S and sameness, S_D and D, arrangement, tiers, maps).
* `sppt_longitudinal()` compares the first and last periods as two
  patterns: robust S with its sameness test (new argument `sameness_tau`),
  and S_D with D, printed in the report and `summary()` and stored as
  `attr(x, "first_last")`.
* New `sppt_clusters()` tests clusters of flagged units against the
  clusters that false positives alone produce on the same map, under the
  one-gap or contiguous rule. Clusters are built within each direction of
  difference by default; `direction = FALSE` reproduces the published
  reference values.
* `sppt()` gains `correction` for multiple testing corrections (methods of
  `p.adjust()`; `"BY"` controls the false discovery rate under any
  dependence). Raw p-values stay in `pvalue` and adjusted ones go in
  `pvalue_adj`; units are flagged on the adjusted values and the sameness
  test is not reported under a correction.
* Large numbers of replicates are computed in blocks, so memory no longer
  grows with `B` (without `calculate_ci`). Results are the same as the
  single-matrix calculation.
* `plot()` methods for `sppt()` and `sppt_longitudinal()` results.
* `sppt_longitudinal()` results now have class `sppt_long_result`, with
  `summary()` and `print()` methods.

## Bug fixes

* In `sppt_longitudinal()`, the p-value average index
  (`max_overlap_proportion_p`, and the S-index average built on it) and
  `prop_sig_p` are now shares of the T - 1 test periods. Version 1.0.0
  divided by all T periods, counting the base period, which is never
  compared with itself, as non-significant. The average index was therefore
  slightly too high, by the share of test periods that differ divided by T:
  with 23 periods, an average of 0.809 becomes about 0.800.

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
* Console reports and `summary()` output are wrapped to fit 70 columns, and
  the two-way verdict is printed in summaries.
* A failed arrangement test now gives a warning rather than a message, so
  missing tiers cannot pass unnoticed.
* `title = FALSE` in `sppt_map()` and `plot()` removes the subtitle as well
  as the title.

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
