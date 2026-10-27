# sppt

Spatial point pattern test for areal data.

`sppt` compares spatial distributions of events aggregated to areal units,
such as crimes counted by census tract.

* `sppt()` compares two patterns unit by unit with bounded Monte Carlo
  bootstrap p-values, and summarises the map with the robust S index, a
  sameness test and a test of the spatial arrangement of the differences.
  Flagged units are ranked by confidence tier.
* `sppt_clusters()` tests whether clusters of flagged units are larger than
  false positives alone would produce.
* `sppt_longitudinal()` compares every period of a time series with the
  first and classifies each unit's trajectory.
* `sppt_sameness()` tests an S index against sameness.
* `sppt_map()` redraws the maps from either analysis with a chosen layout.

## Installation

```r
install.packages("sppt")
```

## Example

```r
library(sppt)

res <- sppt(vancouver_ct, group_col = "id",
            count_col = c("BNER2003", "BNER2025"), seed = 42)
summary(res)
plot(res, "pvalue", legend_position = "right")
sppt_clusters(res)

long <- sppt_longitudinal(vancouver_ct, group_col = "id",
                          count_col = paste0("BNER", 2003:2025), seed = 42)
```

See `vignette("sppt")` and `vignette("sppt_longitudinal")`.

## References

Andresen, M.A. (2009). Testing for similarity in area-based spatial patterns:
A nonparametric Monte Carlo approach. *Applied Geography*, 29(3), 333-345.
<https://doi.org/10.1016/j.apgeog.2008.12.004>

Andresen, M.A. (2016). An area-based nonparametric spatial point pattern
test: The test, its applications, and the future. *Methodological
Innovations*, 9, Article 12. <https://doi.org/10.1177/2059799116630659>

Andresen, M.A. (2026). Testing longitudinal change in spatial distributions
and patterns: An extension of the spatial point pattern test. *Annals of the
American Association of Geographers*, in press.
<https://doi.org/10.1080/24694452.2026.2733304>
