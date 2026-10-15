#' Residential and commercial burglary in Vancouver, 2003-2025
#'
#' @description
#' Annual counts of residential burglary (break and enter, residential) and
#' commercial burglary (break and enter, commercial) in the City of
#' Vancouver, British Columbia, from 2003 to 2025, aggregated to 2021 Census
#' spatial units.
#'
#' * `vancouver_ct`: 127 census tracts.
#' * `vancouver_da`: 1,016 dissemination areas.
#'
#' Census tracts are coarse enough that most period comparisons have good
#' power, which makes `vancouver_ct` the better choice for examples and first
#' analyses. Dissemination areas show the fine-scale behaviour of the test,
#' including the low-power warnings.
#'
#' @format An `sf` object with polygon geometry in NAD83(CSRS) / UTM zone 10N
#'   and the following columns:
#' \describe{
#'   \item{`id`}{Integer. Sequential unit identifier, for use as
#'     `group_col`. It is not a Statistics Canada code.}
#'   \item{`BNER2003`, ..., `BNER2025`}{Integer. Residential burglaries in
#'     each year.}
#'   \item{`BNEC2003`, ..., `BNEC2025`}{Integer. Commercial burglaries in
#'     each year.}
#'   \item{`geometry`}{`sfc_MULTIPOLYGON` boundaries of the unit.}
#' }
#' `vancouver_ct` has 127 rows and `vancouver_da` has 1,016 rows.
#'
#' @source Crime counts: Vancouver Police Department open crime data,
#'   \url{https://geodash.vpd.ca/opendata/}. Boundaries: Statistics Canada,
#'   2021 Census boundary files,
#'   \url{https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/index2021-eng.cfm}.
#'   Events were aggregated to the 2021 units by the package author.
#'
#' @examples
#' vancouver_ct
#' colSums(sf::st_drop_geometry(vancouver_ct)[, c("BNER2003", "BNER2025")])
#' plot(vancouver_ct["BNER2025"])
#'
#' @name vancouver
#' @aliases vancouver_ct vancouver_da
#' @keywords datasets
"vancouver_ct"

#' @rdname vancouver
"vancouver_da"
