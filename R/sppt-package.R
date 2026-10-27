#' sppt: Spatial Point Pattern Test for Areal Data
#'
#' Compares spatial distributions of events aggregated to areal units.
#'
#' * [sppt()] compares two patterns (a base and a test) unit by unit, using
#'   bounded Monte Carlo bootstrap p-values and, optionally, marginal
#'   confidence interval overlap.
#' * [sppt_longitudinal()] compares every period of a time series against the
#'   first period and classifies unit-level trajectories.
#' * [sppt_sameness()] tests a robust S index against sameness.
#' * [sppt_map()] redraws the maps from either analysis with a chosen layout.
#'
#' @references
#' Andresen, M.A. (2009). Testing for similarity in area-based spatial
#' patterns: A nonparametric Monte Carlo approach. *Applied Geography*,
#' 29(3), 333-345. \doi{10.1016/j.apgeog.2008.12.004}
#'
#' Andresen, M.A. (2016). An area-based nonparametric spatial point pattern
#' test: The test, its applications, and the future. *Methodological
#' Innovations*, 9, Article 12. \doi{10.1177/2059799116630659}
#'
#' Andresen, M.A. (2026). Testing longitudinal change in spatial
#' distributions and patterns: An extension of the spatial point pattern
#' test. *Annals of the American Association of Geographers*, in press.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom rlang .data :=
NULL

# Write console output only when the caller asked for it.
say_if <- function(verbose) {
  force(verbose)
  function(...) if (isTRUE(verbose)) cat(...)
}

# Directory check shared by every export option. Nothing is ever written
# unless the user names a directory.
check_export_dir <- function(export_dir, what) {
  if (is.null(export_dir) || !is.character(export_dir) ||
      length(export_dir) != 1L || is.na(export_dir))
    stop("`export_dir` must be supplied when `", what, " = TRUE`.",
         call. = FALSE)
  if (!dir.exists(export_dir))
    dir.create(export_dir, recursive = TRUE)
  export_dir
}

# Write an sf result to disk in one of the supported formats.
write_result <- function(result, dir, stem, format) {
  fmt <- tolower(format)
  fp <- switch(fmt,
    gpkg = , shp = file.path(dir, paste0(stem, ".", fmt)),
    csv = , txt = file.path(dir, paste0(stem, ".", fmt)),
    rds = file.path(dir, paste0(stem, ".rds")),
    stop("Unsupported `export_format` '", format,
         "'. Use one of: gpkg, shp, csv, txt, rds.", call. = FALSE)
  )
  plain <- result
  class(plain) <- setdiff(class(plain), c("sppt_result", "sppt_long_result"))
  if (fmt %in% c("gpkg", "shp")) {
    sf::st_write(plain, fp, delete_dsn = TRUE, quiet = TRUE)
  } else if (fmt %in% c("csv", "txt")) {
    utils::write.csv(sf::st_drop_geometry(plain), fp, row.names = FALSE)
  } else {
    saveRDS(result, fp)
  }
  fp
}
