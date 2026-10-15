# Builds data/vancouver_ct.rda and data/vancouver_da.rda from the two layers
# of Vancouver_BNER_BNEC_CRAN.gpkg. Run from the package root.

library(sf)

gpkg <- file.path("..", "Data", "Vancouver_BNER_BNEC_CRAN.gpkg")

prepare <- function(layer) {
  x <- st_read(gpkg, layer = layer, quiet = TRUE)
  st_geometry(x) <- "geometry"
  counts <- grep("^BNE[RC][0-9]{4}$", names(x), value = TRUE)
  counts <- counts[order(substr(counts, 1, 4), substr(counts, 5, 8))]
  x$id <- seq_len(nrow(x))
  x <- x[, c("id", counts, "geometry")]
  for (cc in counts) x[[cc]] <- as.integer(x[[cc]])
  x
}

vancouver_ct <- prepare("CT_2021")
vancouver_da <- prepare("DA_2021")

usethis::use_data(vancouver_ct, vancouver_da, overwrite = TRUE,
                  compress = "xz")
