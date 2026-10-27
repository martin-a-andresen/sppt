run_sppt <- function(...) {
  sppt(vancouver_ct, group_col = "id", count_col = c("BNER2003", "BNER2025"),
       B = 199, seed = 1, verbose = FALSE, create_maps = FALSE, ...)
}

test_that("p-values are bounded below by 1 / (B + 1)", {
  res <- run_sppt(arrangement = FALSE)
  expect_s3_class(res, "sppt_result")
  expect_s3_class(res, "sf")
  p <- res$pvalue[!is.na(res$pvalue)]
  expect_true(all(p >= 1 / 200 & p <= 1))
})

test_that("results are reproducible with a seed", {
  a <- run_sppt(arrangement = FALSE)
  b <- run_sppt(arrangement = FALSE)
  expect_identical(a$pvalue, b$pvalue)
})

test_that("summary works with and without the arrangement test", {
  expect_output(summary(run_sppt(arrangement = FALSE)), "Robust S")
  expect_output(summary(run_sppt()), "Arrangement")
})

test_that("exporting without a directory is an error", {
  expect_error(run_sppt(arrangement = FALSE, export_results = TRUE),
               "export_dir")
})

test_that("inputs are validated", {
  expect_error(sppt(vancouver_ct, group_col = "id", count_col = "BNER2003"),
               "length 2")
  expect_error(sppt(vancouver_ct, group_col = "nope",
                    count_col = c("BNER2003", "BNER2025")), "not found")
  expect_error(run_sppt(alpha = 1.5), "alpha")
})

test_that("sppt_sameness centres S at 1 - alpha", {
  s <- sppt_sameness(s_index = 0.95, n_sig = 10, n_eligible = 200,
                     alpha = 0.05)
  expect_equal(s$s_expected, 0.95)
  expect_true(s$same)
  expect_false(sppt_sameness(0.70, 60, 200, 0.05)$same)
})

test_that("sppt_map redraws maps with a new layout", {
  res <- run_sppt(arrangement = FALSE)
  p <- sppt_map(res, "pvalue", legend_position = "right", north_arrow = "tl")
  expect_s3_class(p, "ggplot")
  expect_length(sppt_map(res), 2L)
  expect_error(sppt_map(res, "ci"), "not available")
  expect_error(sppt_map(res, north_arrow = "middle"), "north_arrow")
})

test_that("confidence tiers cover exactly the flagged units", {
  res <- run_sppt()
  flagged <- !is.na(res$sig_pvalue_dir) & res$sig_pvalue_dir != 0L
  expect_identical(!is.na(res$confidence_tier), flagged)
  expect_true(all(res$confidence_tier[flagged] %in%
                    c("core", "edge", "outside")))
  s <- attr(res, "sppt_summary")$arrangement
  expect_equal(s$n_core + s$n_edge + s$n_outside, sum(flagged))
})

test_that("cluster outlines carry the sign of their label", {
  res <- run_sppt()
  expect_true(all(res$difference[res$moran_hh] > 0))
  expect_true(all(res$difference[res$moran_ll] < 0))
})

test_that("sppt() has no correction argument", {
  expect_false("correction" %in% names(formals(sppt)))
})

test_that("conf_level defaults to 0.85", {
  expect_identical(formals(sppt)$conf_level, 0.85)
  res <- run_sppt(arrangement = FALSE, calculate_ci = TRUE)
  expect_equal(attr(res, "sppt_summary")$conf_level, 0.85)
})

test_that("plot draws the tier map", {
  res <- run_sppt()
  p <- NULL
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  p <- plot(res, "tiers")
  expect_s3_class(p, "ggplot")
  expect_named(sppt_map(res), c("difference", "pvalue", "tiers"))
})

test_that("sppt_clusters labels clusters consistently", {
  res <- run_sppt()
  cl <- sppt_clusters(res, direction = FALSE, nsim = 99, seed = 1)
  expect_s3_class(cl, "sppt_clusters")
  flagged <- !is.na(res$sig_pvalue_dir) & res$sig_pvalue_dir != 0L
  expect_identical(!is.na(cl$membership), flagged)
  expect_equal(as.integer(table(cl$membership)), cl$clusters$size)
  expect_false(is.unsorted(rev(cl$clusters$size)))
  expect_false(is.unsorted(rev(cl$ceiling$min_size)))
  expect_output(print(cl), "ceiling")
  # Clusters within direction are never larger than undirected ones
  cd <- sppt_clusters(res, nsim = 99, seed = 1)
  expect_true(all(cd$clusters$direction != "mixed"))
  expect_lte(max(cd$clusters$size), max(cl$clusters$size))
})

test_that("console output fits in 70 columns", {
  wide <- function(x) max(nchar(strsplit(paste(x, collapse = "\n"), "\n")[[1]]))
  out <- capture.output(res <- suppressWarnings(sppt(
    vancouver_da, group_col = "id", count_col = c("BNER2021", "BNEC2021"),
    B = 199, seed = 1, create_maps = FALSE)))
  expect_lte(wide(out), 70)
  expect_lte(wide(capture.output(summary(res))), 70)
  expect_lte(wide(capture.output(print(sppt_clusters(res, nsim = 99,
                                                     seed = 1)))), 70)
})
