test_that("longitudinal p-values are bounded and trends are classified", {
  res <- sppt_longitudinal(vancouver_ct, group_col = "id",
                           count_col = paste0("BNER", 2003:2010),
                           B = 99, seed = 1, verbose = FALSE,
                           create_maps = FALSE, create_time_plot = FALSE)
  pm <- attr(res, "pvalue_matrix")[, -1]
  expect_true(all(pm >= 1 / 100 & pm <= 1))
  expect_true(all(res$BNER_Trend_p %in% c(-2L, -1L, 0L, 1L, 2L)))
  expect_true(all(c("BNER_intervals_overlap_p", "BNER_Linear_bt",
                    "BNER2003_L", "BNER2003_U") %in% names(res)))
  # Linear trends: 7 envelope maps and 8 p-value maps
  expect_length(sppt_map(res), 15L)
})

test_that("fix_base runs through the p-value step", {
  res <- suppressWarnings(sppt_longitudinal(
    vancouver_ct, group_col = "id", count_col = paste0("BNER", 2003:2007),
    B = 99, seed = 1, fix_base = TRUE, verbose = FALSE,
    create_maps = FALSE, create_time_plot = FALSE))
  expect_equal(res$BNER2003_L, res$BNER2003_U)
  expect_false(anyNA(attr(res, "pvalue_matrix")[, -1]))
})

test_that("envelope bounds follow the rows of data, not sorted IDs", {
  # Character IDs sort as "1", "10", "11", ..., unlike the row order. A unit
  # with no events must get a zero envelope, whatever its position.
  d <- vancouver_ct[1:20, ]
  d$id <- as.character(d$id)
  d$BNER2003[c(2, 9, 15)] <- 0L
  res <- sppt_longitudinal(d, group_col = "id",
                           count_col = paste0("BNER", 2003:2005),
                           B = 99, seed = 1, verbose = FALSE,
                           create_maps = FALSE, create_time_plot = FALSE)
  zero <- d$BNER2003 == 0
  expect_true(all(res$BNER2003_U[zero] == 0))
  expect_true(all(res$BNER2003_U[!zero] > 0))
})

test_that("fewer than three periods is an error", {
  expect_error(sppt_longitudinal(vancouver_ct, group_col = "id",
                                 count_col = c("BNER2003", "BNER2004")),
               "at least 3")
})

test_that("longitudinal results have a class with summary and plot", {
  res <- sppt_longitudinal(vancouver_ct, group_col = "id",
                           count_col = paste0("BNER", 2003:2006),
                           B = 99, seed = 1, verbose = FALSE,
                           create_maps = FALSE, create_time_plot = FALSE)
  expect_s3_class(res, "sppt_long_result")
  expect_s3_class(res, "sf")
  expect_output(s <- summary(res), "Longitudinal SPPT Summary")
  expect_equal(unname(s$indices$pvalue["all", "overlap"]),
               attr(res, "si_overlap_p"))
  expect_identical(formals(sppt_longitudinal)$conf_level, 0.85)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_s3_class(plot(res), "ggplot")
})
