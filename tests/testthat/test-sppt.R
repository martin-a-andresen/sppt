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
  expect_output(summary(run_sppt(arrangement = FALSE)), "Robust S-index")
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
