test_that("sppt_gini matches the lorenzgini documented values", {
  expect_equal(sppt_gini(c(10, 8, 2, rep(0, 7)), FALSE), 0.78)
  expect_equal(sppt_gini(c(10, 8, 2, rep(0, 7)), TRUE), 0.78)
  expect_equal(sppt_gini(c(5, 2, 1, rep(0, 7)), FALSE), 0.80)
  expect_equal(sppt_gini(c(5, 2, 1, rep(0, 7)), TRUE), 0.75)
  expect_equal(sppt_gini(c(10, rep(0, 9)), FALSE), 0.90)
  expect_equal(sppt_gini(c(10, rep(0, 99)), FALSE), 0.99)
  expect_equal(sppt_gini(c(10, rep(0, 9)), FALSE, unbiased = TRUE), 1)
  expect_equal(sppt_gini(c(10, rep(0, 99)), FALSE, unbiased = TRUE), 1)
  expect_true(is.na(sppt_gini(c(1, -1))))
  expect_true(is.na(sppt_gini(c(0, 0))))
})

test_that("sppt_gini agrees with independent formulas", {
  gmd <- function(x) {
    n <- length(x)
    sum(abs(outer(x, x, "-"))) / (2 * n^2 * mean(x))
  }
  # Gao and Andresen (2025), equation 3
  eq3 <- function(x) {
    n <- length(x); cc <- sum(x); y <- sort(x) / cc
    max(1 / cc, 1 / n) * (2 * sum(seq_len(n) * y) - n - 1) -
      max(n / cc, 1) + 1
  }
  set.seed(1)
  for (x in list(vancouver_ct$BNER2021, rpois(400, 0.3),
                 rnbinom(200, mu = 3, size = 0.4))) {
    expect_equal(sppt_gini(x, FALSE), gmd(x))
    expect_equal(sppt_gini(x), eq3(x))
  }
  set.seed(2)
  ci <- sppt_gini(rpois(300, 2), conf.level = 0.95, R = 499)
  expect_named(ci, c("gini", "lwr.ci", "upr.ci"))
  expect_true(ci[["lwr.ci"]] < ci[["gini"]] && ci[["gini"]] < ci[["upr.ci"]])
})

test_that("concentration metrics count the unit that crosses the line", {
  x <- c(5, 3, 1, 1, 0, 0, 0, 0, 0, 0)     # 10 events in 10 units
  m <- sppt_concentration(x)
  expect_equal(m$pct_units_50, 10)        # 5 of 10 events in 1 unit
  expect_equal(m$pct_units_25, 10)
  y <- c(3, 2, 1, rep(0, 17))             # 6 events, 20 units
  m <- sppt_concentration(y)
  expect_equal(m$pct_units_50, 5)         # 3 >= 3
  expect_equal(m$pct_units_50_gen, 100 / 6)
  expect_equal(m$gini_generalized, sppt_gini(y))
  m2 <- sppt_concentration(vancouver_ct, c("BNER2021", "BNEC2021"))
  expect_equal(m2$variable, c("BNER2021", "BNEC2021"))
  expect_equal(m2$n_events, c(sum(vancouver_ct$BNER2021),
                              sum(vancouver_ct$BNEC2021)))
  expect_error(sppt_concentration(vancouver_ct, "nope"), "not found")
  expect_s3_class(sppt_lorenz(y), "ggplot")
  expect_s3_class(sppt_lorenz(y, rescale = TRUE), "ggplot")
})

test_that("sppt_gini_test reproduces the published code with equal probability", {
  x <- vancouver_da$BNER2021
  N <- sum(x); n <- length(x)
  set.seed(221220)
  g_ref <- replicate(200, sppt_gini(sort(stats::rmultinom(1, N, rep(1, n)),
                                         decreasing = TRUE)))
  set.seed(221220)
  sims <- replicate(200, stats::rmultinom(1, N, rep(1, n)))[, 1, ]
  g <- sppt_gini_test(vancouver_da, "BNER2021", random_by = "units",
                      R = 200, seed = 221220, verbose = FALSE)
  t <- attr(g, "gini_test")
  expect_equal(t$gini, sppt_gini(x))
  expect_equal(t$expected, mean(g_ref))
  expect_equal(c(t$lower, t$upper),
               unname(stats::quantile(g_ref, c(0.025, 0.975))))
  expect_equal(t$marginal_gini, t$gini - t$expected)
  expect_equal(g$BNER2021_lower,
               unname(apply(sims, 1, stats::quantile, probs = 0.025)))
  expect_equal(t$dissimilarity, t$concentration + t$dispersion)
  expect_equal(t$concentration, mean(g$BNER2021_local == 1))
  expect_equal(sum(g$BNER2021_expected), N)
  expect_output(print(g), "Marginal Gini")
})

test_that("sppt_gini_test uses exposure by default and checks inputs", {
  g <- sppt_gini_test(vancouver_ct, c("BNER2021", "BNEC2021"), R = 99,
                      seed = 1, verbose = FALSE)
  a <- as.numeric(sf::st_area(vancouver_ct))
  expect_equal(g$BNER2021_expected, sum(vancouver_ct$BNER2021) * a / sum(a))
  expect_equal(attr(g, "gini_test")$variable, c("BNER2021", "BNEC2021"))
  expect_true(grepl("area", attr(g, "gini_test_meta")$random))
  d <- sf::st_drop_geometry(vancouver_ct)
  expect_error(sppt_gini_test(d, "BNER2021"), "needs an sf layer")
  expect_silent(sppt_gini_test(d, "BNER2021", random_by = "units", R = 99,
                               verbose = FALSE))
  expect_error(sppt_gini_test(vancouver_ct, "BNER2021", R = 5), "at least 19")
  expect_error(sppt_gini_test(vancouver_ct, "BNER2021", level = 2), "level")
  pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_s3_class(plot(g, "BNEC2021"), "ggplot")
})

test_that("sppt_envelope matches the published rank envelope", {
  x <- vancouver_da$BNER2021
  w <- as.numeric(sf::st_area(vancouver_da)); w <- w / sum(w)
  set.seed(8722)
  sim <- vapply(seq_len(200), function(i)
    sort(as.integer(stats::rmultinom(1, sum(x), w)), decreasing = TRUE),
    integer(length(x)))
  up <- apply(sim, 1, stats::quantile, probs = 0.975)
  xs <- sort(x, decreasing = TRUE)
  e <- sppt_envelope(vancouver_da, "BNER2021", group_col = "id", R = 200,
                     seed = 8722, verbose = FALSE)
  expect_equal(e$detail$upper, unname(up))
  expect_equal(e$summary$n_over, sum(xs > up))
  expect_equal(e$summary$pct_events_over, 100 * sum(pmax(xs - up, 0)) / sum(x))
  expect_equal(e$summary$n_over + e$summary$n_under + e$summary$n_within,
               nrow(vancouver_da))
  expect_equal(e$detail$unit[1], vancouver_da$id[which.max(x)])
  expect_output(print(e), "Rank Envelope")
  pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_s3_class(plot(e), "ggplot")
})
