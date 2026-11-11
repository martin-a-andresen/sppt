# Crime concentration: Gini and generalised Gini, Lorenz curve, 50-X and
# generalised 50-X, the Gini against random allocation with its local
# counterpart, and the rank envelope.

#' Gini and generalised Gini coefficient
#'
#' @description
#' The Gini coefficient of event counts across units, and the generalised
#' Gini of Bernasco and Steenbeek (2017) for when there are fewer events than
#' units. The calculation follows the lorenzgini package (Steenbeek and
#' Bernasco 2018), written here from the published formulas.
#'
#' @details
#' With counts \eqn{x_{(1)} \le \dots \le x_{(n)}} over \eqn{n} units and
#' \eqn{N} events, the Gini is
#' \deqn{G = \frac{2 \sum_i i\, x_{(i)}}{n N} - \frac{n + 1}{n}.}
#' All events in one unit give \eqn{G = 1 - 1/n}; the unbiased Gini
#' multiplies by \eqn{n/(n - 1)} so that it gives 1. With fewer events than
#' units, perfect equality is impossible (at most \eqn{N} units can have an
#' event), and the generalised Gini rescales the Gini against that limit:
#' \deqn{G' = \max(1, n/N)(G - 1) + 1,}
#' which equals \eqn{G} when \eqn{N \ge n}. Its unbiased form multiplies by
#' \eqn{m/(m - 1)}, \eqn{m = \min(n, N)}.
#'
#' With `conf.level`, a bootstrap confidence interval is computed by
#' resampling units with replacement ([boot::boot()] and [boot::boot.ci()]),
#' as in lorenzgini.
#'
#' @param x Non-negative counts, one per unit. Returns `NA` if any value is
#'   missing or negative, or if there are no events.
#' @param generalized Logical. Generalised Gini? Default `TRUE`.
#' @param unbiased Logical. Unbiased form? Default `FALSE`, as in Gao and
#'   Andresen (2025). (The lorenzgini default is `TRUE`.)
#' @param conf.level Confidence level for a bootstrap interval, or `NA`
#'   (default) for none.
#' @param R Number of bootstrap replicates. Default 1000. The BCa interval
#'   needs more replicates than units, so `R` is raised (with a message) to
#'   the next multiple of 1000 above the number of units when it is not.
#' @param type Interval type passed to [boot::boot.ci()]. Default `"bca"`.
#'
#' @return The coefficient, or with `conf.level` a named vector `gini`,
#'   `lwr.ci`, `upr.ci`.
#'
#' @references
#' Bernasco, W. and Steenbeek, W. (2017). More places than crimes:
#' Implications for evaluating the law of crime concentration at place.
#' *Journal of Quantitative Criminology*, 33(3), 451-467.
#' \doi{10.1007/s10940-016-9324-7}
#'
#' Steenbeek, W. and Bernasco, W. (2018). lorenzgini: Generalized Gini for
#' sparse data situations. R package version 0.1.2.
#' <https://github.com/wsteenbeek/lorenzgini>
#'
#' @examples
#' sppt_gini(c(10, 8, 2, rep(0, 7)), generalized = FALSE)  # 0.78
#' sppt_gini(c(5, 2, 1, rep(0, 7)), generalized = FALSE)   # 0.80
#' sppt_gini(c(5, 2, 1, rep(0, 7)))                        # 0.75
#' sppt_gini(c(10, rep(0, 9)), generalized = FALSE)        # 0.90
#' sppt_gini(c(10, rep(0, 9)), generalized = FALSE, unbiased = TRUE)  # 1
#' @seealso [sppt_concentration()], [sppt_lorenz()], [sppt_gini_test()]
#' @export
sppt_gini <- function(x, generalized = TRUE, unbiased = FALSE,
                      conf.level = NA, R = 1000, type = "bca") {
  x <- as.numeric(x)
  if (!length(x) || anyNA(x) || any(x < 0) || sum(x) <= 0) return(NA_real_)
  if (is.na(conf.level))
    return(gini_value(x, generalized, unbiased))
  # The BCa interval estimates its adjustment by regression on the
  # replicates, which needs more replicates than units.
  if (identical(type, "bca") && R <= length(x)) {
    R_new <- ceiling((length(x) + 1) / 1000) * 1000
    message("BCa needs more replicates than units: R raised from ", R,
            " to ", R_new, ".")
    R <- R_new
  }
  b <- boot::boot(x, function(v, i) gini_value(v[i], generalized, unbiased),
                  R = R)
  ci <- boot::boot.ci(b, conf = conf.level, type = type)
  lim <- ci[[length(ci)]]
  k <- ncol(lim)
  c(gini = b$t0, lwr.ci = lim[1L, k - 1L], upr.ci = lim[1L, k])
}

# Gini of non-negative counts; NA when undefined.
gini_value <- function(x, generalized = TRUE, unbiased = FALSE) {
  n <- length(x)
  N <- sum(x)
  if (n < 1L || N <= 0) return(NA_real_)
  x <- sort(x)
  g <- 2 * sum(seq_len(n) * x) / (n * N) - (n + 1) / n
  if (generalized) g <- max(1, n / N) * (g - 1) + 1
  if (unbiased) {
    m <- if (generalized) min(n, N) else n
    if (m <= 1) return(NA_real_)
    g <- g * m / (m - 1)
  }
  g
}

# Share of units (percent) holding at least share p of the events: the
# smallest k with the k largest counts summing to at least p * N.
units_for_share <- function(x, p, generalized = FALSE) {
  N <- sum(x)
  if (N <= 0) return(NA_real_)
  cs <- cumsum(sort(x, decreasing = TRUE))
  k <- which(cs >= p * N - 1e-9 * N)[1L]
  100 * k / (if (generalized) min(length(x), N) else length(x))
}

#' Crime concentration metrics
#'
#' @description
#' For each count column: the Gini and generalised Gini coefficients
#' ([sppt_gini()]), and the percentage of units that hold 50 and 25 per
#' cent of the events (Weisburd 2015), also in the generalised form of
#' Andresen and Weisburd (2026), which divides by the number of units that
#' could hold an event, \eqn{\min(n, N)}.
#'
#' @details
#' The percentage of units holding a share of the events counts the
#' smallest number of units, taken in decreasing order of counts, whose
#' events reach at least that share, so the unit that crosses the line is
#' included. All units are counted in the denominator, including those
#' with no events.
#'
#' @param data A data frame or `sf` object with the counts, or a numeric
#'   vector of counts.
#' @param count_col Character. One or more count columns of `data`. Ignored
#'   when `data` is a vector.
#' @param unbiased Logical. Unbiased Gini? Default `FALSE`.
#'
#' @return A data frame with one row per count column: `variable`,
#'   `n_units`, `n_events`, `gini`, `gini_generalized`, `pct_units_50`,
#'   `pct_units_25`, `pct_units_50_gen` and `pct_units_25_gen`.
#'
#' @references
#' Andresen, M.A. and Weisburd, D. (2026). Ockham's razor and the
#' measurement of crime concentrations at places. *Journal of Quantitative
#' Criminology*, 42(3), 569-588.
#'
#' Bernasco, W. and Steenbeek, W. (2017). More places than crimes:
#' Implications for evaluating the law of crime concentration at place.
#' *Journal of Quantitative Criminology*, 33(3), 451-467.
#' \doi{10.1007/s10940-016-9324-7}
#'
#' Weisburd, D. (2015). The law of crime concentration and the criminology
#' of place. *Criminology*, 53(2), 133-157. \doi{10.1111/1745-9125.12070}
#'
#' @examples
#' sppt_concentration(vancouver_ct, c("BNER2021", "BNEC2021"))
#' @seealso [sppt_gini()], [sppt_gini_test()], [sppt_envelope()]
#' @export
sppt_concentration <- function(data, count_col = NULL, unbiased = FALSE) {
  cols <- concentration_counts(data, count_col)
  do.call(rbind, lapply(names(cols), function(nm) {
    x <- cols[[nm]]
    data.frame(variable = nm, n_units = length(x), n_events = sum(x),
               gini = gini_value(x, FALSE, unbiased),
               gini_generalized = gini_value(x, TRUE, unbiased),
               pct_units_50 = units_for_share(x, 0.50),
               pct_units_25 = units_for_share(x, 0.25),
               pct_units_50_gen = units_for_share(x, 0.50, TRUE),
               pct_units_25_gen = units_for_share(x, 0.25, TRUE))
  }))
}

# Named list of validated count vectors.
concentration_counts <- function(data, count_col) {
  if (is.numeric(data) && is.null(dim(data))) {
    x <- list(as.numeric(data))
    names(x) <- if (is.null(count_col)) "x" else count_col[1L]
  } else {
    if (!is.data.frame(data))
      stop("`data` must be a data frame, an sf object or a numeric vector.",
           call. = FALSE)
    if (is.null(count_col) || !length(count_col))
      stop("`count_col` must name one or more count columns.", call. = FALSE)
    miss <- setdiff(count_col, names(data))
    if (length(miss))
      stop("Column(s) not found in `data`: ", paste(miss, collapse = ", "),
           call. = FALSE)
    d <- if (inherits(data, "sf")) sf::st_drop_geometry(data) else data
    x <- lapply(count_col, function(nm) as.numeric(d[[nm]]))
    names(x) <- count_col
  }
  for (nm in names(x)) {
    v <- x[[nm]]
    if (anyNA(v) || any(v < 0))
      stop("`", nm, "` must be non-negative counts with no missing values.",
           call. = FALSE)
    if (sum(v) <= 0) stop("`", nm, "` has no events.", call. = FALSE)
  }
  x
}

#' Lorenz curve
#'
#' @description
#' The Lorenz curve of event counts across units: the cumulative share of
#' events against the cumulative share of units, units in decreasing order
#' of counts. When there are fewer events than units, the dotted line marks
#' maximal equality under that constraint (every event in a different
#' unit), reached at \eqn{c/n} on the horizontal axis (Bernasco and
#' Steenbeek 2017). After the lorenzgini package (Steenbeek and Bernasco
#' 2018).
#'
#' @param x Non-negative counts, one per unit.
#' @param prop Logical. Proportions (default) or percentages on the axes.
#' @param rescale Logical. With fewer events than units, keep only the
#'   \eqn{c} busiest units, so the constrained line of equality is the
#'   diagonal. Default `FALSE`.
#'
#' @return A `ggplot` object.
#' @inherit sppt_gini references
#' @examples
#' sppt_lorenz(vancouver_ct$BNER2021)
#' sppt_lorenz(c(3, 1, 1, rep(0, 150)), rescale = TRUE)
#' @export
sppt_lorenz <- function(x, prop = TRUE, rescale = FALSE) {
  x <- as.numeric(x)
  if (!length(x) || anyNA(x) || any(x < 0) || sum(x) <= 0)
    stop("`x` must be non-negative counts with at least one event.",
         call. = FALSE)
  x <- sort(x, decreasing = TRUE)
  N <- sum(x)
  if (rescale && N < length(x)) x <- x[seq_len(N)]
  n <- length(x)
  s <- if (prop) 1 else 100
  d <- data.frame(units = c(0, seq_len(n) / n) * s,
                  events = c(0, cumsum(x) / N) * s)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$units, y = .data$events)) +
    ggplot2::annotate("segment", x = 0, y = 0, xend = s, yend = s,
                      linetype = "longdash") +
    ggplot2::geom_line() +
    ggplot2::scale_x_continuous(limits = c(0, s)) +
    ggplot2::scale_y_continuous(limits = c(0, s)) +
    ggplot2::labs(
      x = paste("Cumulative", if (prop) "proportion" else "percentage",
                "of units"),
      y = paste("Cumulative", if (prop) "proportion" else "percentage",
                "of events")) +
    ggplot2::coord_equal() +
    ggplot2::theme_light()
  if (N < n) {
    cn <- N / n * s
    p <- p +
      ggplot2::annotate("segment", x = 0, y = 0, xend = cn, yend = s,
                        linetype = "dotted") +
      ggplot2::annotate("segment", x = cn, y = 0, xend = cn, yend = s,
                        linetype = "dotted") +
      ggplot2::annotate("text", x = cn, y = 0, label = "c/n", vjust = -0.5,
                        hjust = -0.2)
  }
  p
}

# Exposure weights for the random allocations of the concentration
# functions; as in sppt(), but a data frame without geometry can only use
# "units" or a column.
concentration_weights <- function(data, random_by, n) {
  if (is.numeric(data) && is.null(dim(data))) {
    if (!identical(random_by, "units"))
      stop("With a vector of counts, `random_by` must be \"units\".",
           call. = FALSE)
    return(list(w = rep(1, n), label = "equal probability per unit",
                by = "units"))
  }
  if (identical(random_by, "exposure") && !inherits(data, "sf"))
    stop("`random_by = \"exposure\"` needs an sf layer (area or length); ",
         "use `random_by = \"units\"` or name an exposure column.",
         call. = FALSE)
  if (!inherits(data, "sf") && !identical(random_by, "units")) {
    if (!random_by %in% names(data))
      stop("No column \"", random_by, "\" in `data`.", call. = FALSE)
    w <- as.numeric(data[[random_by]])
    if (anyNA(w) || any(w < 0) || sum(w) <= 0)
      stop("The exposure must be finite and non-negative, with a positive ",
           "total.", call. = FALSE)
    return(list(w = w, label = paste0("probability in proportion to `",
                                      random_by, "`"), by = random_by))
  }
  if (identical(random_by, "units"))
    return(list(w = rep(1, n), label = "equal probability per unit",
                by = "units"))
  exposure_weights(data, random_by)
}

#' Gini coefficient against random allocation, with its local areas
#'
#' @description
#' Tests whether events are more concentrated than a random allocation
#' would make them, and finds the units that drive the result, following
#' Gao and Andresen (2025).
#'
#' **Global.** The events are allocated at random to the units `R` times,
#' by sampling units with replacement (a multinomial allocation of the
#' observed total), and the generalised Gini is computed for each
#' allocation. The expected value is their mean, and the observed Gini is
#' significant when it lies outside the central `level` share of them. The
#' marginal generalised Gini is the observed minus the expected value,
#' \eqn{\Delta G' = G' - E(G')}. The same allocations give the marginal
#' crime concentration of Chalfin, Kaplan and Cuellar (2021),
#' \eqn{\Delta cc = E(cc) - cc}, where \eqn{cc} is the percentage of units
#' holding 50 per cent of the events.
#'
#' **Local.** For each unit, the central `level` share of its counts across
#' the allocations is its non-rejection envelope. Units above it hold more
#' events than random, units below fewer. The dissimilarity index is the
#' share of all units outside their envelopes, and it splits into a
#' concentration index (above) and a dispersion index (below).
#'
#' @details
#' Gao and Andresen (2025) give every unit equal probability. Here the
#' default, as in [sppt()], is `random_by = "exposure"`: in proportion to
#' area for polygons or length for lines such as street segments. Use
#' `random_by = "units"` for the published method; with `R = 1000` and the
#' seed of the published code, it reproduces that code's results. The
#' indices use all units as the denominator, including units with no
#' events and no expected events.
#'
#' @param data An `sf` object or data frame with the counts, or a numeric
#'   vector of counts (then `random_by` must be `"units"`).
#' @param count_col Character. One or more count columns of `data`.
#' @param random_by `"exposure"` (default), `"units"` or the name of a
#'   non-negative numeric column, as in [sppt()].
#' @param R Number of random allocations. Default 1000.
#' @param level Width of the envelopes. Default 0.95 (2.5th and 97.5th
#'   percentiles).
#' @param unbiased Logical. Unbiased generalised Gini? Default `FALSE`.
#' @param seed Optional random seed.
#' @param verbose Logical. Print the report? Default `TRUE`.
#'
#' @return `data` with, for each count column `X`, columns `X_expected`
#'   (expected count), `X_lower` and `X_upper` (the unit's envelope) and
#'   `X_local` (1 above, -1 below, 0 inside), and class `sppt_gini_test`.
#'   Attribute `"gini_test"` holds one row per count column: `variable`,
#'   `n_units`, `n_events`, `gini` (observed generalised Gini), `expected`,
#'   `lower`, `upper`, `marginal_gini`, `significant`, `p_value` (share of
#'   allocations at least as concentrated, \eqn{(1 + r)/(R + 1)}),
#'   `pct_units_50`, `pct_units_50_random`, `marginal_concentration`,
#'   `dissimilarity`, `concentration` and `dispersion`.
#'
#' @references
#' Chalfin, A., Kaplan, J. and Cuellar, M. (2021). Measuring marginal crime
#' concentration: A new solution to an old problem. *Journal of Research in
#' Crime and Delinquency*, 58(4), 467-504.
#'
#' Gao, Y. and Andresen, M.A. (2025). Generalized Gini coefficient, its
#' statistical significance, and the local areas driving the global
#' result. *Journal of Criminal Justice*, 101, 102524.
#' \doi{10.1016/j.jcrimjus.2025.102524}
#'
#' @examples
#' \donttest{
#' g <- sppt_gini_test(vancouver_ct, "BNER2021", seed = 1)
#' plot(g, "BNER2021")
#' }
#' @seealso [sppt_concentration()], [sppt_envelope()], [sppt()]
#' @export
sppt_gini_test <- function(data, count_col = NULL, random_by = "exposure",
                           R = 1000L, level = 0.95, unbiased = FALSE,
                           seed = NULL, verbose = TRUE) {
  check_level(level)
  R <- check_R(R)
  cols <- concentration_counts(data, count_col)
  n <- length(cols[[1L]])
  expo <- concentration_weights(data, random_by, n)
  if (!is.null(seed)) set.seed(seed)
  lo_p <- (1 - level) / 2
  hi_p <- 1 - lo_p

  out <- if (is.numeric(data) && is.null(dim(data)))
    data.frame(x = cols[[1L]]) else data
  rows <- list()
  for (nm in names(cols)) {
    x <- cols[[nm]]
    N <- sum(x)
    draws <- matrix(stats::rmultinom(R, N, expo$w), nrow = n)
    g_obs <- gini_value(x, TRUE, unbiased)
    g_sim <- apply(draws, 2L, gini_value, generalized = TRUE,
                   unbiased = unbiased)
    q <- stats::quantile(g_sim, c(lo_p, hi_p), names = FALSE)
    cc_obs <- units_for_share(x, 0.5)
    cc_sim <- apply(draws, 2L, units_for_share, p = 0.5)
    lower <- apply(draws, 1L, stats::quantile, probs = lo_p, names = FALSE)
    upper <- apply(draws, 1L, stats::quantile, probs = hi_p, names = FALSE)
    loc <- ifelse(x > upper, 1L, ifelse(x < lower, -1L, 0L))
    out[[paste0(nm, "_expected")]] <- N * expo$w / sum(expo$w)
    out[[paste0(nm, "_lower")]] <- lower
    out[[paste0(nm, "_upper")]] <- upper
    out[[paste0(nm, "_local")]] <- loc
    rows[[nm]] <- data.frame(
      variable = nm, n_units = n, n_events = N, gini = g_obs,
      expected = mean(g_sim), lower = q[1L], upper = q[2L],
      marginal_gini = g_obs - mean(g_sim),
      significant = g_obs > q[2L] || g_obs < q[1L],
      p_value = (1 + sum(g_sim >= g_obs - 1e-12)) / (R + 1),
      pct_units_50 = cc_obs, pct_units_50_random = mean(cc_sim),
      marginal_concentration = mean(cc_sim) - cc_obs,
      dissimilarity = mean(loc != 0L), concentration = mean(loc == 1L),
      dispersion = mean(loc == -1L))
  }
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  attr(out, "gini_test") <- tab
  attr(out, "gini_test_meta") <- list(R = R, level = level,
                                      random = expo$label,
                                      unbiased = unbiased)
  class(out) <- unique(c("sppt_gini_test", class(out)))
  if (verbose) print(out)
  invisible(out)
}

#' @rdname sppt_gini_test
#' @param x A result of `sppt_gini_test()`.
#' @param ... Unused.
#' @export
print.sppt_gini_test <- function(x, ...) {
  tab <- attr(x, "gini_test")
  meta <- attr(x, "gini_test_meta")
  cat("\n", rule_line(), "Gini Test Against Random Allocation\n",
      rule_line(), sep = "")
  cat(field_lines("Random", sprintf("%d allocations, %s", meta$R,
                                    meta$random)))
  for (i in seq_len(nrow(tab))) {
    r <- tab[i, ]
    pc <- function(v) sprintf("%.1f%%", v)
    cat(rule_line("-"))
    cat(field_lines("Variable", sprintf("%s (N = %s, %d units)", r$variable,
                                        format(r$n_events, big.mark = ","),
                                        r$n_units)))
    cat(field_lines("Gini", sprintf(
      "%.4f (generalised%s); expected under random %.4f, %g%% envelope %.4f to %.4f",
      r$gini, if (meta$unbiased) ", unbiased" else "", r$expected,
      100 * meta$level, r$lower, r$upper)))
    cat(field_lines("Marginal Gini", sprintf(
      "%.4f; %s (p %s)", r$marginal_gini,
      if (r$gini > r$upper) "more concentrated than random" else
        if (r$gini < r$lower) "less concentrated than random" else
          "not different from random", fmt_p(r$p_value))))
    cat(field_lines("50% of events", sprintf(
      "in %s of units (random: %s); marginal concentration %.1f",
      pc(r$pct_units_50), pc(r$pct_units_50_random),
      r$marginal_concentration)))
    cat(field_lines("Local", sprintf(
      "dissimilarity %.3f = concentration %.3f + dispersion %.3f (share of all %d units)",
      r$dissimilarity, r$concentration, r$dispersion, r$n_units)))
  }
  cat(rule_line())
  invisible(x)
}

#' @rdname sppt_gini_test
#' @param which Character. The count column to map. Default the first.
#' @param title Logical. Show a title? Default `FALSE`.
#' @param legend_position,north_arrow,scale_bar Layout, as in [sppt_map()].
#' @export
plot.sppt_gini_test <- function(x, which = NULL, title = FALSE,
                                legend_position = "bottom",
                                north_arrow = "tr", scale_bar = "bl", ...) {
  if (!inherits(x, "sf"))
    stop("Mapping needs an sf object.", call. = FALSE)
  tab <- attr(x, "gini_test")
  v <- if (is.null(which)) tab$variable[1L] else which[1L]
  if (!v %in% tab$variable)
    stop("`which` must be one of: ", paste(tab$variable, collapse = ", "),
         call. = FALSE)
  check_location(north_arrow, "north_arrow")
  check_location(scale_bar, "scale_bar")
  loc <- factor(x[[paste0(v, "_local")]], levels = c(1L, -1L, 0L),
                labels = c("Actual > random", "Actual < random",
                           "Not significant"))
  d <- sf::st_as_sf(data.frame(local = loc, geometry = sf::st_geometry(x)))
  p <- ggplot2::ggplot(d) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data$local), colour = "grey50",
                     linewidth = 0.1) +
    ggplot2::scale_fill_manual(
      values = c("Actual > random" = "#d7301f",
                 "Actual < random" = "#6baed6",
                 "Not significant" = "white"),
      drop = FALSE, name = NULL) +
    ggplot2::theme_void() +
    legend_theme(legend_position) +
    north_arrow_layer(north_arrow) +
    scale_bar_layer(scale_bar, 0.3)
  if (isTRUE(title))
    p <- p + ggplot2::labs(title = paste0("Local concentration: ", v))
  print(p)
  invisible(p)
}

#' Rank envelope of concentration against random allocation
#'
#' @description
#' Compares the observed counts, sorted from the busiest unit down, with
#' the same ranks under random allocation: the busiest unit with the
#' busiest unit, the second with the second, and so on. Location plays no
#' part, so this asks whether events are more concentrated than chance,
#' not where.
#'
#' @details
#' The observed total is allocated at random `R` times (a multinomial
#' allocation, by default in proportion to exposure, as in [sppt()]), each
#' allocation is sorted in decreasing order, and the central `level` share
#' of the counts at each rank is the envelope for that rank. Units are
#' above, below or within the envelope at their rank. The excess above the
#' upper bounds, as a percentage of all events, is the share of events
#' concentrated beyond what random allocation produces.
#'
#' @inheritParams sppt_gini_test
#' @param group_col Optional unit identifier column, carried into the
#'   detail table.
#'
#' @return An object of class `sppt_envelope`, a list with `summary` (one
#'   row per count column: `variable`, `n_units`, `n_events`, `n_over`,
#'   `n_under`, `n_within`, `pct_events_over`) and `detail` (one row per
#'   unit and count column: `variable`, `unit`, `rank`, `actual`, `lower`,
#'   `upper`, `excess`, positive above the envelope and negative below).
#'
#' @examples
#' \donttest{
#' e <- sppt_envelope(vancouver_ct, "BNER2021", group_col = "id", seed = 1)
#' plot(e)
#' }
#' @seealso [sppt_gini_test()], [sppt_concentration()]
#' @export
sppt_envelope <- function(data, count_col = NULL, random_by = "exposure",
                          R = 1000L, level = 0.95, group_col = NULL,
                          seed = NULL, verbose = TRUE) {
  check_level(level)
  R <- check_R(R)
  cols <- concentration_counts(data, count_col)
  n <- length(cols[[1L]])
  expo <- concentration_weights(data, random_by, n)
  ids <- if (!is.null(group_col)) {
    if (!is.data.frame(data) || !group_col %in% names(data))
      stop("No column \"", group_col, "\" in `data`.", call. = FALSE)
    data[[group_col]]
  } else seq_len(n)
  if (!is.null(seed)) set.seed(seed)
  w <- expo$w / sum(expo$w)
  lo_p <- (1 - level) / 2
  hi_p <- 1 - lo_p
  summ <- list()
  det <- list()
  for (nm in names(cols)) {
    x <- cols[[nm]]
    N <- sum(x)
    ord <- order(x, decreasing = TRUE)
    xs <- x[ord]
    sim <- matrix(NA_integer_, nrow = n, ncol = R)
    for (i in seq_len(R))
      sim[, i] <- sort(as.integer(stats::rmultinom(1L, N, w)),
                       decreasing = TRUE)
    lower <- apply(sim, 1L, stats::quantile, probs = lo_p, names = FALSE)
    upper <- apply(sim, 1L, stats::quantile, probs = hi_p, names = FALSE)
    excess <- ifelse(xs > upper, xs - upper,
                     ifelse(xs < lower, xs - lower, 0))
    summ[[nm]] <- data.frame(
      variable = nm, n_units = n, n_events = N,
      n_over = sum(xs > upper), n_under = sum(xs < lower),
      n_within = sum(xs >= lower & xs <= upper),
      pct_events_over = 100 * sum(pmax(xs - upper, 0)) / N)
    det[[nm]] <- data.frame(variable = nm, unit = ids[ord],
                            rank = seq_len(n), actual = xs, lower = lower,
                            upper = upper, excess = excess)
  }
  out <- list(summary = `rownames<-`(do.call(rbind, summ), NULL),
              detail = `rownames<-`(do.call(rbind, det), NULL),
              R = R, level = level, random = expo$label)
  class(out) <- "sppt_envelope"
  if (verbose) print(out)
  invisible(out)
}

#' @rdname sppt_envelope
#' @param x A result of `sppt_envelope()`.
#' @param ... Unused.
#' @export
print.sppt_envelope <- function(x, ...) {
  cat("\n", rule_line(), "Rank Envelope Against Random Allocation\n",
      rule_line(), sep = "")
  cat(field_lines("Random", sprintf("%d allocations, %s; %g%% envelope",
                                    x$R, x$random, 100 * x$level)))
  for (i in seq_len(nrow(x$summary))) {
    r <- x$summary[i, ]
    cat(rule_line("-"))
    cat(field_lines("Variable", sprintf("%s (N = %s, %d units)", r$variable,
                                        format(r$n_events, big.mark = ","),
                                        r$n_units)))
    cat(field_lines("By rank", sprintf(
      "%d units above the envelope, %d below, %d within",
      r$n_over, r$n_under, r$n_within)))
    cat(field_lines("Above random", sprintf(
      "%.2f%% of events lie above the upper bounds", r$pct_events_over)))
  }
  cat(rule_line())
  invisible(x)
}

#' @rdname sppt_envelope
#' @param which Character. The count column to plot. Default the first.
#' @param log Logical. Logarithmic count axis (counts plus one)? Default
#'   `TRUE`.
#' @export
plot.sppt_envelope <- function(x, which = NULL, log = TRUE, ...) {
  v <- if (is.null(which)) x$summary$variable[1L] else which[1L]
  d <- x$detail[x$detail$variable == v, , drop = FALSE]
  if (!nrow(d)) stop("No count column \"", v, "\".", call. = FALSE)
  d$status <- factor(ifelse(d$excess > 0, "Above", ifelse(d$excess < 0,
                     "Below", "Within")), levels = c("Above", "Within",
                                                     "Below"))
  f <- if (log) function(v) v + 1 else identity
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$rank)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = f(.data$lower),
                                      ymax = f(.data$upper)),
                         fill = "grey80") +
    ggplot2::geom_point(ggplot2::aes(y = f(.data$actual),
                                     colour = .data$status), size = 0.8) +
    ggplot2::scale_colour_manual(values = c(Above = "#d7301f",
                                            Within = "grey30",
                                            Below = "#2171b5"),
                                 drop = FALSE, name = NULL) +
    ggplot2::labs(x = "Rank (busiest unit first)",
                  y = if (log) "Events + 1 (log scale)" else "Events") +
    ggplot2::theme_light() +
    ggplot2::theme(legend.position = "bottom")
  if (log) p <- p + ggplot2::scale_y_log10()
  print(p)
  invisible(p)
}

check_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("`level` must be a number in (0, 1).", call. = FALSE)
}

check_R <- function(R) {
  if (!is.numeric(R) || length(R) != 1L || is.na(R) || R < 19)
    stop("`R` must be a number of at least 19.", call. = FALSE)
  as.integer(R)
}
