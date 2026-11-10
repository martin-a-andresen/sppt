#' Degree of similarity between two spatial patterns
#'
#' @description
#' \eqn{S_D} is one minus the smallest share of events that would have to be
#' redistributed, within one pattern, for the two patterns to be considered
#' the same: for a retest of the redistributed patterns to give
#' \eqn{S \ge 1 - \alpha} with probability `level`. It is computed in both
#' directions, moving the base pattern with the test pattern fixed and the
#' reverse, and \eqn{S_D} uses the smaller share.
#'
#' `sppt_similarity()` computes it for an existing [sppt()] result. Running
#' [sppt()] with `similarity = TRUE` computes it as part of the test, and
#' `sppt_similarity()` then returns the stored values.
#'
#' @details
#' **What is redistributed.** Each pattern is treated as shares of its own
#' total, at double precision, so a redistribution moves shares, and the
#' number of events it represents is that share times the pattern's total.
#' The pattern not being moved, and the totals, stay as observed.
#'
#' **When a unit stops differing.** A unit is flagged when its p-value,
#' \eqn{(1 + r)/(B + 1)}, is at most \eqn{\alpha}, where \eqn{r} counts the
#' centred bootstrap differences at least as large as the observed one. So a
#' unit is not flagged exactly when its share difference is within a
#' half-width given by its bootstrap distribution. Moving a unit only to that
#' boundary is not enough: on a retest it is flagged about half the time,
#' through Monte Carlo variation, and its interval narrows where the moved
#' pattern loses events. Each unit is therefore moved until the share of
#' bootstrap differences at least as large as its own reaches a margin
#' \eqn{q \ge \alpha}, and the intervals are recomputed at the moved shares
#' until they settle.
#'
#' **The minimum.** For a given margin, the over-represented units must give
#' up \eqn{E} and the under-represented units must gain \eqn{F}; the smallest
#' redistribution moves \eqn{\max(E, F)}, the difference going to, or coming
#' from, units with room inside their intervals.
#'
#' **The retest.** At the moved shares, each unit's chance of being flagged
#' on a retest is a binomial probability, from a separate bootstrap, so the
#' number flagged has a Poisson-binomial distribution and the probability of
#' \eqn{S \ge 1 - \alpha} follows exactly. The margin is the smallest
#' \eqn{q}, found by bisection to within `tol`, at which that probability
#' reaches `level`. A retest here means testing the redistributed patterns
#' again, with the same totals and a new bootstrap.
#'
#' **Properties.** \eqn{S_D} lies in \[0, 1\] and equals 1 when no
#' redistribution is needed. The two directions differ because moving a
#' pattern changes its own sampling variance, and with it the intervals. The
#' redistribution needed depends on power: with few events per unit the
#' intervals are wide and less has to move, so \eqn{S_D} at a fine scale is
#' not comparable with \eqn{S_D} at a coarse one.
#'
#' The calculation always uses the uncorrected flags and
#' \eqn{B} = 4999 or the larger number that `alpha` needs (as with
#' `B = "auto"` without a correction), or `B` itself if a number was given
#' and no correction was used. Against a random distribution only the
#' observed pattern is moved. With `fix_base = TRUE` the base is not
#' bootstrapped, as in the test. It takes roughly 30 to 80 bootstraps per
#' direction: seconds for a hundred units, a few minutes for a thousand.
#' The random-number state is restored afterwards, so the rest of the
#' result is the same with or without it.
#'
#' @param result A result returned by [sppt()] with p-values.
#' @param level Probability that a retest of the redistributed patterns
#'   gives \eqn{S \ge 1 - \alpha}. Default 0.95.
#' @param tol Bisection tolerance on the margin. Default 0.0005.
#' @param seed Optional seed for the bootstraps. Without one the current
#'   random-number stream is used, and then restored.
#'
#' @return A data frame with one row per direction (one row against a
#'   random distribution): `moved` (the pattern redistributed), `share` (the
#'   share of that pattern moved), `events` (the number of events moved),
#'   `margin` (the margin \eqn{q}) and `p_retest` (the probability that a
#'   retest gives \eqn{S \ge 1 - \alpha}). Attributes: `s_index_d`
#'   (\eqn{S_D}, one minus the smaller share), `level`, `alpha`, `B`, and
#'   `shares`, a list with the redistributed shares of each unit for each
#'   direction (proportions, in the order of the units), which can be mapped
#'   against the observed shares to see where events move.
#'
#' @examples
#' \donttest{
#' res <- sppt(vancouver_ct, group_col = "id",
#'             count_col = c("BNER2021", "BNEC2021"), seed = 1,
#'             verbose = FALSE, create_maps = FALSE)
#' sppt_similarity(res, seed = 1)
#' }
#'
#' @seealso [sppt()], [sppt_sameness()]
#' @export
sppt_similarity <- function(result, level = 0.95, tol = 5e-4, seed = NULL) {
  if (!inherits(result, "sppt_result"))
    stop("`result` must be returned by sppt().", call. = FALSE)
  s <- attr(result, "sppt_summary")
  if (!"pvalue" %in% names(result) || is.null(s$alpha))
    stop("`result` has no p-values; run sppt() with ",
         "calculate_pvalues = TRUE.", call. = FALSE)
  if (!is.null(s$similarity) && is.null(seed) &&
      isTRUE(all.equal(s$similarity$level, level)) &&
      isTRUE(all.equal(s$similarity$tol, tol)))
    return(similarity_frame(s$similarity))
  cnt <- function(nm) {
    x <- as.numeric(result[[nm]])
    x[is.na(x)] <- 0
    x
  }
  sim <- similarity_search(cnt(s$base_name), cnt(s$test_name),
                           base_name = s$base_name, test_name = s$test_name,
                           alpha = s$alpha, B = similarity_B(s),
                           level = level, tol = tol,
                           fix_base = isTRUE(s$fix_base),
                           random_side = s$random$side, seed = seed)
  similarity_frame(sim)
}

similarity_frame <- function(sim) {
  out <- sim$directions
  attr(out, "s_index_d") <- sim$s_index_d
  attr(out, "level") <- sim$level
  attr(out, "alpha") <- sim$alpha
  attr(out, "B") <- sim$B
  attr(out, "shares") <- sim$shares
  out
}

# B for the similarity search: the B of the run without a correction,
# otherwise what B = "auto" gives without one.
similarity_B <- function(s) {
  if (is.null(s$correction) || s$correction == "none") s$B
  else choose_B(s$alpha, s$n_eligible, "none")$B
}

# Bootstrap differences for the active units, as in the test. `abs` holds
# |differences centred on their replicate mean| (what the test compares with
# the observed difference); `u` the differences centred on their exact
# expectation, and `se` the standard error of the replicate mean. A retest
# centres on its own replicate mean, which differs from the expectation by
# about N(0, se). That shift moves every replicate relative to the observed
# difference at once, and decides ties: where one pattern has no events in
# a unit, the replicate in which the other also draws none sits exactly at
# the observed difference and counts or not with equal chance. Against a
# random distribution the test compares with the known expectation, so
# there is no shift.
sim_boot <- function(pa, pb, Na, Nb, B, act, fix_base, random_side) {
  n <- length(pa)
  if (!is.null(random_side)) {
    ref <- if (random_side == "base") pa else pb
    N   <- if (random_side == "base") Na else Nb
    x   <- matrix(stats::rmultinom(B, N, ref), nrow = n)[act, , drop = FALSE] / N
    u   <- x - ref[act]
    return(list(abs = abs(u), u = u, se = rep(0, sum(act))))
  }
  b <- matrix(stats::rmultinom(B, Nb, pb), nrow = n)[act, , drop = FALSE] / Nb
  a <- if (fix_base) pa[act] else
    matrix(stats::rmultinom(B, Na, pa), nrow = n)[act, , drop = FALSE] / Na
  d <- b - a
  mu <- rowMeans(d)
  sdv <- sqrt(pmax(rowMeans(d^2) - mu^2, 0))
  list(abs = abs(d - mu), u = d - (pb - pa)[act], se = sdv / sqrt(B))
}

# Each unit's chance of being flagged on a retest: the exceedance count is
# about binomial given the tail share, averaged over the retest's centring
# shift
# (four-point Gauss-Hermite rule; no node at zero, where a tie would count
# with certainty instead of half the time).
sim_flag_prob <- function(bt, dd, r_max, B) {
  nodes <- c(-2.3344142183389773, -0.7419637843027259,
             0.7419637843027259, 2.3344142183389773)
  w <- c(0.0458758547680685, 0.4541241452319315,
         0.4541241452319315, 0.0458758547680685)
  out <- 0
  for (k in seq_along(nodes)) {
    q <- rowMeans(abs(bt$u - nodes[k] * bt$se) >= dd - 1e-12)
    # The tail share is itself estimated from B replicates, which roughly
    # doubles the variance of the retest's count: normal approximation
    z <- (r_max + 0.5 - B * q) / sqrt(2 * B * q * (1 - q))
    out <- out + w[k] * stats::pnorm(z)
  }
  out
}

# Half-width at tail share q: the ceiling(q * B)-th largest value per row.
sim_halfwidth <- function(bt, q) {
  ab <- bt$abs
  k <- max(1L, as.integer(ceiling(q * ncol(ab) - 1e-9)))
  apply(ab, 1L, function(x) -sort(-x, partial = k)[k])
}

# Smallest move of shares x so that |x' - y| <= cc on the active units.
sim_project <- function(x, y, cc) {
  d <- x - y
  t <- x
  over  <- d > cc
  under <- -d > cc
  t[over]  <- y[over] + cc[over]
  t[under] <- pmin(y[under] - cc[under], 1)
  t[under] <- pmax(t[under], 0)
  gap <- sum(x[over] - t[over]) - sum(t[under] - x[under])
  ok  <- !over & !under
  if (gap > 0) {
    room <- pmax(pmin(y + cc, 1) - t, 0) * ok
    if (sum(room) < gap) return(NULL)
    t <- t + gap * room / sum(room)
  } else if (gap < 0) {
    room <- pmax(t - pmax(y - cc, 0), 0) * ok
    if (sum(room) < -gap) return(NULL)
    t <- t + gap * room / sum(room)
  }
  t
}

# P(at most m successes) for independent Bernoulli(p), truncated at m.
sim_pbinom_sum <- function(p, m) {
  f <- c(1, numeric(m))
  for (pi in p[p > 0]) f <- f * (1 - pi) + c(0, f[-(m + 1L)]) * pi
  sum(f)
}

# The search behind S_D. Counts are on the full set of units; base_name and
# test_name label the directions; random_side is NULL, "base" or "test".
similarity_search <- function(counts_a, counts_b, base_name, test_name,
                              alpha, B, level = 0.95, tol = 5e-4,
                              fix_base = FALSE, random_side = NULL,
                              seed = NULL, iters = 6L, q_max = 0.5) {
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    stats::runif(1L)
  old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  if (!is.null(seed)) set.seed(seed)
  s0 <- sample.int(.Machine$integer.max - 1L, 1L)

  Na <- sum(counts_a)
  Nb <- sum(counts_b)
  pa <- counts_a / Na
  pb <- counts_b / Nb
  act <- (counts_a + counts_b) > 0
  n_el <- sum(act)
  m <- floor(alpha * n_el + 1e-9)           # flags allowed for S >= 1 - alpha
  r_max <- floor(alpha * (B + 1) - 1 + 1e-9)  # flagged when r <= r_max

  dirs <- if (is.null(random_side)) c("base", "test") else
    if (random_side == "base") "test" else "base"

  one_direction <- function(dir) {
    x <- if (dir == "test") pb else pa       # pattern moved
    y <- if (dir == "test") pa else pb       # pattern held
    xa <- x[act]
    ya <- y[act]
    shares <- function(t) {                  # full-length shares, moved side
      v <- x
      v[act] <- t
      if (dir == "test") list(pa = pa, pb = v) else list(pa = v, pb = pb)
    }
    boot <- function(t, s) {
      set.seed(s)
      sh <- shares(t)
      sim_boot(sh$pa, sh$pb, Na, Nb, B, act, fix_base, random_side)
    }
    # Intervals do not depend on the moved pattern against a random
    # distribution, nor on the base under fix_base.
    fixed <- !is.null(random_side) || (fix_base && dir == "base")
    ab_fixed <- if (fixed) boot(xa, s0) else NULL
    ab_check <- if (fixed) boot(xa, s0 + 1L) else NULL
    warm <- xa

    evaluate <- function(q) {
      t <- warm
      for (i in seq_len(if (fixed) 1L else iters)) {
        ab <- if (fixed) ab_fixed else boot(t, s0)
        t_new <- sim_project(xa, ya, sim_halfwidth(ab, q))
        if (is.null(t_new)) return(NULL)
        done <- max(abs(t_new - t)) < 1e-10
        t <- t_new
        if (done) break
      }
      ab2 <- if (fixed) ab_check else boot(t, s0 + 1L)
      p_flag <- sim_flag_prob(ab2, abs(t - ya), r_max, B)
      list(q = q, t = t, share = 0.5 * sum(abs(t - xa)),
           p_ok = sim_pbinom_sum(p_flag, m))
    }

    # Bracket the margin, then bisect
    lo <- NULL
    hi <- evaluate(alpha)
    if (!is.null(hi) && hi$p_ok < level) {
      lo <- hi
      hi <- NULL
      step <- 0.005
      repeat {
        q <- min(lo$q + step, q_max)
        warm <- lo$t
        e <- evaluate(q)
        if (!is.null(e) && e$p_ok >= level) { hi <- e; break }
        if (is.null(e) || q >= q_max) break
        lo <- e
        step <- 2 * step
      }
      while (!is.null(hi) && hi$q - lo$q > tol) {
        warm <- lo$t
        e <- evaluate((lo$q + hi$q) / 2)
        if (is.null(e)) break
        if (e$p_ok >= level) hi <- e else lo <- e
      }
    }
    N_moved <- if (dir == "test") Nb else Na
    lab <- if (dir == "test") test_name else base_name
    if (is.null(hi))
      return(list(row = data.frame(moved = lab, share = NA_real_,
                                   events = NA_real_, margin = NA_real_,
                                   p_retest = NA_real_), shares = NULL))
    list(row = data.frame(moved = lab, share = hi$share,
                          events = hi$share * N_moved, margin = hi$q,
                          p_retest = hi$p_ok),
         shares = shares(hi$t)[[if (dir == "test") "pb" else "pa"]])
  }

  res <- lapply(dirs, one_direction)
  directions <- do.call(rbind, lapply(res, `[[`, "row"))
  rownames(directions) <- NULL
  shares <- stats::setNames(lapply(res, `[[`, "shares"), directions$moved)
  best <- suppressWarnings(min(directions$share, na.rm = TRUE))
  list(s_index_d = if (is.finite(best)) 1 - best else NA_real_,
       directions = directions, shares = shares, level = level, tol = tol,
       alpha = alpha, B = B)
}

# Console lines for S_D, shared by sppt(), its summary() method and the
# longitudinal first-against-last comparison.
similarity_lines <- function(sim) {
  d <- sim$directions
  if (is.na(sim$s_index_d))
    return(field_lines("S_D", paste(
      "not found: no redistribution up to the largest margin gives a",
      "retest S >= 1 - alpha with the required probability")))
  d <- d[order(d$share), , drop = FALSE]
  move <- vapply(seq_len(nrow(d)), function(i) sprintf(
    "%.1f%% of %s (%s events)", 100 * d$share[i], d$moved[i],
    format(round(d$events[i]), big.mark = ",")), character(1))
  paste0(
    field_lines("S_D", sprintf("%.4f", sim$s_index_d)),
    field_lines("Redistribution", sprintf(
      "%s%s; then a retest gives S >= %.2f with probability at least %.2f",
      move[1L], if (length(move) > 1L) paste0(", or ", move[2L]) else "",
      1 - sim$alpha, sim$level)))
}
