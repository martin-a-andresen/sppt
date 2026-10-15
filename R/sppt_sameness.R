#' Test a robust S index against sameness
#'
#' @description
#' Tests whether two spatial patterns are more different than two realisations
#' of a single process would be. The null distribution of the S index follows
#' from the calibration of the test itself: with units flagged at `alpha`, a
#' proportion `alpha` of eligible units flag by chance under identity, so the
#' number flagged is \eqn{\mathrm{Binomial}(n, \alpha)} and S is centred at
#' \eqn{1 - \alpha} with standard error \eqn{\sqrt{\alpha(1 - \alpha)/n}},
#' where \eqn{n} is the number of eligible units.
#'
#' Verified against 10,120 simulated combinations (8 crime types, 23 years,
#' 55 spatial scales; 5,000 iterations each). The mean error of the
#' closed-form threshold is +0.0099 at `alpha = 0.05` across the full grid,
#' and +0.0006 where at least 95 percent of units are eligible. Where many
#' units are empty in both patterns the test is mildly conservative, because
#' discreteness in sparse units lowers the realised flagging rate below
#' nominal.
#'
#' [sppt()] calls this function automatically; it is exported so that the
#' test can be applied to S indices computed elsewhere.
#'
#' @param s_index Numeric. Observed robust S index.
#' @param n_sig Integer. Number of units flagged as dissimilar. Used for the
#'   exact binomial tail; may be `NA`, in which case the normal approximation
#'   is used.
#' @param n_eligible Integer. Number of units in the S index denominator, that
#'   is, the total less units empty in both patterns.
#' @param alpha Numeric in (0, 1). The unit-level significance level used in
#'   the spatial point pattern test. Default `0.05`.
#' @param tau Numeric in (0, 1). Significance level for the sameness test
#'   itself. This is a separate choice from `alpha`: `alpha` governs when an
#'   individual unit is flagged, `tau` governs how much of the lower tail of
#'   the null distribution of S is treated as ordinary. Default `0.05`.
#' @param exact Logical or `NULL`. Use the exact binomial tail rather than the
#'   normal approximation. `NULL` (default) selects the exact test when
#'   `n_eligible < 200`, where the two diverge materially and the normal
#'   approximation is anti-conservative.
#'
#' @return A list with elements `s_index`, `n_eligible`, `alpha`, `tau`,
#'   `s_expected`, `se`, `threshold`, `z`, `p_value`, `method`, and `same`
#'   (`TRUE` when the two patterns are no more different than two
#'   realisations of one process).
#'
#' @seealso [sppt()]
#' @examples
#' sppt_sameness(s_index = 0.90, n_sig = 20, n_eligible = 200, alpha = 0.05)
#'
#' @export
sppt_sameness <- function(s_index, n_sig, n_eligible, alpha = 0.05,
                          tau = 0.05, exact = NULL) {

  out <- list(s_index = s_index, n_eligible = n_eligible, alpha = alpha,
              tau = tau, s_expected = NA_real_, se = NA_real_,
              threshold = NA_real_, z = NA_real_, p_value = NA_real_,
              method = NA_character_, same = NA)

  if (length(s_index) != 1L || is.na(s_index) ||
      length(n_eligible) != 1L || is.na(n_eligible) || n_eligible < 1L ||
      length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1)
    return(out)

  if (is.null(exact)) exact <- n_eligible < 200L
  use_exact <- isTRUE(exact) && !is.na(n_sig)

  s_exp <- 1 - alpha
  se    <- sqrt(alpha * (1 - alpha) / n_eligible)
  thr   <- s_exp - stats::qnorm(1 - tau) * se
  z     <- (s_index - s_exp) / se

  # One-sided: is S LOWER than sameness would produce?
  p <- if (use_exact) {
    # P(flagged >= n_sig) under Binomial(n_eligible, alpha)
    stats::pbinom(n_sig - 1L, size = n_eligible, prob = alpha,
                  lower.tail = FALSE)
  } else {
    stats::pnorm(z)
  }

  out$s_expected <- s_exp
  out$se         <- se
  out$threshold  <- thr
  out$z          <- z
  out$p_value    <- p
  out$method     <- if (use_exact) "exact binomial" else "normal"
  out$same       <- p >= tau
  out
}
