#' Degree of similarity between two spatial patterns
#'
#' @description
#' Returns the degree-of-similarity measures of a [sppt()] result: the S
#' index, \eqn{S_D}, the index of dissimilarity \eqn{D}, the share of
#' \eqn{D} held by flagged units, the flagged share differences in each
#' direction, and the share and number of events of each pattern that would
#' have to move to remove the differences at the flagged units.
#'
#' @details
#' With \eqn{p_{ai}} and \eqn{p_{bi}} each unit's share of all events in the
#' base and test patterns (as proportions, whatever `use_percentages` was),
#' \deqn{D = \frac{1}{2}\sum_i |p_{bi} - p_{ai}|, \qquad
#' S_D = 1 - \frac{1}{2}\sum_{i \in F} |p_{bi} - p_{ai}|,}{D = 1/2 * sum over
#' all units of |p_bi - p_ai|,  S_D = 1 - 1/2 * sum over flagged units of
#' |p_bi - p_ai|,}
#' where \eqn{F} is the set of units flagged at `alpha` without a multiple
#' testing correction (units empty in both patterns are never flagged).
#'
#' \eqn{D} is the index of dissimilarity (Duncan and Duncan 1955): the share
#' of either pattern that would have to move to other units for the two to
#' match, the redistribution reading of Hoover (1936). \eqn{S_D} counts only
#' the flagged units. It lies in \[0, 1\], is symmetric in base and test,
#' and is at least \eqn{1 - D}.
#'
#' Over all units the shares where the test pattern is higher and where the
#' base pattern is higher balance, so \eqn{D} is exactly the share to move.
#' Over the flagged units alone they need not balance. Write \eqn{H_b} for
#' the sum of \eqn{p_{bi} - p_{ai}} over flagged units where the test
#' pattern is higher and \eqn{H_a} for the sum of \eqn{p_{ai} - p_{bi}}
#' where the base pattern is higher. Removing the differences at the flagged
#' units moves \eqn{\max(H_a, H_b)} of either pattern, the excess going to
#' or coming from unflagged units, while \eqn{1 - S_D = (H_a + H_b)/2}. So
#' \eqn{1 - S_D} approximates, and never exceeds, the share to move; the two
#' are equal when \eqn{H_a = H_b}, as when every unit is flagged.
#' `share_to_move` and the events to move use \eqn{\max(H_a, H_b)}.
#'
#' The two indices answer different questions. S tests sameness: it is
#' compared with \eqn{1 - \alpha}, what two realisations of one process give,
#' and says nothing about how far apart two different patterns are.
#' \eqn{S_D} measures that distance in events. Where power is limited many
#' genuine differences are not flagged, so \eqn{S_D} overstates similarity;
#' `share_d_flagged`, \eqn{(1 - S_D)/D}, shows how much of the total
#' dissimilarity the flagged units hold. [sppt()] prints a note when it is
#' below 0.6. That cut-off is a judgement call, not a calibrated threshold.
#' \eqn{D} itself is not a test: two realisations of one process give
#' \eqn{D > 0} through sampling variation alone, the more so the fewer events
#' per unit.
#'
#' @param result A result returned by [sppt()] with p-values.
#'
#' @return A one-row data frame with columns `s_index` (robust S),
#'   `s_index_d` (\eqn{S_D}), `dissimilarity_d` (\eqn{D}), `share_d_flagged`,
#'   `n_flagged` (units flagged without a correction), `n_eligible`,
#'   `flagged_test_higher` (\eqn{H_b}), `flagged_base_higher`
#'   (\eqn{H_a}), `share_to_move` (\eqn{\max(H_a, H_b)}), and
#'   `events_to_move_base` and `events_to_move_test` (`share_to_move` times
#'   each pattern's total).
#'
#' @references
#' Duncan, O.D. and Duncan, B. (1955). A methodological analysis of
#' segregation indexes. *American Sociological Review*, 20, 210-217.
#' \doi{10.2307/2088328}
#'
#' Hoover, E.M. (1936). The measurement of industrial localization. *The
#' Review of Economic Statistics*, 18, 162-171.
#'
#' @seealso [sppt()], [sppt_sameness()]
#' @examples
#' res <- sppt(vancouver_ct, group_col = "id",
#'             count_col = c("BNER2021", "BNEC2021"),
#'             B = 999, seed = 1, verbose = FALSE, create_maps = FALSE)
#' sppt_similarity(res)
#'
#' @export
sppt_similarity <- function(result) {
  if (!inherits(result, "sppt_result"))
    stop("`result` must be returned by sppt().", call. = FALSE)
  if (!"pvalue" %in% names(result))
    stop("`result` has no p-values; run sppt() with ",
         "calculate_pvalues = TRUE.", call. = FALSE)
  s <- attr(result, "sppt_summary")
  # Counts may be the expected counts of a random pattern, so not integers
  cnt <- function(nm) {
    x <- as.numeric(result[[nm]])
    x[is.na(x)] <- 0
    x
  }
  counts_a <- cnt(s$base_name)
  counts_b <- cnt(s$test_name)
  # Uncorrected flags. Without a correction these are the flags of the
  # result itself; with one, p <= alpha on the raw p-values.
  flag <- if (is.null(s$correction) || s$correction == "none") {
    !is.na(result$sig_pvalue) & result$sig_pvalue == 1L
  } else {
    !is.na(result$pvalue) & result$pvalue <= s$alpha
  }
  sim <- similarity_d(counts_a, counts_b, flag)
  data.frame(s_index = s$s_index_pv,
             s_index_d = sim$s_index_d,
             dissimilarity_d = sim$dissimilarity_d,
             share_d_flagged = sim$share_d_flagged,
             n_flagged = sim$n_flagged,
             n_eligible = s$n_eligible,
             flagged_test_higher = sim$flagged_test_higher,
             flagged_base_higher = sim$flagged_base_higher,
             share_to_move = sim$share_to_move,
             events_to_move_base = sim$events_to_move[["base"]],
             events_to_move_test = sim$events_to_move[["test"]])
}

# S_D, D and the share of D at flagged units, from counts and flags.
similarity_d <- function(counts_a, counts_b, flag) {
  N_a <- sum(counts_a)
  N_b <- sum(counts_b)
  d    <- counts_b / N_b - counts_a / N_a
  absd <- abs(d)
  flag <- !is.na(flag) & flag
  D    <- 0.5 * sum(absd)
  S_D  <- 1 - 0.5 * sum(absd[flag])
  # Flagged share differences in each direction. They need not balance, so
  # removing the flagged differences moves the larger of the two, which
  # 1 - S_D (their average) approximates from below.
  H_b  <- sum(d[flag & d > 0])
  H_a  <- -sum(d[flag & d < 0])
  move <- max(H_a, H_b)
  list(s_index_d           = S_D,
       dissimilarity_d     = D,
       share_d_flagged     = if (D > 0) (1 - S_D) / D else NA_real_,
       n_flagged           = sum(flag),
       flagged_test_higher = H_b,
       flagged_base_higher = H_a,
       share_to_move       = move,
       events_to_move      = c(base = move * N_a, test = move * N_b))
}

# Console lines for S_D, shared by sppt() and its summary() method.
similarity_lines <- function(sim, base_name, test_name, correction = "none") {
  out <- field_lines("S_D", sprintf(
    "%.4f (D = %.4f; flagged units hold %.0f%% of D%s)",
    sim$s_index_d, sim$dissimilarity_d, 100 * sim$share_d_flagged,
    if (correction != "none") "; uncorrected flags" else ""))
  out <- paste0(out, field_lines("Flagged shares", sprintf(
    "%s higher by %.4f, %s higher by %.4f",
    test_name, sim$flagged_test_higher, base_name, sim$flagged_base_higher)))
  out <- paste0(out, field_lines("Events to move", sprintf(
    "%.1f%% of either pattern: %s %s, or %s %s",
    100 * sim$share_to_move,
    format(round(sim$events_to_move[["base"]]), big.mark = ","), base_name,
    format(round(sim$events_to_move[["test"]]), big.mark = ","), test_name)))
  if (!is.na(sim$share_d_flagged) && sim$share_d_flagged < 0.6)
    out <- paste0(out, field_lines("Note", sprintf(paste(
      "flagged units hold only %.0f%% of the total dissimilarity; S_D",
      "overstates similarity (limited power)"), 100 * sim$share_d_flagged)))
  out
}
