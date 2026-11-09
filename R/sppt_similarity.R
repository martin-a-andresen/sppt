#' Degree of similarity between two spatial patterns
#'
#' @description
#' Returns the degree-of-similarity measures of a [sppt()] result: the S
#' index, \eqn{S_D}, the index of dissimilarity \eqn{D}, the share of
#' \eqn{D} held by flagged units, and the number of events of each pattern
#' that would have to move.
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
#' the redistribution at flagged units, so it is the share of each pattern
#' already in place once only significant differences count. It lies in
#' \[0, 1\], is symmetric in base and test, and is at least \eqn{1 - D}.
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
#'   `events_to_move_base` and `events_to_move_test`
#'   (\eqn{(1 - S_D)} times each pattern's total).
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
  cnt <- function(nm) {
    x <- as.integer(result[[nm]])
    x[is.na(x)] <- 0L
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
             events_to_move_base = sim$events_to_move[["base"]],
             events_to_move_test = sim$events_to_move[["test"]])
}

# S_D, D and the share of D at flagged units, from counts and flags.
similarity_d <- function(counts_a, counts_b, flag) {
  N_a <- sum(counts_a)
  N_b <- sum(counts_b)
  absd <- abs(counts_b / N_b - counts_a / N_a)
  flag <- !is.na(flag) & flag
  D    <- 0.5 * sum(absd)
  S_D  <- 1 - 0.5 * sum(absd[flag])
  list(s_index_d       = S_D,
       dissimilarity_d = D,
       share_d_flagged = if (D > 0) (1 - S_D) / D else NA_real_,
       n_flagged       = sum(flag),
       events_to_move  = c(base = (1 - S_D) * N_a, test = (1 - S_D) * N_b))
}

# Console lines for S_D, shared by sppt() and its summary() method.
similarity_lines <- function(sim, base_name, test_name, correction = "none") {
  out <- field_lines("S_D", sprintf(
    "%.4f (D = %.4f; flagged units hold %.0f%% of D%s)",
    sim$s_index_d, sim$dissimilarity_d, 100 * sim$share_d_flagged,
    if (correction != "none") "; uncorrected flags" else ""))
  out <- paste0(out, field_lines("Events to move", sprintf(
    "%s %s, or %s %s",
    format(round(sim$events_to_move[["base"]]), big.mark = ","), base_name,
    format(round(sim$events_to_move[["test"]]), big.mark = ","), test_name)))
  if (!is.na(sim$share_d_flagged) && sim$share_d_flagged < 0.6)
    out <- paste0(out, field_lines("Note", sprintf(paste(
      "flagged units hold only %.0f%% of the total dissimilarity; S_D",
      "overstates similarity (limited power)"), 100 * sim$share_d_flagged)))
  out
}
