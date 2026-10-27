#' Spatial point pattern test with bootstrap p-values
#'
#' @description
#' Compares two spatial point patterns, a base and a test, across a set of
#' areal units. Each unit is tested with a bounded Monte Carlo bootstrap
#' p-value and, optionally, with marginal confidence interval overlap
#' (Andresen 2009). Both are derived from a single shared set of `B`
#' bootstrap replicates.
#'
#' Global similarity is summarised by the robust S index, the proportion of
#' eligible units (those not empty in both patterns) that are not flagged as
#' different. It is accompanied by two further tests:
#'
#' * a **sameness test** ([sppt_sameness()]), which asks whether the two
#'   patterns are more different than two realisations of a single process;
#' * an **arrangement test**, global Moran's I on the signed evidence surface
#'   \eqn{\mathrm{sign}(d) \times -\log_{10}(p)}, which asks whether the
#'   differences, however many there are, are spatially organised. Local
#'   Moran clusters are outlined on the p-value map, and each flagged unit is
#'   given a confidence tier from its position relative to those clusters.
#'
#' @details
#' **P-values.** For each unit the observed difference in proportions
#' (test minus base) is compared with the centred bootstrap distribution of
#' that difference. The p-value is the bounded Monte Carlo form
#' \deqn{p = \frac{1 + r}{B + 1},}{p = (1 + r) / (B + 1),}
#' where \eqn{r} counts centred bootstrap differences at least as large in
#' absolute value as the observed difference. The smallest attainable
#' p-value is therefore \eqn{1/(B + 1)}. Units with no events in either
#' pattern are excluded from inference and receive `NA`.
#'
#' The bootstrap resamples events at their observed locations and
#' re-aggregates them to the areal units; drawing counts from the
#' multinomial distribution of observed shares is the closed-form equivalent.
#' Locations are never permuted, so the null is sampling variability around
#' the observed intensity surface, not spatial randomness.
#'
#' **Multiple testing.** No correction is applied. Simulation evidence
#' shows that corrections remove power at fine spatial scales without a
#' meaningful benefit in false discovery rate, and that at coarse scales a
#' stricter `alpha` achieves the same as a correction. The sameness and
#' arrangement tests and [sppt_clusters()] address the global question
#' instead. The raw p-values are returned, so a correction can still be
#' applied with [stats::p.adjust()] if required.
#'
#' **Power.** `lambda_min`, the mean number of events per unit in the rarer
#' pattern, governs power. Below 1.5 a low-power warning is issued and below
#' 0.5 a critical warning; at those resolutions, read core and edge units
#' rather than isolated flags, test the clusters with [sppt_clusters()], or
#' aggregate to coarser units.
#'
#' **Confidence tiers.** With `arrangement = TRUE`, each flagged unit is
#' placed in one tier:
#' * `"core"`: inside a local Moran outline of the same direction;
#' * `"edge"`: not core, but sharing a border (queen contiguity) with a core
#'   unit whose difference has the same sign;
#' * `"outside"`: any other flagged unit.
#'
#' Every tier is significant at `alpha`; the tiers add spatial context and
#' rank confidence, they do not certify a difference as real. In simulation
#' (21 crime type pairs, 7 spatial scales, 200 iterations, `alpha = 0.05`,
#' medians across pairs) the share of flagged units that were true
#' differences was 0.34 to 0.88 among all flagged units, 0.85 to 0.98 among
#' core units and 0.81 to 0.97 among core plus edge units. Edge units are
#' least informative when one direction dominates, because a widespread
#' shift places flagged units next to clusters by chance; the console notes
#' this case.
#'
#' Local Moran classifies units around the mean of the signed surface, so a
#' unit is outlined as "test > base" only if its own difference is positive,
#' and as "base > test" only if it is negative.
#'
#' **Relation to published descriptions.** The original test (Andresen 2009)
#' compared confidence intervals only; p-values, the sameness and
#' arrangement tests and the tiers are later additions. The p-value here is
#' the bounded form \eqn{(1 + r)/(B + 1)} rather than \eqn{r/B}.
#'
#' @param data An `sf` object with one row per spatial unit.
#' @param group_col Character. Name of the column that uniquely identifies
#'   each spatial unit (for example `"DAUID"` or `"CTUID"`).
#' @param count_col Character vector of length 2 giving the base and test
#'   count columns, in that order (for example `c("BNER2003", "BNER2025")`).
#' @param B Integer. Number of bootstrap replicates. Default `4999`; an odd
#'   value avoids ties.
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#' @param alpha Numeric in (0, 1). Significance level for flagging a unit.
#'   Default `0.05`. Where control of false discoveries is the priority and
#'   the analysis is at census tracts or an equally coarse scale
#'   (`lambda_min` above roughly 10 to 15), use `0.01`.
#' @param calculate_pvalues Logical. Compute bootstrap p-values. Default
#'   `TRUE`.
#' @param sameness_tau Numeric in (0, 1). Significance level for the sameness
#'   test on the robust S index. Separate from `alpha`, which governs the
#'   flagging of individual units. Default `0.05`.
#' @param arrangement Logical. Run the arrangement test (global and local
#'   Moran's I on the signed evidence surface) and assign confidence tiers.
#'   Requires polygon geometry. Default `TRUE`.
#' @param nsim_moran Integer. Permutations for the global Moran's I test.
#'   Default `999`.
#' @param calculate_ci Logical. Also compute marginal bootstrap confidence
#'   intervals and the CI-based S index (Andresen 2009). Default `FALSE`.
#' @param conf_level Numeric in (0, 1). Confidence level for the interval
#'   overlap criterion. Default `0.85`: two marginal 85% intervals fail to
#'   overlap when the difference exceeds about 2.04 standard errors, which
#'   corresponds to a two-sample test at about `alpha = 0.05` (Andresen
#'   2026). `0.95` reproduces the original, more conservative analysis. Only
#'   used when `calculate_ci = TRUE`.
#' @param fix_base Logical. If `TRUE`, the base pattern is not bootstrapped
#'   and its observed proportions are treated as fixed. This gives an
#'   asymmetric test, intended only for comparison against a known reference
#'   distribution. Default `FALSE`.
#' @param use_percentages Logical. Express proportions as percentages
#'   (0-100) rather than as proportions (0-1). Default `TRUE`.
#' @param verbose Logical. Print the test summary to the console. Default
#'   `TRUE`.
#' @param create_maps Logical. Build the maps and print them. The maps are
#'   always stored on the result (see Value) when this is `TRUE`. Default
#'   `TRUE`.
#' @param map_crs CRS for the maps, passed to [sf::st_transform()]. `NULL`
#'   uses the CRS of `data`.
#' @param export_maps Logical. Save the maps as PNG files in `export_dir`.
#'   Default `FALSE`.
#' @param export_dir Character. Directory for exported maps and results.
#'   Required when `export_maps` or `export_results` is `TRUE`; created if it
#'   does not exist.
#' @param map_dpi Integer. Resolution of exported maps. Default `300`.
#' @param export_results Logical. Write the result to `export_results_dir`.
#'   Default `FALSE`.
#' @param export_format Character. Format for exported results: `"gpkg"`
#'   (default), `"shp"`, `"csv"`, `"txt"` or `"rds"`. Only `"rds"` keeps the
#'   summary and the maps.
#' @param export_results_dir Character or `NULL`. Directory for exported
#'   results. Defaults to `export_dir`.
#'
#' @return An `sf` object of class `sppt_result` containing the columns of
#'   `data` plus:
#' \describe{
#'   \item{`base_prop`, `test_prop`}{Observed share of base and test events
#'     in each unit.}
#'   \item{`difference`}{`test_prop` minus `base_prop`.}
#'   \item{`lambda_min`}{Mean events per unit for the rarer pattern (the same
#'     for every row).}
#'   \item{`low_power_warning`}{`TRUE` when `lambda_min < 1.5`.}
#'   \item{`pvalue`}{Bootstrap p-value; `NA` for units empty in both
#'     patterns.}
#'   \item{`sig_pvalue`}{1 if `pvalue <= alpha`, otherwise 0.}
#'   \item{`sig_pvalue_dir`}{1 = test > base, -1 = base > test, 0 = not
#'     significant.}
#'   \item{`signed_evidence`}{\eqn{\mathrm{sign}(d) \times -\log_{10}(p)}.}
#'   \item{`moran_hh`, `moran_ll`, `moran_cluster`}{Local Moran clusters of
#'     the signed evidence surface (arrangement test only).}
#'   \item{`confidence_tier`}{`"core"`, `"edge"` or `"outside"` for flagged
#'     units, `NA` otherwise (arrangement test only).}
#'   \item{`ci_lower`, `ci_upper`, `ci_lower_test`, `ci_upper_test`}{Marginal
#'     confidence bounds for base and test (CI only).}
#'   \item{`intervals_overlap`}{1 = intervals overlap (similar), 0 = they do
#'     not (CI only).}
#'   \item{`sig_ci_dir`}{Direction of CI-based differences (CI only).}
#' }
#' The attribute `"sppt_summary"` holds the robust S index, the sameness and
#' arrangement tests, the tier counts and the settings; `summary()` prints
#' it. The attribute `"sppt_maps"` holds the maps as a named list of
#' `ggplot` objects (`difference`, `pvalue`, `tiers` with the arrangement
#' test, and `ci` with CI); use [sppt_map()] or `plot()` to redraw them with
#' a different layout.
#'
#' @references
#' Andresen, M.A. (2009). Testing for similarity in area-based spatial
#' patterns: A nonparametric Monte Carlo approach. *Applied Geography*,
#' 29(3), 333-345. \doi{10.1016/j.apgeog.2008.12.004}
#'
#' Andresen, M.A. (2016). An area-based nonparametric spatial point pattern
#' test: The test, its applications, and the future. *Methodological
#' Innovations*, 9, Article 12. \doi{10.1177/2059799116630659}
#'
#' Andresen, M.A. (2026). Testing longitudinal change in spatial
#' distributions and patterns: An extension of the spatial point pattern
#' test. *Annals of the American Association of Geographers*, in press.
#' \doi{10.1080/24694452.2026.2733304}
#'
#' @seealso [sppt_longitudinal()], [sppt_sameness()], [sppt_clusters()],
#'   [sppt_map()], [vancouver_ct]
#' @examples
#' res <- sppt(vancouver_ct, group_col = "id",
#'             count_col = c("BNER2003", "BNER2025"),
#'             B = 199, seed = 42, create_maps = FALSE)
#' summary(res)
#' table(res$sig_pvalue_dir)
#' table(res$confidence_tier, useNA = "ifany")
#'
#' @export
sppt <- function(data,
                 group_col          = "group",
                 count_col,
                 B                  = 4999L,
                 seed               = NULL,
                 alpha              = 0.05,
                 calculate_pvalues  = TRUE,
                 sameness_tau       = 0.05,
                 arrangement        = TRUE,
                 nsim_moran         = 999L,
                 calculate_ci       = FALSE,
                 conf_level         = 0.85,
                 fix_base           = FALSE,
                 use_percentages    = TRUE,
                 verbose            = TRUE,
                 create_maps        = TRUE,
                 map_crs            = NULL,
                 export_maps        = FALSE,
                 export_dir         = NULL,
                 map_dpi            = 300L,
                 export_results     = FALSE,
                 export_format      = "gpkg",
                 export_results_dir = NULL) {

  say <- say_if(verbose)

  # -- 1. Input validation ----------------------------------------------------
  if (!inherits(data, "sf"))
    stop("`data` must be an sf object.", call. = FALSE)

  if (!group_col %in% names(data))
    stop("`group_col` '", group_col, "' not found in data.", call. = FALSE)

  if (anyDuplicated(data[[group_col]]) > 0L)
    stop("`group_col` '", group_col, "' does not uniquely identify rows.",
         call. = FALSE)

  if (length(count_col) != 2L)
    stop("`count_col` must be a character vector of length 2 (base, test).",
         call. = FALSE)

  for (col in count_col) {
    if (!col %in% names(data))
      stop("`count_col` column '", col, "' not found in data.", call. = FALSE)
    vals <- data[[col]]
    if (!is.numeric(vals))
      stop("`count_col` column '", col, "' must be numeric.", call. = FALSE)
    if (any(!is.na(vals) & vals < 0))
      stop("`count_col` column '", col, "' contains negative values.",
           call. = FALSE)
  }

  if (!is.numeric(B) || length(B) != 1L || B < 99)
    stop("`B` must be a single integer >= 99.", call. = FALSE)
  B <- as.integer(B)

  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1)
    stop("`alpha` must be a single numeric value in (0, 1).", call. = FALSE)

  if (!is.numeric(sameness_tau) || length(sameness_tau) != 1L ||
      sameness_tau <= 0 || sameness_tau >= 1)
    stop("`sameness_tau` must be a single numeric value in (0, 1).",
         call. = FALSE)

  if (!is.numeric(nsim_moran) || length(nsim_moran) != 1L || nsim_moran < 99L)
    stop("`nsim_moran` must be a single integer >= 99.", call. = FALSE)
  nsim_moran <- as.integer(nsim_moran)

  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be a single numeric value in (0, 1).",
         call. = FALSE)

  if (export_maps) export_dir <- check_export_dir(export_dir, "export_maps")
  if (export_results) {
    if (is.null(export_results_dir)) export_results_dir <- export_dir
    export_results_dir <- check_export_dir(export_results_dir,
                                           "export_results")
  }

  if (fix_base)
    warning(
      "fix_base = TRUE produces an asymmetric test: the base variable is not\n",
      "bootstrapped, which is inconsistent with the symmetric design principle\n",
      "of SPPT. It also yields a degenerate reference distribution at any unit\n",
      "where the test variable has no events, since the bootstrap replicates\n",
      "for that unit are then identically zero and the p-value takes its\n",
      "smallest attainable value, 1/(B+1), regardless of alpha. The sameness\n",
      "test on the S index is not calibrated for this variant. Not\n",
      "recommended for standard use.",
      call. = FALSE
    )

  # -- 2. Setup ---------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)

  base_name <- count_col[1L]
  test_name <- count_col[2L]
  n_units   <- nrow(data)

  counts_a <- as.integer(data[[base_name]])
  counts_b <- as.integer(data[[test_name]])
  counts_a[is.na(counts_a)] <- 0L
  counts_b[is.na(counts_b)] <- 0L

  N_a <- sum(counts_a)
  N_b <- sum(counts_b)

  # -- 3. lambda_min and decision tier ----------------------------------------
  lambda_min <- min(N_a, N_b) / n_units

  decision_tier <- if (lambda_min >= 17) {
    "Coarse scale (lambda_min >= 17): CI and PV methods both reliable.\n    PV alpha = 0.05; FDR-priority use alpha = 0.01."
  } else if (lambda_min >= 6) {
    "Medium-high scale (6 <= lambda_min < 17): PV alpha = 0.05 recommended;\n    FDR-priority use alpha = 0.01 once lambda_min is above about 10."
  } else if (lambda_min >= 1.5) {
    "Medium scale (1.5 <= lambda_min < 6): PV alpha = 0.05 recommended."
  } else if (lambda_min >= 0.5) {
    "Fine scale (0.5 <= lambda_min < 1.5): Power is limited; interpret\n    core and edge units rather than isolated flags."
  } else {
    "Very fine scale (lambda_min < 0.5): All methods have severely limited\n    power. Consider aggregating to a coarser spatial unit."
  }

  if (lambda_min < 0.5) {
    warning(sprintf(paste0(
      "\n*** CRITICAL POWER WARNING ***\n",
      "lambda_min = %.3f (< 0.5). At this spatial resolution, all methods\n",
      "have severely limited statistical power and FDR exceeds power under\n",
      "all significance thresholds. Results should not be interpreted as\n",
      "reliable spatial pattern estimates. Consider aggregating to a coarser\n",
      "spatial unit where lambda_min >= 1.5."), lambda_min), call. = FALSE)
  } else if (lambda_min < 1.5) {
    warning(sprintf(paste0(
      "\n*** LOW POWER WARNING ***\n",
      "lambda_min = %.3f (< 1.5). At this spatial resolution, statistical\n",
      "power is limited. P-values will be computed but should be interpreted\n",
      "cautiously. What the test detects at this resolution is more spatially\n",
      "coherent than chance, while false positives are scattered: interpret\n",
      "core and edge units (column confidence_tier) rather than isolated\n",
      "flags, test the clusters with sppt_clusters(), or aggregate to a\n",
      "coarser spatial unit where lambda_min >= 1.5."),
      lambda_min), call. = FALSE)
  }

  # -- 4. Console header ------------------------------------------------------
  say("\n")
  say("==================================================================\n")
  say("Spatial Point Pattern Test\n")
  say("==================================================================\n")
  say(sprintf("  Base variable : %s  (N = %d)\n", base_name, N_a))
  say(sprintf("  Test variable : %s  (N = %d)\n", test_name, N_b))
  say(sprintf("  Spatial units : %d\n", n_units))
  say(sprintf("  lambda_min    : %.3f\n", lambda_min))
  say(sprintf("  Bootstrap B   : %d\n", B))
  if (calculate_pvalues)
    say(sprintf("  P-value alpha : %.3f\n", alpha))
  if (calculate_ci)
    say(sprintf("  CI level      : %.1f%%  (marginal overlap criterion)\n",
                conf_level * 100))
  say(sprintf("  Decision tier : %s\n", decision_tier))
  say("------------------------------------------------------------------\n")

  # -- 5. Observed proportions and difference ---------------------------------
  scale_factor <- if (use_percentages) 100 else 1

  prop_a <- if (N_a > 0) counts_a / N_a * scale_factor else rep(0, n_units)
  prop_b <- if (N_b > 0) counts_b / N_b * scale_factor else rep(0, n_units)
  obs_diff <- prop_b - prop_a   # test minus base

  # Units empty in both patterns are excluded from inference
  zero_mask <- counts_a == 0L & counts_b == 0L
  n_zero    <- sum(zero_mask)

  # -- 6. Shared bootstrap ----------------------------------------------------
  # One multinomial draw per variable, B replicates, shared by CI and PV.
  # Drawing counts from the multinomial of observed shares is the closed form
  # of resampling events at their observed locations and re-aggregating.
  if (N_a == 0L || N_b == 0L) {
    warning("One or both count columns sum to zero. No inference possible.",
            call. = FALSE)
    boot_a <- boot_b <- boot_diffs_c <- matrix(0, nrow = n_units, ncol = B)
    valid_bootstrap <- FALSE
  } else {
    p_a <- counts_a / N_a
    p_b <- counts_b / N_b

    boot_a <- if (fix_base) {
      matrix(rep(p_a * scale_factor, B), nrow = n_units, ncol = B)
    } else {
      matrix(stats::rmultinom(B, size = N_a, prob = p_a),
             nrow = n_units) / N_a * scale_factor
    }
    boot_b <- matrix(stats::rmultinom(B, size = N_b, prob = p_b),
                     nrow = n_units) / N_b * scale_factor

    # Centred difference distribution for p-values
    boot_diffs   <- boot_b - boot_a
    boot_diffs_c <- boot_diffs - rowMeans(boot_diffs)
    valid_bootstrap <- TRUE
  }

  # -- 7. Bootstrap p-values --------------------------------------------------
  pvalue     <- rep(NA_real_, n_units)
  sig_pvalue <- rep(NA_integer_, n_units)
  sig_pv_dir <- rep(NA_integer_, n_units)

  if (calculate_pvalues && valid_bootstrap) {
    active <- !zero_mask
    # Bounded Monte Carlo p-value: p = (1 + r) / (B + 1), where r counts the
    # centred bootstrap differences at least as large in absolute value as
    # the observed difference. The unbounded form r / B admits p = 0, which
    # is anti-conservative and unbounded on the -log10 scale.
    pvalue[active] <- (1 + rowSums(
      abs(boot_diffs_c[active, , drop = FALSE]) >= abs(obs_diff[active])
    )) / (B + 1)

    sig_pvalue[active] <- as.integer(pvalue[active] <= alpha)
    sig_pv_dir[active] <- ifelse(sig_pvalue[active] == 1L,
                                 as.integer(sign(obs_diff[active])), 0L)
  }

  # -- 8. Marginal CI overlap (optional) --------------------------------------
  ci_lower <- ci_upper <- ci_lower_test <- ci_upper_test <-
    rep(NA_real_, n_units)
  intervals_overlap <- sig_ci_dir <- rep(NA_integer_, n_units)

  if (calculate_ci && valid_bootstrap) {
    ql <- (1 - conf_level) / 2
    qu <- 1 - ql
    active <- !zero_mask

    ci_lower[active]      <- apply(boot_a[active, , drop = FALSE], 1,
                                   stats::quantile, probs = ql)
    ci_upper[active]      <- apply(boot_a[active, , drop = FALSE], 1,
                                   stats::quantile, probs = qu)
    ci_lower_test[active] <- apply(boot_b[active, , drop = FALSE], 1,
                                   stats::quantile, probs = ql)
    ci_upper_test[active] <- apply(boot_b[active, , drop = FALSE], 1,
                                   stats::quantile, probs = qu)

    # Overlap: 1 = similar, 0 = dissimilar
    intervals_overlap[active] <- as.integer(
      pmax(ci_lower[active], ci_lower_test[active]) <=
        pmin(ci_upper[active], ci_upper_test[active])
    )
    sig_ci_dir[active] <- ifelse(intervals_overlap[active] == 0L,
                                 as.integer(sign(obs_diff[active])), 0L)
  }

  # -- 9. Robust S index and sameness test ------------------------------------
  # The denominator excludes zero-zero units, so the robust index equals the
  # plain index whenever every unit contains at least one event.
  n_eligible <- n_units - n_zero

  s_index_pv <- NA_real_
  n_sig_pv   <- NA_integer_
  sameness   <- NULL
  if (calculate_pvalues && valid_bootstrap && n_eligible > 0L) {
    n_sig_pv   <- sum(sig_pvalue[!zero_mask] == 1L, na.rm = TRUE)
    s_index_pv <- 1 - n_sig_pv / n_eligible
    # Only calibrated under the symmetric bootstrap: fix_base = TRUE breaks
    # the null flagging rate at sparse units.
    if (!fix_base)
      sameness <- sppt_sameness(s_index_pv, n_sig_pv, n_eligible, alpha,
                                tau = sameness_tau)
  }

  # No sameness test for the CI-based index: the null flagging rate of the
  # interval-overlap rule is not 1 - conf_level.
  s_index_ci <- NA_real_
  n_sig_ci   <- NA_integer_
  if (calculate_ci && valid_bootstrap && n_eligible > 0L) {
    n_sig_ci   <- sum(intervals_overlap[!zero_mask] == 0L, na.rm = TRUE)
    s_index_ci <- 1 - n_sig_ci / n_eligible
  }

  # -- 10. Spatial arrangement and confidence tiers ---------------------------
  # Under sameness the per-unit p-values are uniform, so the signed surface
  # sign(difference) x -log10(p) carries no spatial structure. High-high marks
  # clusters of "test > base", low-low clusters of "base > test". Zero-zero
  # units are dropped rather than filled, since a block of identical filled
  # values would itself register as autocorrelation.
  arrange_res <- NULL
  if (isTRUE(arrangement) && calculate_pvalues && valid_bootstrap &&
      n_eligible > 2L) {
    arrange_res <- tryCatch(
      arrangement_test(data, obs_diff, pvalue, zero_mask, sig_pv_dir,
                       nsim_moran),
      error = function(e) {
        message("Arrangement test failed: ", conditionMessage(e))
        NULL
      })
  }

  # Widespread shift: one direction holds at least 80% of flagged units and
  # at least half of those lie outside any cluster. Clusters can then sit on
  # top of a citywide difference, and edge units carry less information.
  wide_shift <- FALSE
  if (!is.null(arrange_res)) {
    n_up <- sum(sig_pv_dir ==  1L, na.rm = TRUE)
    n_dn <- sum(sig_pv_dir == -1L, na.rm = TRUE)
    n_flag_tot <- n_up + n_dn
    if (n_flag_tot > 0L) {
      dom   <- if (n_up >= n_dn) 1L else -1L
      n_dom <- max(n_up, n_dn)
      n_dom_out <- sum(sig_pv_dir == dom & !is.na(arrange_res$tier) &
                         arrange_res$tier == "outside", na.rm = TRUE)
      wide_shift <- n_dom / n_flag_tot >= 0.8 && n_dom_out / n_dom >= 0.5
      arrange_res$dominant_share <- n_dom / n_flag_tot
    }
  }

  # -- 11. Append result columns ----------------------------------------------
  result <- data
  result$base_prop         <- round(prop_a, 6)
  result$test_prop         <- round(prop_b, 6)
  result$difference        <- round(obs_diff, 6)
  result$lambda_min        <- lambda_min
  result$low_power_warning <- lambda_min < 1.5

  if (calculate_pvalues) {
    result$pvalue          <- round(pvalue, 6)
    result$sig_pvalue      <- sig_pvalue
    result$sig_pvalue_dir  <- sig_pv_dir
    result$signed_evidence <- round(
      sign(obs_diff) * -log10(pmax(pvalue, .Machine$double.eps)), 6)
  }
  if (!is.null(arrange_res)) {
    result$moran_hh        <- arrange_res$hh
    result$moran_ll        <- arrange_res$ll
    result$moran_cluster   <- arrange_res$cluster
    result$confidence_tier <- arrange_res$tier
  }
  if (calculate_ci) {
    result$ci_lower          <- round(ci_lower, 6)
    result$ci_upper          <- round(ci_upper, 6)
    result$ci_lower_test     <- round(ci_lower_test, 6)
    result$ci_upper_test     <- round(ci_upper_test, 6)
    result$intervals_overlap <- intervals_overlap
    result$sig_ci_dir        <- sig_ci_dir
  }

  # -- 12. Console results ----------------------------------------------------
  if (calculate_pvalues && !is.na(s_index_pv)) {
    say(sprintf("  Robust S-index  : %.4f  (%d of %d units dissimilar at alpha = %.3f)\n",
                s_index_pv, n_sig_pv, n_eligible, alpha))
    say(sprintf("    %s > %s : %d units\n", test_name, base_name,
                sum(sig_pv_dir[!zero_mask] == 1L, na.rm = TRUE)))
    say(sprintf("    %s > %s : %d units\n", base_name, test_name,
                sum(sig_pv_dir[!zero_mask] == -1L, na.rm = TRUE)))
  }

  if (!is.null(sameness) && !is.na(sameness$p_value)) {
    say(sprintf("  Sameness test   : S = %.4f vs %.4f expected, p %s (%s)\n",
                sameness$s_index, sameness$s_expected,
                fmt_p(sameness$p_value), sameness$method))
    say(sprintf("    Threshold     : S >= %.4f at tau = %.2f\n",
                sameness$threshold, sameness$tau))
    say(sprintf("    Verdict       : %s\n",
                if (isTRUE(sameness$same))
                  "no more different than two realisations of one process"
                else
                  "more different than two realisations of one process"))
    if (n_eligible < n_units * 0.75)
      say("    Note          : many units empty in both patterns; the sameness\n",
          "                    test is conservative here.\n", sep = "")
  } else if (calculate_pvalues && !is.na(s_index_pv) && fix_base) {
    say("  Sameness test   : not reported (not calibrated for fix_base = TRUE).\n")
  }

  if (calculate_ci && !is.na(s_index_ci)) {
    say(sprintf("  Robust CI-index : %.4f  (%d of %d units dissimilar at CI = %.1f%%)\n",
                s_index_ci, n_sig_ci, n_eligible, conf_level * 100))
    say(sprintf("    %s > %s : %d units\n", test_name, base_name,
                sum(sig_ci_dir[!zero_mask] == 1L, na.rm = TRUE)))
    say(sprintf("    %s > %s : %d units\n", base_name, test_name,
                sum(sig_ci_dir[!zero_mask] == -1L, na.rm = TRUE)))
  }

  if (!is.null(arrange_res)) {
    say(sprintf("  Arrangement     : Moran's I = %+.4f on signed evidence, p %s (%d perms)\n",
                arrange_res$moran_I, fmt_p(arrange_res$moran_p),
                arrange_res$nsim))
    say(sprintf("    Local clusters: %d units where %s > %s, %d units where %s > %s\n",
                arrange_res$n_hh, test_name, base_name,
                arrange_res$n_ll, base_name, test_name))
    say(sprintf("  Flagged units   : %d core, %d edge, %d outside\n",
                arrange_res$n_core, arrange_res$n_edge,
                arrange_res$n_outside))
    say("    core = inside a cluster; edge = borders a core unit with the same\n",
        "    sign; outside = all other flagged units (least reliable)\n",
        sep = "")
    if (wide_shift && arrange_res$n_edge > 0L)
      say(sprintf(paste0(
        "    Note: %.0f%% of flagged units share one direction. Edge units\n",
        "    are still significant, but a widespread shift places flagged\n",
        "    units beside clusters by chance; weigh core units more heavily.\n"),
        100 * arrange_res$dominant_share))
  }

  # Two-way verdict: the count question and the arrangement question are
  # separate, and the pair says more than either alone.
  if (!is.null(sameness) && !is.na(sameness$p_value) && !is.null(arrange_res)) {
    count_hi <- !isTRUE(sameness$same)
    arr_hi   <- isTRUE(arrange_res$clustered)
    mark <- function(a, b) if (count_hi == a && arr_hi == b) " <<<" else ""
    say("\n  Two-way verdict\n")
    say("                       | arrangement random    | arrangement clustered\n")
    say("    -----------------------------------------------------------------\n")
    say(sprintf("    count at chance    | same process%-10s| few differences, but concentrated%s\n",
                mark(FALSE, FALSE), mark(FALSE, TRUE)))
    say(sprintf("    count above chance | many, diffuse%-10s| differences real and localised%s\n",
                mark(TRUE, FALSE), mark(TRUE, TRUE)))
    say(sprintf("    This result        : %s, %s%s\n",
                if (count_hi) "count above chance" else "count at chance",
                if (arr_hi) "arrangement clustered" else "arrangement random",
                # Clusters can sit on top of a citywide shift in one
                # direction; "localised" alone would undersell the shift.
                if (count_hi && arr_hi && wide_shift)
                  ",\n                         plus a widespread shift" else ""))
  }

  if (n_zero > 0L)
    say(sprintf("  Note: %d zero-zero units excluded from S-index and inference.\n",
                n_zero))
  say("==================================================================\n\n")

  # -- 13. Summary attribute and class ----------------------------------------
  attr(result, "sppt_summary") <- list(
    base_name        = base_name,
    test_name        = test_name,
    n_units          = n_units,
    n_eligible       = n_eligible,
    n_zero           = n_zero,
    lambda_min       = lambda_min,
    B                = B,
    alpha            = if (calculate_pvalues) alpha else NULL,
    s_index_pv       = s_index_pv,
    n_sig_pv         = n_sig_pv,
    sameness         = sameness,
    arrangement      = if (is.null(arrange_res)) NULL else
      arrange_res[c("moran_I", "moran_p", "nsim", "n_hh", "n_ll",
                    "n_core", "n_edge", "n_outside", "clustered")],
    widespread_shift = wide_shift,
    conf_level       = if (calculate_ci) conf_level else NULL,
    s_index_ci       = s_index_ci,
    n_sig_ci         = n_sig_ci,
    use_percentages  = use_percentages,
    fix_base         = fix_base,
    decision_tier    = decision_tier
  )
  class(result) <- c("sppt_result", class(result))

  # -- 14. Maps ---------------------------------------------------------------
  if (create_maps && valid_bootstrap) {
    maps <- tryCatch(
      build_maps_sppt(result, map_crs = map_crs),
      error = function(e) {
        message("Could not create maps: ", conditionMessage(e))
        NULL
      })
    if (!is.null(maps)) {
      for (m in maps) print(m)
      if (export_maps) {
        files <- c(difference = "map1_difference.png",
                   pvalue     = "map2_pvalue_significance.png",
                   tiers      = "map2b_confidence_tiers.png",
                   ci         = "map3_ci_significance.png")
        for (nm in names(maps))
          ggplot2::ggsave(file.path(export_dir, files[[nm]]), maps[[nm]],
                          width = 8, height = if (nm == "tiers") 7 else 6,
                          dpi = map_dpi)
        say("Maps exported to:", export_dir, "\n")
      }
      attr(result, "sppt_maps") <- maps
    }
  }

  # -- 15. Export results -----------------------------------------------------
  if (export_results) {
    fp <- tryCatch(
      write_result(result, export_results_dir,
                   paste0("sppt_", paste(count_col, collapse = "_vs_")),
                   export_format),
      error = function(e) {
        warning("Export failed: ", conditionMessage(e), call. = FALSE)
        NULL
      })
    if (!is.null(fp)) say("Results exported:", fp, "\n")
  }

  result
}

# -- Internal helpers ---------------------------------------------------------

fmt_p <- function(p) if (p < 0.0001) "< 0.0001" else sprintf("= %.4f", p)

# Global and local Moran's I on the signed evidence surface, and the
# confidence tier of every flagged unit.
arrangement_test <- function(data, obs_diff, pvalue, zero_mask, sig_pv_dir,
                             nsim_moran) {
  n_units <- length(obs_diff)
  nb_all  <- spdep::poly2nb(data, queen = TRUE)
  keep    <- !zero_mask
  nb_sub  <- subset(nb_all, keep)
  lw      <- spdep::nb2listw(nb_sub, style = "W", zero.policy = TRUE)
  # The bounded p-value caps this at -log10(1/(B+1)); the pmax is a guard
  # that should never bind.
  surf <- sign(obs_diff[keep]) *
    -log10(pmax(pvalue[keep], .Machine$double.eps))

  gm <- spdep::moran.mc(surf, lw, nsim = nsim_moran,
                        zero.policy = TRUE, alternative = "greater")

  lm_res <- spdep::localmoran(surf, lw, zero.policy = TRUE)
  lagv   <- spdep::lag.listw(lw, surf, zero.policy = TRUE)
  pcol   <- if ("Pr(z != E(Ii))" %in% colnames(lm_res))
    "Pr(z != E(Ii))" else 5L
  loc_sig <- !is.na(lm_res[, pcol]) & lm_res[, pcol] <= 0.05
  mu_s <- mean(surf)
  mu_l <- mean(lagv)

  # Local Moran sorts units around the MEAN of the surface. When most units
  # lean one way the mean is far from zero, and "below the mean" can include
  # units where the test variable is still higher. Each outlined unit must
  # therefore also have the sign its label states.
  hh_sub <- (surf > mu_s) & (lagv > mu_l) & loc_sig & (surf > 0)
  ll_sub <- (surf < mu_s) & (lagv < mu_l) & loc_sig & (surf < 0)

  hh <- rep(FALSE, n_units); hh[which(keep)[hh_sub]] <- TRUE
  ll <- rep(FALSE, n_units); ll[which(keep)[ll_sub]] <- TRUE

  cl <- rep(NA_character_, n_units)
  cl[hh] <- "test > base"
  cl[ll] <- "base > test"

  # Confidence tiers. Core = flagged and outlined. Edge = flagged, not core,
  # and a queen neighbour of a core unit whose difference has the same sign.
  # Outside = any other flagged unit.
  flag <- !is.na(sig_pv_dir) & sig_pv_dir != 0L
  core <- flag & (hh | ll)
  edge <- vapply(seq_len(n_units), function(i) {
    if (!flag[i] || core[i]) return(FALSE)
    nbv <- nb_all[[i]]
    nbv <- nbv[nbv > 0L]
    any(core[nbv] & sig_pv_dir[nbv] == sig_pv_dir[i])
  }, logical(1))
  tier <- rep(NA_character_, n_units)
  tier[flag] <- "outside"
  tier[edge] <- "edge"
  tier[core] <- "core"

  list(moran_I   = as.numeric(gm$statistic),
       moran_p   = gm$p.value,
       nsim      = nsim_moran,
       n_hh      = sum(hh),
       n_ll      = sum(ll),
       hh        = hh,
       ll        = ll,
       cluster   = cl,
       tier      = tier,
       n_core    = sum(core),
       n_edge    = sum(edge),
       n_outside = sum(flag) - sum(core) - sum(edge),
       clustered = gm$p.value <= 0.05)
}

#' Summarise an SPPT result
#'
#' Prints the robust S index, the sameness and arrangement tests, the
#' confidence tier counts and the settings of a [sppt()] result.
#'
#' @param object,x An object of class `sppt_result`, as returned by [sppt()].
#' @param ... Further arguments passed to the `sf` method.
#' @return `summary()` invisibly returns the summary list stored in the
#'   `"sppt_summary"` attribute. `print()` invisibly returns `x`.
#' @examples
#' res <- sppt(vancouver_ct, group_col = "id",
#'             count_col = c("BNER2003", "BNER2025"),
#'             B = 199, seed = 1, verbose = FALSE, create_maps = FALSE)
#' summary(res)
#' @export
summary.sppt_result <- function(object, ...) {
  s <- attr(object, "sppt_summary")
  if (is.null(s)) return(NextMethod())

  cat("\n")
  cat("==================================================================\n")
  cat("SPPT Summary\n")
  cat("==================================================================\n")
  cat(sprintf("  Base             : %s\n", s$base_name))
  cat(sprintf("  Test             : %s\n", s$test_name))
  cat(sprintf("  Spatial units    : %d  (%d eligible, %d zero-zero excluded)\n",
              s$n_units, s$n_eligible, s$n_zero))
  cat(sprintf("  lambda_min       : %.3f\n", s$lambda_min))
  cat(sprintf("  Bootstrap B      : %d\n", s$B))
  if (!is.null(s$alpha) && !is.na(s$s_index_pv))
    cat(sprintf("  Robust S-index   : %.4f  (alpha = %.3f; %d dissimilar)\n",
                s$s_index_pv, s$alpha, s$n_sig_pv))
  if (!is.null(s$sameness) && !is.na(s$sameness$p_value)) {
    cat(sprintf("  Sameness test    : p %s  (expected S = %.4f, SE = %.4f, %s)\n",
                fmt_p(s$sameness$p_value), s$sameness$s_expected,
                s$sameness$se, s$sameness$method))
    cat(sprintf("  Sameness cut-off : S >= %.4f at tau = %.2f  [%s]\n",
                s$sameness$threshold, s$sameness$tau,
                if (isTRUE(s$sameness$same)) "same" else "different"))
  }
  if (!is.null(s$conf_level) && !is.na(s$s_index_ci))
    cat(sprintf("  Robust CI-index  : %.4f  (CI = %.1f%%; %d dissimilar)\n",
                s$s_index_ci, s$conf_level * 100, s$n_sig_ci))
  if (!is.null(s$arrangement)) {
    cat(sprintf("  Arrangement      : Moran's I = %+.4f on signed evidence, p %s\n",
                s$arrangement$moran_I, fmt_p(s$arrangement$moran_p)))
    cat(sprintf("  Local clusters   : %d (%s > %s), %d (%s > %s)\n",
                s$arrangement$n_hh, s$test_name, s$base_name,
                s$arrangement$n_ll, s$base_name, s$test_name))
    if (!is.null(s$arrangement$n_core))
      cat(sprintf("  Flagged by tier  : %d core, %d edge, %d outside%s\n",
                  s$arrangement$n_core, s$arrangement$n_edge,
                  s$arrangement$n_outside,
                  if (isTRUE(s$widespread_shift)) " (widespread shift)"
                  else ""))
  }
  cat(sprintf("  Decision tier    : %s\n", s$decision_tier))
  cat("==================================================================\n\n")
  invisible(s)
}

#' @rdname summary.sppt_result
#' @export
print.sppt_result <- function(x, ...) {
  summary(x)
  NextMethod()
  invisible(x)
}
