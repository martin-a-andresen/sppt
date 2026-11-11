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
#' **Testing against a random distribution.** With `count_col = "BNER2021"`
#' or `c("random", "BNER2021")`, the column is compared with a random
#' allocation of its own total to the units. By default
#' (`random_by = "exposure"`) each event is placed in a unit with probability
#' in proportion to the unit's exposure: its area for polygons, its length
#' for lines such as street segments. `random_by = "units"` gives every unit
#' equal probability, and the name of a numeric column (population, say)
#' uses that column as the exposure. The random distribution is
#' known exactly, so the test is a Monte Carlo test against it: `B` random
#' allocations give each unit's envelope of shares, and the p-value is the
#' share of allocations that depart from the expected share by at least as
#' much as the observed pattern does, \eqn{p = (1 + r)/(B + 1)}. Every unit
#' with any exposure is eligible, including units with no events, which are
#' informative here. Everything else is reported as usual: the S index and
#' its sameness test (now: no more different than a random allocation),
#' \eqn{S_D} with `similarity = TRUE` (only the observed pattern is moved),
#' the arrangement test, tiers and maps.
#' The expected counts of the random pattern are added as column `random`.
#' Confidence intervals and `fix_base` do not apply.
#'
#' **Sameness and degree of similarity.** The S index tests sameness: two
#' realisations of one process give \eqn{S \approx 1 - \alpha}, and the
#' sameness test asks whether S is lower than that. It does not measure how
#' similar two patterns are. With `similarity = TRUE`, `sppt()` also reports
#' the degree of similarity \eqn{S_D}: one minus the smallest share of
#' events that would have to be redistributed within one pattern for a
#' retest of the two patterns to give \eqn{S \ge 1 - \alpha} with
#' probability `similarity_level`. It is computed moving each pattern in
#' turn, and \eqn{S_D} uses the smaller share; the report gives the share
#' and number of events for both. It takes extra bootstraps (seconds for a
#' hundred units, minutes for a thousand), so it is off by default, and it
#' can also be computed afterwards with [sppt_similarity()], which describes
#' the method. It uses the uncorrected flags, \eqn{p \le \alpha}.
#'
#' **Multiple testing.** No correction is applied by default. Simulation
#' evidence shows that corrections remove much of the power at fine spatial
#' scales, and that at coarse scales a stricter `alpha` achieves much the same
#' as a correction; the sameness and arrangement tests and [sppt_clusters()]
#' address the global question instead. Where a list of units with
#' controlled error is wanted, `correction` applies a method of
#' [stats::p.adjust()]. `"BY"` (Benjamini and Yekutieli) controls the false
#' discovery rate under any dependence between units, which the fixed totals
#' of the two patterns create, so it is the safe choice; `"BH"` assumes
#' independent or positively dependent tests. With a correction, `pvalue`
#' holds the raw p-values, `pvalue_adj` the adjusted ones, units are flagged
#' on the adjusted values, the S index counts those flags, and the sameness
#' test is not reported, because the number of units flagged by chance is no
#' longer binomial with rate `alpha`. The arrangement test and the tiers use
#' the raw p-values for the evidence surface and the corrected flags for the
#' tiers.
#'
#' **Number of replicates.** With `B = "auto"` (default) the number of
#' replicates is the largest of three requirements, rounded up to a number
#' ending in 999: a minimum of 4999; a Monte Carlo standard error of a
#' p-value at `alpha` no more than 10% of `alpha`, which needs
#' \eqn{B \ge 100(1 - \alpha)/\alpha}; and, under a correction, a smallest
#' attainable p-value, \eqn{1/(B + 1)}, small enough for the most significant
#' unit to survive it, \eqn{B + 1 \ge n/\alpha} for most methods and
#' \eqn{B + 1 \ge n H_n/\alpha} for BY, where \eqn{n} is the number of
#' eligible units and \eqn{H_n = 1 + 1/2 + \cdots + 1/n}. So `B = "auto"`
#' gives 4999 at `alpha = 0.05` without a correction, 9999 at `alpha = 0.01`,
#' and 129,999 for BY at `alpha = 0.05` across 824 units. Large numbers of
#' replicates are computed in blocks, so memory does not grow with `B`; the
#' result is the same as the single-matrix calculation. The confidence
#' interval option needs the full matrices, so with `calculate_ci = TRUE`
#' memory grows with \eqn{n \times B}.
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
#'   One element may be `"random"`, or a single column may be given alone,
#'   to test that column against a random distribution of its events (see
#'   Details); a single column is tested as `c("random", column)`, so that
#'   positive differences mark units with more events than random.
#' @param B `"auto"` (default) or an integer of at least 99: the number of
#'   bootstrap replicates. `"auto"` chooses it from `alpha`, the number of
#'   eligible units and `correction` (see Details). An integer is used as
#'   given; an odd value avoids ties.
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
#' @param correction Character. Multiple testing correction, a method of
#'   [stats::p.adjust()]: `"none"` (default), `"BY"`, `"BH"`, `"bonferroni"`,
#'   `"holm"`, `"hochberg"` or `"hommel"`. See Details.
#' @param calculate_ci Logical. Also compute marginal bootstrap confidence
#'   intervals and the CI-based S index (Andresen 2009). Default `FALSE`.
#' @param conf_level Numeric in (0, 1). Confidence level for the interval
#'   overlap criterion. Default `0.85`: two marginal 85% intervals fail to
#'   overlap when the difference exceeds about 2.04 standard errors, which
#'   corresponds to a two-sample test at about `alpha = 0.05` (Andresen
#'   2026). `0.95` reproduces the original, more conservative analysis. Only
#'   used when `calculate_ci = TRUE`. The correspondence assumes intervals of
#'   similar width. Where one pattern has no events in a unit its interval
#'   is \eqn{[0, 0]}, and non-overlap then needs only about 1.44 standard
#'   errors of the other pattern instead of 1.96, so at fine scales with
#'   many such units overlap flags far more units than the p-values; a note
#'   is printed when they are at least 20% of eligible units.
#' @param fix_base Logical. If `TRUE`, the base pattern is not bootstrapped
#'   and its observed proportions are treated as fixed. This gives an
#'   asymmetric test, intended only for comparison against a known reference
#'   distribution. Default `FALSE`.
#' @param random_by How the random distribution spreads events when one
#'   pattern is `"random"`. `"exposure"` (default): in proportion to each
#'   unit's area for polygons or length for lines (street segments, say).
#'   `"units"`: equal probability for every unit. Or the name of a
#'   non-negative numeric column in `data` to use as the exposure.
#' @param use_percentages Logical. Express proportions as percentages
#'   (0-100) rather than as proportions (0-1). Default `TRUE`.
#' @param similarity Logical. Compute the degree of similarity \eqn{S_D}
#'   (see Details and [sppt_similarity()]). Default `FALSE`, because it takes
#'   extra bootstraps.
#' @param similarity_level Numeric in (0, 1). Probability that a retest of
#'   the redistributed patterns gives \eqn{S \ge 1 - \alpha}. Default 0.95.
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
#'   \item{`pvalue_adj`}{Adjusted p-value (only with a `correction`).}
#'   \item{`sig_pvalue`}{1 if the unit is flagged at `alpha` (on the adjusted
#'     p-value under a correction), otherwise 0.}
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
#' The attribute `"sppt_summary"` holds the robust S index, \eqn{S_D} with
#' `similarity = TRUE` (element `similarity`), the sameness and
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
#' @seealso [sppt_longitudinal()], [sppt_sameness()], [sppt_similarity()],
#'   [sppt_clusters()],
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
                 B                  = "auto",
                 seed               = NULL,
                 alpha              = 0.05,
                 calculate_pvalues  = TRUE,
                 sameness_tau       = 0.05,
                 arrangement        = TRUE,
                 nsim_moran         = 999L,
                 correction         = "none",
                 calculate_ci       = FALSE,
                 conf_level         = 0.85,
                 fix_base           = FALSE,
                 random_by          = "exposure",
                 use_percentages    = TRUE,
                 similarity         = FALSE,
                 similarity_level   = 0.95,
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

  # One column, or one column and "random": test the column against a random
  # allocation of its events to the units.
  if (!is.character(count_col) || !length(count_col) %in% 1:2)
    stop("`count_col` must be a character vector of length 2 (base, test), ",
         "or a single column to test against a random distribution.",
         call. = FALSE)
  if (length(count_col) == 1L) count_col <- c("random", count_col)
  is_rand <- tolower(count_col) == "random" & !count_col %in% names(data)
  if (all(is_rand))
    stop("At most one element of `count_col` can be \"random\".", call. = FALSE)
  random_mode <- any(is_rand)
  if (random_mode) {
    count_col[is_rand] <- "random"
    if (fix_base)
      stop("`fix_base` does not apply when one pattern is \"random\": the ",
           "random distribution is already known.", call. = FALSE)
    if (calculate_ci) {
      warning("`calculate_ci` is not available against a random ",
              "distribution; p-values are used.", call. = FALSE)
      calculate_ci <- FALSE
    }
  }

  for (col in count_col[!is_rand]) {
    if (!col %in% names(data))
      stop("`count_col` column '", col, "' not found in data.", call. = FALSE)
    vals <- data[[col]]
    if (!is.numeric(vals))
      stop("`count_col` column '", col, "' must be numeric.", call. = FALSE)
    if (any(!is.na(vals) & vals < 0))
      stop("`count_col` column '", col, "' contains negative values.",
           call. = FALSE)
  }

  B_auto <- is.character(B) && length(B) == 1L && identical(tolower(B), "auto")
  if (!B_auto && (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 99))
    stop("`B` must be \"auto\" or a single integer >= 99.", call. = FALSE)
  if (!B_auto) B <- as.integer(B)

  valid_corrections <- c("none", "BY", "BH", "bonferroni", "holm",
                         "hochberg", "hommel")
  if (!is.character(correction) || length(correction) != 1L ||
      !correction %in% valid_corrections)
    stop("`correction` must be one of: ",
         paste(valid_corrections, collapse = ", "), ".", call. = FALSE)

  if (!is.logical(similarity) || length(similarity) != 1L || is.na(similarity))
    stop("`similarity` must be TRUE or FALSE.", call. = FALSE)
  want_similarity <- similarity
  if (!is.numeric(similarity_level) || length(similarity_level) != 1L ||
      similarity_level <= 0 || similarity_level >= 1)
    stop("`similarity_level` must be a number in (0, 1).", call. = FALSE)
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

  read_counts <- function(nm) {
    x <- as.integer(data[[nm]])
    x[is.na(x)] <- 0L
    x
  }
  if (random_mode) {
    # The random pattern holds the observed total, spread over the units in
    # proportion to their exposure (area, length, a column, or equal). Its
    # "counts" are the expected counts, which need not be whole numbers.
    obs_name   <- count_col[!is_rand]
    counts_obs <- read_counts(obs_name)
    expo       <- exposure_weights(data, random_by)
    w          <- expo$w
    p_ref        <- w / sum(w)
    counts_rand  <- sum(counts_obs) * p_ref
    counts_a <- if (is_rand[1L]) counts_rand else counts_obs
    counts_b <- if (is_rand[2L]) counts_rand else counts_obs
  } else {
    counts_a <- read_counts(base_name)
    counts_b <- read_counts(test_name)
  }

  N_a <- if (random_mode && is_rand[1L]) sum(counts_obs) else sum(counts_a)
  N_b <- if (random_mode && is_rand[2L]) sum(counts_obs) else sum(counts_b)

  # Units empty in both patterns are excluded from inference, and from the
  # number of tests that sets B under a correction.
  zero_mask  <- counts_a == 0L & counts_b == 0L
  n_zero     <- sum(zero_mask)
  n_eligible <- n_units - n_zero

  B_info <- choose_B(alpha, n_eligible, correction)
  if (B_auto) {
    B <- B_info$B
  } else if (calculate_pvalues && correction != "none" && B < B_info$B_correction) {
    warning(sprintf(paste0(
      "With correction = \"%s\" and %d eligible units, no unit can be ",
      "flagged unless B >= %d: the smallest attainable p-value, 1/(B + 1), ",
      "is too large to survive the correction. Use B = \"auto\"."),
      correction, n_eligible, B_info$B_correction), call. = FALSE)
  }
  B_note <- if (!B_auto) "set by user" else if (B_info$binding == "correction")
    sprintf("auto: lets %s flag a unit among %d", correction, n_eligible)
  else if (B_info$binding == "alpha")
    sprintf("auto: Monte Carlo SE of p at alpha = %.2g", sqrt(alpha * (1 - alpha) / B))
  else "auto: default minimum"

  # -- 3. lambda_min and decision tier ----------------------------------------
  lambda_min <- min(N_a, N_b) / n_units

  decision_tier <- if (lambda_min >= 17) {
    "Coarse scale (lambda_min >= 17): CI and PV methods both reliable. PV alpha = 0.05; FDR-priority use alpha = 0.01."
  } else if (lambda_min >= 6) {
    "Medium-high scale (6 <= lambda_min < 17): PV alpha = 0.05 recommended; FDR-priority use alpha = 0.01 once lambda_min is above about 10."
  } else if (lambda_min >= 1.5) {
    "Medium scale (1.5 <= lambda_min < 6): PV alpha = 0.05 recommended."
  } else if (lambda_min >= 0.5) {
    "Fine scale (0.5 <= lambda_min < 1.5): power is limited; interpret core and edge units rather than isolated flags."
  } else {
    "Very fine scale (lambda_min < 0.5): all methods have severely limited power. Consider aggregating to a coarser spatial unit."
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
  # Lines are wrapped to fit 70 columns; see field_lines().
  say("\n", rule_line("="), "Spatial Point Pattern Test\n", rule_line("="),
      sep = "")
  rand_label <- function(N) sprintf(
    "random (N = %d, %s)", as.integer(round(N)), expo$label)
  say(field_lines("Base", if (random_mode && is_rand[1L]) rand_label(N_a)
                  else sprintf("%s (N = %d)", base_name, N_a)))
  say(field_lines("Test", if (random_mode && is_rand[2L]) rand_label(N_b)
                  else sprintf("%s (N = %d)", test_name, N_b)))
  say(field_lines("Units", n_units))
  say(field_lines("lambda_min", sprintf("%.3f", lambda_min)))
  say(field_lines("Bootstrap B", sprintf("%d (%s)", B, B_note)))
  if (calculate_pvalues)
    say(field_lines("alpha", sprintf("%.3f", alpha)))
  if (calculate_pvalues && correction != "none")
    say(field_lines("Correction", correction))
  if (calculate_ci)
    say(field_lines("CI level", sprintf("%.1f%% (marginal overlap)",
                                        conf_level * 100)))
  say(field_lines("Decision tier", decision_tier))
  say(rule_line("-"))

  # -- 5. Observed proportions and difference ---------------------------------
  scale_factor <- if (use_percentages) 100 else 1

  prop_a <- if (N_a > 0) counts_a / N_a * scale_factor else rep(0, n_units)
  prop_b <- if (N_b > 0) counts_b / N_b * scale_factor else rep(0, n_units)
  obs_diff <- prop_b - prop_a   # test minus base

  # -- 6. Shared bootstrap ----------------------------------------------------
  # One multinomial draw per variable, B replicates, shared by CI and PV.
  # Drawing counts from the multinomial of observed shares is the closed form
  # of resampling events at their observed locations and re-aggregating.
  # When the n x B matrices would be large and no intervals are needed, the
  # p-values are computed in blocks of replicates instead (blocked_counts()),
  # which gives the same result without holding the matrices in memory.
  valid_bootstrap <- N_a > 0L && N_b > 0L
  if (!valid_bootstrap)
    warning("One or both count columns sum to zero. No inference possible.",
            call. = FALSE)
  active   <- !zero_mask
  max_cells <- getOption("sppt.max_cells", 2.5e7)
  blocked  <- valid_bootstrap && !calculate_ci &&
    as.numeric(n_units) * B > max_cells
  r_count  <- NULL

  if (valid_bootstrap && random_mode) {
    # Monte Carlo test against the random distribution. That distribution is
    # known, so only the random allocations vary: B allocations of the
    # observed total, and each unit's observed share compared with the
    # envelope of its shares across them. Resampling the observed pattern as
    # well would count its sampling variation twice.
    if (calculate_pvalues)
      r_count <- random_counts(sum(counts_obs), p_ref, B, scale_factor,
                               counts_obs, active)
  } else if (valid_bootstrap) {
    p_a <- counts_a / N_a
    p_b <- counts_b / N_b
    if (blocked) {
      if (calculate_pvalues) {
        say(field_lines("Note", sprintf(
          "computing p-values in blocks (%d units x %d replicates)",
          n_units, B)))
        r_count <- blocked_counts(N_a, N_b, p_a, p_b, B, fix_base,
                                  scale_factor, obs_diff, active)
      }
    } else {
      boot_a <- if (fix_base) {
        matrix(rep(p_a * scale_factor, B), nrow = n_units, ncol = B)
      } else {
        matrix(stats::rmultinom(B, size = N_a, prob = p_a),
               nrow = n_units) / N_a * scale_factor
      }
      boot_b <- matrix(stats::rmultinom(B, size = N_b, prob = p_b),
                       nrow = n_units) / N_b * scale_factor

      # Centred difference distribution for p-values
      boot_diffs_c <- boot_b - boot_a
      boot_diffs_c <- boot_diffs_c - rowMeans(boot_diffs_c)
      if (calculate_pvalues)
        r_count <- rowSums(abs(boot_diffs_c[active, , drop = FALSE]) >=
                             abs(obs_diff[active]))
      rm(boot_diffs_c)
      if (!calculate_ci) rm(boot_a, boot_b)
    }
  }

  # -- 7. Bootstrap p-values --------------------------------------------------
  pvalue     <- rep(NA_real_, n_units)
  pvalue_adj <- rep(NA_real_, n_units)
  sig_pvalue <- rep(NA_integer_, n_units)
  sig_pv_dir <- rep(NA_integer_, n_units)

  if (calculate_pvalues && valid_bootstrap) {
    # Bounded Monte Carlo p-value: p = (1 + r) / (B + 1), where r counts the
    # centred bootstrap differences at least as large in absolute value as
    # the observed difference. The unbounded form r / B admits p = 0, which
    # is anti-conservative and unbounded on the -log10 scale.
    pvalue[active] <- (1 + r_count) / (B + 1)

    # Under a correction units are flagged on the adjusted p-values; the raw
    # p-values are kept for the evidence surface and for reporting.
    p_flag <- pvalue
    if (correction != "none") {
      pvalue_adj[active] <- stats::p.adjust(pvalue[active],
                                            method = correction)
      p_flag <- pvalue_adj
    }
    sig_pvalue[active] <- as.integer(p_flag[active] <= alpha)
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

  s_index_pv <- NA_real_
  n_sig_pv   <- NA_integer_
  sameness   <- NULL
  if (calculate_pvalues && valid_bootstrap && n_eligible > 0L) {
    n_sig_pv   <- sum(sig_pvalue[!zero_mask] == 1L, na.rm = TRUE)
    s_index_pv <- 1 - n_sig_pv / n_eligible
    # Only calibrated for uncorrected p-values under the symmetric bootstrap:
    # a correction moves the null flagging rate away from alpha (under
    # sameness it flags no unit at all in most samples), and fix_base = TRUE
    # breaks it at sparse units.
    if (correction == "none" && !fix_base)
      sameness <- sppt_sameness(s_index_pv, n_sig_pv, n_eligible, alpha,
                                tau = sameness_tau)
  }

  # No sameness test for the CI-based index: the null flagging rate of the
  # interval-overlap rule is not 1 - conf_level.
  s_index_ci <- NA_real_
  n_sig_ci   <- NA_integer_
  ci_one_empty <- NA_real_
  if (calculate_ci && valid_bootstrap && n_eligible > 0L) {
    n_sig_ci   <- sum(intervals_overlap[!zero_mask] == 0L, na.rm = TRUE)
    s_index_ci <- 1 - n_sig_ci / n_eligible
    # Share of eligible units where one pattern has no events: there its
    # bootstrap interval collapses to [0, 0] and overlap is too liberal
    ci_one_empty <- mean(((counts_a == 0) != (counts_b == 0))[!zero_mask])
  }

  # Degree of similarity, S_D, on the uncorrected flags (see Details). The
  # search restores the random-number state, so nothing else changes.
  if (want_similarity) {
    if (calculate_pvalues && valid_bootstrap && n_eligible > 0L) {
      similarity <- similarity_search(
        counts_a, counts_b, base_name = base_name, test_name = test_name,
        alpha = alpha,
        B = if (correction == "none") B else
          choose_B(alpha, n_eligible, "none")$B,
        level = similarity_level, fix_base = fix_base,
        random_side = if (random_mode) (if (is_rand[1L]) "base" else "test"))
    } else {
      similarity <- NULL
      warning("S_D needs p-values from a valid bootstrap; not computed.",
              call. = FALSE)
    }
  } else {
    similarity <- NULL
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
        warning("Arrangement test failed, so no tiers are assigned: ",
                conditionMessage(e), call. = FALSE)
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
  if (random_mode) result[["random"]] <- round(counts_rand, 6)
  result$base_prop         <- round(prop_a, 6)
  result$test_prop         <- round(prop_b, 6)
  result$difference        <- round(obs_diff, 6)
  result$lambda_min        <- lambda_min
  result$low_power_warning <- lambda_min < 1.5

  if (calculate_pvalues) {
    result$pvalue          <- round(pvalue, 6)
    if (correction != "none") result$pvalue_adj <- round(pvalue_adj, 6)
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
    say(field_lines("Robust S", sprintf("%.4f (%d of %d units differ)",
                                        s_index_pv, n_sig_pv, n_eligible)))
    say(field_lines("", sprintf("%d units %s > %s",
                                sum(sig_pv_dir[!zero_mask] == 1L, na.rm = TRUE),
                                test_name, base_name)))
    say(field_lines("", sprintf("%d units %s > %s",
                                sum(sig_pv_dir[!zero_mask] == -1L, na.rm = TRUE),
                                base_name, test_name)))
  }
  if (!is.null(similarity))
    say(similarity_lines(similarity))

  if (!is.null(sameness) && !is.na(sameness$p_value)) {
    say(field_lines("Sameness", sprintf(
      "p %s (%s); S >= %.4f needed at tau = %.2f; %s",
      fmt_p(sameness$p_value), sameness$method, sameness$threshold,
      sameness$tau,
      if (random_mode) {
        if (isTRUE(sameness$same)) "no more different than a random allocation"
        else "more different than a random allocation"
      } else if (isTRUE(sameness$same))
        "no more different than two realisations of one process"
      else "more different than two realisations of one process")))
    if (n_eligible < n_units * 0.75)
      say(field_lines("", paste("Many units are empty in both patterns;",
                                "the sameness test is conservative here.")))
  } else if (calculate_pvalues && !is.na(s_index_pv) && correction != "none") {
    say(field_lines("Sameness", paste(
      "not reported: under a correction the number of units flagged by",
      "chance is no longer binomial with rate alpha")))
  } else if (calculate_pvalues && !is.na(s_index_pv) && fix_base) {
    say(field_lines("Sameness", "not reported (not calibrated for fix_base = TRUE)"))
  }

  if (calculate_ci && !is.na(s_index_ci)) {
    say(field_lines("Robust S (CI)", sprintf(
      "%.4f (%d of %d units differ at CI = %.1f%%)", s_index_ci, n_sig_ci,
      n_eligible, conf_level * 100)))
    if (isTRUE(ci_one_empty >= 0.2))
      say(field_lines("Note", sprintf(paste(
        "in %.0f%% of units one pattern has no events, so its interval is",
        "[0, 0] there and overlap flags differences the p-values do not;",
        "prefer the p-values at this scale"), 100 * ci_one_empty)))
  }

  if (!is.null(arrange_res)) {
    say(field_lines("Arrangement", sprintf(
      "Moran's I = %+.3f on signed evidence, p %s (%d permutations)",
      arrange_res$moran_I, fmt_p(arrange_res$moran_p), arrange_res$nsim)))
    say(field_lines("Clusters", sprintf(
      "%d units %s > %s; %d units %s > %s", arrange_res$n_hh, test_name,
      base_name, arrange_res$n_ll, base_name, test_name)))
    say(field_lines("Tiers", sprintf(
      "%d core, %d edge, %d outside (core = inside a cluster; edge = borders a core unit with the same sign)",
      arrange_res$n_core, arrange_res$n_edge, arrange_res$n_outside)))
    if (wide_shift && arrange_res$n_edge > 0L)
      say(field_lines("Note", sprintf(paste(
        "%.0f%% of flagged units share one direction. Edge units are still",
        "significant, but a widespread shift places flagged units beside",
        "clusters by chance; weigh core units more heavily."),
        100 * arrange_res$dominant_share)))
  }

  # Two-way verdict: the count question and the arrangement question are
  # separate, and the pair says more than either alone.
  if (!is.null(sameness) && !is.na(sameness$p_value) && !is.null(arrange_res)) {
    say("\n", verdict_table(!isTRUE(sameness$same),
                             isTRUE(arrange_res$clustered), wide_shift),
        sep = "")
  }

  if (n_zero > 0L)
    say(field_lines("Note", sprintf(
      "%d units empty in both patterns excluded from the S index and inference.",
      n_zero)))
  say(rule_line("="), "\n", sep = "")

  # -- 13. Summary attribute and class ----------------------------------------
  attr(result, "sppt_summary") <- list(
    base_name        = base_name,
    random           = if (random_mode) list(side = if (is_rand[1L]) "base"
                                             else "test", by = expo$by),
    test_name        = test_name,
    n_units          = n_units,
    n_eligible       = n_eligible,
    n_zero           = n_zero,
    lambda_min       = lambda_min,
    B                = B,
    alpha            = if (calculate_pvalues) alpha else NULL,
    correction       = correction,
    B_note           = B_note,
    s_index_pv       = s_index_pv,
    similarity       = similarity,
    n_sig_pv         = n_sig_pv,
    sameness         = sameness,
    arrangement      = if (is.null(arrange_res)) NULL else
      arrange_res[c("moran_I", "moran_p", "nsim", "n_hh", "n_ll",
                    "n_core", "n_edge", "n_outside", "clustered")],
    widespread_shift = wide_shift,
    conf_level       = if (calculate_ci) conf_level else NULL,
    s_index_ci       = s_index_ci,
    n_sig_ci         = n_sig_ci,
    ci_one_empty     = ci_one_empty,
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

# Number of bootstrap replicates for B = "auto". Three requirements, the
# largest wins, rounded up to a number ending in 999:
#   * a minimum of 4999;
#   * Monte Carlo standard error of a p-value at alpha no more than 10% of
#     alpha, sqrt(alpha (1 - alpha) / B) <= 0.1 alpha, so
#     B >= 100 (1 - alpha) / alpha;
#   * under a correction, the smallest attainable p-value, 1 / (B + 1), small
#     enough for the most significant of n tests to survive it:
#     1 / (B + 1) <= alpha / n, or alpha / (n H_n) for BY, where
#     H_n = 1 + 1/2 + ... + 1/n.
choose_B <- function(alpha, n_tests, correction = "none") {
  round_up <- function(x) {
    k <- max(3, floor(log10(x)) - 1)
    as.integer(ceiling((x + 1) / 10^k) * 10^k - 1)
  }
  n_tests <- max(n_tests, 1L)
  need <- c(minimum = 4999,
            alpha = 100 * (1 - alpha) / alpha,
            correction = if (correction == "none") 0 else
              n_tests * (if (correction == "BY")
                sum(1 / seq_len(n_tests)) else 1) / alpha)
  list(B = round_up(max(need)),
       binding = names(need)[which.max(need)],
       B_correction = if (correction == "none") 0L else
         round_up(need[["correction"]]))
}

# Exceedance counts for the Monte Carlo test against a random distribution:
# B allocations of N events with probabilities p_ref, in blocks. The observed
# and simulated deviations are computed with the same arithmetic so that
# equal counts compare as equal.
# Exposure for the random distribution: unit area for polygons, length for
# lines, equal weights, or a user-supplied column. Returns the weights, a
# console label and the resolved choice.
exposure_weights <- function(data, random_by) {
  if (!is.character(random_by) || length(random_by) != 1L || is.na(random_by))
    stop("`random_by` must be \"exposure\", \"units\" or the name of a ",
         "column in `data`.", call. = FALSE)
  n <- nrow(data)
  if (random_by == "units")
    return(list(w = rep(1, n), label = "equal probability per unit",
                by = "units"))
  if (random_by == "exposure") {
    gt <- unique(as.character(sf::st_geometry_type(data, by_geometry = TRUE)))
    poly <- c("POLYGON", "MULTIPOLYGON")
    line <- c("LINESTRING", "MULTILINESTRING")
    if (all(gt %in% poly)) {
      w <- as.numeric(sf::st_area(data)); what <- "area"
    } else if (all(gt %in% line)) {
      w <- as.numeric(sf::st_length(data)); what <- "length"
    } else {
      stop("`random_by = \"exposure\"` needs all polygons (area) or all ",
           "lines (length); this layer has ", paste(gt, collapse = ", "),
           ". Use `random_by = \"units\"` or name an exposure column.",
           call. = FALSE)
    }
    label <- paste("probability in proportion to", what)
  } else {
    if (!random_by %in% names(data))
      stop("`random_by` must be \"exposure\", \"units\" or the name of a ",
           "column in `data`; no column \"", random_by, "\".", call. = FALSE)
    w <- sf::st_drop_geometry(data)[[random_by]]
    if (!is.numeric(w))
      stop("Exposure column `", random_by, "` must be numeric.", call. = FALSE)
    what  <- random_by
    label <- paste0("probability in proportion to `", random_by, "`")
  }
  if (any(!is.finite(w)) || any(w < 0) || sum(w) <= 0)
    stop("The exposure (", what, ") must be finite and non-negative, with ",
         "a positive total.", call. = FALSE)
  if (any(w == 0))
    warning(sum(w == 0), " unit(s) have zero exposure (", what, "); a random ",
            "allocation never places events there.", call. = FALSE)
  list(w = w, label = label, by = what)
}

random_counts <- function(N, p_ref, B, scale_factor, counts_obs, active) {
  n <- length(p_ref)
  block <- max(1L, as.integer(getOption("sppt.block_cells", 5e6) %/% n))
  sizes <- diff(unique(c(seq.int(0L, B, by = block), B)))
  act <- which(active)
  ref <- p_ref[act] * scale_factor
  obs <- abs(counts_obs[act] / N * scale_factor - ref)
  r <- numeric(length(act))
  for (m in sizes) {
    x <- stats::rmultinom(m, N, p_ref)[act, , drop = FALSE] / N * scale_factor
    r <- r + rowSums(abs(x - ref) >= obs)
  }
  r
}

# P-value exceedance counts computed in blocks of replicates. Two passes over
# the same random numbers: the first finds each unit's mean bootstrap
# difference (the centre), the second counts centred differences at least as
# large as the observed one. The base and test streams are kept apart by
# saving and restoring the generator state, so the replicates, and hence the
# counts, are those of the single-matrix calculation, up to floating-point
# rounding of the centre. The generator is left where that calculation would
# leave it.
blocked_counts <- function(N_a, N_b, p_a, p_b, B, fix_base, scale_factor,
                           obs_diff, active) {
  n <- length(p_a)
  block <- max(1L, as.integer(getOption("sppt.block_cells", 5e6) %/% n))
  sizes <- diff(unique(c(seq.int(0L, B, by = block), B)))
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    stats::runif(1L)
  get_seed <- function() get(".Random.seed", envir = globalenv(),
                             inherits = FALSE)
  put_seed <- function(x) assign(".Random.seed", x, envir = globalenv())

  # Pass 1: row sums of the counts, base stream then test stream
  s0 <- get_seed()
  sum_a <- numeric(n)
  if (!fix_base)
    for (m in sizes) sum_a <- sum_a + rowSums(stats::rmultinom(m, N_a, p_a))
  s1 <- get_seed()
  sum_b <- numeric(n)
  for (m in sizes) sum_b <- sum_b + rowSums(stats::rmultinom(m, N_b, p_b))
  s2 <- get_seed()
  mean_a <- if (fix_base) p_a * scale_factor else sum_a / B / N_a * scale_factor
  centre <- sum_b / B / N_b * scale_factor - mean_a

  # Pass 2: the same draws again, paired block by block
  act <- which(active)
  obs <- abs(obs_diff[act])
  r <- numeric(length(act))
  sa <- s0
  sb <- s1
  for (m in sizes) {
    if (fix_base) {
      a <- p_a[act] * scale_factor
    } else {
      put_seed(sa)
      a <- stats::rmultinom(m, N_a, p_a)[act, , drop = FALSE] / N_a *
        scale_factor
      sa <- get_seed()
    }
    put_seed(sb)
    b <- stats::rmultinom(m, N_b, p_b)[act, , drop = FALSE] / N_b *
      scale_factor
    sb <- get_seed()
    r <- r + rowSums(abs(b - a - centre[act]) >= obs)
  }
  put_seed(s2)
  r
}

fmt_p <- function(p) if (p < 0.0001) "< 0.0001" else sprintf("= %.4f", p)

# The four cells of the two-way verdict: count above chance (sameness
# rejected) by arrangement clustered (global Moran's I significant).
verdict_label <- function(count_hi, arr_hi, wide_shift = FALSE) {
  lab <- if (!count_hi && !arr_hi) "same process"
  else if (!count_hi && arr_hi) "few differences, but concentrated"
  else if (count_hi && !arr_hi) "many differences, diffuse"
  else "differences real and localised"
  if (count_hi && arr_hi && wide_shift) lab <- paste(lab, "plus a widespread shift")
  lab
}

verdict_table <- function(count_hi, arr_hi, wide_shift = FALSE) {
  mark <- function(a, b) if (count_hi == a && arr_hi == b) " <<<" else ""
  paste0(
    "  Two-way verdict    | arrangement random | arrangement clustered\n",
    sprintf("  count at chance    | same process%s%s| few, concentrated%s\n",
            mark(FALSE, FALSE), strrep(" ", 7 - nchar(mark(FALSE, FALSE))),
            mark(FALSE, TRUE)),
    sprintf("  count above chance | many, diffuse%s%s| real, localised%s\n",
            mark(TRUE, FALSE), strrep(" ", 6 - nchar(mark(TRUE, FALSE))),
            mark(TRUE, TRUE)),
    field_lines("This result", verdict_label(count_hi, arr_hi, wide_shift)))
}

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

  cat("\n", rule_line("="), "SPPT Summary\n", rule_line("="), sep = "")
  cat(field_lines("Base", s$base_name))
  cat(field_lines("Test", s$test_name))
  cat(field_lines("Units", sprintf("%d (%d eligible, %d empty in both)",
                                   s$n_units, s$n_eligible, s$n_zero)))
  cat(field_lines("lambda_min", sprintf("%.3f", s$lambda_min)))
  cat(field_lines("Bootstrap B", if (is.null(s$B_note)) s$B else
    sprintf("%d (%s)", s$B, s$B_note)))
  if (!is.null(s$correction) && s$correction != "none")
    cat(field_lines("Correction", s$correction))
  if (!is.null(s$alpha) && !is.na(s$s_index_pv))
    cat(field_lines("Robust S", sprintf("%.4f (alpha = %.3f; %d differ)",
                                        s$s_index_pv, s$alpha, s$n_sig_pv)))
  if (!is.null(s$similarity))
    cat(similarity_lines(s$similarity))
  if (!is.null(s$sameness) && !is.na(s$sameness$p_value))
    cat(field_lines("Sameness", sprintf(
      "p %s (%s); S >= %.4f needed at tau = %.2f [%s]",
      fmt_p(s$sameness$p_value), s$sameness$method, s$sameness$threshold,
      s$sameness$tau, if (isTRUE(s$sameness$same)) "same" else "different")))
  if (!is.null(s$conf_level) && !is.na(s$s_index_ci)) {
    cat(field_lines("Robust S (CI)", sprintf("%.4f (CI = %.1f%%; %d differ)",
                                             s$s_index_ci, s$conf_level * 100,
                                             s$n_sig_ci)))
    if (isTRUE(s$ci_one_empty >= 0.2))
      cat(field_lines("Note", sprintf(paste(
        "one pattern has no events in %.0f%% of units; CI overlap is",
        "liberal there, prefer the p-values"), 100 * s$ci_one_empty)))
  }
  if (!is.null(s$arrangement)) {
    cat(field_lines("Arrangement", sprintf(
      "Moran's I = %+.3f on signed evidence, p %s",
      s$arrangement$moran_I, fmt_p(s$arrangement$moran_p))))
    cat(field_lines("Clusters", sprintf(
      "%d units %s > %s; %d units %s > %s", s$arrangement$n_hh, s$test_name,
      s$base_name, s$arrangement$n_ll, s$base_name, s$test_name)))
    if (!is.null(s$arrangement$n_core))
      cat(field_lines("Tiers", sprintf(
        "%d core, %d edge, %d outside%s", s$arrangement$n_core,
        s$arrangement$n_edge, s$arrangement$n_outside,
        if (isTRUE(s$widespread_shift)) " (widespread shift)" else "")))
  }
  if (!is.null(s$sameness) && !is.na(s$sameness$p_value) &&
      !is.null(s$arrangement))
    cat(field_lines("Verdict", verdict_label(
      !isTRUE(s$sameness$same), isTRUE(s$arrangement$clustered),
      isTRUE(s$widespread_shift))))
  cat(field_lines("Decision tier", s$decision_tier))
  cat(rule_line("="), "\n", sep = "")
  invisible(s)
}

#' @rdname summary.sppt_result
#' @export
print.sppt_result <- function(x, ...) {
  summary(x)
  NextMethod()
  invisible(x)
}
