#' Longitudinal spatial point pattern test
#'
#' @description
#' Extends the spatial point pattern test to a time series of three or more
#' periods (Andresen 2026). Every period is compared with the first period
#' (the base), and each unit's trajectory is classified.
#'
#' Two parallel sets of results are produced:
#'
#' * **Envelope results.** Bootstrap confidence intervals are computed for
#'   every unit and period. A unit is stable when its intervals overlap in all
#'   periods; the largest set of mutually overlapping periods gives a
#'   partial-stability measure. Trends are fitted only for units that are not
#'   stable.
#' * **P-value results** (`calculate_pvalues = TRUE`). Each period is tested
#'   against the base with the bounded Monte Carlo bootstrap p-value
#'   \eqn{p = (1 + r)/(B + 1)} used by [sppt()]. A unit is stable when no
#'   period is significant at `alpha`; trends are fitted only for units with
#'   at least one significant period.
#'
#' Trends are ordinary least squares fits of each unit's share (or count)
#' on time, linear or quadratic, tested at `trend_alpha`. Units are also
#' grouped into Jenks natural-breaks classes of their base-period count, and
#' trend outcomes are tabulated by class. Units that had events and then
#' recorded none in each of the last `min_zero_streak` periods are flagged
#' as sustained zeros.
#'
#' @inheritParams sppt
#' @param count_col Character vector of three or more count columns, one per
#'   period, in time order (for example `paste0("BNER", 2003:2025)`). The
#'   first is the base period. Column names should end in the period number,
#'   which is used to label the time axis.
#' @param prefix Character or `NULL`. Prefix for the output column names, so
#'   that results for several crime types can sit in one object. `NULL`
#'   (default) derives it from `count_col`, for example `"BNER_"`.
#' @param conf_level Numeric in (0, 1). Confidence level of the bootstrap
#'   envelopes. Default `0.85`: two marginal 85% intervals fail to overlap
#'   when the difference exceeds about 2.04 standard errors, which
#'   corresponds to a two-sample test at about `alpha = 0.05` (Andresen
#'   2026).
#' @param fix_base Logical. If `TRUE`, the base period is not bootstrapped and
#'   its observed shares are treated as fixed (the approach of earlier
#'   longitudinal applications). Not recommended for standard use. Default
#'   `FALSE`.
#' @param min_zero_streak Integer. Number of final periods with no events
#'   that, after earlier events, mark a unit as a sustained zero. Default `3`.
#' @param trend_type Character. `"linear"` (default) or `"quadratic"`. A
#'   quadratic trend needs four or more periods.
#' @param trend_alpha Numeric in (0, 1). Significance level for the trend
#'   regressions; independent of `alpha`. Default `0.05`.
#' @param jenks `"optimal"` (default) to choose the number of Jenks classes
#'   by the goodness of variance fit (GVF) stopping rule, or an integer of at
#'   least 2 to fix it.
#' @param jenks_gvf_threshold Numeric in (0, 1). With `jenks = "optimal"`,
#'   classes are added until the GVF improves by less than this. Default
#'   `0.05`.
#' @param create_maps Logical. Build the maps and print them. Default `TRUE`.
#' @param create_time_plot Logical. Plot total events per period with the
#'   summed envelope bounds. Default `TRUE`.
#' @param export_dir Character. Directory for every export option. Required
#'   when any `export_*` argument is `TRUE`; created if it does not exist.
#' @param export_maps Logical. Save the maps in `export_dir`. Default `FALSE`.
#' @param map_format Character. File type for exported maps: `"png"`
#'   (default), `"tif"`, `"eps"`, `"svg"` or `"pdf"`.
#' @param map_dpi Integer. Resolution of exported maps and the time plot.
#'   Default `300`; use 600 for journal TIFF submission.
#' @param map_width_in,map_height_in Numeric. Size of exported maps in
#'   inches. Defaults `7` and `5.5`.
#' @param include_title Logical. Keep map titles in exported maps. `FALSE`
#'   removes them, as journals usually require. Default `TRUE`.
#' @param figure_prefix Character or `NULL`. If given, exported maps are named
#'   `<figure_prefix>1`, `<figure_prefix>2`, and so on in map order (for
#'   example `"SmithFigure"`) instead of descriptive names.
#' @param export_time_plot Logical. Save the time plot as a PNG in
#'   `export_dir`. Default `FALSE`.
#' @param export_tables Logical. Save the summary indices and the Jenks trend
#'   tables as a CSV in `export_dir`. Default `FALSE`.
#' @param export_results Logical. Write the result to `export_dir`. Default
#'   `FALSE`.
#'
#' @return An `sf` object of class `sppt_long_result` containing the
#'   columns of `data` plus, for each
#'   period, `<count_col>_L` and `<count_col>_U` (envelope bounds) and, with
#'   p-values, `pvalue_<count_col>` (`NA` for the base period). Unit-level
#'   results are named with `prefix`:
#' \describe{
#'   \item{`intervals_overlap`}{1 if all envelopes overlap (stable).}
#'   \item{`max_overlap_proportion`}{Largest share of periods whose envelopes
#'     mutually overlap.}
#'   \item{`sustained_zero`}{1 if the unit declined to zero and stayed there.}
#'   \item{`Trend`}{Trend class. Linear: 2 / -2 significant increase /
#'     decrease, 1 / -1 non-significant, 0 stable. Quadratic: 1 / -1
#'     significant linear only, 2 / 3 significant increase with negative /
#'     positive curvature, -2 / -3 significant decrease with positive /
#'     negative curvature, 0 otherwise.}
#'   \item{`Trend_All`, `Trend_Linear_Sig`}{Linear direction, and significant
#'     linear direction (1, -1, 0).}
#'   \item{`Linear_R2`, `Linear_b0`, `Linear_bt`, `Linear_bt_sig`}{Linear fit;
#'     `Linear_bt_sig` is the slope for significant units only.}
#'   \item{`Quad_R2`, `Quad_b0`, `Quad_bt`, `Quad_bt2`, `Quad_Sig`}{Quadratic
#'     fit and the sign of a significant quadratic term (quadratic only).}
#'   \item{`*_p` versions}{The same quantities gated on p-value significance:
#'     `intervals_overlap_p`, `max_overlap_proportion_p` (share of
#'     non-significant periods), `prop_sig_p`, `sustained_zero_p`, `Trend_p`,
#'     `Trend_p_All`, `Trend_p_Linear_Sig`, `Linear_p_R2`, `Linear_p_b0`,
#'     `Linear_p_bt`, `Linear_p_bt_sig`, `Quad_p_R2`, `Quad_p_b0`,
#'     `Quad_p_bt`, `Quad_p_bt2`, `Quad_p_Sig`.}
#'   \item{`jenks_level`, `jenks_min`, `jenks_max`, `baseline_interval`}{Jenks
#'     class of the base-period count, its bounds and label.}
#' }
#' The S indices, Jenks breaks, trend tables and p-value matrix are stored as
#' attributes (for example `attr(x, "s_index_overlap")`,
#' `attr(x, "interval_table_linear")`, `attr(x, "pvalue_matrix")`). With
#' `create_maps = TRUE` the maps are stored as a named list of `ggplot`
#' objects in `attr(x, "sppt_maps")`; use [sppt_map()] or `plot()` to redraw
#' them with a different layout. `summary()` prints the S indices and trend
#' counts.
#'
#' **Relation to the published description.** Andresen (2026) gives the
#' p-value as \eqn{r/B}; this implementation uses the bounded form
#' \eqn{(1 + r)/(B + 1)}, as in [sppt()], so no p-value is exactly zero. The
#' bundled `vancouver_da` data have 1,016 dissemination areas, so results
#' differ slightly from the 1,019 units analysed there.
#'
#' @references
#' Andresen, M.A. (2026). Testing longitudinal change in spatial
#' distributions and patterns: An extension of the spatial point pattern
#' test. *Annals of the American Association of Geographers*, in press.
#' \doi{10.1080/24694452.2026.2733304}
#'
#' @seealso [sppt()], [sppt_map()], [vancouver_ct]
#' @examples
#' \donttest{
#' res <- sppt_longitudinal(vancouver_ct, group_col = "id",
#'                          count_col = paste0("BNER", 2003:2025),
#'                          B = 199, seed = 42, verbose = FALSE,
#'                          create_maps = FALSE, create_time_plot = FALSE)
#' table(res$BNER_Trend_p)
#' attr(res, "si_overlap_p")
#' }
#'
#' @export
sppt_longitudinal <- function(data,
                              group_col           = "group",
                              count_col,
                              B                   = 4999L,
                              seed                = NULL,
                              alpha               = 0.05,
                              calculate_pvalues   = TRUE,
                              conf_level          = 0.85,
                              fix_base            = FALSE,
                              use_percentages     = TRUE,
                              prefix              = NULL,
                              min_zero_streak     = 3L,
                              trend_type          = "linear",
                              trend_alpha         = 0.05,
                              jenks               = "optimal",
                              jenks_gvf_threshold = 0.05,
                              verbose             = TRUE,
                              create_maps         = TRUE,
                              create_time_plot    = TRUE,
                              export_dir          = NULL,
                              export_maps         = FALSE,
                              map_format          = "png",
                              map_dpi             = 300L,
                              map_width_in        = 7,
                              map_height_in       = 5.5,
                              include_title       = TRUE,
                              figure_prefix       = NULL,
                              export_time_plot    = FALSE,
                              export_tables       = FALSE,
                              export_results      = FALSE,
                              export_format       = "gpkg") {

  say <- say_if(verbose)

  # ---------------------------------------------------------------------------
  # Validate parameters
  # ---------------------------------------------------------------------------
  if (!inherits(data, "sf"))
    stop("`data` must be an sf object.", call. = FALSE)
  if (!group_col %in% names(data))
    stop("`group_col` '", group_col, "' not found in data.", call. = FALSE)
  if (anyDuplicated(data[[group_col]]) > 0L)
    stop("`group_col` '", group_col, "' does not uniquely identify rows.",
         call. = FALSE)
  if (length(count_col) < 3L)
    stop("`count_col` must name at least 3 periods; use sppt() to compare ",
         "two.", call. = FALSE)
  for (cc in count_col) {
    if (!cc %in% names(data))
      stop("`count_col` column '", cc, "' not found in data.", call. = FALSE)
    if (!is.numeric(data[[cc]]))
      stop("`count_col` column '", cc, "' must be numeric.", call. = FALSE)
    if (any(!is.na(data[[cc]]) & data[[cc]] < 0))
      stop("`count_col` column '", cc, "' contains negative values.",
           call. = FALSE)
  }
  if (!is.numeric(B) || length(B) != 1L || B < 99)
    stop("`B` must be a single integer >= 99.", call. = FALSE)
  B <- as.integer(B)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1)
    stop("`alpha` must be a single numeric value in (0, 1).", call. = FALSE)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be a single numeric value in (0, 1).",
         call. = FALSE)
  if (!is.numeric(trend_alpha) || length(trend_alpha) != 1L ||
      trend_alpha <= 0 || trend_alpha >= 1)
    stop("`trend_alpha` must be a single numeric value in (0, 1).",
         call. = FALSE)
  trend_type <- match.arg(trend_type, c("linear", "quadratic"))
  if (!identical(jenks, "optimal") &&
      !(is.numeric(jenks) && length(jenks) == 1L && jenks >= 2 &&
        jenks == round(jenks)))
    stop("`jenks` must be \"optimal\" or an integer >= 2.", call. = FALSE)
  if (!is.numeric(jenks_gvf_threshold) || jenks_gvf_threshold <= 0 ||
      jenks_gvf_threshold >= 1)
    stop("`jenks_gvf_threshold` must be numeric in (0, 1).", call. = FALSE)
  map_format <- match.arg(map_format, c("png", "tif", "eps", "svg", "pdf"))

  if (export_maps || export_time_plot || export_tables || export_results)
    export_dir <- check_export_dir(
      export_dir, c("export_maps", "export_time_plot", "export_tables",
                    "export_results")[c(export_maps, export_time_plot,
                                        export_tables, export_results)][1])

  if (fix_base)
    warning(
      "fix_base = TRUE produces an asymmetric test: the base period is not ",
      "bootstrapped.\nNot recommended for standard use.",
      call. = FALSE
    )

  # Quadratic requires >= 4 periods
  if (trend_type == "quadratic" && length(count_col) < 4L) {
    say("Note: Defaulting to linear trend - quadratic requires 4+ time periods.\n")
    say("      You have", length(count_col), "time periods.\n\n")
    trend_type <- "linear"
  }

  # Prefix for output names, e.g. "BNER2003" -> "BNER_"
  if (is.null(prefix)) {
    var_stem <- gsub("[0-9]+$", "", count_col[1])
    prefix   <- if (nchar(var_stem) > 0) paste0(var_stem, "_") else ""
  } else if (nchar(prefix) > 0 && !endsWith(prefix, "_")) {
    prefix <- paste0(prefix, "_")
  }
  pv <- function(nm) paste0(prefix, nm)

  if (!is.null(seed)) set.seed(seed)

  result <- data
  N  <- nrow(result)
  T_ <- length(count_col)
  lower_cols <- paste0(count_col, "_L")
  upper_cols <- paste0(count_col, "_U")

  alpha_ci   <- 1 - conf_level
  lower_prob <- alpha_ci / 2
  upper_prob <- 1 - alpha_ci / 2

  # ---------------------------------------------------------------------------
  # Bootstrap envelopes for each period. The N x B matrices are kept for the
  # p-values, so no second bootstrap is needed. Rows follow the row order of
  # `data`.
  # ---------------------------------------------------------------------------
  boot_list <- vector("list", T_)

  for (i in seq_len(T_)) {
    obs_counts_i <- result[[count_col[i]]]
    obs_counts_i[is.na(obs_counts_i)] <- 0
    n_i <- sum(obs_counts_i)

    if (fix_base && i == 1L) {
      # Base not bootstrapped: envelope collapses to the observed values
      obs_i <- if (use_percentages) {
        if (n_i > 0) obs_counts_i / n_i * 100 else rep(0, N)
      } else obs_counts_i
      result[[lower_cols[i]]] <- obs_i
      result[[upper_cols[i]]] <- obs_i
      boot_list[[i]] <- matrix(obs_i, nrow = N, ncol = B)
    } else if (n_i == 0) {
      result[[lower_cols[i]]] <- 0
      result[[upper_cols[i]]] <- 0
      boot_list[[i]] <- matrix(0, nrow = N, ncol = B)
    } else {
      if (verbose)
        message("Running ", B, " bootstrap samples on ", n_i,
                " events for ", count_col[i], "...")
      if (!is.null(seed)) set.seed(seed + i - 1)

      # Expand events, build a sparse one-hot matrix, multiply by the
      # multinomial weight matrix
      event_expand_i <- rep(seq_len(N), times = obs_counts_i)
      onehot_i <- Matrix::sparseMatrix(i = seq_len(n_i), j = event_expand_i,
                                       x = 1, dims = c(n_i, N))
      W_i <- stats::rmultinom(B, size = n_i, prob = rep(1 / n_i, n_i))
      boot_counts_i <- as.matrix(Matrix::t(onehot_i) %*% W_i)

      group_values_i <- if (use_percentages) {
        col_sums_i <- colSums(boot_counts_i)
        col_sums_i[col_sums_i == 0] <- 1
        sweep(boot_counts_i, 2, col_sums_i, "/") * 100
      } else {
        boot_counts_i
      }
      boot_list[[i]] <- group_values_i

      result[[lower_cols[i]]] <- apply(group_values_i, 1, stats::quantile,
                                       probs = lower_prob)
      result[[upper_cols[i]]] <- apply(group_values_i, 1, stats::quantile,
                                       probs = upper_prob)
    }
  }

  # ---------------------------------------------------------------------------
  # Jenks classification of the base period
  # ---------------------------------------------------------------------------
  first_period_counts <- data[[count_col[1]]]
  jenks_k <- if (identical(jenks, "optimal")) {
    compute_jenks_k(first_period_counts, k_max = 20,
                    gvf_thresh = jenks_gvf_threshold)
  } else {
    as.integer(jenks)
  }
  say(sprintf("\nJenks classification: k = %d (based on %s)\n", jenks_k,
              if (identical(jenks, "optimal"))
                sprintf("optimal GVF threshold = %.0f%%",
                        jenks_gvf_threshold * 100)
              else "user-specified k"))

  # Breaks on non-zero values, with zero as its own lower bound
  x_nz <- first_period_counts[!is.na(first_period_counts) &
                                first_period_counts > 0]
  jenks_brks <- tryCatch({
    ji <- suppressWarnings(
      classInt::classIntervals(x_nz, n = jenks_k, style = "jenks"))
    c(0, ji$brks)
  }, error = function(e) {
    warning("Jenks classification failed; falling back to equal-count ",
            "breaks. Error: ", conditionMessage(e), call. = FALSE)
    stats::quantile(first_period_counts,
                    probs = seq(0, 1, length.out = jenks_k + 1), na.rm = TRUE)
  })
  jenks_brks <- unique(jenks_brks)

  result[[pv("jenks_level")]] <- as.integer(cut(
    result[[count_col[1]]], breaks = jenks_brks,
    labels = seq_len(length(jenks_brks) - 1), include.lowest = TRUE))
  result[[pv("jenks_min")]] <- jenks_brks[result[[pv("jenks_level")]]]
  result[[pv("jenks_max")]] <- jenks_brks[result[[pv("jenks_level")]] + 1L]
  say(sprintf("Jenks breaks: %s\n\n",
              paste(round(jenks_brks, 1), collapse = " | ")))

  # ---------------------------------------------------------------------------
  # Time trend plot (uses the envelope bounds above)
  # ---------------------------------------------------------------------------
  var_prefix   <- gsub("[0-9]+$", "", count_col[1])
  time_periods <- suppressWarnings(as.numeric(gsub(var_prefix, "", count_col)))
  if (anyNA(time_periods)) time_periods <- seq_len(T_)
  totals_vec <- vapply(count_col, function(v) sum(result[[v]], na.rm = TRUE),
                       numeric(1))

  if (create_time_plot || export_time_plot) {
    if (use_percentages) {
      lower_bounds <- vapply(seq_len(T_), function(i)
        sum(result[[lower_cols[i]]] / 100 * totals_vec[i], na.rm = TRUE),
        numeric(1))
      upper_bounds <- vapply(seq_len(T_), function(i)
        sum(result[[upper_cols[i]]] / 100 * totals_vec[i], na.rm = TRUE),
        numeric(1))
    } else {
      lower_bounds <- vapply(lower_cols, function(v)
        sum(result[[v]], na.rm = TRUE), numeric(1))
      upper_bounds <- vapply(upper_cols, function(v)
        sum(result[[v]], na.rm = TRUE), numeric(1))
    }
    draw_time_plot <- function() {
      old_par <- graphics::par(mar = c(7, 4, 4, 2))
      on.exit(graphics::par(old_par), add = TRUE)
      y_range <- range(c(totals_vec, lower_bounds, upper_bounds),
                       na.rm = TRUE)
      y_limits <- c(y_range[1], y_range[2] + diff(y_range) * 0.15)
      graphics::plot(time_periods, totals_vec, type = "b",
                     xlab = "Time Period", ylab = paste("Total", var_prefix),
                     main = paste(var_prefix, "Trend Over Time"),
                     pch = 19, col = "blue", xaxt = "n", ylim = y_limits,
                     cex = 1.2, lwd = 2)
      graphics::lines(time_periods, lower_bounds, col = "red", lty = 2,
                      lwd = 1.5)
      graphics::lines(time_periods, upper_bounds, col = "red", lty = 2,
                      lwd = 1.5)
      graphics::axis(1, at = time_periods, labels = time_periods, las = 2,
                     cex.axis = 0.9)
      graphics::grid()
      graphics::text(time_periods, totals_vec,
                     labels = format(totals_vec, big.mark = ","),
                     pos = 3, cex = 0.7, col = "darkblue")
      graphics::legend("topright",
                       legend = c("Observed",
                                  paste0(conf_level * 100, "% CI Bounds")),
                       col = c("blue", "red"), lty = c(1, 2),
                       lwd = c(2, 1.5), pch = c(19, NA), bty = "n",
                       cex = 0.9)
    }
    if (create_time_plot) draw_time_plot()
    if (export_time_plot) {
      tp_path <- file.path(export_dir,
                           paste0("time_trend_", var_prefix, ".png"))
      grDevices::png(filename = tp_path, width = 8 * map_dpi,
                     height = 6 * map_dpi, res = map_dpi)
      tryCatch(draw_time_plot(), finally = grDevices::dev.off())
      say("Time trend plot exported to:", tp_path, "\n\n")
    }
  }

  # ---------------------------------------------------------------------------
  # Envelope overlap and sustained zeros
  # ---------------------------------------------------------------------------
  say("Calculating overlap statistics...\n")
  L_mat <- as.matrix(sf::st_drop_geometry(result)[, lower_cols, drop = FALSE])
  U_mat <- as.matrix(sf::st_drop_geometry(result)[, upper_cols, drop = FALSE])
  result[[pv("intervals_overlap")]] <- vapply(seq_len(N), function(r)
    as.integer(max(L_mat[r, ], na.rm = TRUE) <= min(U_mat[r, ], na.rm = TRUE)),
    integer(1))
  result[[pv("max_overlap_proportion")]] <- vapply(seq_len(N), function(r) {
    ms <- find_max_overlap_size(L_mat[r, ], U_mat[r, ])
    if (ms <= 1) 0 else ms / T_
  }, numeric(1))

  say("Calculating sustained zeros...\n")
  C_mat <- as.matrix(sf::st_drop_geometry(result)[, count_col, drop = FALSE])
  if (T_ >= min_zero_streak) {
    last_idx <- (T_ - min_zero_streak + 1):T_
    result[[pv("sustained_zero")]] <- vapply(seq_len(N), function(r)
      as.integer(sum(C_mat[r, ] > 0, na.rm = TRUE) > 0 &&
                   all(C_mat[r, last_idx] == 0)), integer(1))
  } else {
    result[[pv("sustained_zero")]] <- 0L
  }

  # ---------------------------------------------------------------------------
  # Trends (envelope-gated)
  # ---------------------------------------------------------------------------
  say("Calculating trends...\n")
  Y_mat <- if (use_percentages) {
    sweep(C_mat, 2, ifelse(totals_vec > 0, totals_vec, 1), "/") * 100 *
      rep(totals_vec > 0, each = N)
  } else {
    C_mat
  }
  trend_cols <- c("Trend", "Trend_All", "Trend_Linear_Sig", "Linear_R2",
                  "Linear_b0", "Linear_bt", "Quad_R2", "Quad_b0", "Quad_bt",
                  "Quad_bt2", "Quad_Sig")
  env <- fit_unit_trends(Y_mat, result[[pv("intervals_overlap")]] == 1,
                         trend_type, trend_alpha)
  for (nm in trend_cols) result[[pv(nm)]] <- env[[nm]]
  result[[pv("Linear_bt_sig")]] <- ifelse(
    result[[pv("Trend_Linear_Sig")]] == 0, NA, result[[pv("Linear_bt")]])

  # ---------------------------------------------------------------------------
  # Bootstrap p-values (each period vs the base)
  # ---------------------------------------------------------------------------
  pvalue_matrix <- NULL
  if (calculate_pvalues) {

    # lambda_min diagnostic for each period pair
    pu_base <- mean(result[[count_col[1]]], na.rm = TRUE)
    lm_period <- count_col[-1]
    lm_vals <- vapply(count_col[-1], function(v)
      min(pu_base, mean(result[[v]], na.rm = TRUE)), numeric(1))
    lm_crit <- lm_vals < 0.5
    lm_low  <- lm_vals < 1.5

    say("\n")
    say("========================================================================\n")
    say("SIMULATION-BASED P-VALUE CALCULATION\n")
    say("========================================================================\n")
    say(sprintf("Method : Nonparametric centred bootstrap difference (vs %s)\n",
                count_col[1]))
    say(sprintf("Alpha  : %.4g\n", alpha))
    say(sprintf("Units  : %d  |  Periods: %d  |  B: %d\n\n", N, T_, B))
    say(sprintf("  %-15s  %10s  %8s  %8s\n", "Period", "lambda_min", "<1.5",
                "<0.5"))
    say("  ", paste(rep("-", 47), collapse = ""), "\n", sep = "")
    for (i in seq_along(lm_period))
      say(sprintf("  %-15s  %10.3f  %8s  %8s\n", lm_period[i], lm_vals[i],
                  if (lm_low[i]) "[!]" else "",
                  if (lm_crit[i]) "[!!]" else ""))
    say("  [!] = low power (lambda_min < 1.5); [!!] = critical (lambda_min < 0.5)\n")

    n_low  <- sum(lm_low)
    n_crit <- sum(lm_crit)
    if (n_crit > 0 || n_low > 0) {
      msg <- character(0)
      if (n_crit > 0)
        msg <- c(msg, sprintf(
          "%d of %d period pairs (%.1f%%) have lambda_min < 0.5 [critical: severely limited power]",
          n_crit, T_ - 1, 100 * n_crit / (T_ - 1)))
      if (n_low > n_crit)
        msg <- c(msg, sprintf(
          "%d of %d period pairs (%.1f%%) have lambda_min < 1.5 [low power: interpret with caution]",
          n_low, T_ - 1, 100 * n_low / (T_ - 1)))
      warning(paste(c("\n*** POWER WARNING ***", msg,
                      "Consider aggregating to a coarser spatial unit where lambda_min >= 1.5."),
                    collapse = "\n"), call. = FALSE)
    }
    say("\n")

    pvalue_matrix <- matrix(NA_real_, nrow = N, ncol = T_)
    colnames(pvalue_matrix) <- count_col

    say("Computing p-values...\n")
    boot_base    <- boot_list[[1]]
    cnt_base     <- C_mat[, 1]
    tot_base     <- totals_vec[1]
    obs_pct_base <- if (tot_base > 0) cnt_base / tot_base * 100 else rep(0, N)

    for (t in 2:T_) {
      cnt_x2     <- C_mat[, t]
      tot_x2     <- totals_vec[t]
      obs_pct_x2 <- if (tot_x2 > 0) cnt_x2 / tot_x2 * 100 else rep(0, N)
      obs_diff   <- if (use_percentages) obs_pct_x2 - obs_pct_base
                    else cnt_x2 - cnt_base
      boot_diffs         <- boot_list[[t]] - boot_base
      boot_diffs_centred <- boot_diffs - rowMeans(boot_diffs)
      # Bounded Monte Carlo p-value, p = (1 + r) / (B + 1), as in sppt()
      pvalue_matrix[, t] <- (1 + rowSums(abs(boot_diffs_centred) >=
                                           abs(obs_diff))) / (B + 1)
    }

    for (t in seq_len(T_))
      result[[paste0("pvalue_", count_col[t])]] <- pvalue_matrix[, t]

    say("========================================================================\n")
    say(sprintf("  %-15s  %10s  %12s\n", "Period", "N_events",
                "Sig (p<=alpha)"))
    say("  ", paste(rep("-", 42), collapse = ""), "\n", sep = "")
    for (t in seq_len(T_)) {
      n_sig <- if (t == 1) NA else sum(pvalue_matrix[, t] <= alpha,
                                       na.rm = TRUE)
      say(sprintf("  %-15s  %10s  %12s\n", count_col[t],
                  format(totals_vec[t], big.mark = ","),
                  if (is.na(n_sig)) "(baseline)" else as.character(n_sig)))
    }
    say("========================================================================\n\n")

    # P-value analogues of the envelope variables
    is_sig <- matrix(0L, nrow = N, ncol = T_)
    for (t in 2:T_)
      is_sig[, t] <- as.integer(!is.na(pvalue_matrix[, t]) &
                                  pvalue_matrix[, t] <= alpha)
    result[[pv("intervals_overlap_p")]]      <- as.integer(rowSums(is_sig) == 0)
    result[[pv("max_overlap_proportion_p")]] <- 1 - rowMeans(is_sig)
    result[[pv("prop_sig_p")]]               <- rowMeans(is_sig)
    result[[pv("sustained_zero_p")]]         <- result[[pv("sustained_zero")]]

    say("Calculating p-value trends...\n")
    pvt <- fit_unit_trends(Y_mat, result[[pv("intervals_overlap_p")]] == 1,
                           trend_type, trend_alpha)
    p_names <- c(Trend = "Trend_p", Trend_All = "Trend_p_All",
                 Trend_Linear_Sig = "Trend_p_Linear_Sig",
                 Linear_R2 = "Linear_p_R2", Linear_b0 = "Linear_p_b0",
                 Linear_bt = "Linear_p_bt", Quad_R2 = "Quad_p_R2",
                 Quad_b0 = "Quad_p_b0", Quad_bt = "Quad_p_bt",
                 Quad_bt2 = "Quad_p_bt2", Quad_Sig = "Quad_p_Sig")
    for (nm in names(p_names)) result[[pv(p_names[[nm]])]] <- pvt[[nm]]
    result[[pv("Linear_p_bt_sig")]] <- ifelse(
      result[[pv("Trend_p_Linear_Sig")]] == 0, NA, result[[pv("Linear_p_bt")]])
  }

  # ---------------------------------------------------------------------------
  # Summary indices
  # ---------------------------------------------------------------------------
  say("Computing summary indices...\n")
  nonzero_mask <- rowSums(C_mat > 0, na.rm = TRUE) > 0
  nonzero_obs  <- sum(nonzero_mask)

  # Shares of units with each trend code, over all units and over units with
  # at least one event
  trend_codes <- if (trend_type == "quadratic") {
    c(increasing = 1, decreasing = -1, pos_lin_neg_quad = 2,
      pos_lin_pos_quad = 3, neg_lin_pos_quad = -2, neg_lin_neg_quad = -3,
      not_significant = 0)
  } else {
    c(increasing = 2, decreasing = -2, nonsig_increasing = 1,
      nonsig_decreasing = -1)
  }
  index_set <- function(overlap, average, sustained, trend) {
    either <- overlap == 1 | sustained == 1
    share  <- function(x, keep) {
      n <- sum(keep)
      if (n == 0) NA_real_ else sum(x[keep], na.rm = TRUE) / n
    }
    mk <- function(keep) c(
      overlap   = share(overlap, keep),
      average   = if (sum(keep) == 0) NA_real_ else
        mean(average[keep], na.rm = TRUE),
      sustained = share(sustained, keep),
      vapply(trend_codes, function(code) share(trend == code, keep),
             numeric(1)),
      combined  = share(either, keep))
    list(all = mk(rep(TRUE, N)), robust = mk(nonzero_mask))
  }
  si_env <- index_set(result[[pv("intervals_overlap")]],
                      result[[pv("max_overlap_proportion")]],
                      result[[pv("sustained_zero")]], result[[pv("Trend")]])
  si_p <- if (calculate_pvalues)
    index_set(result[[pv("intervals_overlap_p")]],
              result[[pv("max_overlap_proportion_p")]],
              result[[pv("sustained_zero_p")]], result[[pv("Trend_p")]])

  # ---------------------------------------------------------------------------
  # Jenks trend tables
  # ---------------------------------------------------------------------------
  first_period <- count_col[1]
  jenks_range_labels <- sprintf("[%.0f, %.0f]",
                                jenks_brks[-length(jenks_brks)],
                                jenks_brks[-1])
  result[[pv("baseline_interval")]] <- factor(
    result[[pv("jenks_level")]], levels = seq_len(length(jenks_brks) - 1L),
    labels = jenks_range_labels)

  interval_tables <- function(lin_all, lin_sig, trend) {
    df <- sf::st_drop_geometry(result)
    g  <- df[[pv("baseline_interval")]]
    keep <- !is.na(g)
    by_level <- function(f) {
      out <- tapply(seq_len(N)[keep], g[keep], f)
      out[is.na(out)] <- 0L
      as.vector(out)
    }
    lev <- factor(levels(g), levels = levels(g))
    present <- levels(g) %in% unique(as.character(g[keep]))
    base_cnt <- df[[first_period]]
    lin <- data.frame(
      baseline_interval = lev,
      Min_Count      = by_level(function(ix) min(base_cnt[ix], na.rm = TRUE)),
      Max_Count      = by_level(function(ix) max(base_cnt[ix], na.rm = TRUE)),
      Increasing_All = by_level(function(ix) sum(df[[lin_all]][ix] == 1,
                                                 na.rm = TRUE)),
      Decreasing_All = by_level(function(ix) sum(df[[lin_all]][ix] == -1,
                                                 na.rm = TRUE)),
      Sig_Increasing = by_level(function(ix) sum(df[[lin_sig]][ix] == 1,
                                                 na.rm = TRUE)),
      Sig_Decreasing = by_level(function(ix) sum(df[[lin_sig]][ix] == -1,
                                                 na.rm = TRUE))
    )[present, , drop = FALSE]
    names(lin)[1] <- pv("baseline_interval")
    quadt <- NULL
    if (trend_type == "quadratic") {
      cnt <- function(code) by_level(function(ix)
        sum(df[[trend]][ix] == code, na.rm = TRUE))
      quadt <- data.frame(
        baseline_interval = lev,
        Min_Count       = lin$Min_Count[match(levels(g), levels(g)[present])],
        Max_Count       = lin$Max_Count[match(levels(g), levels(g)[present])],
        Not_Significant = cnt(0), Sig_Increasing = cnt(1),
        Sig_Decreasing  = cnt(-1), SigInc_NegQuad = cnt(2),
        SigInc_PosQuad  = cnt(3), SigDec_PosQuad = cnt(-2),
        SigDec_NegQuad  = cnt(-3)
      )[present, , drop = FALSE]
      names(quadt)[1] <- pv("baseline_interval")
      rownames(quadt) <- NULL
    }
    rownames(lin) <- NULL
    list(linear = lin, quad = quadt)
  }
  it_env <- interval_tables(pv("Trend_All"), pv("Trend_Linear_Sig"),
                            pv("Trend"))
  it_p <- if (calculate_pvalues)
    interval_tables(pv("Trend_p_All"), pv("Trend_p_Linear_Sig"),
                    pv("Trend_p"))

  # ---------------------------------------------------------------------------
  # Console report
  # ---------------------------------------------------------------------------
  if (verbose)
    print_long_report(
      say = say, time_periods = time_periods, totals_vec = totals_vec,
      T_ = T_, trend_type = trend_type, use_percentages = use_percentages,
      conf_level = conf_level, alpha = alpha, trend_alpha = trend_alpha,
      si_env = si_env, si_p = si_p, result = result, pv = pv, N = N,
      nonzero_obs = nonzero_obs, first_period = first_period,
      jenks_k = jenks_k, jenks_gvf_threshold = jenks_gvf_threshold,
      it_env = it_env, it_p = it_p, calculate_pvalues = calculate_pvalues)

  # ---------------------------------------------------------------------------
  # Export tables
  # ---------------------------------------------------------------------------
  if (export_tables) {
    fp <- file.path(export_dir,
                    paste0("sppt_longitudinal_tables_", var_prefix, ".csv"))
    tryCatch({
      write_long_tables(fp, time_periods, totals_vec, T_, trend_type,
                        use_percentages, conf_level, alpha, trend_alpha,
                        si_env, si_p, first_period, it_env, it_p,
                        calculate_pvalues)
      say("Tables exported to:", fp, "\n\n")
    }, error = function(e)
      warning("Failed to export tables: ", conditionMessage(e), call. = FALSE))
  }

  # ---------------------------------------------------------------------------
  # Attributes
  # ---------------------------------------------------------------------------
  for (nm in names(si_env$all))
    attr(result, paste0("s_index_", nm)) <- unname(si_env$all[[nm]])
  for (nm in names(si_env$robust))
    attr(result, paste0("robust_s_index_", nm)) <- unname(si_env$robust[[nm]])
  attr(result, "jenks_k")               <- jenks_k
  attr(result, "jenks_breaks")          <- jenks_brks
  attr(result, "jenks_gvf_threshold")   <- jenks_gvf_threshold
  attr(result, "use_percentages")       <- use_percentages
  attr(result, "min_zero_streak")       <- min_zero_streak
  attr(result, "trend_type")            <- trend_type
  attr(result, "interval_table_linear") <- it_env$linear
  attr(result, "interval_table_quad")   <- it_env$quad
  if (calculate_pvalues) {
    attr(result, "pvalue_matrix") <- pvalue_matrix
    attr(result, "alpha")         <- alpha
    for (nm in names(si_p$all))
      attr(result, paste0("si_", nm, "_p")) <- unname(si_p$all[[nm]])
    for (nm in names(si_p$robust))
      attr(result, paste0("robust_si_", nm, "_p")) <-
        unname(si_p$robust[[nm]])
    attr(result, "interval_table_p_linear") <- it_p$linear
    attr(result, "interval_table_p_quad")   <- it_p$quad
  }
  attr(result, "sppt_long_meta") <- list(
    prefix = prefix, trend_type = trend_type, alpha = alpha,
    calculate_pvalues = calculate_pvalues, count_col = count_col,
    conf_level = conf_level, trend_alpha = trend_alpha, B = B,
    n_units = N, nonzero_units = nonzero_obs, jenks_k = jenks_k,
    totals = stats::setNames(as.numeric(totals_vec), count_col))
  class(result) <- c("sppt_long_result", class(result))

  # ---------------------------------------------------------------------------
  # Maps
  # ---------------------------------------------------------------------------
  if (create_maps || export_maps) {
    say("Creating maps...\n")
    maps <- tryCatch(build_maps_long(result), error = function(e) {
      message("Could not create maps: ", conditionMessage(e))
      NULL
    })
    if (!is.null(maps)) {
      if (create_maps) for (m in maps) print(m)
      if (export_maps) {
        for (k in seq_along(maps)) {
          nm <- names(maps)[k]
          stem <- if (!is.null(figure_prefix)) {
            paste0(figure_prefix, k)
          } else if (nm == "p_trend") {
            if (trend_type == "quadratic") "pmap9_quadratic_trends_p"
            else "pmap9_linear_trends_p"
          } else {
            long_map_files[[nm]]
          }
          p_out <- drop_title(maps[[nm]], include_title)
          ggplot2::ggsave(
            filename = file.path(export_dir, paste0(stem, ".", map_format)),
            plot = p_out, width = map_width_in, height = map_height_in,
            dpi = map_dpi, units = "in")
        }
        say("Maps exported to:", export_dir, "\n")
      }
      say(sprintf("%d maps created.\n\n", length(maps)))
      attr(result, "sppt_maps") <- maps
    }
  }

  # ---------------------------------------------------------------------------
  # Export results
  # ---------------------------------------------------------------------------
  if (export_results) {
    fp <- tryCatch(
      write_result(result, export_dir,
                   paste0("sppt_longitudinal_output_", var_prefix),
                   export_format),
      error = function(e) {
        warning("Failed to export results: ", conditionMessage(e),
                call. = FALSE)
        NULL
      })
    if (!is.null(fp)) say("Results exported:", fp, "\n")
  }

  result
}

# -- Internal helpers ---------------------------------------------------------

# Optimal number of Jenks classes by the GVF stopping rule.
compute_jenks_k <- function(x, k_max = 20, gvf_thresh = 0.05) {
  x_nz <- x[!is.na(x) & x > 0]
  if (length(x_nz) < 3) return(2L)
  sdam <- sum((x_nz - mean(x_nz))^2)
  if (sdam == 0) return(2L)
  prev_gvf <- 0
  best_k   <- 2L
  for (k in 2:min(k_max, length(x_nz) - 1)) {
    gvf_k <- tryCatch({
      brks      <- suppressWarnings(
        classInt::classIntervals(x_nz, n = k, style = "jenks"))
      classes   <- classInt::findCols(brks)
      cls_means <- tapply(x_nz, classes, mean)
      sdcm      <- sum((x_nz - cls_means[classes])^2)
      (sdam - sdcm) / sdam
    }, error = function(e) NA_real_)
    if (is.na(gvf_k)) break
    if (k > 2 && (gvf_k - prev_gvf) < gvf_thresh) {
      best_k <- k - 1L
      break
    }
    best_k   <- k
    prev_gvf <- gvf_k
  }
  best_k
}

# Size of the largest set of mutually overlapping intervals (greedy clique).
find_max_overlap_size <- function(lowers, uppers) {
  n <- length(lowers)
  overlap_matrix <- matrix(FALSE, n, n)
  for (ii in seq_len(n)) for (jj in seq_len(n)) {
    overlap_matrix[ii, jj] <- if (ii == jj) TRUE else
      max(lowers[ii], lowers[jj]) <= min(uppers[ii], uppers[jj])
  }
  degrees <- rowSums(overlap_matrix)
  clique <- which.max(degrees)
  candidates <- setdiff(seq_len(n), clique)
  while (length(candidates) > 0) {
    found <- FALSE
    for (v in candidates) {
      if (all(overlap_matrix[v, clique])) {
        clique <- c(clique, v)
        candidates <- setdiff(candidates, v)
        found <- TRUE
        break
      }
    }
    if (!found) break
  }
  length(clique)
}

# OLS trend summaries for every unit. Units where `stable` is TRUE get the
# stable defaults (0 for classes, NA for coefficients).
fit_unit_trends <- function(Y, stable, trend_type, alpha_trend) {
  n <- nrow(Y)
  out <- list(
    Trend = integer(n), Trend_All = integer(n), Trend_Linear_Sig = integer(n),
    Linear_R2 = rep(NA_real_, n), Linear_b0 = rep(NA_real_, n),
    Linear_bt = rep(NA_real_, n), Quad_R2 = rep(NA_real_, n),
    Quad_b0 = rep(NA_real_, n), Quad_bt = rep(NA_real_, n),
    Quad_bt2 = rep(NA_real_, n), Quad_Sig = integer(n))
  quad <- trend_type == "quadratic"
  time <- seq_len(ncol(Y))
  tsq  <- time^2

  for (r in which(!stable)) {
    y <- as.numeric(Y[r, ])
    fit_l <- tryCatch(suppressWarnings(stats::lm(y ~ time)),
                      error = function(e) NULL)
    fit_q <- if (quad) tryCatch(suppressWarnings(stats::lm(y ~ time + tsq)),
                                error = function(e) NULL)
    sl <- if (is.null(fit_l)) NA_real_ else unname(stats::coef(fit_l)[2])
    p_lin <- function() suppressWarnings(summary(fit_l))$coefficients[2, 4]

    out$Trend[r] <- tryCatch({
      if (quad) {
        lin_b  <- stats::coef(fit_q)[2]
        quad_b <- stats::coef(fit_q)[3]
        cf     <- suppressWarnings(summary(fit_q))$coefficients
        lin_sig  <- cf[2, 4] <= alpha_trend
        quad_sig <- cf[3, 4] <= alpha_trend
        if      (lin_b >  0 && lin_sig && !quad_sig)               1L
        else if (lin_b <  0 && lin_sig && !quad_sig)              -1L
        else if (lin_b >  0 && lin_sig && quad_sig && quad_b < 0)  2L
        else if (lin_b >  0 && lin_sig && quad_sig && quad_b > 0)  3L
        else if (lin_b <  0 && lin_sig && quad_sig && quad_b > 0) -2L
        else if (lin_b <  0 && lin_sig && quad_sig && quad_b < 0) -3L
        else 0L
      } else {
        pl <- p_lin()
        if (sl > 0) { if (pl <= alpha_trend) 2L else 1L } else
          { if (pl <= alpha_trend) -2L else -1L }
      }
    }, error = function(e) 0L)

    out$Trend_All[r] <- tryCatch(
      if (sl > 0) 1L else if (sl < 0) -1L else 0L,
      error = function(e) 0L)

    out$Trend_Linear_Sig[r] <- tryCatch({
      pl <- p_lin()
      if      (sl > 0 && pl <= alpha_trend)  1L
      else if (sl < 0 && pl <= alpha_trend) -1L
      else 0L
    }, error = function(e) 0L)

    if (!is.null(fit_l)) {
      out$Linear_R2[r] <- tryCatch(suppressWarnings(summary(fit_l))$r.squared,
                                   error = function(e) NA_real_)
      out$Linear_b0[r] <- unname(stats::coef(fit_l)[1])
      out$Linear_bt[r] <- sl
    }

    if (quad && !is.null(fit_q)) {
      out$Quad_R2[r] <- tryCatch(suppressWarnings(summary(fit_q))$r.squared,
                                 error = function(e) NA_real_)
      cq <- stats::coef(fit_q)
      out$Quad_b0[r]  <- unname(cq[1])
      out$Quad_bt[r]  <- unname(cq[2])
      out$Quad_bt2[r] <- unname(cq[3])
      out$Quad_Sig[r] <- tryCatch({
        qp <- suppressWarnings(summary(fit_q))$coefficients[3, 4]
        if      (cq[3] > 0 && qp <= alpha_trend)  1L
        else if (cq[3] < 0 && qp <= alpha_trend) -1L
        else 0L
      }, error = function(e) 0L)
    }
  }
  out
}

# Labels for the S index tables, keyed by the names used in index_set().
long_index_labels <- function(trend_type, p = FALSE) {
  if (!p) {
    base <- c(overlap = "S-Index (Complete Overlap)",
              average = "S-Index Average",
              sustained = "S-Index (Sustained Zeros)",
              increasing = "S-Index (Sig Increasing Trend)",
              decreasing = "S-Index (Sig Decreasing Trend)")
    extra <- if (trend_type == "quadratic") c(
      pos_lin_neg_quad = "S-Index (Pos Lin, Neg Quad)",
      pos_lin_pos_quad = "S-Index (Pos Lin, Pos Quad)",
      neg_lin_pos_quad = "S-Index (Neg Lin, Pos Quad)",
      neg_lin_neg_quad = "S-Index (Neg Lin, Neg Quad)",
      not_significant  = "S-Index (Not Significant)")
    else c(nonsig_increasing = "S-Index (Nonsig Increasing)",
           nonsig_decreasing = "S-Index (Nonsig Decreasing)")
    c(base, extra, combined = "Combined S-Index")
  } else {
    base <- c(overlap = "SI (All Periods Non-Sig)",
              average = "SI Average (Prop Non-Sig Periods)",
              sustained = "SI (Sustained Zeros)",
              increasing = "SI (Sig Increasing Trend_p)",
              decreasing = "SI (Sig Decreasing Trend_p)")
    extra <- if (trend_type == "quadratic") c(
      pos_lin_neg_quad = "SI (Pos Lin, Neg Quad)",
      pos_lin_pos_quad = "SI (Pos Lin, Pos Quad)",
      neg_lin_pos_quad = "SI (Neg Lin, Pos Quad)",
      neg_lin_neg_quad = "SI (Neg Lin, Neg Quad)",
      not_significant  = "SI (Not Significant)")
    else c(nonsig_increasing = "SI (Nonsig Increasing Trend_p)",
           nonsig_decreasing = "SI (Nonsig Decreasing Trend_p)")
    c(base, extra, combined = "Combined SI")
  }
}

# Observation counts reported alongside the indices.
long_obs_counts <- function(result, pv, trend_type, p, N, nonzero_obs) {
  ov <- if (p) "intervals_overlap_p" else "intervals_overlap"
  sz <- if (p) "sustained_zero_p" else "sustained_zero"
  tr <- if (p) "Trend_p" else "Trend"
  qs <- if (p) "Quad_p_Sig" else "Quad_Sig"
  n_code <- function(code) sum(result[[pv(tr)]] == code, na.rm = TRUE)
  out <- c("Total observations" = N,
           "Observations with complete overlap" =
             sum(result[[pv(ov)]], na.rm = TRUE),
           "Observations with sustained zeros" =
             sum(result[[pv(sz)]], na.rm = TRUE))
  if (trend_type == "quadratic") {
    out <- c(out,
             "Observations with sig increasing trend" = n_code(1),
             "Observations with sig decreasing trend" = n_code(-1),
             "Observations with Pos Lin, Neg Quad" = n_code(2),
             "Observations with Pos Lin, Pos Quad" = n_code(3),
             "Observations with Neg Lin, Pos Quad" = n_code(-2),
             "Observations with Neg Lin, Neg Quad" = n_code(-3),
             "Observations not significant" = n_code(0),
             "Obs with sig convex quad" =
               sum(result[[pv(qs)]] == 1, na.rm = TRUE),
             "Obs with sig concave quad" =
               sum(result[[pv(qs)]] == -1, na.rm = TRUE),
             "Obs not sig quad" = sum(result[[pv(qs)]] == 0, na.rm = TRUE))
  } else {
    out <- c(out,
             "Observations with sig increasing trend" = n_code(2),
             "Observations with sig decreasing trend" = n_code(-2),
             "Observations with nonsig incr trend" = n_code(1),
             "Observations with nonsig decr trend" = n_code(-1))
  }
  if (p) names(out)[2] <- "Observations all periods non-sig"
  c(out, "Observations with non-zero counts" = nonzero_obs)
}

print_long_report <- function(say, time_periods, totals_vec, T_, trend_type,
                              use_percentages, conf_level, alpha,
                              trend_alpha, si_env, si_p, result, pv, N,
                              nonzero_obs, first_period, jenks_k,
                              jenks_gvf_threshold, it_env, it_p,
                              calculate_pvalues) {
  say("\n========================================\n")
  say("TOTAL COUNTS BY TIME PERIOD\n")
  say("========================================\n")
  say(sprintf("%-10s %15s\n", "Period", "Total Count"))
  say("----------------------------------------\n")
  for (i in seq_along(time_periods))
    say(sprintf("%-10s %15s\n", format(time_periods[i]),
                format(totals_vec[i], big.mark = ",")))
  say("========================================\n")
  say("Longitudinal Spatial Pattern Statistics\n")
  say(sprintf("Time Periods: %d\n", T_))
  say(sprintf("Trend Type: %s\n", tools::toTitleCase(trend_type)))
  say(if (use_percentages) "Using: Percentages\n" else "Using: Counts\n")
  say(sprintf("Conf Level (CI): %.1f%%\n", conf_level * 100))
  say(sprintf("Alpha (p-value): %.4g\n", alpha))
  say(sprintf("Trend Alpha (OLS): %.4g\n", trend_alpha))
  say("========================================\n")

  show_indices <- function(si, labels, robust_prefix) {
    for (nm in names(labels))
      say(sprintf("  %-34s %s\n", paste0(labels[[nm]], ":"),
                  format(round(si$all[[nm]], 4))))
    say("----------------------------------------\n")
    say(robust_prefix, "\n")
    for (nm in names(labels))
      say(sprintf("  %-34s %s\n", paste0("Robust ", labels[[nm]], ":"),
                  format(round(si$robust[[nm]], 4))))
  }
  say("ENVELOPE-BASED INDICES (all locations):\n")
  show_indices(si_env, long_index_labels(trend_type),
               "ENVELOPE-BASED ROBUST INDICES (non-zero):")
  if (calculate_pvalues) {
    say("----------------------------------------\n")
    say(paste0("P-VALUE-BASED INDICES [alpha=", alpha, "] (all locations):\n"))
    show_indices(si_p, long_index_labels(trend_type, p = TRUE),
                 paste0("P-VALUE-BASED INDICES [alpha=", alpha,
                        "] (non-zero):"))
  }

  show_counts <- function(x) for (nm in names(x))
    say(sprintf("%-41s %s\n", paste0(nm, ":"), format(x[[nm]])))
  say("----------------------------------------\n")
  say("ENVELOPE OBSERVATION COUNTS:\n")
  show_counts(long_obs_counts(result, pv, trend_type, FALSE, N, nonzero_obs))
  if (calculate_pvalues) {
    say("----------------------------------------\n")
    say(paste0("P-VALUE OBSERVATION COUNTS [alpha=", alpha, "]:\n"))
    show_counts(long_obs_counts(result, pv, trend_type, TRUE, N, nonzero_obs))
  }
  say("========================================\n\n")

  show_tables <- function(it, header) {
    say(header, "\n")
    say("(Based on", first_period, ")\n")
    say(sprintf("Jenks k = %d | GVF threshold = %.0f%%\n", jenks_k,
                jenks_gvf_threshold * 100))
    say("\n--- LINEAR TREND RESULTS ---\n")
    say(sprintf("%-14s %10s %10s %12s %12s %12s %12s\n", "Jenks Level",
                "Min", "Max", "Incr All", "Decr All", "Sig Incr",
                "Sig Decr"))
    say(paste(rep("-", 90), collapse = ""), "\n")
    lt <- it$linear
    for (i in seq_len(nrow(lt)))
      say(sprintf("%-14s %10.0f %10.0f %12d %12d %12d %12d\n",
                  as.character(lt[[1]][i]), lt$Min_Count[i],
                  lt$Max_Count[i], lt$Increasing_All[i],
                  lt$Decreasing_All[i], lt$Sig_Increasing[i],
                  lt$Sig_Decreasing[i]))
    if (!is.null(it$quad)) {
      qt <- it$quad
      say("\n--- QUADRATIC TREND RESULTS ---\n")
      say(sprintf("%-14s %10s %10s %12s %12s %12s %12s %12s %12s %12s\n",
                  "Jenks Level", "Min", "Max", "Not Sig", "Sig Incr",
                  "Sig Decr", "+L -Q", "+L +Q", "-L +Q", "-L -Q"))
      say(paste(rep("-", 130), collapse = ""), "\n")
      for (i in seq_len(nrow(qt)))
        say(sprintf("%-14s %10.0f %10.0f %12d %12d %12d %12d %12d %12d %12d\n",
                    as.character(qt[[1]][i]), qt$Min_Count[i],
                    qt$Max_Count[i], qt$Not_Significant[i],
                    qt$Sig_Increasing[i], qt$Sig_Decreasing[i],
                    qt$SigInc_NegQuad[i], qt$SigInc_PosQuad[i],
                    qt$SigDec_PosQuad[i], qt$SigDec_NegQuad[i]))
    }
    say("========================================\n\n")
  }
  show_tables(it_env, "TREND BREAKDOWN BY BASELINE JENKS LEVEL")
  if (calculate_pvalues)
    show_tables(it_p, paste0("TREND_P BREAKDOWN BY BASELINE JENKS LEVEL [alpha=",
                             alpha, "]"))
}

write_long_tables <- function(fp, time_periods, totals_vec, T_, trend_type,
                              use_percentages, conf_level, alpha,
                              trend_alpha, si_env, si_p, first_period,
                              it_env, it_p, calculate_pvalues) {
  con <- file(fp, "w")
  on.exit(close(con), add = TRUE)
  wl  <- function(...) writeLines(paste0(...), con)
  wcsv <- function(x) utils::write.csv(x, con, row.names = FALSE)
  si_frame <- function(si, labels) {
    data.frame(
      Metric = c(unname(labels), "",
                 paste("Robust", unname(labels))),
      Value  = c(round(unname(si$all[names(labels)]), 4), NA,
                 round(unname(si$robust[names(labels)]), 4)))
  }

  wl("TOTAL COUNTS BY TIME PERIOD"); wl("")
  wcsv(data.frame(Period = time_periods, Total_Count = unname(totals_vec)))
  wl(""); wl("")
  wl("LONGITUDINAL SPATIAL PATTERN STATISTICS (ENVELOPE)")
  wl("Time Periods: ", T_)
  wl("Trend Type: ", tools::toTitleCase(trend_type))
  wl("Using: ", if (use_percentages) "Percentages" else "Counts")
  wl("Conf Level (CI): ", conf_level * 100, "%")
  wl("Alpha (p-value): ", alpha)
  wl("Trend Alpha (OLS): ", trend_alpha)
  wl("")
  wcsv(si_frame(si_env, long_index_labels(trend_type)))
  wl(""); wl("")
  wl("TREND BREAKDOWN BY BASELINE JENKS LEVEL")
  wl("(Based on ", first_period, ")"); wl("")
  wl("--- LINEAR TREND RESULTS ---"); wl("")
  wcsv(it_env$linear)
  if (!is.null(it_env$quad)) {
    wl(""); wl("--- QUADRATIC TREND RESULTS ---"); wl("")
    wcsv(it_env$quad)
  }
  if (calculate_pvalues) {
    wl(""); wl("")
    wl("P-VALUE-BASED STATISTICS [alpha=", alpha, "]"); wl("")
    wcsv(si_frame(si_p, long_index_labels(trend_type, p = TRUE)))
    wl(""); wl("")
    wl("TREND_P BREAKDOWN BY BASELINE JENKS LEVEL [alpha=", alpha, "]")
    wl("(Based on ", first_period, ")"); wl("")
    wl("--- LINEAR TREND RESULTS ---"); wl("")
    wcsv(it_p$linear)
    if (!is.null(it_p$quad)) {
      wl(""); wl("--- QUADRATIC TREND RESULTS ---"); wl("")
      wcsv(it_p$quad)
    }
  }
  invisible(fp)
}

#' Summarise a longitudinal SPPT result
#'
#' Prints the settings, the S indices and the trend classification counts
#' of a [sppt_longitudinal()] result.
#'
#' @param object,x An object of class `sppt_long_result`, as returned by
#'   [sppt_longitudinal()].
#' @param ... Further arguments passed to the `sf` method.
#' @return `summary()` invisibly returns a list with the settings, the
#'   S indices (all units and units with at least one event) and the trend
#'   counts. `print()` invisibly returns `x`.
#' @examples
#' \donttest{
#' res <- sppt_longitudinal(vancouver_ct, group_col = "id",
#'                          count_col = paste0("BNER", 2003:2010),
#'                          B = 199, seed = 42, verbose = FALSE,
#'                          create_maps = FALSE, create_time_plot = FALSE)
#' summary(res)
#' }
#' @export
summary.sppt_long_result <- function(object, ...) {
  meta <- attr(object, "sppt_long_meta")
  if (is.null(meta)) return(NextMethod())
  col <- function(nm) paste0(meta$prefix, nm)
  idx <- function(stem) {
    nm <- c("overlap", "average", "sustained", "combined")
    stats::setNames(vapply(nm, function(k) {
      v <- attr(object, paste0(stem, k, if (grepl("^robust_si_|^si_", stem)) "_p" else ""))
      if (is.null(v)) NA_real_ else v
    }, numeric(1)), nm)
  }
  env_all <- idx("s_index_")
  env_rob <- idx("robust_s_index_")
  p_all <- if (isTRUE(meta$calculate_pvalues)) idx("si_") else NULL
  p_rob <- if (isTRUE(meta$calculate_pvalues)) idx("robust_si_") else NULL
  trend_env <- table(factor(object[[col("Trend")]]))
  trend_p <- if (isTRUE(meta$calculate_pvalues))
    table(factor(object[[col("Trend_p")]])) else NULL

  cc <- meta$count_col
  cat("\n", rule_line("="), "Longitudinal SPPT Summary\n", rule_line("="),
      sep = "")
  cat(field_lines("Periods", sprintf("%d (%s to %s; base %s)", length(cc),
                                     cc[1], cc[length(cc)], cc[1])))
  cat(field_lines("Units", sprintf("%d (%d with at least one event)",
                                   meta$n_units, meta$nonzero_units)))
  cat(field_lines("Bootstrap B", meta$B))
  cat(field_lines("Envelope level", sprintf("%.1f%%", meta$conf_level * 100)))
  if (isTRUE(meta$calculate_pvalues))
    cat(field_lines("alpha", sprintf("%.3f", meta$alpha)))
  cat(field_lines("Trend", sprintf("%s, alpha = %.3f", meta$trend_type,
                                   meta$trend_alpha)))
  cat(rule_line("-"))
  cat("                      all units   units with events\n")
  row <- function(lbl, a, b)
    cat(sprintf("    %-16s %11.4f %19.4f\n", lbl, a, b))
  if (!is.null(p_all)) {
    cat("  P-value indices\n")
    row("stable (S)", p_all["overlap"], p_rob["overlap"])
    row("average", p_all["average"], p_rob["average"])
    row("sustained zero", p_all["sustained"], p_rob["sustained"])
    row("combined", p_all["combined"], p_rob["combined"])
  }
  cat("  Envelope indices\n")
  row("stable (S)", env_all["overlap"], env_rob["overlap"])
  row("average", env_all["average"], env_rob["average"])
  row("sustained zero", env_all["sustained"], env_rob["sustained"])
  row("combined", env_all["combined"], env_rob["combined"])
  cat(rule_line("-"))
  tab_text <- function(tb)
    paste(sprintf("%s: %d", names(tb), as.integer(tb)), collapse = ", ")
  if (!is.null(trend_p)) cat(field_lines("Trend_p codes", tab_text(trend_p)))
  cat(field_lines("Trend codes", tab_text(trend_env)))
  cat(field_lines("", "see ?sppt_longitudinal for the trend codes"))
  cat(rule_line("="), "\n", sep = "")

  invisible(list(meta = meta,
                 indices = list(envelope = rbind(all = env_all,
                                                 robust = env_rob),
                                pvalue = if (!is.null(p_all))
                                  rbind(all = p_all, robust = p_rob)),
                 trend = list(envelope = trend_env, pvalue = trend_p)))
}

#' @rdname summary.sppt_long_result
#' @export
print.sppt_long_result <- function(x, ...) {
  summary(x)
  NextMethod()
  invisible(x)
}
