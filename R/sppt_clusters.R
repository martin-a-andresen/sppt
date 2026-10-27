#' Test whether clusters of flagged units exceed chance
#'
#' @description
#' At fine spatial scales unit-level power is low and the false discovery
#' rate is high, but chance cannot produce large clusters of false
#' positives. `sppt_clusters()` groups the flagged units of a [sppt()] result
#' into clusters and compares their sizes with the clusters that false
#' positives alone produce on the same map.
#'
#' The chance ceiling is simulated on the study area's own adjacency
#' structure: in each of `nsim` draws a share `alpha` of the eligible units,
#' the null flagging rate of the test, is flagged at random, and the
#' resulting clusters are recorded. For each number of
#' clusters \eqn{k = 1, 2, \ldots}, the ceiling is the smallest size such
#' that chance produces \eqn{k} clusters at least that large in no more than
#' a proportion `level` of draws. The observed map is then tested cluster by
#' cluster: the \eqn{k}-th largest observed cluster against the
#' \eqn{k}-th largest null cluster.
#'
#' @details
#' **Cluster rules.** Under `rule = "contiguous"` flagged units join a
#' cluster only if they share a border or a corner (queen contiguity). Under
#' `rule = "one_gap"` they also join across a single intervening unit of any
#' kind, which suits fine units, where real spatial processes rarely produce
#' solid blocks. One-gap clusters are larger, so their ceilings are higher.
#'
#' **Direction.** With `direction = TRUE` (default) only flagged units with
#' the same sign of difference join, and the null assigns each chance flag a
#' random sign. Each cluster is then a directional pattern ("test > base" or
#' "base > test"); without this, the one-gap rule can merge neighbouring
#' areas of opposite difference into a single uninformative cluster. With
#' `direction = FALSE` units join regardless of sign, as in the simulation
#' that produced the published reference values; directional ceilings are
#' lower because chance flags of one sign are half as dense.
#'
#' **Number of chance flags.** With `null_count = "fixed"` (default) each
#' draw flags exactly `round(alpha * n)` of the `n` eligible units, the
#' design of the simulation behind the published reference values, which it
#' reproduces. With `null_count = "binomial"` each eligible unit is flagged
#' independently with probability `alpha`, so the number of chance flags
#' varies between draws. This is closer to the sampling behaviour of the
#' test and gives somewhat higher ceilings, most noticeably where units are
#' few (for Vancouver census tracts, a single-cluster ceiling of 7 rather
#' than 5 under the one-gap rule).
#'
#' **Reading the result.** A cluster above its ceiling is a pattern that
#' false positives are unlikely to produce, and can be read as such even
#' where its individual units have unreliable p-values. The tests for
#' different \eqn{k} are not adjusted for each other; requiring any one of
#' them to hold, as the published guidance does, is a union of rules each at
#' level `level`. A map with no cluster above its ceiling is not evidence of
#' no difference: power may be too low.
#'
#' **Reference values.** For Vancouver, `alpha = 0.05`, the one-gap rule
#' and undirected clusters (`direction = FALSE`), the published
#' single-cluster ceilings are about
#' 11 units at dissemination areas, 9 at 15 ha hexagons, 7 at 30 ha
#' hexagons and 5 at census tracts. Ceilings depend on the number of units
#' and on how they connect, which is why they are simulated here for each
#' map rather than taken from a table.
#'
#' @param result A result returned by [sppt()] with p-values.
#' @param rule `"one_gap"` (default) or `"contiguous"`. See Details.
#' @param direction Logical. Build clusters within each direction of
#'   difference only. Default `TRUE`.
#' @param level Numeric in (0, 1). Probability that chance reaches the
#'   ceiling. Default `0.05`.
#' @param null_count `"fixed"` (default) or `"binomial"`: how many units
#'   each null draw flags. See Details.
#' @param max_clusters Integer. Largest number of clusters for which a
#'   ceiling is computed and tested. Default `6`.
#' @param nsim Integer. Number of null draws. Default `999`.
#' @param seed Integer or `NULL`. Random seed for the null draws.
#'
#' @return An object of class `sppt_clusters`, a list with:
#' \describe{
#'   \item{`clusters`}{Data frame with one row per observed cluster, largest
#'     first: `cluster`, `size`, `direction` (`"test > base"`,
#'     `"base > test"` or `"mixed"`) and `exceeds`, whether the size reaches
#'     the single-cluster ceiling.}
#'   \item{`membership`}{Integer vector, one element per row of `result`,
#'     giving the cluster of each flagged unit (`NA` for other units).}
#'   \item{`ceiling`}{Data frame of the chance ceiling by number of clusters
#'     required.}
#'   \item{`tests`}{Data frame comparing the \eqn{k}-th largest observed
#'     cluster with its ceiling, with a Monte Carlo p-value
#'     \eqn{(1 + r)/(nsim + 1)}.}
#'   \item{`settings`}{List of the settings used.}
#' }
#'
#' @references
#' Andresen, M.A. (2025). Properties of the spatial point pattern test:
#' Unit-level performance, error geography, and cluster-level inference
#' across scales. Manuscript under review.
#'
#' @seealso [sppt()]
#' @examples
#' res <- sppt(vancouver_ct, group_col = "id",
#'             count_col = c("BNER2003", "BNER2025"),
#'             B = 199, seed = 42, verbose = FALSE, create_maps = FALSE)
#' cl <- sppt_clusters(res, nsim = 199, seed = 1)
#' cl
#' cl$ceiling
#'
#' @export
sppt_clusters <- function(result, rule = c("one_gap", "contiguous"),
                          direction = TRUE, level = 0.05,
                          null_count = c("fixed", "binomial"),
                          max_clusters = 6L, nsim = 999L, seed = NULL) {
  if (!inherits(result, "sppt_result"))
    stop("`result` must be returned by sppt().", call. = FALSE)
  if (!"sig_pvalue_dir" %in% names(result))
    stop("`result` has no p-values; run sppt() with ",
         "calculate_pvalues = TRUE.", call. = FALSE)
  rule <- match.arg(rule)
  null_count <- match.arg(null_count)
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("`level` must be a single numeric value in (0, 1).", call. = FALSE)
  if (!is.numeric(nsim) || length(nsim) != 1L || nsim < 99)
    stop("`nsim` must be a single integer >= 99.", call. = FALSE)
  if (!is.numeric(max_clusters) || length(max_clusters) != 1L ||
      max_clusters < 1)
    stop("`max_clusters` must be a positive integer.", call. = FALSE)
  nsim <- as.integer(nsim)
  max_clusters <- as.integer(max_clusters)

  s <- attr(result, "sppt_summary")
  alpha <- s$alpha

  dir_obs  <- result$sig_pvalue_dir
  eligible <- !is.na(dir_obs)
  n        <- length(dir_obs)
  n_elig   <- sum(eligible)

  # Pairs of units that may join a cluster: queen neighbours, plus
  # neighbours of neighbours under the one-gap rule. Built once on the full
  # geometry, so a gap may be any unit, including one empty in both
  # patterns.
  edges <- cluster_edges(result, rule)

  # Observed clusters
  flag_obs <- eligible & dir_obs != 0L
  obs <- label_clusters(flag_obs, if (direction) dir_obs else NULL,
                        edges, n)

  # Null clusters: a share alpha of eligible units flagged at random (a
  # fixed number, or each unit independently), with a random sign when
  # clusters are directional.
  if (!is.null(seed)) set.seed(seed)
  elig_idx <- which(eligible)
  n_fp <- round(alpha * n_elig)
  null_sizes <- matrix(0L, nrow = nsim, ncol = max_clusters)
  for (b in seq_len(nsim)) {
    f <- logical(n)
    if (null_count == "fixed") {
      f[elig_idx[sample.int(n_elig, n_fp)]] <- TRUE
    } else {
      f[elig_idx[stats::runif(n_elig) < alpha]] <- TRUE
    }
    sg <- if (direction) {
      v <- integer(n)
      v[f] <- sample(c(-1L, 1L), sum(f), replace = TRUE)
      v
    } else NULL
    sz <- label_clusters(f, sg, edges, n)$sizes
    k <- min(length(sz), max_clusters)
    if (k > 0L) null_sizes[b, seq_len(k)] <- sz[seq_len(k)]
  }

  # Ceiling: smallest size s with P(k-th largest chance cluster >= s) <= level
  ceiling_k <- vapply(seq_len(max_clusters), function(k) {
    x <- null_sizes[, k]
    for (sz in seq_len(max(x) + 1L))
      if (mean(x >= sz) <= level) return(as.integer(sz))
    as.integer(max(x) + 1L)
  }, integer(1))
  ceiling_df <- data.frame(n_clusters = seq_len(max_clusters),
                           min_size = ceiling_k)

  # Observed clusters, largest first
  n_obs <- length(obs$sizes)
  cl_dir <- vapply(seq_len(n_obs), function(j) {
    d <- unique(dir_obs[which(obs$membership == j)])
    if (length(d) > 1L) "mixed" else if (d == 1L) "test > base"
    else "base > test"
  }, character(1))
  clusters <- data.frame(cluster = seq_len(n_obs), size = obs$sizes,
                         direction = cl_dir,
                         exceeds = obs$sizes >= ceiling_k[1L],
                         stringsAsFactors = FALSE)

  k_test <- seq_len(min(n_obs, max_clusters))
  tests <- data.frame(
    k = k_test,
    observed = obs$sizes[k_test],
    ceiling = ceiling_k[k_test],
    exceeds = obs$sizes[k_test] >= ceiling_k[k_test],
    p_value = vapply(k_test, function(k)
      (1 + sum(null_sizes[, k] >= obs$sizes[k])) / (nsim + 1), numeric(1)))

  structure(list(
    clusters   = clusters,
    membership = obs$membership,
    ceiling    = ceiling_df,
    tests      = tests,
    settings   = list(rule = rule, direction = direction, alpha = alpha,
                      null_count = null_count,
                      level = level, nsim = nsim, n_units = n,
                      n_eligible = n_elig, n_flagged = sum(flag_obs),
                      base_name = s$base_name, test_name = s$test_name)),
    class = "sppt_clusters")
}

#' @rdname sppt_clusters
#' @param x An object of class `sppt_clusters`.
#' @param ... Ignored.
#' @export
print.sppt_clusters <- function(x, ...) {
  st <- x$settings
  cat("\n", rule_line("="), "SPPT Cluster Test\n", rule_line("="), sep = "")
  cat(field_lines("Rule", paste0(
    if (st$rule == "one_gap") "one gap" else "contiguous",
    if (st$direction) ", within direction" else ", either direction")))
  cat(field_lines("Flagged units", sprintf("%d of %d eligible (alpha = %.3f)",
                                           st$n_flagged, st$n_eligible,
                                           st$alpha)))
  cat(field_lines("Null draws", sprintf(
    "%d (%s); ceiling at level %.2f", st$nsim,
    if (st$null_count == "fixed") "fixed share flagged" else "binomial flags",
    st$level)))
  cat(rule_line("-"))
  if (nrow(x$tests) == 0L) {
    cat("  No flagged units form a cluster.\n")
  } else {
    cat("  k   k-th largest   ceiling   p-value\n")
    for (i in seq_len(nrow(x$tests)))
      cat(sprintf("  %-3d %12d %9d   %7.4f%s\n", x$tests$k[i],
                  x$tests$observed[i], x$tests$ceiling[i],
                  x$tests$p_value[i],
                  if (x$tests$exceeds[i]) "  *" else ""))
    cat(field_lines("", paste(
      "* the k largest clusters all reach the ceiling, which chance",
      "does in no more than the stated share of draws")))
    big <- x$clusters[x$clusters$exceeds, , drop = FALSE]
    if (nrow(big) > 0L)
      cat(field_lines("Above ceiling", paste(
        sprintf("%d units (%s)", big$size, big$direction), collapse = "; ")))
  }
  cat(rule_line("="), "\n", sep = "")
  invisible(x)
}

# -- Internal helpers ---------------------------------------------------------

# Unit pairs (i < j) that may join a cluster under the given rule.
cluster_edges <- function(x, rule) {
  nb <- spdep::poly2nb(x, queen = TRUE)
  n  <- length(nb)
  from <- rep(seq_len(n), lengths(nb))
  to   <- unlist(nb, use.names = FALSE)
  ok   <- to > 0L
  A <- Matrix::sparseMatrix(i = from[ok], j = to[ok], x = 1, dims = c(n, n))
  A <- ((A + Matrix::t(A)) > 0) * 1
  R <- if (rule == "one_gap") ((A + A %*% A) > 0) * 1 else A
  R <- Matrix::triu(R, k = 1)
  R <- methods::as(R, "TsparseMatrix")
  list(i = R@i + 1L, j = R@j + 1L)
}

# Connected components of the flagged units (union-find over the edges
# whose two ends are both flagged, and of the same sign if `sgn` is given).
# Returns cluster sizes, largest first, and the cluster of each unit.
label_clusters <- function(flag, sgn, edges, n) {
  membership <- rep(NA_integer_, n)
  if (!any(flag)) return(list(sizes = integer(0), membership = membership))
  ok <- flag[edges$i] & flag[edges$j]
  if (!is.null(sgn)) ok <- ok & sgn[edges$i] == sgn[edges$j]
  ei <- edges$i[ok]
  ej <- edges$j[ok]
  parent <- seq_len(n)
  find <- function(v) {
    while (parent[v] != v) {
      parent[v] <<- parent[parent[v]]
      v <- parent[v]
    }
    v
  }
  for (k in seq_along(ei)) {
    a <- find(ei[k])
    b <- find(ej[k])
    if (a != b) parent[a] <- b
  }
  units <- which(flag)
  roots <- vapply(units, find, integer(1))
  sizes <- table(roots)
  ord   <- order(-as.integer(sizes), as.integer(names(sizes)))
  roots_sorted <- as.integer(names(sizes))[ord]
  membership[units] <- match(roots, roots_sorted)
  list(sizes = as.integer(sizes)[ord], membership = membership)
}
