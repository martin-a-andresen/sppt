#' Redraw SPPT maps with a different layout
#'
#' @description
#' [sppt()] and [sppt_longitudinal()] draw their maps with a fixed layout
#' that suits Vancouver: north arrow top right, scale bar along the bottom,
#' and legends placed where the city's shape leaves room. `sppt_map()`
#' rebuilds any of those maps from a result object with the legend, north
#' arrow and scale bar placed wherever suits your study area.
#'
#' Every map is an ordinary `ggplot` object, so anything not covered by the
#' arguments here can be changed by adding layers or a theme, for example
#' `sppt_map(res, "pvalue") + ggplot2::theme(legend.text =
#' ggplot2::element_text(size = 12))`.
#'
#' @param result A result returned by [sppt()] or [sppt_longitudinal()].
#' @param which Character vector of map names. `NULL` (default) returns every
#'   map available for the result. For [sppt()] results: `"difference"`,
#'   `"pvalue"`, `"tiers"` (confidence tiers; needs the arrangement test) and
#'   `"ci"`. For [sppt_longitudinal()] results:
#'   `"overlap"`, `"s_index_average"`, `"sustained_zero"`, `"trend_all"`,
#'   `"trend_linear_sig"`, `"linear_r2"`, `"linear_slope"`, `"quad_r2"`,
#'   `"quad_trend"`, and the p-value versions `"p_overlap"`,
#'   `"p_sig_proportion"`, `"p_sustained_zero"`, `"p_trend_all"`,
#'   `"p_trend_linear_sig"`, `"p_linear_r2"`, `"p_linear_slope"`,
#'   `"p_quad_r2"`, `"p_trend"`. The quadratic maps exist only when
#'   `trend_type = "quadratic"`.
#' @param legend_position Where to put the legend: `"right"`, `"left"`,
#'   `"top"`, `"bottom"`, `"none"`, or a numeric vector `c(x, y)` giving the
#'   top-left corner of the legend in panel coordinates (0 to 1; negative `x`
#'   pushes it into the left margin). `NULL` (default) keeps each map's
#'   built-in placement.
#' @param north_arrow Location of the north arrow: `"tr"`, `"tl"`, `"br"`,
#'   `"bl"`, or `"none"` to omit it. `NULL` (default) uses `"tr"`.
#' @param scale_bar Location of the scale bar, with the same options as
#'   `north_arrow`. `NULL` (default) uses `"bl"` for [sppt()] maps and `"br"`
#'   for [sppt_longitudinal()] maps.
#' @param crs CRS to project the maps to, passed to [sf::st_transform()].
#'   `NULL` keeps the CRS of the result.
#' @param title Logical. Keep the map titles. `FALSE` removes them, as
#'   journals usually require. Default `TRUE`.
#'
#' @return A `ggplot` object when `which` names a single map, otherwise a
#'   named list of `ggplot` objects.
#'
#' @seealso [sppt()], [sppt_longitudinal()]
#' @examples
#' res <- sppt(vancouver_ct, group_col = "id",
#'             count_col = c("BNER2003", "BNER2025"),
#'             B = 199, seed = 42, verbose = FALSE, create_maps = FALSE)
#' p <- sppt_map(res, "pvalue", legend_position = "right",
#'               north_arrow = "tl", scale_bar = "br")
#' p
#'
#' @export
sppt_map <- function(result, which = NULL, legend_position = NULL,
                     north_arrow = NULL, scale_bar = NULL, crs = NULL,
                     title = TRUE) {
  check_location(north_arrow, "north_arrow")
  check_location(scale_bar, "scale_bar")
  if (!is.null(legend_position) && !is.numeric(legend_position) &&
      !legend_position %in% c("right", "left", "top", "bottom", "none"))
    stop("`legend_position` must be \"right\", \"left\", \"top\", ",
         "\"bottom\", \"none\" or a numeric vector c(x, y).", call. = FALSE)

  layout <- list(legend_position = legend_position,
                 north_arrow = north_arrow, scale_bar = scale_bar,
                 title = title)

  maps <- if (inherits(result, "sppt_result")) {
    build_maps_sppt(result, map_crs = crs, layout = layout)
  } else if (!is.null(attr(result, "sppt_long_meta"))) {
    build_maps_long(result, map_crs = crs, layout = layout)
  } else {
    stop("`result` must be returned by sppt() or sppt_longitudinal().",
         call. = FALSE)
  }

  if (is.null(which)) return(maps)
  unknown <- setdiff(which, names(maps))
  if (length(unknown) > 0L)
    stop("Map(s) not available for this result: ",
         paste(unknown, collapse = ", "), ". Available: ",
         paste(names(maps), collapse = ", "), ".", call. = FALSE)
  if (length(which) == 1L) maps[[which]] else maps[which]
}

#' Plot an SPPT result
#'
#' Draws one of the maps stored with a [sppt()] or [sppt_longitudinal()]
#' result. A thin wrapper around [sppt_map()], so every layout argument of
#' that function can be passed here.
#'
#' @param x A result returned by [sppt()] or [sppt_longitudinal()].
#' @param which Name of the map to draw. Defaults to `"pvalue"` for [sppt()]
#'   results and `"p_trend"` for [sppt_longitudinal()] results (`"trend_all"`
#'   when p-values were not computed). See [sppt_map()] for every name.
#' @param ... Layout arguments passed to [sppt_map()], such as
#'   `legend_position`, `north_arrow`, `scale_bar`, `crs` and `title`.
#' @return The `ggplot` object, invisibly, after printing it.
#' @seealso [sppt_map()]
#' @examples
#' res <- sppt(vancouver_ct, group_col = "id",
#'             count_col = c("BNER2003", "BNER2025"),
#'             B = 199, seed = 42, verbose = FALSE, create_maps = FALSE)
#' plot(res)
#' plot(res, "tiers", legend_position = "right")
#' @export
plot.sppt_result <- function(x, which = "pvalue", ...) {
  p <- sppt_map(x, which = which[1L], ...)
  print(p)
  invisible(p)
}

#' @rdname plot.sppt_result
#' @export
plot.sppt_long_result <- function(x, which = NULL, ...) {
  if (is.null(which)) {
    meta <- attr(x, "sppt_long_meta")
    which <- if (isTRUE(meta$calculate_pvalues)) "p_trend" else "trend_all"
  }
  p <- sppt_map(x, which = which[1L], ...)
  print(p)
  invisible(p)
}

# -- Layout helpers -----------------------------------------------------------

check_location <- function(x, arg) {
  if (!is.null(x) && !(is.character(x) && length(x) == 1L &&
                       x %in% c("tr", "tl", "br", "bl", "none")))
    stop("`", arg, "` must be one of \"tr\", \"tl\", \"br\", \"bl\" or ",
         "\"none\".", call. = FALSE)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Legend placement. A numeric position puts the legend inside the panel with
# its top-left corner at that point; a negative x needs extra left margin.
legend_theme <- function(position) {
  if (is.numeric(position)) {
    th <- ggplot2::theme(
      legend.position        = "inside",
      legend.position.inside = position[1:2],
      legend.justification   = c(0, 1),
      legend.background      = ggplot2::element_rect(fill = "white",
                                                     colour = "black")
    )
    if (position[1] < 0)
      th <- th + ggplot2::theme(
        plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 70,
                                      unit = "pt"))
    th
  } else {
    ggplot2::theme(legend.position = position)
  }
}

north_arrow_layer <- function(location, minimal = FALSE) {
  if (identical(location, "none")) return(NULL)
  if (minimal) {
    ggspatial::annotation_north_arrow(
      location = location, which_north = "true",
      pad_x = ggplot2::unit(0.5, "cm"), pad_y = ggplot2::unit(0.2, "cm"),
      height = ggplot2::unit(1, "cm"), width = ggplot2::unit(0.5, "cm"),
      style = ggspatial::north_arrow_minimal())
  } else {
    ggspatial::annotation_north_arrow(
      location = location, which_north = "true",
      height = ggplot2::unit(1, "cm"), width = ggplot2::unit(0.5, "cm"))
  }
}

scale_bar_layer <- function(location, width_hint, pad = NULL) {
  if (identical(location, "none")) return(NULL)
  if (is.null(pad)) {
    ggspatial::annotation_scale(location = location, width_hint = width_hint)
  } else {
    ggspatial::annotation_scale(location = location, width_hint = width_hint,
                                pad_x = ggplot2::unit(pad[1], "cm"),
                                pad_y = ggplot2::unit(pad[2], "cm"))
  }
}

drop_title <- function(p, keep) {
  if (isTRUE(keep)) p else p + ggplot2::labs(title = NULL)
}

strip_sppt_class <- function(x) {
  class(x) <- setdiff(class(x), c("sppt_result", "sppt_long_result"))
  x
}

# Dissolve the outlined units of one direction into a single boundary, so a
# cluster is drawn by its outer edge rather than by every unit border.
dissolve_units <- function(x, idx) {
  if (!length(idx)) return(NULL)
  g <- tryCatch(
    sf::st_union(sf::st_make_valid(sf::st_geometry(x[idx, ]))),
    error = function(e) sf::st_geometry(x[idx, ]))
  sf::st_sf(geometry = g)
}

# -- Maps for sppt() ----------------------------------------------------------

build_maps_sppt <- function(result, map_crs = NULL, layout = list()) {
  s <- attr(result, "sppt_summary")
  base_name <- s$base_name
  test_name <- s$test_name

  rm <- strip_sppt_class(result)
  if (!is.null(map_crs)) rm <- sf::st_transform(rm, crs = map_crs)

  # Legend sits below the map on two rows by default. Placing it inside the
  # panel depends on the shape of the study area, and a single wide row of
  # keys stretches the plot device and squeezes the map.
  base_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.direction = "horizontal",
      legend.box       = "horizontal",
      legend.key.size  = ggplot2::unit(0.4, "cm"),
      legend.margin    = ggplot2::margin(t = 2, b = 2),
      legend.title     = ggplot2::element_text(face = "bold", size = 9),
      legend.text      = ggplot2::element_text(size = 8),
      panel.grid       = ggplot2::element_blank(),
      axis.text        = ggplot2::element_blank(),
      axis.ticks       = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(size = 13, face = "bold",
                                               hjust = 0.5),
      plot.subtitle    = ggplot2::element_text(size = 8, hjust = 0.5,
                                               lineheight = 1.1)
    ) +
    legend_theme(layout$legend_position %||% "bottom")
  legend_2row <- ggplot2::guides(
    fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE,
                                 title.position = "top"))
  carto <- list(north_arrow_layer(layout$north_arrow %||% "tr"),
                scale_bar_layer(layout$scale_bar %||% "bl", 0.35))

  vs <- paste0("Base: ", base_name, "  |  Test: ", test_name)
  lv <- c(paste0(test_name, " > ", base_name),
          paste0(base_name, " > ", test_name),
          "Not significant", "Excluded (0,0)")
  sig_colors <- stats::setNames(c("black", "gray60", "white", "gray90"), lv)
  label_dir <- function(d) factor(
    ifelse(is.na(d), lv[4],
           ifelse(d == 1L, lv[1], ifelse(d == -1L, lv[2], lv[3]))),
    levels = lv)

  maps <- list()

  # Proportional difference
  maps$difference <- ggplot2::ggplot(rm) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[["difference"]]),
                     colour = "gray50", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      name = if (isTRUE(s$use_percentages)) "Difference\n(%)" else "Difference",
      na.value = "gray80") +
    carto + base_theme +
    ggplot2::labs(title = "Proportional Difference", subtitle = vs)

  # P-value significance, with a self-contained subtitle
  if (!is.null(s$alpha) && "sig_pvalue_dir" %in% names(rm)) {
    rm$sig_label <- label_dir(rm$sig_pvalue_dir)
    arr <- s$arrangement
    sub <- vs
    if (!is.na(s$s_index_pv)) {
      line2 <- sprintf("Robust S = %.3f", s$s_index_pv)
      if (!is.null(s$sameness) && !is.na(s$sameness$p_value))
        line2 <- paste0(line2, ", p ", fmt_p(s$sameness$p_value))
      if (!is.null(arr))
        line2 <- paste0(line2, sprintf("   |   Moran's I = %+.3f, p %s",
                                       arr$moran_I, fmt_p(arr$moran_p)))
      sub <- paste0(sub, "\n", line2)
      if (!is.null(arr) && arr$n_hh > 0L)
        sub <- paste0(sub, sprintf("\nSolid outline: %d units, %s > %s",
                                   arr$n_hh, test_name, base_name))
      if (!is.null(arr) && arr$n_ll > 0L)
        sub <- paste0(sub, sprintf("\nDashed outline: %d units, %s > %s",
                                   arr$n_ll, base_name, test_name))
    }

    # Unit borders are light and thin so the cluster boundaries read at fine
    # scales, where many small units would otherwise form a mesh.
    unit_col <- "gray70"
    unit_lwd <- 0.1
    cl_up <- if ("moran_hh" %in% names(rm))
      dissolve_units(rm, which(rm$moran_hh)) else NULL
    cl_dn <- if ("moran_ll" %in% names(rm))
      dissolve_units(rm, which(rm$moran_ll)) else NULL
    add_clusters <- function(p) {
      # Solid = test > base, dashed = base > test
      if (!is.null(cl_up))
        p <- p + ggplot2::geom_sf(data = cl_up, inherit.aes = FALSE,
                                  fill = NA, colour = "red",
                                  linewidth = 0.8, linetype = "solid")
      if (!is.null(cl_dn))
        p <- p + ggplot2::geom_sf(data = cl_dn, inherit.aes = FALSE,
                                  fill = NA, colour = "red",
                                  linewidth = 0.8, linetype = "22")
      p
    }

    p <- ggplot2::ggplot(rm) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data[["sig_label"]]),
                       colour = unit_col, linewidth = unit_lwd) +
      ggplot2::scale_fill_manual(values = sig_colors, name = "Pattern",
                                 na.value = "gray80", drop = FALSE)
    maps$pvalue <- add_clusters(p) + carto + base_theme + legend_2row +
      ggplot2::labs(title = sprintf("P-value Significance  (alpha = %.3f)",
                                    s$alpha),
                    subtitle = sub)

    # Confidence tiers: fill carries tier and direction, and each legend key
    # carries its count, so an empty category reads as (0).
    if ("confidence_tier" %in% names(rm)) {
      up <- lv[1]
      dn <- lv[2]
      tl <- c(paste("Core:", up), paste("Edge:", up), paste("Outside:", up),
              "Not significant",
              paste("Core:", dn), paste("Edge:", dn), paste("Outside:", dn),
              lv[4])
      dr  <- rm$sig_pvalue_dir
      lab <- ifelse(is.na(dr), lv[4],
             ifelse(dr == 0L, "Not significant",
                    paste0(tools::toTitleCase(rm$confidence_tier), ": ",
                           ifelse(dr == 1L, up, dn))))
      rm$tier_label <- factor(lab, levels = tl)
      tier_n <- table(rm$tier_label)
      tier_keys <- stats::setNames(
        paste0(tl, " (", as.integer(tier_n[tl]), ")"), tl)
      tier_colors <- stats::setNames(
        c("#b2182b", "#ef8a62", "#fddbc7", "white",
          "#2166ac", "#67a9cf", "#d1e5f0", "gray90"), tl)
      arr <- s$arrangement
      pt <- ggplot2::ggplot(rm) +
        ggplot2::geom_sf(ggplot2::aes(fill = .data[["tier_label"]]),
                         colour = unit_col, linewidth = unit_lwd) +
        ggplot2::scale_fill_manual(values = tier_colors,
                                   labels = tier_keys,
                                   name = "Flagged units",
                                   na.value = "gray80", drop = FALSE)
      maps$tiers <- add_clusters(pt) + carto + base_theme +
        ggplot2::guides(fill = ggplot2::guide_legend(
          ncol = 2, byrow = FALSE, title.position = "top")) +
        ggplot2::labs(
          title = sprintf("Confidence Tiers  (alpha = %.3f)", s$alpha),
          subtitle = paste0(
            vs,
            sprintf("\nCore %d, edge %d, outside %d", arr$n_core,
                    arr$n_edge, arr$n_outside),
            "\nCore = inside a cluster; edge = borders a core unit, same sign"))
    }
  }

  # CI significance
  if (!is.null(s$conf_level) && "sig_ci_dir" %in% names(rm)) {
    rm$ci_label <- label_dir(rm$sig_ci_dir)
    maps$ci <- ggplot2::ggplot(rm) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data[["ci_label"]]),
                       colour = "gray50", linewidth = 0.2) +
      ggplot2::scale_fill_manual(values = sig_colors, name = "Pattern",
                                 na.value = "gray80", drop = FALSE) +
      carto + base_theme + legend_2row +
      ggplot2::labs(title = sprintf("CI Significance  (conf_level = %.1f%%)",
                                    s$conf_level * 100),
                    subtitle = vs)
  }

  # The coordinate system goes last: a later geom_sf() layer would otherwise
  # replace it with the default and bring the graticule back.
  maps <- lapply(maps, function(p) p + ggplot2::coord_sf(datum = NA))
  lapply(maps, drop_title, keep = layout$title %||% TRUE)
}

# -- Maps for sppt_longitudinal() ---------------------------------------------

# File names used when exporting, in figure order.
long_map_files <- c(
  overlap            = "map1_complete_overlap",
  s_index_average    = "map2_s_index_average",
  sustained_zero     = "map3_sustained_zeros",
  trend_all          = "map4_all_trends",
  trend_linear_sig   = "map5_significant_linear_trends",
  linear_r2          = "map6_linear_r2",
  linear_slope       = "map7_linear_slope",
  quad_r2            = "map8_quadratic_r2",
  quad_trend         = "map9_quadratic_trends",
  p_overlap          = "pmap1_complete_nonsig",
  p_sig_proportion   = "pmap2_sig_proportion",
  p_sustained_zero   = "pmap3_sustained_zeros_p",
  p_trend_all        = "pmap4_all_trends_p",
  p_trend_linear_sig = "pmap5_significant_linear_trends_p",
  p_linear_r2        = "pmap6_linear_r2_p",
  p_linear_slope     = "pmap7_linear_slope_p",
  p_quad_r2          = "pmap8_quadratic_r2_p",
  p_trend            = "pmap9_trends_p"
)

build_maps_long <- function(result, map_crs = NULL, layout = list()) {
  meta <- attr(result, "sppt_long_meta")
  col  <- function(nm) paste0(meta$prefix, nm)
  quad <- identical(meta$trend_type, "quadratic")

  rm <- result
  if (!is.null(map_crs)) rm <- sf::st_transform(rm, crs = map_crs)

  # Built-in placement: categorical legends inside the panel at the top left,
  # continuous colour bars on the right. Suits Vancouver.
  left  <- layout$legend_position %||% c(-0.10, 0.98)
  right <- layout$legend_position %||% "right"

  base_theme <- ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "white", colour = NA),
    plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
    panel.grid       = ggplot2::element_blank(),
    panel.border     = ggplot2::element_blank(),
    axis.text        = ggplot2::element_blank(),
    axis.ticks       = ggplot2::element_blank(),
    axis.title       = ggplot2::element_blank(),
    legend.text      = ggplot2::element_text(size = 10),
    plot.title       = ggplot2::element_text(size = 12, hjust = 0.5,
                                             face = "bold")
  )
  carto <- list(north_arrow_layer(layout$north_arrow %||% "tr", minimal = TRUE),
                scale_bar_layer(layout$scale_bar %||% "br", 0.15,
                                pad = c(0.65, 0.2)))
  finish <- function(p, pos) p + carto + base_theme + legend_theme(pos)

  bar <- ggplot2::guides(fill = ggplot2::guide_colourbar(
    barheight = ggplot2::unit(4, "cm"), barwidth = ggplot2::unit(0.5, "cm")))
  grey_r2 <- ggplot2::scale_fill_gradientn(
    colours = c("white", "gray20", "gray40", "gray60", "gray80", "black"),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1),
    name = NULL, na.value = "white")
  slope_scale <- function(lim) ggplot2::scale_fill_gradientn(
    colours = c("#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0",
                "white", "#fddbc7", "#f4a582", "#d6604d", "#b2182b",
                "#67001f"),
    values = scales::rescale(c(-1, -0.8, -0.6, -0.4, -0.2, 0,
                               0.2, 0.4, 0.6, 0.8, 1)),
    limits = c(-lim, lim), name = NULL, na.value = "white")
  slope_limit <- function(x) {
    m <- suppressWarnings(max(abs(x), na.rm = TRUE))
    if (!is.finite(m) || m == 0) 1 else m
  }
  quad_scale <- ggplot2::scale_fill_manual(
    values = c("-3" = "#053061", "-2" = "#2166ac", "-1" = "#4393c3",
               "0" = "white", "1" = "#d6604d", "2" = "#b2182b",
               "3" = "#67001f"),
    labels = c("-3" = "Sig Dec + Neg Quad", "-2" = "Sig Dec + Pos Quad",
               "-1" = "Sig Decreasing", "0" = "Not significant",
               "1" = "Sig Increasing", "2" = "Sig Inc + Neg Quad",
               "3" = "Sig Inc + Pos Quad"),
    name = NULL, drop = TRUE)
  fill_map <- function(var, discrete = TRUE, levels = NULL) {
    v <- rm[[var]]
    if (discrete) v <- if (is.null(levels)) factor(v) else factor(v, levels)
    rm$.fill <- v
    ggplot2::ggplot() +
      ggplot2::geom_sf(data = rm, ggplot2::aes(fill = .data[[".fill"]]),
                       colour = "gray50", linewidth = 0.3)
  }

  maps <- list()

  # -- Envelope maps --
  maps$overlap <- finish(
    fill_map(col("intervals_overlap")) +
      ggplot2::scale_fill_manual(
        values = c("0" = "black", "1" = "white"),
        labels = c("0" = "No complete overlap", "1" = "Complete overlap"),
        name = NULL) +
      ggplot2::labs(title = "SPPT Output: No Overlap (Instability)"), left)

  maps$s_index_average <- finish(
    fill_map(col("max_overlap_proportion"), discrete = FALSE) +
      ggplot2::scale_fill_gradientn(
        colours = c("black", "gray30", "gray50", "gray70", "gray90", "white"),
        breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1),
        name = NULL) +
      bar + ggplot2::labs(title = "S-Index Average: Degree of Instability"),
    right)

  maps$sustained_zero <- finish(
    fill_map(col("sustained_zero")) +
      ggplot2::scale_fill_manual(
        values = c("0" = "white", "1" = "black"),
        labels = c("0" = "Active or always zero",
                   "1" = "Sustained zero decline"),
        name = NULL) +
      ggplot2::labs(title = "Sustained Zero Events (Decline Only)"), left)

  maps$trend_all <- finish(
    fill_map(col("Trend_All")) +
      ggplot2::scale_fill_manual(
        values = c("-1" = "gray80", "0" = "white", "1" = "black"),
        labels = c("-1" = "Decreasing", "0" = "No trend",
                   "1" = "Increasing"),
        name = NULL) +
      ggplot2::labs(title = "Linear Trend Direction (All Observations)"),
    left)

  maps$trend_linear_sig <- finish(
    fill_map(col("Trend_Linear_Sig")) +
      ggplot2::scale_fill_manual(
        values = c("-1" = "gray70", "0" = "white", "1" = "black"),
        labels = c("-1" = "Sig Decreasing", "0" = "Not significant",
                   "1" = "Sig Increasing"),
        name = NULL) +
      ggplot2::labs(title = "Significant Linear Trends"), left)

  maps$linear_r2 <- finish(
    fill_map(col("Linear_R2"), discrete = FALSE) + grey_r2 + bar +
      ggplot2::labs(title = "Linear Trend R-squared"), right)

  maps$linear_slope <- finish(
    fill_map(col("Linear_bt_sig"), discrete = FALSE) +
      slope_scale(slope_limit(rm[[col("Linear_bt_sig")]])) + bar +
      ggplot2::labs(title = "Linear Trend Slope (Significant Only)"), right)

  if (quad) {
    maps$quad_r2 <- finish(
      fill_map(col("Quad_R2"), discrete = FALSE) + grey_r2 + bar +
        ggplot2::labs(title = "Quadratic Trend R-squared"), right)
    maps$quad_trend <- finish(
      fill_map(col("Trend")) + quad_scale +
        ggplot2::labs(title = "Quadratic Trend Analysis: Significant Trends"),
      right)
  }

  # -- P-value maps --
  if (isTRUE(meta$calculate_pvalues)) {
    pv_sub <- paste0("\n[alpha=", meta$alpha, "]")

    maps$p_overlap <- finish(
      fill_map(col("intervals_overlap_p")) +
        ggplot2::scale_fill_manual(
          values = c("0" = "black", "1" = "white"),
          labels = c("0" = "At least one sig period",
                     "1" = "All periods non-significant"),
          name = NULL) +
        ggplot2::labs(title = paste0("P-Value: All Periods Non-Significant",
                                     pv_sub)), left)

    maps$p_sig_proportion <- finish(
      fill_map(col("prop_sig_p"), discrete = FALSE) +
        ggplot2::scale_fill_gradientn(
          colours = c("white", "gray90", "gray70", "gray50", "gray30",
                      "black"),
          breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1),
          name = NULL) +
        bar + ggplot2::labs(title = paste0(
          "P-Value: Proportion of Significant Periods", pv_sub)), right)

    maps$p_sustained_zero <- finish(
      fill_map(col("sustained_zero_p")) +
        ggplot2::scale_fill_manual(
          values = c("0" = "white", "1" = "black"),
          labels = c("0" = "Active or always zero",
                     "1" = "Sustained zero decline"),
          name = NULL) +
        ggplot2::labs(title = paste0(
          "P-Value: Sustained Zero Events (Count-Based)", pv_sub)), left)

    maps$p_trend_all <- finish(
      fill_map(col("Trend_p_All"), levels = c(-1L, 0L, 1L)) +
        ggplot2::scale_fill_manual(
          values = c("-1" = "gray80", "0" = "white", "1" = "black"),
          labels = c("-1" = "Decreasing", "0" = "No trend",
                     "1" = "Increasing"),
          name = NULL, drop = TRUE) +
        ggplot2::labs(title = paste0(
          "P-Value: Linear Trend Direction (All Obs)", pv_sub)), left)

    maps$p_trend_linear_sig <- finish(
      fill_map(col("Trend_p_Linear_Sig"), levels = c(-1L, 0L, 1L)) +
        ggplot2::scale_fill_manual(
          values = c("-1" = "gray70", "0" = "white", "1" = "black"),
          labels = c("-1" = "Sig Decreasing", "0" = "Not significant",
                     "1" = "Sig Increasing"),
          name = NULL, drop = TRUE) +
        ggplot2::labs(title = paste0("P-Value: Significant Linear Trends",
                                     pv_sub)), left)

    maps$p_linear_r2 <- finish(
      fill_map(col("Linear_p_R2"), discrete = FALSE) + grey_r2 + bar +
        ggplot2::labs(title = paste0("P-Value: Linear Trend R-squared",
                                     pv_sub)), right)

    maps$p_linear_slope <- finish(
      fill_map(col("Linear_p_bt_sig"), discrete = FALSE) +
        slope_scale(slope_limit(rm[[col("Linear_p_bt_sig")]])) + bar +
        ggplot2::labs(title = paste0(
          "P-Value: Linear Trend Slope (Significant Only)", pv_sub)), right)

    if (quad) {
      maps$p_quad_r2 <- finish(
        fill_map(col("Quad_p_R2"), discrete = FALSE) + grey_r2 + bar +
          ggplot2::labs(title = paste0("P-Value: Quadratic Trend R-squared",
                                       pv_sub)), right)
      maps$p_trend <- finish(
        fill_map(col("Trend_p")) + quad_scale +
          ggplot2::labs(title = paste0(
            "P-Value: Quadratic Trend Analysis: Significant Trends",
            pv_sub)), right)
    } else {
      maps$p_trend <- finish(
        fill_map(col("Trend_p")) +
          ggplot2::scale_fill_manual(
            values = c("-2" = "gray50", "-1" = "gray80", "0" = "white",
                       "1" = "gray30", "2" = "black"),
            labels = c("-2" = "Sig Decreasing", "-1" = "Nonsig Decreasing",
                       "0" = "Not significant", "1" = "Nonsig Increasing",
                       "2" = "Sig Increasing"),
            name = NULL, drop = TRUE) +
          ggplot2::labs(title = paste0("P-Value: Trend Classification (Linear)",
                                       pv_sub)), left)
    }
  }

  lapply(maps, drop_title, keep = layout$title %||% TRUE)
}
