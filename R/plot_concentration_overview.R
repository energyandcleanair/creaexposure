#' Plot an overview of all available concentration maps
#'
#' Scans the concentration folder structure, loads each available .tif,
#' and produces a facetted thumbnail plot.
#'
#' @param pollutants Character vector. Pollutants to include. NULL = all.
#' @param sources Character vector. Sources to include. NULL = all.
#' @param max_year_per_source Integer. Maximum number of years to show per
#'   source (evenly spaced). Use NULL to show all. Default 3.
#' @param filename Character. Output filename. If NULL, returns the plot
#'   without saving.
#' @param width Numeric. Plot width in inches.
#' @param height Numeric. Plot height in inches. If NULL, auto-computed.
#'
#' @return A ggplot object (invisibly if saved to file).
#' @export
plot_concentration_overview <- function(pollutants = NULL,
                                        sources = NULL,
                                        max_year_per_source = 3,
                                        filename = NULL,
                                        width = 16,
                                        height = NULL) {

  registry <- .concentration_sources
  if (!is.null(pollutants)) {
    registry <- registry[intersect(names(registry), tolower(pollutants))]
  }

  # Collect all available maps
  entries <- list()

  for (poll in names(registry)) {
    poll_sources <- registry[[poll]]
    if (!is.null(sources)) {
      poll_sources <- poll_sources[intersect(names(poll_sources), sources)]
    }

    for (src in names(poll_sources)) {
      src_config <- poll_sources[[src]]

      for (ver in names(src_config$versions)) {
        ver_config <- src_config$versions[[ver]]

        # Determine available years
        years <- tryCatch(
          get_concentration_available_years(poll, source = src, version = ver),
          error = function(e) integer(0)
        )

        # For sources without years (e.g. geoschem O3)
        has_years <- length(years) > 0

        if (!has_years && is.null(ver_config$year_regex) &&
            is.null(ver_config$fixed_years)) {
          # Yearless source — check if file exists
          variants <- if (!is.null(ver_config$default_variant)) {
            ver_config$default_variant
          } else {
            NULL
          }
          path <- tryCatch({
            fn <- ver_config$file_template(year = NULL, variant = variants)
            .concentration_path(poll, src, ver, fn)
          }, error = function(e) NULL)

          if (!is.null(path) && file.exists(path)) {
            label <- paste0(toupper(poll), " | ", src, " ", ver)
            if (!is.null(variants)) {
              label <- paste0(label, " (", variants, ")")
            }
            entries[[length(entries) + 1]] <- list(
              pollutant = poll, source = src, version = ver,
              year = NA, variant = variants, path = path, label = label
            )
          }

          # Check additional variants (e.g. sm8h for geoschem)
          if (!is.null(ver_config$default_variant)) {
            # Try known variants from file_template
            other_variants <- tryCatch({
              # Test sm8h for geoschem
              candidates <- c("m3m", "sm8h", "no_ssdust")
              candidates[candidates != ver_config$default_variant]
            }, error = function(e) character(0))

            for (v in other_variants) {
              path_v <- tryCatch({
                fn <- ver_config$file_template(year = NULL, variant = v)
                .concentration_path(poll, src, ver, fn)
              }, error = function(e) NULL)

              if (!is.null(path_v) && file.exists(path_v)) {
                label <- paste0(toupper(poll), " | ", src, " ", ver, " (", v, ")")
                entries[[length(entries) + 1]] <- list(
                  pollutant = poll, source = src, version = ver,
                  year = NA, variant = v, path = path_v, label = label
                )
              }
            }
          }
          next
        }

        # Subsample years if needed
        if (!is.null(max_year_per_source) && length(years) > max_year_per_source) {
          idx <- round(seq(1, length(years), length.out = max_year_per_source))
          years <- years[idx]
        }

        for (yr in years) {
          variant <- ver_config$default_variant
          fn <- ver_config$file_template(yr, variant)
          path <- .concentration_path(poll, src, ver, fn)

          if (file.exists(path)) {
            label <- paste0(toupper(poll), " | ", src, " ", ver, " ", yr)
            entries[[length(entries) + 1]] <- list(
              pollutant = poll, source = src, version = ver,
              year = yr, variant = variant, path = path, label = label
            )
          }
        }
      }
    }
  }

  if (length(entries) == 0) {
    stop("No concentration maps found. Run the migration script first.")
  }

  message(glue::glue("Found {length(entries)} concentration maps. Loading..."))

  # Load rasters and convert to data frames for ggplot
  all_dfs <- list()

  for (entry in entries) {
    r <- tryCatch({
      r <- terra::rast(entry$path)
      # Downsample for thumbnails (max ~200 cells per side)
      max_cells <- 200
      if (max(dim(r)[1:2]) > max_cells) {
        fact <- ceiling(max(dim(r)[1:2]) / max_cells)
        r <- terra::aggregate(r, fact = fact, fun = "mean", na.rm = TRUE)
      }
      r
    }, error = function(e) {
      message("  Failed to load: ", entry$path, " — ", conditionMessage(e))
      NULL
    })

    if (is.null(r)) next

    df <- as.data.frame(r, xy = TRUE)
    names(df)[3] <- "value"
    df <- df[!is.na(df$value), ]
    df$label <- entry$label
    df$pollutant <- toupper(entry$pollutant)
    all_dfs[[length(all_dfs) + 1]] <- df
  }

  plot_df <- do.call(rbind, all_dfs)

  # Normalize values per panel to 0–1 so each thumbnail uses full color range
  # Append min/max to label for reference
  plot_df <- do.call(rbind, lapply(split(plot_df, plot_df$label), function(d) {
    q01 <- stats::quantile(d$value, 0.01, na.rm = TRUE)
    q99 <- stats::quantile(d$value, 0.99, na.rm = TRUE)
    rng <- q99 - q01
    if (rng > 0) {
      d$value_norm <- pmin(pmax((d$value - q01) / rng, 0), 1)
    } else {
      d$value_norm <- 0.5
    }
    d$label_full <- paste0(d$label, "\n[", round(q01, 1), " - ", round(q99, 1), "]")
    d
  }))

  n_panels <- length(unique(plot_df$label_full))
  ncol <- min(n_panels, 4)
  nrow <- ceiling(n_panels / ncol)
  if (is.null(height)) height <- 3.5 * nrow

  plt <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y, fill = value_norm)) +
    ggplot2::geom_raster() +
    ggplot2::facet_wrap(~label_full, scales = "free", ncol = ncol) +
    ggplot2::scale_fill_viridis_c(option = "rocket", direction = -1,
                                   name = "Normalized\n(range per panel)",
                                   breaks = c(0, 0.5, 1),
                                   labels = c("low", "mid", "high")) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 7, face = "bold"),
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(2, "cm"),
      legend.key.height = ggplot2::unit(0.3, "cm"),
      panel.spacing = ggplot2::unit(0.3, "lines")
    )

  if (!is.null(filename)) {
    message(glue::glue("Saving to: {filename}"))
    ggplot2::ggsave(filename, plt, width = width, height = height, dpi = 150,
                     limitsize = FALSE)
    return(invisible(plt))
  }

  return(plt)
}
