# ------------------------------------------------------------
# Phyloseq helper
# ------------------------------------------------------------

#' Filter taxa in a phyloseq object by prevalence
#'
#' Keeps only taxa present in at least a specified fraction of samples.
#' Prevalence is defined as the fraction of samples with abundance strictly
#' greater than zero.
#'
#' @param ps A \code{phyloseq} object.
#' @param prevalence_threshold Numeric scalar between 0 and 1 giving the minimum
#'   prevalence required for a taxon to be retained.
#'
#' @return A filtered \code{phyloseq} object containing only taxa with
#'   prevalence greater than or equal to \code{prevalence_threshold}.
#'
#' @examples
#' # ps_filt <- filter_taxa_by_prevalence(ps, prevalence_threshold = 0.1)
filter_taxa_by_prevalence <- function(ps, prevalence_threshold) {
  # Check that the input is a phyloseq object.
  if (!inherits(ps, "phyloseq")) {
    stop("ps must be a phyloseq object.", call. = FALSE)
  }
  
  # Check that the threshold is a single numeric value in [0, 1].
  if (!is.numeric(prevalence_threshold) || length(prevalence_threshold) != 1 ||
      prevalence_threshold < 0 || prevalence_threshold > 1) {
    stop("prevalence_threshold must be a numeric value between 0 and 1.", call. = FALSE)
  }
  
  # Extract the OTU table as a matrix and ensure taxa are rows.
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) {
    otu <- t(otu)
  }
  
  # Compute prevalence as the fraction of samples with non-zero abundance.
  prevalence <- rowSums(otu > 0, na.rm = TRUE) / ncol(otu)
  
  # Identify taxa that pass the prevalence threshold.
  keep_taxa <- names(prevalence)[prevalence >= prevalence_threshold]
  
  # Subset the phyloseq object to retained taxa only.
  ps_filtered <- prune_taxa(keep_taxa, ps)
  
  ps_filtered
}


#' Add a small in-panel legend to a facet-grid upper-bound plot
#'
#' Returns ggplot layers that can be added with `+` to an existing plot.
#' The legend is drawn only inside the facet identified by `legend_test`
#' and `legend_stratum`.
#'
#' @param legend_test Character scalar. Method/facet row where the legend should appear.
#' @param legend_stratum Character scalar. Abundance stratum/facet column where the legend should appear.
#' @param x Numeric vector of length 2 giving the start and end x positions of legend lines.
#' @param y_top Numeric scalar. y position for the first legend line.
#' @param y_gap Numeric scalar. Vertical gap between the first and second legend lines.
#' @param text_gap Numeric scalar. Horizontal gap between line end and label text.
#' @param threshold_label Character scalar. Label for the blue threshold line.
#' @param observed_label Character scalar. Label for the black observed-envelope line.
#' @param threshold_color Character scalar. Color for the threshold line.
#' @param observed_color Character scalar. Color for the observed-envelope line.
#' @param threshold_linetype Character scalar/integer. Linetype for the threshold line.
#' @param observed_linetype Character scalar/integer. Linetype for the observed-envelope line.
#' @param line_width Numeric scalar. Line width for legend segments.
#' @param text_size Numeric scalar. Text size for legend labels.
#'
#' @return A list of ggplot layers that can be added with `+`.
#'
add_inpanel_threshold_legend <- function(
    legend_test,
    legend_stratum,
    test_levels = NULL,
    x_pos = c(0.12, 0.36),
    y_top = 0.78,
    y_gap = 0.10,
    text_gap = 0.04,
    threshold_label = "80% detection barrier",
    observed_label = "Observed effect q95",
    threshold_color = "#2C6BEA",
    observed_color = "black",
    threshold_linetype = "31",
    observed_linetype = "solid",
    line_width = 0.9,
    text_size = 2.5
) {
  stopifnot(length(x_pos) == 2)
  
  seg_df <- tibble::tibble(
    test = factor(c(legend_test, legend_test), levels = test_levels),
    col_label = "Abundance stratum",
    abund_stratum = factor(c(legend_stratum, legend_stratum),
                           levels = c("Low", "Medium", "High")),
    x = c(x_pos[1], x_pos[1]),
    xend = c(x_pos[2], x_pos[2]),
    y = c(y_top, y_top - y_gap),
    yend = c(y_top, y_top - y_gap),
    line_group = c("threshold", "observed")
  )
  
  txt_df <- tibble::tibble(
    test = factor(c(legend_test, legend_test), levels = test_levels),
    col_label = "Abundance stratum",
    abund_stratum = factor(c(legend_stratum, legend_stratum),
                           levels = c("Low", "Medium", "High")),
    x = c(x_pos[2] + text_gap, x_pos[2] + text_gap),
    y = c(y_top, y_top - y_gap),
    label = c(threshold_label, observed_label)
  )
  
  list(
    ggplot2::geom_segment(
      data = seg_df %>% dplyr::filter(line_group == "threshold"),
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = threshold_color,
      linetype = threshold_linetype,
      linewidth = line_width
    ),
    ggplot2::geom_segment(
      data = seg_df %>% dplyr::filter(line_group == "observed"),
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = observed_color,
      linetype = observed_linetype,
      linewidth = line_width
    ),
    ggplot2::geom_label(
      data = txt_df,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0.5,
      size = text_size,
      colour = "black",
      fill = "white",
      linewidth = 0.2,
      label.r = grid::unit(0.08, "lines"),
      label.padding = grid::unit(0.12, "lines")
    )
  )
}

#' Add a small in-panel legend to a facet-grid p-value plot
#'
#' Returns ggplot layers that can be added with `+` to an existing plot.
#' The legend is drawn only inside the facet identified by `legend_test`
#' and `legend_stratum`.
#'
#' @param legend_test Character scalar. Method/facet row where the legend should appear.
#' @param legend_stratum Character scalar. Abundance stratum/facet column where the legend should appear.
#' @param x Numeric vector of length 2 giving the start and end x positions of legend lines.
#' @param y_top Numeric scalar. y position for the first legend line.
#' @param y_gap Numeric scalar. Vertical gap between the first and second legend lines.
#' @param text_gap Numeric scalar. Horizontal gap between line end and label text.
#' @param threshold_label Character scalar. Label for the blue threshold line.
#' @param observed_label Character scalar. Label for the black observed-envelope line.
#' @param threshold_color Character scalar. Color for the threshold line.
#' @param observed_color Character scalar. Color for the observed-envelope line.
#' @param threshold_linetype Character scalar/integer. Linetype for the threshold line.
#' @param observed_linetype Character scalar/integer. Linetype for the observed-envelope line.
#' @param line_width Numeric scalar. Line width for legend segments.
#' @param text_size Numeric scalar. Text size for legend labels.
#'
#' @return A list of ggplot layers that can be added with `+`.
#'
add_inpanel_pval_legend <- function(
    legend_test,
    legend_stratum,
    test_levels = NULL,
    x_pos = c(0.12, 0.36),
    y_top = 0.78,
    y_gap = 0.10,
    text_gap = 0.04,
    threshold_label = "80% detection barrier",
    observed_label = "Observed effect q95",
    threshold_color = "#2C6BEA",
    observed_color = "black",
    threshold_linetype = "31",
    observed_linetype = "solid",
    line_width = 0.9,
    text_size = 2.5
) {
  stopifnot(length(x_pos) == 2)
  
  seg_df <- tibble::tibble(
    test = factor(c(legend_test, legend_test), levels = test_levels),
    abund_stratum = factor(c(legend_stratum, legend_stratum),
                           levels = c("Low", "Medium", "High")),
    x = c(x_pos[1], x_pos[1]),
    xend = c(x_pos[2], x_pos[2]),
    y = c(y_top, y_top - y_gap),
    yend = c(y_top, y_top - y_gap),
    line_group = c("threshold", "observed")
  )
  
  txt_df <- tibble::tibble(
    test = factor(c(legend_test, legend_test), levels = test_levels),
    abund_stratum = factor(c(legend_stratum, legend_stratum),
                           levels = c("Low", "Medium", "High")),
    x = c(x_pos[2] + text_gap, x_pos[2] + text_gap),
    y = c(y_top, y_top - y_gap),
    label = c(threshold_label, observed_label)
  )
  
  list(
    ggplot2::geom_segment(
      data = seg_df %>% dplyr::filter(line_group == "threshold"),
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = threshold_color,
      linetype = threshold_linetype,
      linewidth = line_width
    ),
    ggplot2::geom_segment(
      data = seg_df %>% dplyr::filter(line_group == "observed"),
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = observed_color,
      linetype = observed_linetype,
      linewidth = line_width
    ),
    ggplot2::geom_label(
      data = txt_df,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0.5,
      size = text_size,
      colour = "black",
      fill = "white",
      linewidth = 0.2,
      label.r = grid::unit(0.08, "lines"),
      label.padding = grid::unit(0.12, "lines")
    )
  )
}