# ============================================================
# Baseline abundance utilities
# ============================================================

#' Compute baseline abundance summaries from a SIMBA HDF5 file
#'
#' Reads the original filtered feature table from the HDF5 simulation file,
#' converts it to relative abundance, and computes the median nonzero abundance
#' for each feature.
#'
#' @param sim.location Character scalar giving the HDF5 file path.
#'
#' @return A tibble with columns:
#'   \itemize{
#'     \item \code{feature}
#'     \item \code{abundance_nonzero_median_baseline}
#'     \item \code{log10_abundance_nonzero_median_baseline}
#'   }
compute_baseline_abundance_from_h5 <- function(sim.location) {
  # Read the baseline filtered feature table and feature names.
  feat <- rhdf5::h5read(sim.location, "original_data/filt_features")
  feat_names <- rhdf5::h5read(sim.location, "original_data/filt_feature_names")
  
  rownames(feat) <- feat_names
  
  # Convert to relative abundance so abundance strata have a common meaning.
  feat <- sweep(feat, 2, colSums(feat), FUN = "/")
  
  # Median abundance among nonzero entries only.
  nonzero_median <- apply(feat, 1, function(x) {
    x <- x[x > 0]
    if (length(x) == 0) return(NA_real_)
    stats::median(x, na.rm = TRUE)
  })
  
  tibble::tibble(
    feature = rownames(feat),
    abundance_nonzero_median_baseline = as.numeric(nonzero_median),
    log10_abundance_nonzero_median_baseline = log10(as.numeric(nonzero_median))
  )
}

#' Add common abundance strata to a data frame
#'
#' Joins a feature-level abundance summary table and assigns each feature to a
#' shared abundance stratum, typically low/medium/high tertiles.
#'
#' @param df Data frame containing a \code{feature} column.
#' @param abundance_df Data frame containing \code{feature} and the abundance
#'   column to stratify on.
#' @param abundance_col Character scalar naming the abundance column in
#'   \code{abundance_df} after joining.
#' @param labels Character vector of labels for the strata.
#'
#' @return A copy of \code{df} with joined abundance columns and a new factor
#'   column \code{abund_stratum}.
add_common_abundance_strata <- function(df,
                                        abundance_df,
                                        abundance_col = "log10_abundance_nonzero_median_baseline",
                                        labels = c("Low", "Medium", "High")) {
  # If df already contains the abundance column, drop it before joining.
  # abundance_df is treated as the source of truth.
  df <- df %>%
    dplyr::select(-dplyr::any_of(abundance_col))
  
  abundance_df <- abundance_df %>%
    dplyr::select(feature, dplyr::all_of(abundance_col))
  
  # Join feature-level abundance summaries.
  df2 <- df %>%
    dplyr::left_join(abundance_df, by = "feature")
  
  # Define common cut points across all rows in the joined table.
  qs <- stats::quantile(df2[[abundance_col]], probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  qs <- unique(qs)
  
  # If the data are too degenerate to form tertiles, return NA strata.
  if (length(qs) < 4) {
    df2$abund_stratum <- NA_character_
    return(df2)
  }
  
  df2 %>%
    dplyr::mutate(
      abund_stratum = cut(
        .data[[abundance_col]],
        breaks = qs,
        include.lowest = TRUE,
        labels = labels
      )
    )
}


# ============================================================
# Case/control detectability analysis and plotting utilities
# ============================================================
#
# This file provides a streamlined set of helpers to:
#   1. derive feature-level prevalence and baseline abundance
#   2. standardize real DA results across methods
#   3. estimate abundance-stratified local detectability thresholds
#   4. summarize observed real effect-size envelopes
#   5. create facet-grid plots for:
#        - effect-size upper bounds
#        - raw p-value distributions
#
# Expected inputs:
#   - master_calib: long calibration table returned from binding
#     method_results[[method]]$calib_df across methods
#   - real_results_list: named list with one entry per method;
#     each method entry should contain named numeric vectors for
#     p-values and effect sizes
#   - feat_real: taxa x samples relative abundance matrix
#   - prevalence: named numeric vector of taxon prevalence in feat_real
#
# Final entry points:
#   - run_casecontrol_upper_bound_abundance()
#   - plot_real_pvalue_distribution_abundance()
#
# ============================================================


# ------------------------------------------------------------
# 1. Feature-level prevalence and abundance statistics
# ------------------------------------------------------------

#' Compute prevalence and baseline abundance from a relative-abundance matrix
#'
#' Computes taxon prevalence, median nonzero relative abundance, and its
#' log10-transform from a taxa x samples relative-abundance matrix.
#'
#' @param feat_real Numeric matrix with taxa in rows and samples in columns.
#'   Values should be relative abundances.
#'
#' @return A tibble with columns:
#'   \itemize{
#'     \item \code{feature}: taxon name
#'     \item \code{prevalence}: fraction of samples with abundance > 0
#'     \item \code{abundance_nonzero_median}: median abundance among nonzero samples
#'     \item \code{log10_abundance_nonzero_median}: log10 median nonzero abundance
#'   }
#'
compute_real_feature_stats <- function(feat_real) {
  stopifnot(is.data.frame(feat_real))
  stopifnot(!is.null(rownames(feat_real)))

  prevalence <- rowMeans(feat_real > 0, na.rm = TRUE)
  
  # Median among nonzero entries only. This separates prevalence from abundance level.
  nonzero_median <- apply(feat_real, 1, function(x) {
    x <- x[x > 0]
    if (length(x) == 0) {
      return(NA_real_)
    }
    stats::median(x, na.rm = TRUE)
  })
  
  tibble::tibble(
    feature = rownames(feat_real),
    prevalence = as.numeric(prevalence),
    abundance_nonzero_median = as.numeric(nonzero_median),
    log10_abundance_nonzero_median = log10(as.numeric(nonzero_median))
  )
}


#' Add abundance tertiles to a data frame
#'
#' Adds low/medium/high abundance strata based on tertiles of a chosen abundance
#' column. This is used both for real-data feature summaries and for calibration
#' tables.
#'
#' @param df Data frame containing an abundance column.
#' @param abundance_col Character scalar naming the abundance column.
#' @param labels Character vector of labels for the tertiles.
#'
#' @return Input data frame with a new factor column \code{abund_stratum}.
#'
add_abundance_tertiles <- function(df,
                                   abundance_col,
                                   labels = c("Low", "Medium", "High")) {
  stopifnot(is.data.frame(df))
  stopifnot(abundance_col %in% colnames(df))
  
  qs <- stats::quantile(df[[abundance_col]], probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  qs <- unique(qs)
  
  # If the abundance distribution is too degenerate, we cannot form tertiles.
  if (length(qs) < 4) {
    df$abund_stratum <- NA_character_
    return(df)
  }
  
  df$abund_stratum <- cut(
    df[[abundance_col]],
    breaks = qs,
    include.lowest = TRUE,
    labels = labels
  )
  
  df
}


#' Ensure that a calibration table contains abundance strata
#'
#' If \code{master_calib} already contains \code{abund_stratum}, it is returned
#' unchanged except for factor ordering. Otherwise, the function looks for a
#' baseline abundance column and constructs common tertiles.
#'
#' @param master_calib Calibration data frame.
#' @param abundance_col_candidates Character vector of candidate abundance columns.
#'
#' @return Calibration data frame with factor column \code{abund_stratum}.
#'
ensure_calibration_abundance_strata <- function(
    master_calib,
    abundance_col_candidates = c(
      "log10_abundance_nonzero_median_baseline",
      "log10_abundance_nonzero_median"
    )
) {
  stopifnot(is.data.frame(master_calib))
  
  if ("abund_stratum" %in% colnames(master_calib)) {
    master_calib$abund_stratum <- factor(
      master_calib$abund_stratum,
      levels = c("Low", "Medium", "High")
    )
    return(master_calib)
  }
  
  abundance_col <- intersect(abundance_col_candidates, colnames(master_calib))
  if (length(abundance_col) == 0) {
    stop(
      "master_calib must contain either 'abund_stratum' or one of these abundance columns: ",
      paste(shQuote(abundance_col_candidates), collapse = ", ")
    )
  }
  
  out <- add_abundance_tertiles(master_calib, abundance_col = abundance_col[1])
  out$abund_stratum <- factor(out$abund_stratum, levels = c("Low", "Medium", "High"))
  out
}


# ------------------------------------------------------------
# 2. Real-method results standardization
# ------------------------------------------------------------

#' Resolve a named numeric vector from a method result list entry
#'
#' Helper used internally to support either \code{pval}/\code{eff} or
#' \code{p.val}/\code{eff.size} naming conventions.
#'
#' @param x List-like object for one method.
#' @param candidates Character vector of candidate component names.
#' @param method_name Character scalar used only in error messages.
#' @param component_label Human-readable component label for error messages.
#'
#' @return Named numeric vector.
#'
resolve_named_numeric_component <- function(x,
                                            candidates,
                                            method_name,
                                            component_label) {
  nm <- intersect(candidates, names(x))
  if (length(nm) == 0) {
    stop(
      "Could not find ", component_label, " for method '", method_name,
      "'. Tried: ", paste(shQuote(candidates), collapse = ", ")
    )
  }
  
  out <- x[[nm[1]]]
  
  # Allow one-column matrices/data.frames as input.
  if (is.matrix(out) || is.data.frame(out)) {
    if (ncol(out) != 1) {
      stop(
        component_label, " for method '", method_name,
        "' must be a named numeric vector or one-column matrix/data.frame."
      )
    }
    out <- out[, 1]
  }
  
  if (!is.numeric(out)) {
    stop(component_label, " for method '", method_name, "' is not numeric.")
  }
  if (is.null(names(out))) {
    stop(component_label, " for method '", method_name, "' must be named by feature.")
  }
  
  out
}


#' Build a standardized real-data result table for one method
#'
#' Constructs a method-specific table containing raw p-values, adjusted p-values,
#' effect sizes, prevalence, and an indicator of nominal/significant detection.
#'
#' @param pval Named numeric vector of raw p-values.
#' @param eff Named numeric vector of effect sizes.
#' @param prevalence Named numeric vector of taxon prevalence.
#' @param test Character scalar naming the method.
#' @param adjust Multiple-testing adjustment method passed to \code{p.adjust()}.
#'   Use \code{"none"} to retain raw p-values.
#'
#' @return A tibble with one row per taxon.
#'
make_real_effect_df <- function(pval,
                                eff,
                                prevalence,
                                test,
                                alpha = 0.05,
                                adjust = "none") {
  stopifnot(is.numeric(pval), !is.null(names(pval)))
  stopifnot(is.numeric(eff), !is.null(names(eff)))
  stopifnot(is.numeric(prevalence), !is.null(names(prevalence)))
  
  feat <- Reduce(intersect, list(names(pval), names(eff), names(prevalence)))
  
  if (length(feat) == 0) {
    stop("No overlapping feature names found for method '", test, "'.")
  }
  
  p_raw <- as.numeric(pval[feat])
  
  p_adj <- if (identical(adjust, "none")) {
    p_raw
  } else {
    stats::p.adjust(p_raw, method = adjust)
  }
  
  tibble::tibble(
    feature = feat,
    test = test,
    p_raw = p_raw,
    p_adj = p_adj,
    effect_size = as.numeric(eff[feat]),
    abs_effect_size = abs(as.numeric(eff[feat])),
    prevalence = as.numeric(prevalence[feat]),
    detected = p_adj < alpha
  ) %>%
    filter(!is.na(p_raw))
}


#' Build a standardized real-data result table across methods
#'
#' Converts a named list of method results into one combined long-format table.
#'
#' @param real_results_list Named list with one entry per method.
#'   Each method entry must contain named numeric vectors for p-values and
#'   effect sizes. Supported component names are:
#'   \itemize{
#'     \item p-values: \code{pval} or \code{p.val}
#'     \item effect sizes: \code{eff} or \code{eff.size}
#'   }
#' @param prevalence Named numeric vector of taxon prevalence.
#' @param methods Character vector of methods to include.
#' @param adjust Multiple-testing adjustment method passed to \code{p.adjust()}.
#'
#' @return A tibble with one row per taxon per method.
#'
build_real_effect_df <- function(real_results_list,
                                 prevalence,
                                 methods = names(real_results_list),
                                 adjust = "none",
                                 alpha = 0.05) {
  purrr::map_dfr(methods, function(m) {
    method_obj <- real_results_list[[m]]
    
    if (is.null(method_obj)) {
      stop("Method '", m, "' not found in real_results_list.")
    }
    
    pval <- resolve_named_numeric_component(
      x = method_obj,
      candidates = c("pval", "p.val"),
      method_name = m,
      component_label = "p-values"
    )
    
    eff <- resolve_named_numeric_component(
      x = method_obj,
      candidates = c("eff", "eff.size"),
      method_name = m,
      component_label = "effect sizes"
    )
    
    make_real_effect_df(
      pval = pval,
      eff = eff,
      prevalence = prevalence,
      test = m,
      adjust = adjust,
      alpha = alpha
    )
  })
}


# ------------------------------------------------------------
# 3. Local prevalence-matched summaries for calibration
# ------------------------------------------------------------

#' Fit local detectability thresholds stratified by abundance
#'
#' For each method and abundance stratum, fits a local logistic detection model
#' in moving prevalence windows and inverts it to estimate the effect size
#' threshold corresponding to target detection probabilities.
#'
#' @param calib_df Calibration table. Must contain columns:
#'   \code{test}, \code{is_marker}, \code{detected}, \code{abs_effect_size},
#'   \code{prevalence_full}, and \code{abund_stratum}.
#' @param probs Numeric vector of target detection probabilities.
#' @param prev_grid Numeric vector of prevalence grid values at which local
#'   thresholds are estimated.
#' @param prevalence_col Character scalar naming the prevalence column in
#'   \code{calib_df}.
#' @param bandwidth Numeric half-width of the local prevalence window.
#' @param abundance_stratum_col Character scalar naming the abundance stratum column.
#' @param min_n Minimum number of taxa required in a local window.
#' @param min_detected Minimum number of detected taxa required in a local window.
#' @param grid_length Number of grid points used when inverting the detection curve.
#'
#' @return A tibble with columns including \code{test}, \code{abund_stratum},
#'   \code{prev_mid}, \code{target_prob}, and \code{threshold}.
#'
fit_detection_thresholds_local_abundance <- function(
    calib_df,
    probs = c(0.5, 0.8, 0.95),
    prev_grid = seq(0.10, 0.95, by = 0.05),
    prevalence_col = "prevalence_full",
    bandwidth = 0.08,
    abundance_stratum_col = "abund_stratum",
    min_n = 100,
    min_detected = 20,
    grid_length = 500
) {
  calib_df %>%
    dplyr::filter(
      is_marker,
      !is.na(abs_effect_size),
      !is.na(detected),
      !is.na(.data[[prevalence_col]]),
      !is.na(.data[[abundance_stratum_col]])
    ) %>%
    dplyr::group_by(test, .data[[abundance_stratum_col]]) %>%
    dplyr::group_modify(function(.x, .y) {
      purrr::map_dfr(prev_grid, function(p0) {
        # Local prevalence window around p0
        sub <- .x %>%
          dplyr::filter(abs(.data[[prevalence_col]] - p0) <= bandwidth) %>%
          dplyr::filter(is.finite(abs_effect_size), !is.na(detected))
        
        n_all <- nrow(sub)
        n_det <- sum(sub$detected, na.rm = TRUE)
        
        if (n_all < min_n || n_det < min_detected) {
          return(tibble::tibble(
            prev_mid = p0,
            target_prob = probs,
            threshold = NA_real_,
            n = n_all,
            n_detected = n_det
          ))
        }
        
        # Spline if enough unique effect-size values; otherwise fall back to linear.
        u_eff <- sort(unique(sub$abs_effect_size[is.finite(sub$abs_effect_size)]))
        n_uniq <- length(u_eff)
        
        form <- if (n_uniq >= 6) {
          stats::as.formula(as.numeric(detected) ~ splines::ns(abs_effect_size, df = 3))
        } else if (n_uniq >= 4) {
          stats::as.formula(as.numeric(detected) ~ abs_effect_size)
        } else {
          NULL
        }
        
        if (is.null(form)) {
          return(tibble::tibble(
            prev_mid = p0,
            target_prob = probs,
            threshold = NA_real_,
            n = n_all,
            n_detected = n_det
          ))
        }
        
        fit <- tryCatch(
          stats::glm(form, data = sub, family = stats::binomial()),
          error = function(e) NULL
        )
        
        if (is.null(fit)) {
          return(tibble::tibble(
            prev_mid = p0,
            target_prob = probs,
            threshold = NA_real_,
            n = n_all,
            n_detected = n_det
          ))
        }
        
        x_grid <- seq(0, max(sub$abs_effect_size, na.rm = TRUE), length.out = grid_length)
        
        pred <- stats::predict(
          fit,
          newdata = data.frame(abs_effect_size = x_grid),
          type = "response"
        )
        
        purrr::map_dfr(probs, function(p_target) {
          idx <- which(pred >= p_target)
          tibble::tibble(
            prev_mid = p0,
            target_prob = p_target,
            threshold = if (length(idx) == 0) NA_real_ else x_grid[min(idx)],
            n = n_all,
            n_detected = n_det
          )
        })
      })
    }) %>%
    dplyr::ungroup() %>%
    dplyr::rename(abund_stratum = !!abundance_stratum_col)
}

#' Fit local detectability thresholds stratified by abundance using empirical monotone calibration
#'
#' This is a drop-in empirical replacement for fit_detection_thresholds_local_abundance().
#' Instead of a logistic model, it uses isotonic regression to estimate a monotone
#' detection probability curve as a function of absolute effect size.
#'
#' @param calib_df Calibration table. Must contain columns:
#'   \code{test}, \code{is_marker}, \code{detected}, \code{abs_effect_size},
#'   \code{prevalence_full}, and \code{abund_stratum}.
#' @param probs Numeric vector of target detection probabilities.
#' @param prev_grid Numeric vector of prevalence grid values at which local
#'   thresholds are estimated.
#' @param prevalence_col Character scalar naming the prevalence column in
#'   \code{calib_df}.
#' @param bandwidth Numeric half-width of the local prevalence window.
#' @param abundance_stratum_col Character scalar naming the abundance stratum column.
#' @param min_n Minimum number of taxa required in a local window.
#' @param min_detected Minimum number of detected taxa required in a local window.
#' @param min_unique_effects Minimum number of unique effect sizes required in a
#'   local window to fit a meaningful monotone empirical curve.
#'
#' @return A tibble with columns including \code{test}, \code{abund_stratum},
#'   \code{prev_mid}, \code{target_prob}, and \code{threshold}.
#'
fit_detection_thresholds_local_abundance_empirical <- function(
    calib_df,
    probs = c(0.5, 0.8, 0.95),
    prev_grid = seq(0.10, 0.95, by = 0.05),
    prevalence_col = "prevalence_full",
    bandwidth = 0.08,
    abundance_stratum_col = "abund_stratum",
    min_n = 100,
    min_detected = 20,
    min_unique_effects = 6
) {
  calib_df %>%
    dplyr::filter(
      is_marker,
      !is.na(abs_effect_size),
      !is.na(detected),
      !is.na(.data[[prevalence_col]]),
      !is.na(.data[[abundance_stratum_col]])
    ) %>%
    dplyr::group_by(test, .data[[abundance_stratum_col]]) %>%
    dplyr::group_modify(function(.x, .y) {
      purrr::map_dfr(prev_grid, function(p0) {
        # Local prevalence window
        sub <- .x %>%
          dplyr::filter(abs(.data[[prevalence_col]] - p0) <= bandwidth) %>%
          dplyr::filter(is.finite(abs_effect_size), !is.na(detected))
        
        n_all <- nrow(sub)
        n_det <- sum(sub$detected, na.rm = TRUE)
        
        if (n_all < min_n || n_det < min_detected) {
          return(tibble::tibble(
            prev_mid = p0,
            target_prob = probs,
            threshold = NA_real_,
            n = n_all,
            n_detected = n_det,
            fit_type = "empirical_isotonic"
          ))
        }
        
        # Sort by effect size
        sub <- sub %>%
          dplyr::arrange(abs_effect_size)
        
        x <- sub$abs_effect_size
        y <- as.numeric(sub$detected)
        
        # Need enough unique effect-size support
        n_uniq <- length(unique(x))
        if (n_uniq < min_unique_effects) {
          return(tibble::tibble(
            prev_mid = p0,
            target_prob = probs,
            threshold = NA_real_,
            n = n_all,
            n_detected = n_det,
            fit_type = "empirical_isotonic"
          ))
        }
        
        # Empirical monotone fit of detection probability vs effect size
        iso <- tryCatch(
          stats::isoreg(x = x, y = y),
          error = function(e) NULL
        )
        
        if (is.null(iso)) {
          return(tibble::tibble(
            prev_mid = p0,
            target_prob = probs,
            threshold = NA_real_,
            n = n_all,
            n_detected = n_det,
            fit_type = "empirical_isotonic"
          ))
        }
        
        # Fitted values are stepwise and monotone nondecreasing
        fitted_prob <- pmin(pmax(iso$yf, 0), 1)
        
        purrr::map_dfr(probs, function(p_target) {
          idx <- which(fitted_prob >= p_target)
          tibble::tibble(
            prev_mid = p0,
            target_prob = p_target,
            threshold = if (length(idx) == 0) NA_real_ else x[min(idx)],
            n = n_all,
            n_detected = n_det,
            fit_type = "empirical_isotonic"
          )
        })
      })
    }) %>%
    dplyr::ungroup() %>%
    dplyr::rename(abund_stratum = !!abundance_stratum_col)
}

#' Summarize the local observed effect-size envelope stratified by abundance
#'
#' Computes local quantiles of the observed real effect-size distribution in
#' moving prevalence windows for each method and abundance stratum.
#'
#' @param real_df Real-data result table containing \code{test}, \code{prevalence},
#'   \code{abs_effect_size}, and \code{abund_stratum}.
#' @param prev_grid Numeric vector of prevalence grid values.
#' @param bandwidth Numeric half-width of the local prevalence window.
#' @param abundance_stratum_col Character scalar naming the abundance stratum column.
#' @param min_n Minimum number of taxa required in a local window.
#'
#' @return A tibble with local q50/q90/q95/q99 summaries.
#'
summarise_real_effect_envelope_local_abundance <- function(
    real_df,
    prev_grid = seq(0.10, 0.95, by = 0.05),
    bandwidth = 0.08,
    abundance_stratum_col = "abund_stratum",
    min_n = 30
) {
  real_df %>%
    dplyr::filter(
      !is.na(prevalence),
      !is.na(abs_effect_size),
      !is.na(.data[[abundance_stratum_col]])
    ) %>%
    dplyr::group_by(test, .data[[abundance_stratum_col]]) %>%
    dplyr::group_modify(function(.x, .y) {
      purrr::map_dfr(prev_grid, function(p0) {
        sub <- .x %>%
          dplyr::filter(abs(prevalence - p0) <= bandwidth)
        
        # Fallback to nearest min_n taxa if the local window is sparse.
        if (nrow(sub) < min_n && nrow(.x) >= min_n) {
          sub <- .x %>%
            dplyr::mutate(dist_prev = abs(prevalence - p0)) %>%
            dplyr::arrange(dist_prev) %>%
            dplyr::slice_head(n = min_n)
        }
        
        n_taxa <- nrow(sub)
        
        if (n_taxa == 0) {
          return(tibble::tibble(
            prev_mid = p0,
            n_taxa = 0,
            q50 = NA_real_,
            q90 = NA_real_,
            q95 = NA_real_,
            q99 = NA_real_
          ))
        }
        
        tibble::tibble(
          prev_mid = p0,
          n_taxa = n_taxa,
          q50 = stats::quantile(sub$abs_effect_size, 0.50, na.rm = TRUE),
          q90 = stats::quantile(sub$abs_effect_size, 0.90, na.rm = TRUE),
          q95 = stats::quantile(sub$abs_effect_size, 0.95, na.rm = TRUE),
          q99 = stats::quantile(sub$abs_effect_size, 0.99, na.rm = TRUE)
        )
      })
    }) %>%
    dplyr::ungroup() %>%
    dplyr::rename(abund_stratum = !!abundance_stratum_col)
}


#' Summarize the local p-value envelope stratified by abundance
#'
#' Computes a local summary of \code{-log10(raw p-value)} across prevalence for
#' each method and abundance stratum. This is intended as a support plot to assess
#' whether taxa with apparently larger effect sizes are also associated with
#' stronger nominal statistical evidence.
#'
#' @param df Real-data result table containing \code{test}, \code{abund_stratum},
#'   \code{prevalence}, and \code{logp}.
#' @param prev_grid Numeric vector of prevalence grid values.
#' @param bandwidth Numeric half-width of the local prevalence window.
#' @param min_n Minimum number of taxa required in a local window.
#' @param prob Quantile used for the upper-envelope summary.
#'
#' @return A tibble with columns \code{prev_mid}, \code{q95_logp}, and \code{median_logp}.
#'
summarise_local_pvalue_curve <- function(df,
                                         prev_grid = seq(0.10, 0.95, by = 0.05),
                                         bandwidth = 0.08,
                                         min_n = 30,
                                         prob = 0.95) {
  df %>%
    dplyr::filter(!is.na(prevalence), !is.na(logp), !is.na(abund_stratum)) %>%
    dplyr::group_by(test, abund_stratum) %>%
    dplyr::group_modify(function(.x, .y) {
      purrr::map_dfr(prev_grid, function(p0) {
        sub <- .x %>%
          dplyr::filter(abs(prevalence - p0) <= bandwidth)
        
        # Fallback to nearest min_n taxa if the local window is sparse.
        if (nrow(sub) < min_n && nrow(.x) >= min_n) {
          sub <- .x %>%
            dplyr::mutate(dist_prev = abs(prevalence - p0)) %>%
            dplyr::arrange(dist_prev) %>%
            dplyr::slice_head(n = min_n)
        }
        
        if (nrow(sub) == 0) {
          return(tibble::tibble(
            prev_mid = p0,
            n_taxa = 0,
            q95_logp = NA_real_,
            median_logp = NA_real_
          ))
        }
        
        tibble::tibble(
          prev_mid = p0,
          n_taxa = nrow(sub),
          q95_logp = stats::quantile(sub$logp, prob = prob, na.rm = TRUE),
          median_logp = stats::median(sub$logp, na.rm = TRUE)
        )
      })
    }) %>%
    dplyr::ungroup()
}

# ============================================================
# 4. Model-based abundance-stratified detectability boundaries
# ============================================================

# ------------------------------------------------------------
# Model fitting and boundary prediction helpers
# ------------------------------------------------------------

#' Fit additive detectability models stratified by method
#'
#' Fits one logistic regression model per method to predict marker detectability
#' from absolute effect size, prevalence, and continuous baseline abundance.
#'
#' The intended use is to generate a smooth model-based detectability boundary
#' for abundance-stratified plotting, while keeping abundance as a continuous
#' predictor in the model itself.
#'
#' @param calib_df Calibration table containing at least:
#'   \code{test}, \code{is_marker}, \code{detected}, \code{abs_effect_size},
#'   a prevalence column, and a continuous abundance column.
#' @param methods Character vector of methods to include. Defaults to all methods
#'   present in \code{calib_df}.
#' @param prevalence_col Character scalar naming the prevalence column.
#' @param abundance_col Character scalar naming the continuous abundance column.
#' @param min_rows Minimum number of rows required to fit a model for a method.
#' @param formula Optional model formula. Defaults to additive M3:
#'   \code{detected_bin ~ effect + prevalence + abundance}.
#'
#' @return A named list of fitted \code{glm} objects, one per method.
fit_detectability_boundary_models_abundance <- function(
    calib_df,
    methods = sort(unique(calib_df$test)),
    prevalence_col = "prevalence_full",
    abundance_col = "log10_abundance_nonzero_median_baseline",
    min_rows = 50,
    formula = detected_bin ~ effect + prevalence + abundance
) {
  library("mgcv")
  
  req <- c("test", "is_marker", "detected", "abs_effect_size", prevalence_col, abundance_col)
  miss <- setdiff(req, colnames(calib_df))
  if (length(miss) > 0) {
    stop("calib_df is missing required columns: ",
         paste(shQuote(miss), collapse = ", "))
  }
  
  dat <- calib_df %>%
    dplyr::filter(
      test %in% methods,
      is_marker,
      !is.na(detected),
      !is.na(abs_effect_size),
      !is.na(.data[[prevalence_col]]),
      !is.na(.data[[abundance_col]])
    ) %>%
    dplyr::mutate(
      detected_bin = as.integer(detected),
      effect = as.numeric(abs_effect_size),
      prevalence = as.numeric(.data[[prevalence_col]]),
      abundance = as.numeric(.data[[abundance_col]])
    )
  
  fits <- purrr::map(
    methods,
    function(m) {
      dfm <- dat %>% dplyr::filter(test == m)
      
      if (nrow(dfm) < min_rows) {
        return(NULL)
      }
      if (length(unique(dfm$detected_bin)) < 2) {
        return(NULL)
      }
      
      #stats::glm(formula, data = dfm, family = stats::binomial())
      mgcv::gam(formula, data = dfm, family = stats::binomial())
    }
  )
  
  names(fits) <- methods
  fits[!vapply(fits, is.null, logical(1))]
}


#' Compute representative abundance values for plotting strata
#'
#' Uses a feature-level abundance table to derive one representative abundance
#' value per abundance stratum, typically the median abundance within each
#' stratum. These values are then plugged into the continuous detectability
#' model to generate one smooth boundary per stratum.
#'
#' @param abundance_df Feature-level abundance table containing
#'   \code{abund_stratum} and the chosen abundance column.
#' @param abundance_col Character scalar naming the continuous abundance column.
#' @param abundance_stratum_col Character scalar naming the abundance-stratum
#'   column.
#' @param summary_fun Function used to summarize abundance within each stratum.
#'   Defaults to \code{stats::median}.
#'
#' @return A tibble with columns \code{abund_stratum} and
#'   \code{abundance_value}.
compute_abundance_reference_values <- function(
    abundance_df,
    abundance_col = "log10_abundance_nonzero_median_baseline",
    abundance_stratum_col = "abund_stratum",
    summary_fun = stats::median
) {
  req <- c(abundance_col, abundance_stratum_col)
  miss <- setdiff(req, colnames(abundance_df))
  if (length(miss) > 0) {
    stop("abundance_df is missing required columns: ",
         paste(shQuote(miss), collapse = ", "))
  }
  
  abundance_df %>%
    dplyr::filter(!is.na(.data[[abundance_col]]), !is.na(.data[[abundance_stratum_col]])) %>%
    dplyr::group_by(.data[[abundance_stratum_col]]) %>%
    dplyr::summarise(
      abundance_value = summary_fun(.data[[abundance_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(abund_stratum = !!abundance_stratum_col)
}


#' Predict a GAM-based detectability boundary across prevalence and abundance strata
#'
#' This is a drop-in replacement for the previous parametric boundary predictor.
#' For each fitted method-specific GAM, it predicts detection probability across
#' a prevalence grid at representative abundance values for each stratum, then
#' numerically inverts the fitted surface to find the smallest effect size at
#' which the predicted detection probability reaches the target value.
#'
#' Compared with the earlier GLM-based version, this allows the 80% contour to
#' bend more flexibly and better follow the simulated detected/missed cloud,
#' which is particularly useful in smaller relapse settings.
#'
#' @param model_list Named list of fitted \code{mgcv::gam} objects.
#' @param calib_df Calibration table used to determine a reasonable effect-size
#'   search range per method.
#' @param abundance_ref_df Tibble with columns \code{abund_stratum} and
#'   \code{abundance_value}, typically generated from fixed abundance strata in
#'   the original feature table.
#' @param target_prob Target detection probability to invert, for example 0.8.
#' @param prev_grid Numeric vector of prevalence values at which the boundary
#'   should be predicted.
#' @param effect_col Character scalar naming the effect-size column in
#'   \code{calib_df}.
#' @param grid_length Number of effect-size grid points used for numerical
#'   inversion.
#' @param effect_upper_quantile Upper quantile of observed simulated effect sizes
#'   used to define the search range per method.
#'
#' @return A tibble with columns:
#'   \itemize{
#'     \item \code{test}
#'     \item \code{abund_stratum}
#'     \item \code{prev_mid}
#'     \item \code{target_prob}
#'     \item \code{threshold}
#'     \item \code{fit_type}
#'   }
predict_detectability_boundary_abundance <- function(
    model_list,
    calib_df,
    abundance_ref_df,
    target_prob = 0.8,
    prev_grid = seq(0.10, 0.95, by = 0.05),
    effect_col = "abs_effect_size",
    grid_length = 500,
    effect_upper_quantile = 0.995
) {
  req_calib <- c("test", effect_col)
  miss_calib <- setdiff(req_calib, colnames(calib_df))
  if (length(miss_calib) > 0) {
    stop("calib_df is missing required columns: ",
         paste(shQuote(miss_calib), collapse = ", "))
  }
  
  req_ref <- c("abund_stratum", "abundance_value")
  miss_ref <- setdiff(req_ref, colnames(abundance_ref_df))
  if (length(miss_ref) > 0) {
    stop("abundance_ref_df is missing required columns: ",
         paste(shQuote(miss_ref), collapse = ", "))
  }
  
  purrr::imap_dfr(model_list, function(fit, method_name) {
    # Use the simulated effect-size range for this method to define a search grid.
    dfm <- calib_df %>%
      dplyr::filter(
        test == method_name,
        !is.na(.data[[effect_col]]),
        is.finite(.data[[effect_col]])
      )
    
    if (nrow(dfm) == 0) {
      return(NULL)
    }
    
    eff_max <- stats::quantile(dfm[[effect_col]], probs = effect_upper_quantile, na.rm = TRUE)
    eff_max <- max(eff_max, max(dfm[[effect_col]], na.rm = TRUE), 1e-8)
    
    eff_grid <- seq(0, eff_max, length.out = grid_length)
    
    purrr::map_dfr(seq_len(nrow(abundance_ref_df)), function(i) {
      abund_stratum_i <- abundance_ref_df$abund_stratum[i]
      abundance_i <- abundance_ref_df$abundance_value[i]
      
      purrr::map_dfr(prev_grid, function(p0) {
        newdata <- data.frame(
          effect = eff_grid,
          prevalence = p0,
          abundance = abundance_i
        )
        
        # Predict detectability over the effect-size grid and find the first
        # point at which the fitted surface reaches the target probability.
        pred <- stats::predict(fit, newdata = newdata, type = "response")
        idx <- which(pred >= target_prob)
        
        tibble::tibble(
          test = method_name,
          abund_stratum = abund_stratum_i,
          prev_mid = p0,
          target_prob = target_prob,
          threshold = if (length(idx) == 0) NA_real_ else eff_grid[min(idx)],
          #fit_type = "model_glm_additive"
          fit_type = "model_gam"
        )
      })
    })
  })
}

# ------------------------------------------------------------
# 4. Plotting functions
# ------------------------------------------------------------

#' Plot abundance-stratified real effects against a model-based detectability boundary
#'
#' This is a drop-in replacement for the previous
#' \code{plot_effect_size_upper_bounds_abundance()} function. It expects the
#' output of the redesigned \code{run_casecontrol_upper_bound_abundance()} and
#' keeps the same overall visual structure:
#'
#' \itemize{
#'   \item grey points: all real taxa
#'   \item red points: significant real taxa
#'   \item blue dashed line: model-based detectability boundary
#'   \item black line: local real-data upper envelope
#' }
#'
#' @param res_bound List returned by
#'   \code{run_casecontrol_upper_bound_abundance()}.
#' @param methods Optional character vector of methods to plot. Defaults to all
#'   methods present in \code{res_bound$real_df}.
#' @param envelope_stat Character scalar naming which envelope to plot.
#'   Supported values are \code{"q50"}, \code{"q90"}, \code{"q95"}, and
#'   \code{"q99"}.
#' @param point_alpha Alpha for grey real-data points.
#' @param point_size Point size for real-data points.
#' @param detected_point_size Point size for red detected points.
#' @param threshold_linewidth Line width for the blue boundary.
#' @param envelope_linewidth Line width for the black envelope.
#' @param facet_scales Passed to \code{facet_grid()}.
#' @param threshold_label Label used in the legend for the blue boundary.
#' @param envelope_label Label used in the legend for the black envelope.
#'
#' @return A ggplot object.
plot_effect_size_upper_bounds_abundance <- function(
    res_bound,
    methods = NULL,
    point_alpha = 0.5,
    point_size = 0.8,
    detected_point_size = 2,
    threshold_linewidth = 0.7,
    facet_scales = "free_y",
    significant_label = "Recurrence-associated species (meta_mOTU_v3_12389)",
    insignificant_label = "Species failing the association test",
    threshold_label = "80% detectability boundary (stratum)"
) {
  stopifnot(is.list(res_bound))
  stopifnot(all(c("real_df", "envelope_df", "thresh_df") %in% names(res_bound)))
  
  real_df <- res_bound$real_df
  thresh_df <- res_bound$thresh_df
  
  if (is.null(methods)) {
    methods <- sort(unique(real_df$test))
  }
  
  real_df <- real_df %>% 
    dplyr::filter(test %in% methods) %>%
    mutate(col_label = "Abundance stratum")
  thresh_df <- thresh_df %>%
    dplyr::filter(test %in% methods) %>%
    mutate(col_label = "Abundance stratum")

  p <- ggplot2::ggplot() +
    # All real taxa
    ggplot2::geom_point(
      data = real_df,
      ggplot2::aes(x = prevalence, y = abs_effect_size, colour = insignificant_label),
      alpha = point_alpha,
      size = point_size
    ) +
    # Significant real taxa
    ggplot2::geom_point(
      data = real_df %>% dplyr::filter(detected),
      ggplot2::aes(x = prevalence, y = abs_effect_size, colour = significant_label),
      size = detected_point_size
    ) +
    # Model-based boundary
    ggplot2::geom_line(
      data = thresh_df,
      ggplot2::aes(x = prev_mid, y = threshold, linetype = threshold_label),
      colour = "#8E1F63",
      linewidth = threshold_linewidth,
      na.rm = TRUE
    ) +
    ggh4x::facet_nested(test ~ col_label + abund_stratum, scales = facet_scales) +
    ggplot2::scale_linetype_manual(
      values = stats::setNames(
        c("solid"),
        c(threshold_label)
      ),
      name = NULL
    ) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(
        c("#9E9E9E", "#BE5302"),
        c(insignificant_label, significant_label)
      ),
      name = NULL
    ) +
    ggplot2::labs(
      x = "Taxon prevalence",
      y = "Absolute effect size"
    ) +
    guides(
      colour = guide_legend(
        nrow = 2,
        byrow = TRUE,
        override.aes = list(
          size = c(3, 2),
          alpha = c(1, 1)
        )
      ),
      linetype = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          linewidth = c(1.1)
        )
      ),
      shape = guide_legend(nrow = 2, byrow = TRUE)
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white", colour = "black"),
      strip.text = ggplot2::element_text(size = 11, colour = "black"),
      axis.text = ggplot2::element_text(size = 10, colour = "black"),
      axis.title = ggplot2::element_text(size = 13, colour = "black"),
      panel.spacing = grid::unit(0.8, "lines"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      panel.grid.minor = ggplot2::element_blank()
    )
  p
}


#' Plot abundance-stratified raw p-value distributions
#'
#' Creates a facet-grid plot with methods in rows and abundance strata in columns.
#' Grey points show individual taxa and the dashed line marks raw \eqn{p = 0.05}.
#' The black line summarizes the local p-value envelope across prevalence.
#'
#' @param real_results_list Named list of real-method results.
#' @param feat_real Taxa x samples relative abundance matrix.
#' @param methods Character vector of methods to include and their display order.
#' @param prevalence Optional named numeric vector of taxon prevalence. If \code{NULL},
#'   prevalence is recomputed from \code{feat_real}.
#' @param adjust Multiple-testing adjustment used only to define the \code{detected}
#'   flag. Use \code{"none"} to highlight nominally significant taxa.
#' @param prev_grid Numeric vector of prevalence grid values.
#' @param bandwidth Numeric half-width of the local prevalence window.
#' @param min_n Minimum number of taxa required in a local window.
#' @param summary_stat Which local summary line to draw: \code{"q95"} or \code{"median"}.
#' @param point_alpha Point transparency.
#' @param point_size Point size.
#'
#' @return A ggplot object.
#' 
#' @export
#'
plot_real_pvalue_distribution_abundance <- function(real_results_list,
                                                    feat_real,
                                                    methods,
                                                    prevalence = NULL,
                                                    adjust = "fdr",
                                                    alpha = 0.05,
                                                    prev_grid = seq(0.10, 0.95, by = 0.05),
                                                    bandwidth = 0.08,
                                                    min_n = 30,
                                                    summary_stat = c("q95", "median"),
                                                    show_envelope = FALSE,
                                                    point_alpha = 0.5,
                                                    point_size = 0.8) {
  summary_stat <- match.arg(summary_stat)
  
  # Use supplied prevalence if provided; otherwise compute it from feat_real.
  if (is.null(prevalence)) {
    feature_stats <- compute_real_feature_stats(feat_real)
    prevalence <- stats::setNames(feature_stats$prevalence, feature_stats$feature)
  }
  
  real_df <- build_real_effect_df(
    real_results_list = real_results_list,
    prevalence = prevalence,
    methods = methods,
    adjust = adjust,
    alpha = alpha
  )
  
  feature_stats <- compute_real_feature_stats(feat_real)
  
  plot_df <- real_df %>%
    dplyr::left_join(feature_stats, by = c("feature", "prevalence")) %>%
    add_abundance_tertiles(abundance_col = "log10_abundance_nonzero_median") %>%
    dplyr::mutate(
      test = factor(test, levels = methods),
      col_label = "Abundance stratum",
      abund_stratum = factor(abund_stratum, levels = c("Low", "Medium", "High")),
      p_raw = pmax(p_raw, .Machine$double.xmin),
      logp = -log10(p_raw)
    )
  
  line_df <- summarise_local_pvalue_curve(
    df = plot_df,
    prev_grid = prev_grid,
    bandwidth = bandwidth,
    min_n = min_n
  ) %>%
    dplyr::mutate(
      test = factor(test, levels = methods),
      col_label = "Abundance stratum",
      abund_stratum = factor(abund_stratum, levels = c("Low", "Medium", "High"))
    )
  
  yline_col <- if (summary_stat == "q95") "q95_logp" else "median_logp"
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = prevalence, y = logp)) +
    ggplot2::geom_point(
      colour = "#9E9E9E",
      alpha = point_alpha,
      size = point_size
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(plot_df, detected),
      colour = "#BE5302",
      #alpha = 0.85,
      size = point_size + 0.5
    ) 
  if (show_envelope) {
    p <- p +
      ggplot2::geom_line(
        data = dplyr::filter(line_df, !is.na(.data[[yline_col]])),
        ggplot2::aes(x = prev_mid, y = .data[[yline_col]]),
        inherit.aes = FALSE,
        colour = "black",
        linewidth = 0.9
      )
  }
  p <- p +
    ggh4x::facet_nested(test ~ col_label + abund_stratum) +
    ggplot2::labs(
      x = "Taxon prevalence",
      y = expression(-log[10]("raw p-value"))
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white", colour = "black"),
      strip.text = ggplot2::element_text(size = 11, colour = "black"),
      axis.text = ggplot2::element_text(size = 10, colour = "black"),
      axis.title = ggplot2::element_text(size = 13, colour = "black"),
      panel.spacing = grid::unit(0.8, "lines")
    )
}

#' Plot abundance-stratified simulation clouds for detectability diagnostics
#'
#' Builds a supplementary-style simulation plot showing implanted marker taxa
#' as points in prevalence-effect space, colored by whether they were detected
#' by the DA method. This figure is useful for visually checking whether a
#' model-based detectability boundary separates detected from missed simulated
#' markers in a sensible way.
#'
#' The function expects a calibration table such as the one obtained by
#' \code{bind_rows(purrr::map(method_results, "calib_df"))}. It filters to
#' implanted markers, adds abundance strata using a fixed feature-level
#' abundance table, and plots the simulation cloud faceted by method and
#' abundance stratum.
#'
#' @param master_calib A data frame or tibble containing simulation calibration
#'   results across methods and replicates.
#' @param abundance_df Feature-level abundance table used to assign fixed
#'   abundance strata. Must contain \code{feature} and the chosen abundance
#'   column. This is typically derived from the real relative-abundance matrix.
#' @param methods Optional character vector of methods to plot. If NULL, all
#'   methods present in \code{master_calib} are used.
#' @param thresh_df Optional data frame containing model-based boundary lines to
#'   overlay. If supplied, it should contain at least \code{test},
#'   \code{abund_stratum}, \code{prev_mid}, and \code{threshold}.
#' @param prevalence_col Character scalar naming the prevalence column in
#'   \code{master_calib}.
#' @param effect_col Character scalar naming the absolute effect-size column in
#'   \code{master_calib}.
#' @param abundance_col Character scalar naming the continuous abundance column
#'   used to define abundance strata.
#' @param abundance_labels Character vector giving the desired abundance-stratum
#'   order.
#' @param point_alpha_missed Alpha for missed simulated markers.
#' @param point_alpha_detected Alpha for detected simulated markers.
#' @param point_size_missed Point size for missed simulated markers.
#' @param point_size_detected Point size for detected simulated markers.
#' @param threshold_linewidth Line width for the optional model-based boundary.
#' @param facet_scales Passed to \code{facet_grid()}.
#' @param threshold_label Legend label for the boundary line.
#' @param missed_label Legend label for missed simulated markers.
#' @param detected_label Legend label for detected simulated markers.
#'
#' @return A ggplot object.
plot_simulation_detectability_cloud_abundance <- function(
    master_calib,
    feat_real,
    methods = NULL,
    thresh_df = NULL,
    prevalence_col = "prevalence_full",
    effect_col = "abs_effect_size",
    abundance_col = "log10_abundance_nonzero_median_baseline",
    abundance_labels = c("Low", "Medium", "High"),
    point_alpha_missed = 0.28,
    point_alpha_detected = 0.38,
    point_size_missed = 0.8,
    point_size_detected = 0.7,
    threshold_linewidth = 0.7,
    facet_scales = "free_y",
    threshold_label = "Model-based 80% detectability boundary",
    missed_label = "Missed simulated markers",
    detected_label = "Detected simulated markers"
) {
  req_calib <- c("test", "feature", "is_marker", "detected", prevalence_col, effect_col)
  miss_calib <- setdiff(req_calib, colnames(master_calib))
  if (length(miss_calib) > 0) {
    stop("master_calib is missing required columns: ",
         paste(shQuote(miss_calib), collapse = ", "))
  }
  
  abundance_df <- compute_real_feature_stats(
    sweep(feat_real, 2, colSums(feat_real), FUN = "/")
  ) %>%
    dplyr::transmute(
      feature = feature,
      log10_abundance_nonzero_median_baseline = log10_abundance_nonzero_median
    )
  
  req_abund <- c("feature", abundance_col)
  miss_abund <- setdiff(req_abund, colnames(abundance_df))
  if (length(miss_abund) > 0) {
    stop("abundance_df is missing required columns: ",
         paste(shQuote(miss_abund), collapse = ", "))
  }
  
  if (is.null(methods)) {
    methods <- sort(unique(master_calib$test))
  }
  
  # Keep only implanted markers with complete coordinates for plotting.
  sim_df <- master_calib %>%
    dplyr::filter(
      test %in% methods,
      is_marker,
      !is.na(detected),
      !is.na(.data[[prevalence_col]]),
      !is.na(.data[[effect_col]])
    ) %>%
    dplyr::mutate(
      prevalence = as.numeric(.data[[prevalence_col]]),
      abs_effect_size = as.numeric(.data[[effect_col]])
    )
  
  # Add fixed abundance strata using the original feature-level abundance table.
  sim_df <- add_common_abundance_strata(
    df = sim_df,
    abundance_df = abundance_df,
    abundance_col = abundance_col,
    labels = abundance_labels
  )
  bad_features <- sim_df$feature[is.na(sim_df$abund_stratum)]
  if (length(bad_features) > 0) {
    print(filter(sim_df, is.na(abund_stratum)))
    stop(
      "Some features did not receive an abundance stratum. Examples: ",
      paste(utils::head(unique(bad_features), 10), collapse = ", ")
    )
  }
  
  sim_df$abund_stratum <- factor(sim_df$abund_stratum, levels = abundance_labels)
  sim_df$detected_status <- ifelse(sim_df$detected, detected_label, missed_label)
  sim_df$col_label <- "Abundance stratum"
  
  
  # Restrict and standardize the boundary table if provided.
  if (!is.null(thresh_df)) {
    req_thresh <- c("test", "abund_stratum", "prev_mid", "threshold")
    miss_thresh <- setdiff(req_thresh, colnames(thresh_df))
    if (length(miss_thresh) > 0) {
      stop("thresh_df is missing required columns: ",
           paste(shQuote(miss_thresh), collapse = ", "))
    }
    
    thresh_df <- thresh_df %>%
      dplyr::filter(test %in% methods) %>%
      dplyr::mutate(
        abund_stratum = factor(abund_stratum, levels = abundance_labels),
        col_label = "Abundance stratum"
      )
  }
  
  p <- ggplot2::ggplot() +
    # Missed simulated markers first, so detected points sit on top.
    ggplot2::geom_point(
      data = sim_df %>% dplyr::filter(!detected),
      ggplot2::aes(
        x = prevalence,
        y = abs_effect_size,
        colour = missed_label
      ),
      alpha = point_alpha_missed,
      size = point_size_missed
    ) +
    ggplot2::geom_point(
      data = sim_df %>% dplyr::filter(detected),
      ggplot2::aes(
        x = prevalence,
        y = abs_effect_size,
        colour = detected_label
      ),
      alpha = point_alpha_detected,
      size = point_size_detected
    )
  
  # Optionally overlay the model-based detectability boundary.
  if (!is.null(thresh_df)) {
    p <- p +
      ggplot2::geom_line(
        data = thresh_df,
        ggplot2::aes(
          x = prev_mid,
          y = threshold,
          linetype = threshold_label
        ),
        colour = "#7A1F5C",
        linewidth = threshold_linewidth,
        na.rm = TRUE
      ) +
      ggplot2::scale_linetype_manual(
        values = stats::setNames("solid", threshold_label),
        name = NULL
      )
  }
  
  p +
    #ggplot2::facet_grid(test ~ abund_stratum, scales = facet_scales) +
    ggh4x::facet_nested(test ~ col_label + abund_stratum, scales = facet_scales) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(
        c("#BFBFBF", "#4C78A8"),
        c(missed_label, detected_label)
      ),
      name = NULL
    ) +
    ggplot2::labs(
      x = "Taxon prevalence",
      y = "Absolute effect size"
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white", colour = "black"),
      strip.text = ggplot2::element_text(size = 11, colour = "black"),
      axis.text = ggplot2::element_text(size = 10, colour = "black"),
      axis.title = ggplot2::element_text(size = 13, colour = "black"),
      panel.spacing = grid::unit(0.8, "lines"),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
}

#' Plot unified detectability panels with simulation cloud, real positives, and two boundaries
#'
#' Creates a unified figure that overlays simulated detected/missed markers with
#' real positive taxa and two model-based detectability boundaries:
#' (i) the abundance-stratum boundary and
#' (ii) an exact-abundance boundary for each real positive taxon.
#'
#' This is especially useful for sparse real analyses such as relapse/recurrence,
#' where a single real positive taxon can be interpreted more clearly in the
#' context of the simulated detectability cloud.
#'
#' @param res_bound Output from \code{run_casecontrol_upper_bound_abundance()}.
#'   Must contain at least \code{real_df}, \code{thresh_df},
#'   \code{boundary_models}, and \code{abundance_ref_df}.
#' @param master_calib Calibration table used to build the simulation cloud.
#' @param methods Optional character vector of methods to plot. Defaults to all
#'   methods present in \code{res_bound$real_df}.
#' @param target_prob Detection probability contour to use for the exact-abundance
#'   boundary. Defaults to 0.8.
#' @param prevalence_col Character scalar naming the prevalence column in
#'   \code{master_calib}.
#' @param effect_col Character scalar naming the absolute effect-size column in
#'   \code{master_calib}.
#' @param abundance_col Character scalar naming the continuous abundance column in
#'   both \code{master_calib} and \code{res_bound$real_df}.
#' @param abundance_labels Character vector giving the abundance-stratum order.
#' @param effect_upper_quantile Upper quantile of observed simulated effect sizes
#'   used to define the inversion range for the exact-abundance boundary.
#' @param grid_length Number of effect-size grid points used for numerical
#'   inversion.
#' @param point_alpha_missed Alpha for missed simulated markers.
#' @param point_alpha_detected Alpha for detected simulated markers.
#' @param point_size_missed Point size for missed simulated markers.
#' @param point_size_detected Point size for detected simulated markers.
#' @param real_point_size Point size for real positive taxa.
#' @param stratum_linewidth Line width for the stratum-level boundary.
#' @param exact_linewidth Line width for the exact-abundance boundary.
#' @param facet_scales Passed to \code{facet_grid()}.
#' @param missed_label Legend label for missed simulated markers.
#' @param detected_label Legend label for detected simulated markers.
#' @param stratum_label Legend label for the stratum-level boundary.
#' @param exact_label Legend label for the exact-abundance boundary.
#' @param real_label Legend label for real positive taxa.
#'
#' @return A ggplot object.
plot_unified_detectability_abundance <- function(
    res_bound,
    master_calib,
    methods = NULL,
    target_prob = 0.8,
    prevalence_col = "prevalence_full",
    effect_col = "abs_effect_size",
    abundance_col = "log10_abundance_nonzero_median_baseline",
    abundance_labels = c("Low", "Medium", "High"),
    effect_upper_quantile = 0.995,
    grid_length = 500,
    point_alpha_missed = 0.45,
    point_alpha_detected = 0.35,
    point_size_missed = 0.7,
    point_size_detected = 0.7,
    real_point_size = 2.0,
    stratum_linewidth = 1.1,
    exact_linewidth = 0.9,
    facet_scales = "free_y",
    missed_label = "Undetected implanted markers (FN)",
    detected_label = "Detected implanted markers (TP)",
    stratum_label = "80% detectability boundary (stratum)",
    exact_label = "80% detectability boundary (meta_mOTU_v3_12389)",
    real_label = "Recurrence-associated species (meta_mOTU_v3_12389)"
) {
  stopifnot(is.list(res_bound))
  stopifnot(all(c("real_df", "thresh_df", "boundary_models", "abundance_ref_df") %in% names(res_bound)))
  
  real_df <- res_bound$real_df
  thresh_df <- res_bound$thresh_df
  boundary_models <- res_bound$boundary_models
  abundance_ref_df <- res_bound$abundance_ref_df
  
  if (is.null(methods)) {
    methods <- sort(unique(real_df$test))
  }
  
  req_calib <- c("test", "feature", "is_marker", "detected", prevalence_col, effect_col)
  miss_calib <- setdiff(req_calib, colnames(master_calib))
  if (length(miss_calib) > 0) {
    stop("master_calib is missing required columns: ",
         paste(shQuote(miss_calib), collapse = ", "))
  }
  
  req_real <- c("feature", "test", "prevalence", "abs_effect_size", "detected",
                "abund_stratum", abundance_col)
  miss_real <- setdiff(req_real, colnames(real_df))
  if (length(miss_real) > 0) {
    stop("res_bound$real_df is missing required columns: ",
         paste(shQuote(miss_real), collapse = ", "))
  }
  
  # ------------------------------------------------------------
  # 1. Prepare simulation cloud
  # ------------------------------------------------------------
  sim_df <- master_calib %>%
    dplyr::filter(
      test %in% methods,
      is_marker,
      !is.na(detected),
      !is.na(.data[[prevalence_col]]),
      !is.na(.data[[effect_col]])
    ) %>%
    dplyr::mutate(
      prevalence = as.numeric(.data[[prevalence_col]]),
      abs_effect_size = as.numeric(.data[[effect_col]])
    )
  
  # Use abundance information already present in master_calib if available,
  # otherwise borrow the feature-level abundance reference from real_df.
  if (!(abundance_col %in% colnames(sim_df))) {
    abundance_lookup <- real_df %>%
      dplyr::distinct(feature, .data[[abundance_col]])
    sim_df <- sim_df %>%
      dplyr::left_join(abundance_lookup, by = "feature")
  }
  
  sim_df <- add_common_abundance_strata(
    df = sim_df,
    abundance_df = real_df %>% dplyr::distinct(feature, .data[[abundance_col]]),
    abundance_col = abundance_col,
    labels = abundance_labels
  )
  
  sim_df$abund_stratum <- factor(sim_df$abund_stratum, levels = abundance_labels)
  sim_df$col_label <- "Abundance stratum"
  
  # ------------------------------------------------------------
  # 2. Prepare real positive taxa
  # ------------------------------------------------------------
  real_pos_df <- real_df %>%
    dplyr::filter(
      test %in% methods,
      detected
    ) %>%
    dplyr::mutate(
      abund_stratum = factor(abund_stratum, levels = abundance_labels),
      col_label = "Abundance stratum"
    )
  
  # ------------------------------------------------------------
  # 3. Keep stratum-level boundary
  # ------------------------------------------------------------
  stratum_thresh_df <- thresh_df %>%
    dplyr::filter(test %in% methods) %>%
    dplyr::mutate(
      abund_stratum = factor(abund_stratum, levels = abundance_labels),
      col_label = "Abundance stratum"
    )
  
  # ------------------------------------------------------------
  # 4. Build exact-abundance boundary for each real positive taxon
  # ------------------------------------------------------------
  exact_thresh_df <- purrr::pmap_dfr(
    list(
      feature = real_pos_df$feature,
      test = real_pos_df$test,
      abundance_value = real_pos_df[[abundance_col]],
      abund_stratum = real_pos_df$abund_stratum
    ),
    function(feature, test, abundance_value, abund_stratum) {
      fit <- boundary_models[[test]]
      if (is.null(fit) || is.na(abundance_value) || !is.finite(abundance_value)) {
        return(NULL)
      }
      
      dfm <- master_calib %>%
        dplyr::filter(
          test == test,
          !is.na(.data[[effect_col]]),
          is.finite(.data[[effect_col]])
        )
      
      if (nrow(dfm) == 0) {
        return(NULL)
      }
      
      eff_max <- stats::quantile(dfm[[effect_col]], probs = effect_upper_quantile, na.rm = TRUE)
      eff_max <- max(eff_max, max(dfm[[effect_col]], na.rm = TRUE), 1e-8)
      eff_grid <- seq(0, eff_max, length.out = grid_length)
      
      purrr::map_dfr(seq_len(nrow(stratum_thresh_df %>% dplyr::filter(test == test, abund_stratum == abund_stratum))), function(.i) {
        NULL
      })
      
      prev_grid <- sort(unique(stratum_thresh_df$prev_mid[stratum_thresh_df$test == test &
                                                            stratum_thresh_df$abund_stratum == abund_stratum]))
      
      if (length(prev_grid) == 0) {
        return(NULL)
      }
      
      purrr::map_dfr(prev_grid, function(p0) {
        newdata <- data.frame(
          effect = eff_grid,
          prevalence = p0,
          abundance = abundance_value
        )
        
        pred <- stats::predict(fit, newdata = newdata, type = "response")
        idx <- which(pred >= target_prob)
        
        tibble::tibble(
          feature = feature,
          test = test,
          abund_stratum = abund_stratum,
          col_label = "Abundance stratum",
          prev_mid = p0,
          target_prob = target_prob,
          threshold = if (length(idx) == 0) NA_real_ else eff_grid[min(idx)]
        )
      })
    }
  )
  
  exact_thresh_df$abund_stratum <- factor(exact_thresh_df$abund_stratum, levels = abundance_labels)
  
  # ------------------------------------------------------------
  # 5. Plot
  # ------------------------------------------------------------
  ggplot2::ggplot() +
    # Missed simulated markers first
    ggplot2::geom_point(
      data = sim_df %>% dplyr::filter(!detected),
      ggplot2::aes(
        x = prevalence,
        y = abs_effect_size,
        colour = missed_label
      ),
      alpha = point_alpha_missed,
      size = point_size_missed
    ) +
    # Detected simulated markers second
    ggplot2::geom_point(
      data = sim_df %>% dplyr::filter(detected),
      ggplot2::aes(
        x = prevalence,
        y = abs_effect_size,
        colour = detected_label
      ),
      alpha = point_alpha_detected,
      size = point_size_detected
    ) +
    # Stratum-level boundary
    ggplot2::geom_line(
      data = stratum_thresh_df,
      ggplot2::aes(
        x = prev_mid,
        y = threshold,
        linetype = stratum_label
      ),
      colour = "#8E1F63",
      linewidth = stratum_linewidth,
      na.rm = TRUE
    ) +
    # Exact-abundance boundary
    ggplot2::geom_line(
      data = exact_thresh_df,
      ggplot2::aes(
        x = prev_mid,
        y = threshold,
        group = interaction(test, abund_stratum, feature),
        linetype = exact_label
      ),
      colour = "#CC6D00",
      linewidth = exact_linewidth,
      na.rm = TRUE
    ) +
    # Real positive taxa
    ggplot2::geom_point(
      data = real_pos_df,
      ggplot2::aes(
        x = prevalence,
        y = abs_effect_size,
        colour = real_label
      ),
      size = real_point_size
    ) +
    ggh4x::facet_nested(test ~ col_label + abund_stratum, scales = facet_scales) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(
        c("#C9C9C9", "#B3CDE3", "#BE5302"),
        c(missed_label, detected_label, real_label)
      ),
      name = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = stats::setNames(
        c("solid", "22"),
        c(stratum_label, exact_label)
      ),
      name = NULL
    ) +
    ggplot2::labs(
      x = "Taxon prevalence",
      y = "Absolute effect size"
    ) +
    guides(
      linetype = guide_legend(
        order = 1,
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          linewidth = c(1.1, 0.9)
          )
        ),
      colour = guide_legend(
        order = 2,
        nrow = 2,
        byrow = TRUE,
        override.aes = list(
          size = c(2, 3, 2),
          alpha = c(1, 1, 1)
        )
      )
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white", colour = "black"),
      strip.text = ggplot2::element_text(size = 11, colour = "black"),
      axis.text = ggplot2::element_text(size = 10, colour = "black"),
      axis.title = ggplot2::element_text(size = 13, colour = "black"),
      panel.spacing = grid::unit(0.8, "lines"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.spacing.y = grid::unit(0.2, "lines"),
      legend.box.spacing = grid::unit(0.6, "lines"),
      legend.margin = ggplot2::margin(t = 5, r = 0, b = 0, l = 0),
      legend.key.height = grid::unit(0.45, "lines"),
      panel.grid.minor = ggplot2::element_blank()
    )
  
}

# ------------------------------------------------------------
# 5. Top-level wrapper for abundance-stratified upper-bound analysis
# ------------------------------------------------------------

#' Build abundance-stratified detectability layers using a model-based boundary
#'
#' This is a drop-in replacement for the previous
#' \code{run_casecontrol_upper_bound_abundance()} workflow. It keeps the same
#' overall outputs for plotting but replaces the local simulation threshold with
#' a method-specific model-derived detectability boundary.
#'
#' The returned object contains:
#' \itemize{
#'   \item \code{real_df}: standardized real-data effect table with abundance strata
#'   \item \code{envelope_df}: local observed effect-size envelope in real data
#'   \item \code{thresh_df}: model-based detectability boundary
#'   \item \code{boundary_models}: fitted detectability models by method
#'   \item \code{abundance_ref_df}: representative abundance values per stratum
#' }
#'
#' @param master_calib Long calibration table built from
#'   \code{bind_rows(purrr::map(method_results, "calib_df"))}.
#' @param feat_real Taxa x samples relative-abundance matrix.
#' @param methods Character vector of methods to include.
#' @param real_results_list Named list of real-data DA results by method.
#' @param prevalence Named numeric vector of taxon prevalence in the real data.
#'   If NULL, prevalence is computed from \code{feat_real}.
#' @param target_prob Target detection probability for the model-derived
#'   boundary. Defaults to 0.8.
#' @param adjust Multiple-testing adjustment applied to real results before
#'   flagging real detections.
#' @param envelope_stat Character scalar naming which real-data envelope column
#'   should be emphasized later in plotting. Defaults to \code{"q95"}.
#' @param prev_grid Numeric vector of prevalence grid values used both for the
#'   model-based boundary and the local real-data envelope.
#' @param envelope_bandwidth Half-width of the local prevalence window used for
#'   the real-data envelope.
#' @param alpha Significance threshold used for the real-data detection flag.
#' @param abundance_col Character scalar naming the continuous abundance column
#'   used in the simulation detectability model.
#' @param abundance_col_real_name Character scalar naming the abundance column
#'   to create from \code{feat_real}. This is renamed to \code{abundance_col}
#'   so that the same fixed abundance tertiles are used in both real and
#'   simulated data.
#' @param abundance_labels Character vector of labels for abundance strata.
#' @param min_rows_model Minimum number of simulated marker rows required to fit
#'   a detectability model for a method.
#' @param effect_upper_quantile Upper simulated effect-size quantile used to
#'   define the numerical inversion range for the model-based boundary.
#'
#' @return A list with components \code{real_df}, \code{envelope_df},
#'   \code{thresh_df}, \code{boundary_models}, and \code{abundance_ref_df}.
run_casecontrol_upper_bound_abundance <- function(
    master_calib,
    feat_real,
    methods,
    real_results_list,
    prevalence = NULL,
    target_prob = 0.8,
    adjust = "none",
    envelope_stat = "q95",
    prev_grid = seq(0.10, 0.95, by = 0.05),
    envelope_bandwidth = 0.08,
    alpha = 0.05,
    abundance_col = "log10_abundance_nonzero_median_baseline",
    abundance_col_real_name = "log10_abundance_nonzero_median",
    abundance_labels = c("Low", "Medium", "High"),
    min_rows_model = 50,
    effect_upper_quantile = 0.995
) {
  stopifnot(is.data.frame(feat_real))
  stopifnot(!is.null(rownames(feat_real)))
  
  # Compute real-data feature-level prevalence and abundance from the relative
  # abundance matrix. This provides the fixed abundance strata used for display.
  real_feat_stats <- compute_real_feature_stats(feat_real)
  
  # Align the real-data abundance column name with the simulation-side name so
  # we can use one common stratum definition across real and simulated data.
  abundance_df <- real_feat_stats %>%
    dplyr::transmute(
      feature = feature,
      !!abundance_col := .data[[abundance_col_real_name]],
      prevalence = prevalence
    )
  
  # If prevalence was not supplied, derive it from the real feature table.
  if (is.null(prevalence)) {
    prevalence <- stats::setNames(real_feat_stats$prevalence, real_feat_stats$feature)
  }
  
  # Build a standardized real-data table across methods.
  real_df <- build_real_effect_df(
    real_results_list = real_results_list,
    prevalence = prevalence,
    methods = methods,
    adjust = adjust,
    alpha = alpha
  )
  
  # Add fixed abundance strata to the real-data table using the original
  # feature-level abundance distribution.
  real_df <- add_common_abundance_strata(
    df = real_df,
    abundance_df = abundance_df,
    abundance_col = abundance_col,
    labels = abundance_labels
  )
  
  real_df$abund_stratum <- factor(real_df$abund_stratum, levels = abundance_labels)
  
  # Add the same fixed abundance strata to the calibration table so that model
  # fitting and plotting use the same low/medium/high stratification.
  calib_df <- add_common_abundance_strata(
    df = master_calib,
    abundance_df = abundance_df,
    abundance_col = abundance_col,
    labels = abundance_labels
  )
  
  calib_df$abund_stratum <- factor(calib_df$abund_stratum, levels = abundance_labels)
  
  # Fit one additive detectability model per method using continuous abundance.
  boundary_models <- fit_detectability_boundary_models_abundance(
    calib_df = calib_df,
    methods = methods,
    prevalence_col = "prevalence_full",
    abundance_col = abundance_col,
    min_rows = min_rows_model,
    #formula = detected_bin ~ effect + prevalence + abundance + effect:prevalence + effect:abundance
    formula = detected_bin ~ te(effect, prevalence, k = c(3, 4)) + s(abundance, k = 4, bs = "ts")
  )
  
  # Use the median abundance within each display stratum as the representative
  # continuous abundance value for drawing one smooth boundary per panel.
  abundance_ref_df <- compute_abundance_reference_values(
    abundance_df = add_common_abundance_strata(
      df = abundance_df,
      abundance_df = abundance_df,
      abundance_col = abundance_col,
      labels = abundance_labels
    ),
    abundance_col = abundance_col,
    abundance_stratum_col = "abund_stratum",
    summary_fun = stats::median
  )
  
  abundance_ref_df$abund_stratum <- factor(abundance_ref_df$abund_stratum, levels = abundance_labels)
  
  # Predict the smooth model-based detectability boundary.
  thresh_df <- predict_detectability_boundary_abundance(
    model_list = boundary_models,
    calib_df = calib_df,
    abundance_ref_df = abundance_ref_df,
    target_prob = target_prob,
    prev_grid = prev_grid,
    effect_col = "abs_effect_size",
    effect_upper_quantile = effect_upper_quantile
  )
  
  thresh_df$abund_stratum <- factor(thresh_df$abund_stratum, levels = abundance_labels)
  
  # Keep the real-data local envelope exactly as before.
  envelope_df <- summarise_real_effect_envelope_local_abundance(
    real_df = real_df,
    prev_grid = prev_grid,
    bandwidth = envelope_bandwidth,
    abundance_stratum_col = "abund_stratum"
  )
  
  envelope_df$abund_stratum <- factor(envelope_df$abund_stratum, levels = abundance_labels)
  envelope_df$envelope_stat <- envelope_stat
  
  list(
    real_df = real_df,
    envelope_df = envelope_df,
    thresh_df = thresh_df,
    boundary_models = boundary_models,
    abundance_ref_df = abundance_ref_df
  )
}
