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
      #"log10_abundance_nonzero_median_full",
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


# ------------------------------------------------------------
# 4. Plotting functions
# ------------------------------------------------------------

#' Plot abundance-stratified effect-size upper bounds
#'
#' Creates a facet-grid plot with methods in rows and abundance strata in columns.
#' Grey points show all real taxa, red points show detected taxa (if any), the
#' blue line shows the simulation-derived detection threshold, and the black line
#' shows the observed real effect-size envelope.
#'
#' @param real_df Real-data result table.
#' @param thresh_df Threshold table from \code{fit_detection_thresholds_local_abundance()}.
#' @param envelope_df Envelope table from \code{summarise_real_effect_envelope_local_abundance()}.
#' @param methods Character vector of methods to include and their display order.
#' @param target_prob Detection probability to plot (for example 0.8).
#' @param envelope_stat Which envelope statistic to draw: one of \code{"q95"},
#'   \code{"q90"}, \code{"q99"}, or \code{"q50"}.
#'
#' @return A ggplot object.
#'
plot_effect_size_upper_bounds_abundance <- function(real_df,
                                                    thresh_df,
                                                    envelope_df,
                                                    methods,
                                                    target_prob = 0.8,
                                                    envelope_stat = "q95") {
  thresh_use <- thresh_df %>%
    dplyr::filter(target_prob == !!target_prob, !is.na(threshold)) %>%
    dplyr::mutate(
      #test = factor(test, levels = methods),
      col_label = "Abundance stratum",
      abund_stratum = factor(abund_stratum, levels = c("Low", "Medium", "High"))
    )
  
  env_use <- envelope_df %>%
    dplyr::filter(!is.na(.data[[envelope_stat]])) %>%
    dplyr::mutate(
      test = factor(test, levels = methods),
      col_label = "Abundance stratum",
      abund_stratum = factor(abund_stratum, levels = c("Low", "Medium", "High"))
    )
  
  real_use <- real_df %>%
    dplyr::filter(test %in% methods) %>%
    dplyr::mutate(
      test = factor(test, levels = methods),
      col_label = "Abundance stratum",
      abund_stratum = factor(abund_stratum, levels = c("Low", "Medium", "High"))
    )
  
  ggplot2::ggplot() +
    # All real taxa
    ggplot2::geom_point(
      data = real_use,
      ggplot2::aes(x = prevalence, y = abs_effect_size),
      colour = "grey45",
      alpha = 0.25,
      size = 1.0
    ) +
    # Real taxa that passed the chosen significance threshold, if any
    ggplot2::geom_point(
      data = dplyr::filter(real_use, detected),
      ggplot2::aes(x = prevalence, y = abs_effect_size),
      colour = "red3",
      alpha = 0.85,
      size = 1.5
    ) +
    # Simulation-derived threshold
    ggplot2::geom_line(
      data = thresh_use,
      ggplot2::aes(x = prev_mid, y = threshold),
      colour = "#2C6BEA",
      linetype = "31",
      linewidth = 1.0
    ) +
    ggplot2::geom_point(
      data = thresh_use,
      ggplot2::aes(x = prev_mid, y = threshold),
      colour = "#2C6BEA",
      size = 0.9
    ) +
    # Observed real-data envelope
    ggplot2::geom_line(
      data = env_use,
      ggplot2::aes(x = prev_mid, y = .data[[envelope_stat]]),
      colour = "black",
      linetype = "solid",
      linewidth = 1.1
    ) +
    ggplot2::geom_point(
      data = env_use,
      ggplot2::aes(x = prev_mid, y = .data[[envelope_stat]]),
      colour = "black",
      size = 1.2
    ) +
    ggh4x::facet_nested(test ~ col_label + abund_stratum, scales = "free_y") + 
    ggplot2::theme(strip.placement = "outside") +
    ggplot2::labs(
      # title = paste0(round(target_prob * 100), "% detection thresholds stratified by baseline abundance"),
      # subtitle = paste0(
      #   "Blue: simulation-derived threshold; black: observed real effect-size envelope (",
      #   envelope_stat, ")"
      # ),
      x = "Taxon prevalence",
      y = "Absolute estimated effect size"
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white", colour = "black"),
      strip.text = ggplot2::element_text(size = 11, colour = "black"),
      axis.text = ggplot2::element_text(size = 10, colour = "black"),
      axis.title = ggplot2::element_text(size = 13, colour = "black"),
      plot.title = ggplot2::element_text(size = 15, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 11),
      panel.spacing = grid::unit(0.8, "lines")
    )
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
                                                    point_alpha = 0.25,
                                                    point_size = 1.2) {
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
  
  ggplot2::ggplot(plot_df, ggplot2::aes(x = prevalence, y = logp)) +
    ggplot2::geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::geom_point(
      colour = "grey45",
      alpha = point_alpha,
      size = point_size
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(plot_df, detected),
      colour = "red3",
      alpha = 0.85,
      size = point_size + 0.3
    ) +
    ggplot2::geom_line(
      data = dplyr::filter(line_df, !is.na(.data[[yline_col]])),
      ggplot2::aes(x = prev_mid, y = .data[[yline_col]]),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.9
    ) +
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


# ------------------------------------------------------------
# 5. Top-level wrapper for abundance-stratified upper-bound analysis
# ------------------------------------------------------------

#' Run abundance-stratified upper-bound analysis for a case/control result
#'
#' This wrapper:
#'   \enumerate{
#'     \item standardizes real-method outputs
#'     \item ensures that the calibration table contains abundance strata
#'     \item computes local abundance-stratified detectability thresholds
#'     \item computes the local observed real effect-size envelope
#'     \item returns both intermediate tables and a ready-to-plot figure
#'   }
#'
#' @param master_calib Long-format calibration table created from the per-method
#'   calibration outputs. Must contain, at minimum, \code{test}, \code{feature},
#'   \code{is_marker}, \code{detected}, \code{abs_effect_size}, and a prevalence
#'   column such as \code{prevalence_full}. It must also contain either
#'   \code{abund_stratum} or a baseline abundance column that can be turned into
#'   common tertiles.
#' @param feat_real Taxa x samples relative abundance matrix.
#' @param methods Character vector of methods to include and their display order.
#' @param real_results_list Named list of real-method results.
#' @param prevalence Optional named numeric vector of taxon prevalence. If
#'   \code{NULL}, prevalence is recomputed from \code{feat_real}.
#' @param target_prob Detection probability to plot in the final upper-bound figure.
#' @param probs Numeric vector of target detection probabilities to estimate.
#' @param prevalence_col Character scalar naming the prevalence column in
#'   \code{master_calib}.
#' @param prev_grid Numeric vector of prevalence grid values.
#' @param bandwidth Numeric half-width of the local prevalence window.
#' @param min_n Minimum number of taxa required in a local window.
#' @param min_detected Minimum number of detected taxa required in a local window.
#' @param adjust Multiple-testing adjustment method used only to define the
#'   \code{detected} flag in the real-data table. Use \code{"none"} for nominal
#'   significance.
#' @param envelope_stat Which local effect-size envelope statistic to draw.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{calib_df}: calibration table used for fitting
#'     \item \code{threshold_df}: local detectability thresholds
#'     \item \code{real_df}: standardized real-data result table
#'     \item \code{envelope_df}: local observed real effect-size envelope
#'     \item \code{plot}: ggplot object
#'   }
#'
#' @export
run_casecontrol_upper_bound_abundance <- function(master_calib,
                                                  feat_real,
                                                  methods,
                                                  real_results_list,
                                                  prevalence = NULL,
                                                  target_prob = 0.8,
                                                  probs = c(0.5, 0.8, 0.95),
                                                  prevalence_col = "prevalence_full",
                                                  prev_grid = seq(0.10, 0.95, by = 0.05),
                                                  bandwidth = 0.08,
                                                  min_n = 50,
                                                  min_detected = 10,
                                                  adjust = "fdr",
                                                  alpha = 0.05,
                                                  envelope_stat = "q95") {

  # Use supplied prevalence if available; otherwise compute from feat_real.
  if (is.null(prevalence)) {
    feature_stats <- compute_real_feature_stats(feat_real)
    prevalence <- stats::setNames(feature_stats$prevalence, feature_stats$feature)
  }
  
  # Standardize the real data results across methods.
  real_df <- build_real_effect_df(
    real_results_list = real_results_list,
    prevalence = prevalence,
    methods = methods,
    adjust = adjust,
    alpha = alpha
  )

  # Join abundance statistics from the real relative-abundance matrix.
  feature_stats <- compute_real_feature_stats(feat_real)
  real_df <- real_df %>%
    dplyr::left_join(feature_stats, by = c("feature", "prevalence")) %>%
    add_abundance_tertiles(abundance_col = "log10_abundance_nonzero_median") %>%
    dplyr::mutate(
      test = factor(test, levels = methods),
      abund_stratum = factor(abund_stratum, levels = c("Low", "Medium", "High"))
    )
  
  # Make sure the calibration table uses common abundance tertiles.
  calib_df <- master_calib %>%
    dplyr::filter(test %in% methods) %>%
    ensure_calibration_abundance_strata() %>%
    dplyr::mutate(
      test = factor(test, levels = methods),
      abund_stratum = factor(abund_stratum, levels = c("Low", "Medium", "High"))
    )

  # Estimate simulation-derived local thresholds.
  threshold_df <- fit_detection_thresholds_local_abundance_empirical(
    calib_df = calib_df,
    probs = probs,
    prev_grid = prev_grid,
    prevalence_col = prevalence_col,
    bandwidth = bandwidth,
    abundance_stratum_col = "abund_stratum",
    min_n = min_n,
    min_detected = min_detected
  ) %>%
    dplyr::mutate(
      test = factor(test, levels = methods),
      abund_stratum = factor(abund_stratum, levels = c("Low", "Medium", "High"))
    )
  
  # Summarize the observed real-data effect envelope.
  envelope_df <- summarise_real_effect_envelope_local_abundance(
    real_df = real_df,
    prev_grid = prev_grid,
    bandwidth = bandwidth,
    abundance_stratum_col = "abund_stratum",
    min_n = 30
  ) 
  envelope_df <- envelope_df %>%
    dplyr::mutate(
      test = factor(test, levels = methods),
      abund_stratum = factor(abund_stratum, levels = c("Low", "Medium", "High"))
    )
  
  p <- plot_effect_size_upper_bounds_abundance(
    real_df = real_df,
    thresh_df = threshold_df,
    envelope_df = envelope_df,
    methods = methods,
    target_prob = target_prob,
    envelope_stat = envelope_stat
  )
  
  list(
    calib_df = calib_df,
    threshold_df = threshold_df,
    real_df = real_df,
    envelope_df = envelope_df,
    plot = p
  )
}