# ============================================================
# SIMBA group-name parsing helper
# ============================================================

#' Parse SIMBA group names into simulation metadata
#'
#' Extracts the abundance-scale index, prevalence-scale index, repetition index
#' and corresponding \code{ab.scale} value from group names of the form
#' \code{"ab2_prev1_rep7"}.
#'
#' @param group_name Character scalar giving the SIMBA group name.
#' @param ab.scale Numeric vector of abundance-scale values used in the benchmark.
#'
#' @return A one-row tibble with columns:
#'   \itemize{
#'     \item \code{group}
#'     \item \code{ab_index}
#'     \item \code{prev_index}
#'     \item \code{rep_idx}
#'     \item \code{ab_scale}
#'   }
#'
#' @export
#' 
#' @examples
#' # parse_group_info("ab2_prev1_rep7", ab.scale = c(1, 2, 3))
parse_group_info <- function(group_name, ab.scale) {
  # Match group names like "ab2_prev1_rep7".
  m <- stringr::str_match(group_name, "^ab([0-9]+)_prev([0-9]+)_rep([0-9]+)$")
  if (any(is.na(m))) {
    stop("Could not parse group name: ", group_name)
  }
  
  # Convert captured substrings to integer indices.
  ab_idx   <- as.integer(m[, 2])
  prev_idx <- as.integer(m[, 3])
  rep_idx  <- as.integer(m[, 4])
  
  # Check that the abundance index is valid for the supplied ab.scale vector.
  if (ab_idx < 1 || ab_idx > length(ab.scale)) {
    stop("ab index out of range for group: ", group_name)
  }
  
  tibble::tibble(
    group = group_name,
    ab_index = ab_idx,
    prev_index = prev_idx,
    rep_idx = rep_idx,
    ab_scale = ab.scale[ab_idx]
  )
}


# ============================================================
# Calibration-data extraction from apply.test()
# ============================================================

#' Extract long-format calibration data from a patched SIMBA apply.test() result
#'
#' Converts the nested output of a patched \code{apply.test()} call into a
#' long-format tibble containing p-values, adjusted p-values, effect sizes,
#' marker status and prevalence summaries for each simulated replicate.
#'
#' This function assumes that:
#' \itemize{
#'   \item \code{apply.test()} returns a list with components \code{"p.val"} and
#'         optionally \code{"eff.size"}
#'   \item \code{sim.location} is an HDF5 file containing the simulated group
#'         stored at \code{group}
#'   \item the global objects \code{ab.scale} and \code{alpha} exist, as used in
#'         the surrounding benchmarking workflow
#' }
#'
#' @param apply_res Output from a patched \code{apply.test()} call.
#' @param sim.location Character scalar giving the HDF5 file path.
#' @param group Character scalar naming the simulated group, for example
#'   \code{"ab2_prev1_rep7"}.
#' @param subset_size Integer scalar giving the subset size used in
#'   \code{apply.test()}.
#' @param adjust Character scalar giving the p-value adjustment method passed to
#'   \code{p.adjust()}. Use \code{"pass"} or \code{"PASS"} to skip adjustment.
#'
#' @return A tibble with one row per taxon per replicate column, containing:
#'   p-values, adjusted p-values, effect sizes, marker status, prevalence
#'   summaries, and parsed simulation metadata.
extract_apply_test_long <- function(apply_res,
                                    sim.location,
                                    group,
                                    subset_size,
                                    adjust = "BH",
                                    compute_rep_prevalence = FALSE) {
  # Construct the subset label used inside SIMBA objects.
  subset_name <- paste0("subset_", subset_size)
  
  # Unwrap p-value results. Depending on the patched return structure,
  # an extra list level keyed by group name may be present.
  pval_res <- apply_res[["p.val"]]
  if (is.list(pval_res) && length(pval_res) == 1 && group %in% names(pval_res)) {
    pval_res <- pval_res[[group]]
  }
  pval_mat <- pval_res[[subset_name]]
  
  # Unwrap effect-size results if available.
  eff_res <- apply_res[["eff.size"]]
  if (is.list(eff_res) && length(eff_res) == 1 && group %in% names(eff_res)) {
    eff_res <- eff_res[[group]]
  }
  eff_mat <- NULL
  if (!is.null(eff_res) && length(eff_res) > 0 && subset_name %in% names(eff_res)) {
    eff_mat <- eff_res[[subset_name]]
  }
  
  # Sanity checks on p-value/effect-size matrices.
  if (!is.matrix(pval_mat)) {
    stop("pval_mat must be a matrix")
  }
  if (!is.null(eff_mat)) {
    stopifnot(is.matrix(eff_mat))
    stopifnot(identical(rownames(eff_mat), rownames(pval_mat)))
    stopifnot(identical(colnames(eff_mat), colnames(pval_mat)))
  }
  
  # Read the simulated group from HDF5 and rebuild the feature matrix.
  grp <- rhdf5::h5read(sim.location, group)
  
  feat.sim <- grp$features
  rownames(feat.sim) <- grp$feature_names
  colnames(feat.sim) <- grp$sample_names
  
  # Marker indices/names used when the signal was implanted.
  marker_idx <- grp$marker_idx
  
  # Prevalence over the full simulated group.
  prevalence_full <- rowMeans(feat.sim > 0)
  
  # Pull replicate sample selections if available.
  # These are only needed when replicate-specific prevalence should be computed.
  idx_info <- NULL
  if ("test_idx" %in% names(grp) && subset_name %in% names(grp$test_idx)) {
    idx_info <- grp$test_idx[[subset_name]]
  }
  
  # Parse group-level simulation metadata.
  group_info <- parse_group_info(group_name = group, ab.scale = ab.scale)
  
  # Convert each replicate column to long format.
  out <- purrr::map_dfr(seq_len(ncol(pval_mat)), function(j) {
    rep_col <- colnames(pval_mat)[j]
    
    p <- pval_mat[, j]
    if (!adjust %in% c("pass", "PASS")) {
      padj <- stats::p.adjust(p, method = adjust)
    } else {
      padj <- p
    }
    
    # Effect sizes may be absent for some methods.
    eff <- if (!is.null(eff_mat)) eff_mat[, j] else rep(NA_real_, length(p))
    names(eff) <- names(p) <- rownames(pval_mat)
    
    # Detection status is based on the global alpha used in the workflow.
    detected <- padj < alpha
    is_marker <- rownames(pval_mat) %in% marker_idx
    
    tibble::tibble(
      group = group,
      rep = rep_col,
      feature = rownames(pval_mat),
      pval = as.numeric(p),
      padj = as.numeric(padj),
      is_marker = as.logical(is_marker),
      effect_size = as.numeric(eff),
      abs_effect_size = abs(as.numeric(eff)),
      prevalence_full = as.numeric(prevalence_full[rownames(pval_mat)]),
      ab_index = group_info$ab_index,
      prev_index = group_info$prev_index,
      rep_idx = group_info$rep_idx,
      ab_scale = group_info$ab_scale
    )
  })
  
  out
}

# ============================================================
# Parallelized implementation of SIMBA on multiple methods
# ============================================================

#' Run one method on one simulated group and optionally extract calibration data
#'
#' Executes a single \code{test_name} on a single simulated group, evaluates
#' benchmark performance via \code{eval.test()}, and optionally extracts a
#' long-format calibration table from the raw \code{apply.test()} output.
#'
#' @param g Character scalar naming the simulated group.
#' @param test_name Character scalar naming the differential abundance method.
#' @param sim.implant Character scalar giving the HDF5 file path.
#' @param subset_size Integer scalar passed to \code{apply.test()}.
#' @param ab.scale Numeric vector of abundance-scale values used in the benchmark.
#' @param alpha Numeric scalar giving the significance threshold used in
#'   \code{eval.test()}.
#' @param covariates Optional character vector of metadata columns passed through
#'   the patched SIMBA workflow.
#' @param norm Character scalar naming the normalization option passed to
#'   \code{apply.test()}.
#' @param simba.loc Optional file path to the local SIMBA fork; if supplied,
#'   \code{devtools::load_all()} is called inside the worker.
#' @param return_calib Logical; if \code{TRUE}, also return a calibration table
#'   produced by \code{extract_apply_test_long()}.
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{eval}: one-row or few-row summary from \code{eval.test()}
#'     \item \code{calib}: long-format calibration tibble, or \code{NULL}
#'   }
run_one_group_eval <- function(g, test_name, sim.implant, subset_size,
                               ab.scale, alpha = 0.05, covariates = NULL,
                               norm = "pass", simba.loc = NULL,
                               return_calib = TRUE) {
  # Load the local SIMBA fork inside the worker if requested.
  if (!is.null(simba.loc)) {
    devtools::load_all(simba.loc, quiet = TRUE)
  }
  
  subset_name <- paste0("subset_", subset_size)
  
  # Run the selected DA method on the specified simulated group.
  res <- apply.test(
    sim.location = sim.implant,
    group = g,
    subset = subset_size,
    norm = norm,
    test = test_name,
    covariates = covariates
  )
  
  # Unwrap p-values and effect sizes from the patched SIMBA return structure.
  pval_res <- res[["p.val"]]
  eff_res  <- res[["eff.size"]]
  
  if (is.list(pval_res) && length(pval_res) == 1 && g %in% names(pval_res)) {
    pval_res <- pval_res[[g]]
  }
  if (is.list(eff_res) && length(eff_res) == 1 && g %in% names(eff_res)) {
    eff_res <- eff_res[[g]]
  }
  
  pval_mat <- pval_res[[subset_name]]
  eff_mat <- NULL
  if (!is.null(eff_res) && length(eff_res) > 0 && subset_name %in% names(eff_res)) {
    eff_mat <- eff_res[[subset_name]]
  }
  
  # Evaluate benchmark performance using the known implanted markers.
  ev <- eval.test(
    sim.location = sim.implant,
    group = g,
    res.mat = pval_mat,
    eff.mat = eff_mat,
    alpha = alpha
  ) %>%
    dplyr::bind_cols(parse_group_info(g, ab.scale)) %>%
    dplyr::mutate(
      test = test_name,
      power = R,
      empirical_fdr = FDR,
      precision = PR,
      recall = R
    )
  
  # Optionally extract long-format calibration data for downstream plotting.
  calib <- NULL
  if (isTRUE(return_calib)) {
    calib <- extract_apply_test_long(
      apply_res = res,
      sim.location = sim.implant,
      group = g,
      subset_size = subset_size,
      adjust = "BH",
      compute_rep_prevalence = FALSE
    ) %>%
      dplyr::mutate(
        test = test_name,
        norm = norm
      )
  }
  
  list(
    eval = ev,
    calib = calib
  )
}


#' Run one differential abundance method across multiple simulated groups in parallel
#'
#' Parallelizes over simulated groups for a single method, while keeping methods
#' sequential at the calling level. Returns both the evaluation summaries and,
#' optionally, the combined calibration table.
#'
#' @param test_name Character scalar naming the differential abundance method.
#' @param groups.sub Character vector of simulated group names to evaluate.
#' @param sim.implant Character scalar giving the HDF5 file path.
#' @param subset_size Integer scalar passed to \code{apply.test()}.
#' @param ab.scale Numeric vector of abundance-scale values used in the benchmark.
#' @param alpha Numeric scalar giving the significance threshold used in
#'   \code{eval.test()}.
#' @param norm Character scalar naming the normalization option passed to
#'   \code{apply.test()}.
#' @param covariates Optional character vector of metadata columns passed through
#'   the patched SIMBA workflow.
#' @param seed Seed specification passed to \code{furrr::furrr_options()}.
#' @param show_progress Logical; if \code{TRUE}, show progress via \code{furrr}.
#' @param simba.loc Optional file path to the local SIMBA fork; if supplied,
#'   \code{devtools::load_all()} is called inside each worker.
#' @param return_calib Logical; if \code{TRUE}, also return the combined
#'   calibration table.
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{evals}: combined evaluation summaries across groups
#'     \item \code{calib_df}: combined calibration table, or \code{NULL}
#'   }
#'
#' @export
run_method_eval_parallel <- function(test_name, groups.sub, sim.implant, subset_size,
                                     ab.scale, alpha = 0.05, norm = "pass",
                                     covariates = NULL, seed = TRUE,
                                     show_progress = TRUE, simba.loc = NULL,
                                     return_calib = FALSE) {
  
  run_one_group_eval_local <- run_one_group_eval
  
  # Parallelize over simulated groups for the current method.
  res_list <- furrr::future_map(
    groups.sub,
    ~ run_one_group_eval_local(
      g = .x,
      test_name = test_name,
      sim.implant = sim.implant,
      subset_size = subset_size,
      ab.scale = ab.scale,
      alpha = alpha,
      covariates = covariates,
      norm = norm,
      simba.loc = simba.loc,
      return_calib = return_calib
    ),
    .options = furrr::furrr_options(seed = seed, globals = TRUE),
    .progress = show_progress
  )
  
  # Bind benchmark summaries across groups.
  evals <- dplyr::bind_rows(purrr::map(res_list, "eval"))
  
  # Bind calibration data across groups if requested.
  calib_df <- NULL
  if (isTRUE(return_calib)) {
    baseline_abund_sim <- compute_baseline_abundance_from_h5(sim.implant)
    calib_df <- dplyr::bind_rows(purrr::map(res_list, "calib"))
    calib_df <- add_common_abundance_strata(
      df = calib_df,
      abundance_df = baseline_abund_sim
    )
  }
  
  list(
    evals = evals,
    calib_df = calib_df
  )
}
