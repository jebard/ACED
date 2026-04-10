# instaPrism_strategy.R
# ACeD wrapper for InstaPrism Bayesian deconvolution
#
# InstaPrism is a fast approximation of BayesPrism that uses an EM algorithm
# to estimate cell-type proportions from bulk RNA-seq using a single-cell
# reference. It is the recommended replacement for BayesPrism in large
# cohort analyses due to substantially reduced compute time with equivalent
# accuracy.
#
# Citation (required):
#   Mengying Hu et al. (2023). InstaPrism: an R package for fast implementation
#   of BayesPrism. Bioinformatics Advances, vbad100.
#   https://doi.org/10.1093/bioadv/vbad100
#
# Installation:
#   devtools::install_github("humengying0907/InstaPrism")
#
# Key differences from MuSiC and GEDIT3:
#   - Uses per-cell NORMALIZED expression (log1p data slot), not raw counts
#   - Reference is built from individual cells, not cluster-averaged profiles
#   - Runs Bayesian inference independently per bulk sample — works correctly
#     for single-sample objects without any workaround
#   - Returns both cell-type and cell-state level proportions; ACED uses
#     cell-type level (seurat_clusters) to match all other strategies
#
# @author jbard


# ─────────────────────────────────────────────────────────────────────────────
# MAIN ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

#' Run InstaPrism deconvolution within the ACeD framework
#'
#' Prepares the single-cell reference and pseudobulk query from a Seurat object,
#' runs InstaPrism deconvolution, and returns estimated cell-type proportions
#' in the standard ACeD format (samples x clusters data frame).
#'
#' Both cell.type.labels and cell.state.labels are set to seurat_clusters
#' because ACED operates at a single level of clustering resolution. The
#' cell-type level result (@Post.ini.ct@theta) is returned, which in this
#' context corresponds directly to the Seurat cluster assignments.
#'
#' @param ref_obj A preprocessed Seurat object with seurat_clusters and
#'   orig.ident set. The RNA assay must contain a normalized data layer/slot.
#' @return A data frame (samples x clusters) of estimated cell-type proportions,
#'   with rownames matching orig.ident values and colnames matching cluster labels.
#' @export
run_instaprism <- function(ref_obj) {

  if (!requireNamespace("InstaPrism", quietly = TRUE)) {
    stop("Package 'InstaPrism' is required for strategy='instaprism'.\n",
         "Install with: devtools::install_github('humengying0907/InstaPrism')")
  }

  message("Please cite InstaPrism: Hu et al. (2023) Bioinformatics Advances vbad100 ",
          "https://doi.org/10.1093/bioadv/vbad100")
  message("Preparing InstaPrism reference from single-cell object")

  ref_phi <- instaprism_prep_reference(ref_obj)

  message("Preparing bulk query for InstaPrism")
  bulk_mat <- instaprism_prep_query(ref_obj)

  message("Running InstaPrism deconvolution")
  deconv_res <- InstaPrism::InstaPrism(
    bulk_Expr = bulk_mat,
    refPhi_cs = ref_phi
  )

  # Extract cell-type level proportions: theta is K x M (clusters x samples)
  # Transpose to M x K (samples x clusters) to match ACeD convention
  result <- as.data.frame(t(deconv_res@Post.ini.ct@theta))

  # Restore original sample names as rownames.
  # InstaPrism uses column names of bulk_Expr as sample identifiers.
  # These should already match orig.ident values, but we enforce it explicitly.
  expected_samples <- rownames(result)
  actual_samples   <- unique(ref_obj$orig.ident)
  if (length(expected_samples) == length(actual_samples)) {
    rownames(result) <- actual_samples
  }

  # Strip any 'g' prefix that Seurat v5 may have prepended to numeric cluster names
  if (all(grepl("^g[0-9]", colnames(result)))) {
    colnames(result) <- sub("^g", "", colnames(result))
  }

  return(result)
}


# ─────────────────────────────────────────────────────────────────────────────
# REFERENCE PREPARATION
# ─────────────────────────────────────────────────────────────────────────────

#' Prepare InstaPrism reference object from a Seurat object
#'
#' Extracts the normalized expression matrix and cluster labels from the Seurat
#' object and passes them to InstaPrism::refPrepare(). Both cell.type.labels
#' and cell.state.labels are set to seurat_clusters because ACeD works at a
#' single resolution level — the caller (ACED / ACED_GS / ACED_BRENT) updates
#' seurat_clusters before each call to evaluate_deconvolution().
#'
#' InstaPrism requires the NORMALIZED data slot (log1p-transformed), not raw
#' counts. Using raw counts produces inaccurate reference profiles because
#' BayesPrism/InstaPrism's EM algorithm assumes the input is on a log scale.
#'
#' @param ref_obj Seurat object.
#' @return An InstaPrism refPhi object ready for InstaPrism().
#' @keywords internal
instaprism_prep_reference <- function(ref_obj) {

  ver <- tryCatch(utils::packageVersion("Seurat"), error = function(e) "4.0.0")

  # Extract normalized expression: genes x cells
  # Use the 'data' layer/slot (log1p-normalized), not 'counts'
  if (ver >= "5.0.0") {
    sc_expr <- as.matrix(Seurat::GetAssayData(ref_obj, assay = "RNA",
                                              layer = "data"))
  } else {
    sc_expr <- as.matrix(Seurat::GetAssayData(ref_obj, assay = "RNA",
                                              slot  = "data"))
  }

  # Cell labels: both type and state use seurat_clusters
  # In a typical BayesPrism/InstaPrism workflow, cell.type.labels is coarse
  # (e.g. T cell, B cell) and cell.state.labels is fine (e.g. CD4+ effector).
  # Within ACED, seurat_clusters IS the resolution being optimised, so it
  # serves as both levels — the cell-type output is what ACED evaluates.
  cluster_labels <- as.character(ref_obj$seurat_clusters)

  refPhi_obj <- InstaPrism::refPrepare(
    sc_Expr          = sc_expr,
    cell.type.labels = cluster_labels,
    cell.state.labels = cluster_labels
  )

  return(refPhi_obj)
}


# ─────────────────────────────────────────────────────────────────────────────
# QUERY (BULK) PREPARATION
# ─────────────────────────────────────────────────────────────────────────────

#' Prepare bulk pseudobulk query matrix for InstaPrism
#'
#' InstaPrism expects bulk_Expr as a genes x samples matrix.
#' We pseudobulk by sample identity (orig.ident) using AverageExpression on
#' the normalized data slot to match the scale of the reference expression.
#'
#' Single-sample objects are handled by computing rowMeans directly via
#' GetAssayData, bypassing the Seurat v5 grouping bug (same pattern as other
#' ACeD strategies). InstaPrism handles single-sample input natively — no
#' duplication workaround is needed.
#'
#' @param ref_obj Seurat object.
#' @return A genes x samples plain numeric matrix.
#' @keywords internal
instaprism_prep_query <- function(ref_obj) {

  ver       <- tryCatch(utils::packageVersion("Seurat"), error = function(e) "4.0.0")
  n_samples <- length(unique(ref_obj$orig.ident))

  if (n_samples == 1) {
    message("instaprism_prep_query: single-sample object detected -- ",
            "computing average expression directly.")
    sample_name <- as.character(unique(ref_obj$orig.ident))
    if (ver >= "5.0.0") {
      expr_data <- as.matrix(Seurat::GetAssayData(ref_obj, assay = "RNA",
                                                  layer = "data"))
    } else {
      expr_data <- as.matrix(Seurat::GetAssayData(ref_obj, assay = "RNA",
                                                  slot  = "data"))
    }
    # Build G x 1 matrix with the correct sample name
    bulk_mat <- matrix(rowMeans(expr_data, na.rm = TRUE), ncol = 1,
                       dimnames = list(rownames(expr_data), sample_name))

  } else {
    # Multi-sample: AverageExpression by orig.ident returns G x M
    # Use the normalized data slot to match the reference expression scale
    if (ver >= "5.0.0") {
      avg <- Seurat::AverageExpression(ref_obj, assay = "RNA",
                                       layer    = "data",
                                       group.by = "orig.ident")
    } else {
      avg <- Seurat::AverageExpression(ref_obj, assay = "RNA",
                                       slot     = "data",
                                       group.by = "orig.ident")
    }
    bulk_mat <- as.matrix(avg$RNA)  # G x M, already correct orientation
  }

  # Guard: replace NAs
  if (any(is.na(bulk_mat))) {
    bulk_mat[is.na(bulk_mat)] <- 0
    message("instaprism_prep_query: NA values replaced with 0.")
  }

  return(bulk_mat)
}
