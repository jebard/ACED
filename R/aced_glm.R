# ACeD Internal Deconvolution Strategies
# Provides three lightweight deconvolution methods that operate entirely within R,
# requiring no external Python tools. These complement the primary GEDIT3 and MuSiC
# strategies and are integrated into the ACED() / ACED_GS() / ACED_BRENT() workflows
# by passing strategy="glm", strategy="lasso", or strategy="spillover".
#
# @author jbard

# ─────────────────────────────────────────────────────────────────────────────
# SHARED HELPERS
# ─────────────────────────────────────────────────────────────────────────────

#' Build pseudobulk matrices from a Seurat object
#'
#' Returns a list with two matrices:
#'   $query     : samples × genes  (average expression per sample via orig.ident)
#'   $reference : genes × clusters (average expression per cluster via seurat_clusters)
#'
#' Both use AverageExpression() on the counts slot, which is the correct input
#' for linear-model-based deconvolution.
#'
#' @param ref_obj A preprocessed Seurat object with seurat_clusters and orig.ident set.
#' @return Named list with $query and $reference matrices.
#' @keywords internal
build_pseudobulk_matrices <- function(ref_obj) {
  # Query: M samples × G genes
  query_mat <- t(Seurat::AverageExpression(
    ref_obj, assay = "RNA", slot = "counts", group.by = "orig.ident"
  )$RNA)

  # Reference: G genes × K clusters
  ref_mat <- Seurat::AverageExpression(
    ref_obj, assay = "RNA", slot = "counts", group.by = "seurat_clusters"
  )$RNA

  # Align genes — intersection guards against any edge-case mismatch
  shared_genes <- intersect(colnames(query_mat), rownames(ref_mat))
  if (length(shared_genes) == 0) {
    stop("No shared genes between query and reference matrices. Check Seurat object integrity.")
  }
  if (length(shared_genes) < nrow(ref_mat)) {
    message(paste0("build_pseudobulk_matrices: using ", length(shared_genes),
                   " of ", nrow(ref_mat), " genes shared between query and reference."))
  }

  list(
    query     = query_mat[, shared_genes, drop = FALSE],   # M × G
    reference = ref_mat[shared_genes, , drop = FALSE]       # G × K
  )
}


#' Normalize a proportion vector: zero negatives then sum-normalize.
#'
#' Returns a zero vector (rather than NaN) if all values are <= 0,
#' which prevents silent NaN propagation when LASSO zeros everything out.
#'
#' @param x Numeric vector of raw coefficients.
#' @return Numeric vector summing to 1 (or all-zeros if input is degenerate).
#' @keywords internal
safe_normalize <- function(x) {
  x[x < 0] <- 0
  total <- sum(x, na.rm = TRUE)
  if (total == 0 || is.na(total)) {
    warning("safe_normalize: all coefficients are zero or NA — returning zero vector. ",
            "This may indicate the reference signatures do not explain this sample.")
    return(rep(0, length(x)))
  }
  x / total
}


# ─────────────────────────────────────────────────────────────────────────────
# aced_glm — pseudoinverse (Moore-Penrose) deconvolution
# ─────────────────────────────────────────────────────────────────────────────

#' ACeD pseudoinverse (GLM) deconvolution
#'
#' Deconvolves each pseudobulk sample against the full cell-type reference matrix
#' simultaneously using the Moore-Penrose pseudoinverse. This solves the
#' least-squares problem:
#'
#'   beta = pinv(Reference) %*% bulk_sample
#'
#' where Reference is G×K (genes × clusters) and bulk_sample is G×1.
#' All K coefficients are estimated jointly, so inter-cluster correlations
#' are accounted for — unlike fitting each cluster separately, which
#' produces biased estimates due to shared signal.
#'
#' @param ref_obj A preprocessed Seurat object with seurat_clusters and orig.ident.
#' @return A data frame (samples × clusters) of estimated cell-type proportions,
#'         with rownames = sample IDs and colnames = cluster IDs.
#' @export
aced_glm <- function(ref_obj) {
  message("ACeD: running pseudoinverse (GLM) deconvolution")

  mats <- build_pseudobulk_matrices(ref_obj)
  query_mat <- mats$query      # M × G
  ref_mat   <- mats$reference  # G × K

  num_samples  <- nrow(query_mat)  # M
  num_clusters <- ncol(ref_mat)    # K

  proportions <- matrix(
    0,
    nrow = num_samples,
    ncol = num_clusters,
    dimnames = list(rownames(query_mat), colnames(ref_mat))
  )

  # Compute the pseudoinverse of the reference matrix once (G×K → K×G).
  # For each sample, multiply pinv(Reference) [K×G] by bulk [G×1] → K×1.
  ref_pinv <- MASS::ginv(ref_mat)  # K × G

  for (i in seq_len(num_samples)) {
    bulk_sample <- as.vector(query_mat[i, ])  # G×1

    # Joint least-squares solution across all K clusters simultaneously
    beta <- ref_pinv %*% bulk_sample  # K×1

    proportions[i, ] <- safe_normalize(as.vector(beta))
  }

  return(as.data.frame(proportions))
}


# ─────────────────────────────────────────────────────────────────────────────
# aced_lasso — LASSO-regularized deconvolution
# ─────────────────────────────────────────────────────────────────────────────

#' ACeD LASSO deconvolution
#'
#' Deconvolves each pseudobulk sample against the reference cell-type matrix
#' using LASSO-regularized regression (cv.glmnet, alpha=1). The regularization
#' penalizes small/noisy coefficients toward zero, making this approach more
#' robust than the pseudoinverse when cell-type signatures are highly correlated.
#'
#' The regression is structured as:
#'   bulk_sample (G×1) ~ Reference (G×K) %*% beta (K×1)
#' so all K clusters compete simultaneously — this is the correct joint
#' formulation for deconvolution.
#'
#' Note: set.seed() should be called once in the calling loop (e.g. before
#' the resolution sweep in ACED()) rather than inside this function, so that
#' cross-validation fold assignments vary naturally across resolutions.
#'
#' @param ref_obj A preprocessed Seurat object with seurat_clusters and orig.ident.
#' @return A data frame (samples × clusters) of estimated cell-type proportions,
#'         with rownames = sample IDs and colnames = cluster IDs.
#' @export
aced_lasso <- function(ref_obj) {
  message("ACeD: running LASSO deconvolution")

  mats <- build_pseudobulk_matrices(ref_obj)
  query_mat <- mats$query      # M × G
  ref_mat   <- mats$reference  # G × K

  num_samples  <- nrow(query_mat)
  num_clusters <- ncol(ref_mat)

  proportions <- matrix(
    0,
    nrow = num_samples,
    ncol = num_clusters,
    dimnames = list(rownames(query_mat), colnames(ref_mat))
  )

  for (i in seq_len(num_samples)) {
    bulk_sample <- as.vector(query_mat[i, ])  # G×1

    # LASSO: regress bulk_sample ~ Reference
    # x = G×K reference, y = G×1 bulk — all K clusters compete simultaneously
    lasso_model <- glmnet::cv.glmnet(
      x     = ref_mat,         # G × K predictor matrix (genes as observations)
      y     = bulk_sample,     # G × 1 response
      alpha = 1,               # LASSO penalty
      lower.limits = 0         # enforce non-negativity directly — cleaner than
      # post-hoc zeroing and gives LASSO proper constraints
    )

    # Extract coefficients at lambda that minimises cross-validated error.
    # [-1] drops the intercept; we want only the K cluster coefficients.
    coefs <- as.vector(coef(lasso_model, s = "lambda.min"))[-1]  # K×1

    proportions[i, ] <- safe_normalize(coefs)
  }

  return(as.data.frame(proportions))
}


# ─────────────────────────────────────────────────────────────────────────────
# aced_lasso_spillover — cluster-level spillover analysis
# ─────────────────────────────────────────────────────────────────────────────

#' ACeD LASSO spillover analysis
#'
#' Constructs an artificial "bulk sample" by aggregating the raw counts of one
#' or more specified clusters, then deconvolves it against the full reference
#' using LASSO. A perfectly resolved cluster should return ~100% for itself;
#' a cluster with correlated neighbours will show proportional spillover.
#'
#' This is a diagnostic tool: run it on each cluster at a given resolution to
#' identify which clusters bleed into which others, guiding resolution selection.
#'
#' Design rationale:
#'   - AggregateExpression (sum of counts) is used for the pseudo-bulk "query",
#'     because a real bulk RNA-seq sample represents summed cell expression.
#'   - AverageExpression (mean counts) is used for the reference signatures,
#'     because deconvolution reference profiles should be per-cell averages
#'     for comparability across clusters of different sizes.
#'   The intentional mismatch in scale is addressed by LASSO's internal
#'   normalisation during coefficient estimation.
#'
#' @param ref_obj  A preprocessed Seurat object with seurat_clusters and orig.ident.
#' @param cluster  Character vector of one or more cluster name(s) to use as the
#'                 pseudo-bulk query. Must match names in levels(ref_obj$seurat_clusters).
#' @return Named numeric vector of length K (one entry per cluster), summing to 1.
#'         Names correspond to cluster IDs. Clusters with zero spillover have value 0.
#' @export
aced_lasso_spillover <- function(ref_obj, cluster = NULL) {
  if (is.null(cluster)) {
    stop("'cluster' argument is required: provide a cluster name or character vector of names. ",
         "Valid values: ", paste(levels(ref_obj$seurat_clusters), collapse = ", "))
  }

  valid_clusters <- levels(ref_obj$seurat_clusters)
  bad <- setdiff(cluster, valid_clusters)
  if (length(bad) > 0) {
    stop("Unknown cluster(s): ", paste(bad, collapse = ", "),
         ". Valid clusters: ", paste(valid_clusters, collapse = ", "))
  }

  message("ACeD spillover: building pseudo-bulk from cluster(s): ", paste(cluster, collapse = ", "))

  # --- Query: aggregate counts for the target cluster(s) ---
  # AggregateExpression gives the sum of raw counts per cluster (G×K matrix).
  # We take only the target cluster column(s), then collapse to a single G×1 vector
  # by summing across selected clusters (preserving total count depth).
  agg <- Seurat::AggregateExpression(
    ref_obj, assay = "RNA", slot = "counts", group.by = "seurat_clusters"
  )$RNA                                   # G × K

  if (length(cluster) == 1) {
    pseudo_bulk <- as.vector(agg[, cluster])           # G × 1
  } else {
    pseudo_bulk <- as.vector(rowSums(agg[, cluster]))  # sum across target clusters → G × 1
  }

  # --- Reference: average expression per cluster ---
  ref_mat <- Seurat::AverageExpression(
    ref_obj, assay = "RNA", slot = "counts", group.by = "seurat_clusters"
  )$RNA  # G × K

  # Align genes
  shared_genes <- intersect(names(pseudo_bulk), rownames(ref_mat))
  if (length(shared_genes) == 0) {
    stop("No shared genes between pseudo-bulk query and reference. Check Seurat object.")
  }
  pseudo_bulk <- pseudo_bulk[shared_genes]
  ref_mat     <- ref_mat[shared_genes, , drop = FALSE]

  message("ACeD spillover: fitting LASSO across ", ncol(ref_mat), " clusters using ",
          length(shared_genes), " genes")

  # --- LASSO: decompose pseudo_bulk against reference ---
  lasso_model <- glmnet::cv.glmnet(
    x            = ref_mat,     # G × K
    y            = pseudo_bulk, # G × 1
    alpha        = 1,
    lower.limits = 0            # non-negativity constraint
  )

  coefs <- as.vector(coef(lasso_model, s = "lambda.min"))[-1]  # K×1, drop intercept
  names(coefs) <- colnames(ref_mat)

  result <- safe_normalize(coefs)
  names(result) <- colnames(ref_mat)

  message("ACeD spillover: top contributors:")
  top <- sort(result, decreasing = TRUE)
  for (nm in names(top)[top > 0.01]) {
    message("  ", nm, ": ", round(top[nm] * 100, 1), "%")
  }

  return(result)
}
