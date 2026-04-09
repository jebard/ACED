# ACeD Internal Deconvolution Strategies
# Provides two lightweight deconvolution methods that operate entirely within R,
# requiring no external Python tools. These complement the primary GEDIT3 and MuSiC
# strategies and are integrated into the ACED() / ACED_GS() / ACED_BRENT() workflows
# by passing strategy="glm" (NNLS) or strategy="lasso" (LASSO).
#
# aced_lasso_spillover() is a standalone post-hoc diagnostic — call it directly
# after an ACED run; it cannot be used as a strategy= argument.
#
# Dependencies:
#   strategy="glm"   -> install.packages("nnls")
#   strategy="lasso" -> install.packages("glmnet")  [usually already installed]
#
# @author jbard

# ─────────────────────────────────────────────────────────────────────────────
# SEURAT v4 / v5 COMPATIBILITY WRAPPERS
# ─────────────────────────────────────────────────────────────────────────────

# In Seurat v5, the 'slot' parameter was deprecated in favour of 'layer', and
# the return value of AverageExpression()$RNA may be a sparse dgCMatrix rather
# than a dense base matrix.  t.default() cannot handle sparse matrices, so we
# coerce explicitly with as.matrix() and route the slot/layer argument based on
# the installed Seurat version.

.seurat_avg_expr <- function(obj, assay = "RNA", group.by) {
  ver <- tryCatch(utils::packageVersion("Seurat"), error = function(e) "4.0.0")
  if (ver >= "5.0.0") {
    result <- Seurat::AverageExpression(obj, assay = assay, layer = "counts",
                                        group.by = group.by)
  } else {
    result <- Seurat::AverageExpression(obj, assay = assay, slot = "counts",
                                        group.by = group.by)
  }
  as.matrix(result[[assay]])   # always return a plain dense matrix (genes x groups)
}

.seurat_agg_expr <- function(obj, assay = "RNA", group.by) {
  ver <- tryCatch(utils::packageVersion("Seurat"), error = function(e) "4.0.0")
  if (ver >= "5.0.0") {
    result <- Seurat::AggregateExpression(obj, assay = assay, layer = "counts",
                                          group.by = group.by)
  } else {
    result <- Seurat::AggregateExpression(obj, assay = assay, slot = "counts",
                                          group.by = group.by)
  }
  as.matrix(result[[assay]])   # genes x groups, dense
}


# ─────────────────────────────────────────────────────────────────────────────
# SHARED HELPERS
# ─────────────────────────────────────────────────────────────────────────────

#' Build pseudobulk matrices from a Seurat object (Seurat v4/v5 compatible)
#'
#' Returns a list:
#'   $query     : M samples x G genes  (AverageExpression by orig.ident, transposed)
#'   $reference : G genes x K clusters (AverageExpression by seurat_clusters)
#'
#' @param ref_obj Preprocessed Seurat object with seurat_clusters and orig.ident.
#' @return Named list with $query (M x G) and $reference (G x K).
#' @keywords internal
build_pseudobulk_matrices <- function(ref_obj) {

  # Reference: G x K — average expression per cluster
  ref_mat <- .seurat_avg_expr(ref_obj, group.by = "seurat_clusters")  # G x K

  # Query: G x M from AverageExpression, transposed to M x G
  query_mat <- t(.seurat_avg_expr(ref_obj, group.by = "orig.ident"))   # M x G

  # Align genes (columns of query, rows of reference)
  shared_genes <- intersect(colnames(query_mat), rownames(ref_mat))
  if (length(shared_genes) == 0) {
    stop("build_pseudobulk_matrices: no shared genes between query and reference. ",
         "Check that the Seurat object has been processed correctly.")
  }
  if (length(shared_genes) < nrow(ref_mat)) {
    message("build_pseudobulk_matrices: using ", length(shared_genes),
            " of ", nrow(ref_mat), " genes present in both matrices.")
  }

  list(
    query     = query_mat[, shared_genes, drop = FALSE],   # M x G
    reference = ref_mat[shared_genes,   , drop = FALSE]    # G x K
  )
}


#' Normalize a coefficient vector to cell-type proportions
#'
#' Zeros any negative values, then divides by the total. Returns a zero vector
#' (not NaN) if every coefficient is <= 0, so downstream metrics receive valid
#' numeric input even when a deconvolution solution is degenerate.
#'
#' @param x Numeric vector of raw coefficients, length K.
#' @return Numeric vector of proportions summing to 1 (or all-zeros if degenerate).
#' @keywords internal
safe_normalize <- function(x) {
  x[is.na(x) | x < 0] <- 0
  total <- sum(x)
  if (total == 0) {
    warning("safe_normalize: all coefficients are zero — returning zero proportion vector. ",
            "This resolution may produce degenerate deconvolution results.")
    return(rep(0, length(x)))
  }
  x / total
}


# ─────────────────────────────────────────────────────────────────────────────
# aced_glm — Non-negative least squares (NNLS) deconvolution
# ─────────────────────────────────────────────────────────────────────────────
#
# WHY NOT PSEUDOINVERSE?
# The Moore-Penrose pseudoinverse (ginv) fails on real scRNA-seq data because:
#
#   1. Raw count data is unscaled: AverageExpression returns mean counts per
#      gene per cluster. Clusters with more total RNA (higher total counts)
#      dominate the pseudoinverse solution regardless of composition, because
#      the minimum-norm solution assigns all weight to the direction of maximum
#      variance — which is total expression level, not cell-type identity.
#
#   2. High collinearity: with ~20,000 genes and only K=3-20 clusters, the
#      reference matrix is G >> K overdetermined. Thousands of housekeeping
#      genes are shared across all clusters and produce near-identical columns.
#      The resulting large condition number means tiny numerical differences
#      between cluster signatures produce huge swings in the unconstrained
#      coefficient vector, collapsing to a single dominant cluster.
#
# THE FIX — Non-negative least squares (NNLS):
#   NNLS solves:  min ||Reference %*% beta - bulk||^2  subject to beta >= 0
#
#   The non-negativity constraint regularises the solution: it cannot use
#   large negative and positive cancelling terms to achieve minimum norm,
#   so collinear columns cannot destabilise the result. NNLS is the algorithm
#   underlying CIBERSORT, TIMER, and most standard deconvolution tools.
#
#   We additionally L2-normalise each reference column before fitting (divide
#   by its Euclidean norm) so that no cluster dominates simply by having higher
#   total counts. Coefficients are rescaled back after fitting to be
#   interpretable on the original proportion scale, then sum-normalised to 1.

#' ACeD NNLS deconvolution
#'
#' Deconvolves each pseudobulk sample against the reference using non-negative
#' least squares (NNLS). Reference columns are L2-normalised before fitting to
#' remove the effect of differences in total counts between clusters.
#' This is the same statistical framework used by CIBERSORT and TIMER.
#'
#' Requires the 'nnls' package (install.packages("nnls")).
#'
#' @param ref_obj Preprocessed Seurat object with seurat_clusters and orig.ident.
#' @return Data frame (samples x clusters) of estimated cell-type proportions.
#' @export
aced_glm <- function(ref_obj) {
  message("ACeD: running NNLS deconvolution")

  if (!requireNamespace("nnls", quietly = TRUE)) {
    stop("Package 'nnls' is required for strategy='glm'. ",
         "Install it with: install.packages('nnls')")
  }

  mats     <- build_pseudobulk_matrices(ref_obj)
  query    <- mats$query      # M x G
  ref_mat  <- mats$reference  # G x K

  num_samples  <- nrow(query)
  num_clusters <- ncol(ref_mat)

  # L2-normalise each reference column (cluster signature) so that
  # no cluster dominates the NNLS fit due to higher total count magnitude.
  # We store the norms and divide back out after fitting so that coefficients
  # retain their proportion interpretation before sum-normalisation.
  col_norms <- sqrt(colSums(ref_mat^2))
  col_norms[col_norms == 0] <- 1  # guard against all-zero columns
  ref_normalised <- sweep(ref_mat, 2, col_norms, FUN = "/")  # G x K, unit columns

  proportions <- matrix(
    0, nrow = num_samples, ncol = num_clusters,
    dimnames = list(rownames(query), colnames(ref_mat))
  )

  for (i in seq_len(num_samples)) {
    bulk <- query[i, ]  # G x 1

    # Also L2-normalise the bulk sample for numerical comparability
    bulk_norm <- bulk / max(sqrt(sum(bulk^2)), 1e-10)

    # NNLS: min ||ref_normalised %*% beta - bulk_norm||^2, beta >= 0
    result <- nnls::nnls(A = ref_normalised, b = bulk_norm)
    beta   <- result$x  # K x 1, non-negative

    proportions[i, ] <- safe_normalize(beta)
  }

  return(as.data.frame(proportions))
}


# ─────────────────────────────────────────────────────────────────────────────
# aced_lasso — LASSO-regularized deconvolution
# ─────────────────────────────────────────────────────────────────────────────

#' ACeD LASSO deconvolution
#'
#' Deconvolves each pseudobulk sample using LASSO regression:
#'   bulk_sample (G x 1) ~ Reference (G x K) %*% beta (K x 1)
#' All K clusters compete simultaneously. LASSO regularisation penalises
#' noisy/absent cell types toward zero; lower.limits=0 enforces non-negativity
#' as a hard constraint inside the solver (more principled than post-hoc zeroing).
#'
#' Note: do not call set.seed() inside this function. The caller (evaluate_deconvolution
#' or the resolution loop) should set the seed once before iterating.
#'
#' @param ref_obj Preprocessed Seurat object with seurat_clusters and orig.ident.
#' @return Data frame (samples x clusters) of estimated cell-type proportions.
#' @export
aced_lasso <- function(ref_obj) {
  message("ACeD: running LASSO deconvolution")

  mats     <- build_pseudobulk_matrices(ref_obj)
  query    <- mats$query      # M x G
  ref_mat  <- mats$reference  # G x K

  num_samples  <- nrow(query)
  num_clusters <- ncol(ref_mat)

  proportions <- matrix(
    0, nrow = num_samples, ncol = num_clusters,
    dimnames = list(rownames(query), colnames(ref_mat))
  )

  for (i in seq_len(num_samples)) {
    bulk <- query[i, ]  # G x 1 — current sample's average expression

    # LASSO: predict bulk from reference (G x K predictor, G x 1 response)
    # lower.limits=0 enforces non-negativity — no negative proportions
    lasso_model <- glmnet::cv.glmnet(
      x             = ref_mat,
      y             = bulk,
      alpha         = 1,
      lower.limits  = 0
    )

    # Drop intercept (index 1); keep the K cluster coefficients
    coefs <- as.vector(coef(lasso_model, s = "lambda.min"))[-1]
    names(coefs) <- colnames(ref_mat)

    proportions[i, ] <- safe_normalize(coefs)
  }

  return(as.data.frame(proportions))
}


# ─────────────────────────────────────────────────────────────────────────────
# aced_lasso_spillover — cluster-level spillover diagnostic (standalone only)
# ─────────────────────────────────────────────────────────────────────────────

#' ACeD LASSO spillover diagnostic
#'
#' A DIAGNOSTIC TOOL — not a deconvolution strategy. Use this after running ACED()
#' to inspect how well a specific cluster is resolved at the optimal resolution.
#'
#' Aggregates raw counts for the specified cluster(s) into a pseudo-bulk sample
#' (simulating what that cluster looks like as a bulk RNA-seq library), then
#' deconvolves it against the full reference. A well-resolved cluster should
#' return ~100% for itself. Spillover to neighbours indicates the cluster's
#' signature overlaps with related types at this resolution.
#'
#' Design note: AggregateExpression (sum of counts) is used for the pseudo-bulk
#' query because a real bulk sample represents summed cell signal.
#' AverageExpression (mean counts) is used for the reference so signatures are
#' comparable per-cell across clusters of different sizes. The scale difference
#' is absorbed by LASSO's internal coefficient estimation.
#'
#' @param ref_obj  Preprocessed Seurat object with seurat_clusters and orig.ident.
#' @param cluster  Character vector of one or more cluster names to use as the
#'                 pseudo-bulk query. Must match levels(ref_obj$seurat_clusters).
#' @return Named numeric vector of proportions (length K, sums to 1).
#'         Each element is the estimated contribution of that cluster.
#' @export
aced_lasso_spillover <- function(ref_obj, cluster = NULL) {

  # —— Input validation ——
  if (is.null(cluster)) {
    stop("'cluster' is required. Provide one or more cluster names.\n",
         "Valid clusters: ", paste(levels(ref_obj$seurat_clusters), collapse = ", "))
  }
  invalid <- setdiff(as.character(cluster), levels(ref_obj$seurat_clusters))
  if (length(invalid) > 0) {
    stop("Unknown cluster(s): ", paste(invalid, collapse = ", "), "\n",
         "Valid clusters: ", paste(levels(ref_obj$seurat_clusters), collapse = ", "))
  }

  message("ACeD spillover: building pseudo-bulk from cluster(s): ",
          paste(cluster, collapse = ", "))

  # —— Query: aggregate raw counts for target cluster(s) ——
  # AggregateExpression returns G x K matrix of summed counts
  agg_mat <- .seurat_agg_expr(ref_obj, group.by = "seurat_clusters")  # G x K

  cluster <- as.character(cluster)
  if (length(cluster) == 1) {
    pseudo_bulk <- agg_mat[, cluster]            # G x 1 vector
  } else {
    pseudo_bulk <- rowSums(agg_mat[, cluster])   # sum across target clusters → G
  }

  # —— Reference: average expression per cluster ——
  ref_mat <- .seurat_avg_expr(ref_obj, group.by = "seurat_clusters")  # G x K

  # —— Align genes ——
  shared_genes <- intersect(names(pseudo_bulk), rownames(ref_mat))
  if (length(shared_genes) == 0) {
    stop("No shared genes between pseudo-bulk query and reference matrix.")
  }
  pseudo_bulk <- pseudo_bulk[shared_genes]
  ref_mat     <- ref_mat[shared_genes, , drop = FALSE]

  message("ACeD spillover: fitting LASSO — ", ncol(ref_mat), " clusters, ",
          length(shared_genes), " genes")

  # —— LASSO ——
  lasso_model <- glmnet::cv.glmnet(
    x            = ref_mat,
    y            = pseudo_bulk,
    alpha        = 1,
    lower.limits = 0
  )

  coefs <- as.vector(coef(lasso_model, s = "lambda.min"))[-1]
  names(coefs) <- colnames(ref_mat)

  result <- safe_normalize(coefs)
  names(result) <- colnames(ref_mat)

  # —— Report ——
  message("ACeD spillover: results (>1%):")
  for (nm in names(sort(result, decreasing = TRUE))) {
    if (result[nm] > 0.01) {
      message("  Cluster ", nm, ": ", round(result[nm] * 100, 1), "%")
    }
  }

  return(result)
}
