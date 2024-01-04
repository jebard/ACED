
aced_glm <- function(ref_obj){
print("Defaulting to ACeD pseudoinverse strategy to perform linear regression-based deconvolution")
# Define the number of pseudobulk samples and reference cell types
num_pseudobulk_samples <- length(ref_obj$orig.ident)
num_reference_cell_types <- length(ref_obj$seurat_clusters)

# Simulate pseudobulk single-cell data
set.seed(123)  # For reproducibility

## Reference Pseudobulk
pseudobulk_data <- AverageExpression(ref_obj,assay="RNA",slot = "count",group.by="orig.ident")$RNA
## Reference Query by cluster
reference_data <- AverageExpression(ref_obj,assays = "RNA",slot = "count",group.by="seurat_clusters")$RNA

num_pseudobulk_samples <-ncol(pseudobulk_data)
num_reference_cell_types <-ncol(reference_data)

# Initialize a matrix to store proportions
proportions <- matrix(0, ncol = num_reference_cell_types, nrow = num_pseudobulk_samples)
rownames(proportions) <- colnames(pseudobulk_data)
colnames(proportions) <- colnames(reference_data)
# Loop through each pseudobulk sample
for (i in 1:num_pseudobulk_samples) {
  # Extract the expression data for the current pseudobulk sample
  pseudobulk_sample <- pseudobulk_data[, i]

  # Initialize a vector to store proportions for the current sample
  sample_proportions <- numeric(num_reference_cell_types)

  # Loop through each reference cell type
  for (j in 1:num_reference_cell_types) {
    # Extract the expression data for the current reference cell type
    reference_cell_type <- reference_data[, j]

    # Create a design matrix
    design_matrix <- cbind(Intercept = 1, reference_cell_type)


    # Use pseudoinverse to calculate coefficients
    coefficients <- ginv(design_matrix) %*% pseudobulk_sample
    # Extract and threshold the coefficients
    coef_i <- coefficients[2]

    coef_i[coef_i < 0] <- 0 ## if the estimated proportion is less than zero, set to zero.

    # Store the coefficients in the vector
    sample_proportions[j] <- coef_i
  }

  # Normalize proportions so they sum to 1 for the current sample
  sample_proportions <- sample_proportions / sum(sample_proportions, na.rm = TRUE)

  # Store the proportions in the matrix (column-wise)
  proportions[i, ] <- sample_proportions
}

return(proportions)
}


aced_lasso <- function(ref_obj){
  print("Defaulting to ACeD Lasso Deconvolution")
  # Define the number of pseudobulk samples and reference cell types
  num_pseudobulk_samples <- length(ref_obj$orig.ident)
  num_reference_cell_types <- length(ref_obj$seurat_clusters)

  # Simulate pseudobulk single-cell data
  set.seed(123)  # For reproducibility

  ## Reference Pseudobulk
  pseudobulk_data <- AverageExpression(ref_obj,assay="RNA",slot = "count",group.by="orig.ident")$RNA
  ## Reference Query by cluster
  reference_data <- AverageExpression(ref_obj,assays = "RNA",slot = "count",group.by="seurat_clusters")$RNA

  num_pseudobulk_samples <-ncol(pseudobulk_data)
  num_reference_cell_types <-ncol(reference_data)

  # Loop through each pseudobulk sample

  # Perform LASSO regression
  lasso_model <- cv.glmnet(x = reference_data, y = as.vector(pseudobulk_data[,1]), alpha = 1)  # alpha = 1 specifies LASSO regularization

  # Extract coefficients from the model for the lambda that minimizes cross-validated error
  coefficients <- coef(lasso_model, s = "lambda.min")

  # Extract non-zero coefficients (those selected by LASSO)
  nonzero_coefficients <- coefficients[-1, , drop = FALSE]

  # Sum of non-zero coefficients
  sum_nonzero <- sum(abs(nonzero_coefficients))

  # Calculate proportions by normalizing non-zero coefficients
  proportions <- nonzero_coefficients / sum_nonzero



  #coef_i[coef_i < 0] <- 0 ## if the estimated proportion is less than zero, set to zero.

  # Normalize proportions so they sum to 1 for the current sample
  #coef_i <- coef_i / sum(coef_i, na.rm = TRUE)


  return(proportions)

}

