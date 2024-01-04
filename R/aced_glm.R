
aced_glm <- function(ref_obj){

# Define the number of pseudobulk samples and reference cell types
num_pseudobulk_samples <- length(ref_obj$orig.ident)
num_reference_cell_types <- length(ref_obj$seurat_clusters)

# Simulate pseudobulk single-cell data
set.seed(123)  # For reproducibility
#pseudobulk_data <- matrix(rpois(100 * num_pseudobulk_samples, 5), ncol = num_pseudobulk_samples)
pseudobulk_data <- AverageExpression(ref_obj,assay="RNA",slot = "count",group.by="orig.ident")$RNA
# Simulate reference single-cell data
#reference_data <- matrix(rpois(100 * num_reference_cell_types, 10), ncol = num_reference_cell_types)
reference_data <- AverageExpression(ref_obj,assays = "RNA",slot = "count",group.by="seurat_clusters")$RNA

num_pseudobulk_samples <-ncol(pseudobulk_data)
num_reference_cell_types <-ncol(reference_data)

# Create a sample mixture matrix (you can provide your own)
mixture_matrix <- matrix(runif(num_pseudobulk_samples * num_reference_cell_types), ncol = num_reference_cell_types)

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
    coef_i[coef_i < 0] <- 0

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
