# Deconvolution Strategy Pattern
# Ingests a request to run cellular deconvolution, returns a normalized result regardless of strategy
# @author jbard

#' run_deconvolution
#'
#' @decription This function takes in a reference single-cell object, bulk or pseudobulk query, and
#' requires a chosen deconvolution strategy. It seeks to call the external tool, and return
#' a standardized data table of sample-to-cell type predictions, regardless of implementation
#' @param strategy
#' @param ref_obj
#' @param query_obj
#' @return matrix of sample-to-cell type predicted values.
#' @export
#' @examples
evaluate_deconvolution <- function(ref_obj, query_obj, strategy){
  ### first, launch the deconvolution analysis and get the results
  estimated_proportions <- NULL
  if(strategy=="bayesprism"){
    estimated_proportions <- run_bayesprism()
  } else if (strategy=="music") {
    estimated_proportions <- run_music()
  } else if (strategy=="gedit"){
    estimated_proportions <- run_gedit(ref_obj,query_obj)
  } else {
    warning("Invalid deconvolution strategy specified")
  }
  ### next we gather the various metrics and report back
  actual_proportion <- get_cluster_proportions(ref_obj)
  cluster_cell_counts <- get_cluster_cell_counts(ref_obj)
  random_proportions <- get_random_proportions(ref_obj)

  ## calculate out the important variables to return
  MAE <- calculate_mean_absolute_error(actual_prop = actual_proportion,estimated_proportion = estimated_proportions)
  MAE_RANDOM <- calculate_mean_absolute_error(actual_prop = actual_proportion,estimated_proportion = random_proportions)
  RSE <- calculate_relative_squared_error(actual_prop = actual_proportion,estimated_proportion = estimated_proportions)
  SMAPE <- calculate_smape(actual_prop = actual_proportion,estimated_proportion = estimated_proportions)
  RMSE <- calculate_rmse(actual_prop = actual_proportion,estimated_proportion = estimated_proportions)
  AVP_Z <- count_actual_zero(actual_proportion)
  EVP_Z <- count_predicted_zero(estimated_proportions)
  AE <- calculate_absolute_error(actual_proportion,estimated_proportions)
  AE_CC <- calculate_cell_absolute_error(actual_proportion,estimated_proportions,cluster_cell_counts)
  LM <- calculate_linear_regression(actual_proportion,estimated_proportions)
  ACE <- calculate_absolute_cell_error(ref_obj,actual_proportion,estimated_proportions)
  ACE_Random <- calculate_absolute_cell_error(ref_obj,actual_proportion,random_proportions)
  message("Deconvolution results in: ",MAE,",",RSE,",",SMAPE,",",RMSE,",",AVP_Z,",",EVP_Z,",",AE,",",AE_CC,",",LM,",",ACE,ACE_Random,MAE_RANDOM)
  return(c(MAE,RSE,SMAPE,RMSE,AVP_Z,EVP_Z,AE,AE_CC,LM,ACE,ACE_Random,MAE_RANDOM))
}

get_random_proportions <- function(ref_obj){
  actual_proportion <- get_cluster_proportions(ref_obj)
  # first get the number of samples in the object
  samples <- length(unique(ref_obj$orig.ident))
  # next get the number of clusters in the object
  clusts <- length(levels(ref_obj$seurat_clusters))
  m <- matrix(rnorm(samples * clusts,mean(actual_proportion),sd = sd(actual_proportion)), nrow=samples)
  prob <- exp(m)/rowSums(exp(m))
  return(prob)
}

calculate_linear_regression <- function(actual_prop,estimated_proportion){
  avp <- as.vector(actual_prop)
  evp <- as.vector(as.matrix(estimated_proportion))
  lin_reg <- lm(avp~evp)
  return(sum(abs(lin_reg$residuals)))
}

calculate_absolute_error <- function(actual_prop,estimated_proportion){
  avp <- as.vector(actual_prop)
  evp <- as.vector(as.matrix(estimated_proportion))
  return(sum(abs(avp - evp)))
}

calculate_cell_absolute_error <- function(actual_prop,estimated_proportion,cluster_cell_counts){
  avp <- as.vector(actual_prop)
  evp <- as.vector(as.matrix(estimated_proportion))
  return(sum(abs(avp - evp)*cluster_cell_counts))
}

get_cluster_proportions <- function(ref_obj){
  proportions <- prop.table(table(ref_obj$orig.ident,ref_obj$seurat_clusters),margin=1)
  return(proportions)
}

get_cluster_cell_counts <- function(ref_obj){
  return(as.vector(table(ref_obj$orig.ident,ref_obj$seurat_clusters)))
}

calculate_mean_absolute_error <- function(actual_prop,estimated_proportion){
  avp <- as.vector(actual_prop)
  evp <- as.vector(as.matrix(estimated_proportion))
  return(mean(abs(avp - evp)))
}

calculate_relative_squared_error <- function(actual_prop,estimated_proportion){
  avp <- as.vector(actual_prop)
  evp <- as.vector(as.matrix(estimated_proportion))
  return(Metrics::rse(avp,evp))
}

calculate_smape <- function(actual_prop,estimated_proportion){
  avp <- as.vector(actual_prop)
  evp <- as.vector(as.matrix(estimated_proportion))
  return(Metrics::smape(avp + .0000000001,evp))
}

calculate_rmse <- function(actual_prop,estimated_proportion){
  avp <- as.vector(actual_prop)
  evp <- as.vector(as.matrix(estimated_proportion))
  return(Metrics::rmse(avp,evp))
}

count_actual_zero <- function(actual_prop){
  return(sum(actual_prop==0))
}

count_predicted_zero <- function(estimated_proportion){
  return(sum(estimated_proportion==0))
}

calculate_absolute_cell_error <- function(ref_obj=ref_obj,actual_prop=actual_prop,
                                          estimated_proportion=estimated_proportion){
  cells_per_sample <- as.vector(table(ref_obj$orig.ident))
  a.mat <- unclass(actual_prop)
  b.mat <- as.matrix(estimated_proportion)
  a <- as.vector(a.mat * cells_per_sample)
  b <- as.vector(b.mat * cells_per_sample)
  message("ACE:", length(a))
  message("ACE:", length(b))
  return(sum(abs(a-b)))
}

clear_resolutions <- function(){
  for (res in seq(from=0.01, to=3.00,by=.01)){
    resolution_string <- paste0("integrated_snn_res.",res)
    tryCatch(
      expr ={
      combined.seurat.sct[[resolution_string]] <- NULL
      }, error = function(e){
      })
  }
}

calculate_cluster_tree <- function(){
  combined.seurat.sct <- BuildClusterTree(combined.seurat.sct)
  data.tree <- Tool(object = combined.seurat.sct, slot = "BuildClusterTree")
  return(ape::cophenetic.phylo(data.tree))
}
