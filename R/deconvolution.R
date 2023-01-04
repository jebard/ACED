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

  ### next we gather the TRUE or Expected proportions of the reference object
  actual_proportion <- get_cluster_proportions(ref_obj)
  MAE <- calculate_mean_absolute_error(actual_prop = actual_proportion,estimated_proportion = estimated_proportions )
  message("Deconvolution results in MAE: ",MAE)
}

get_cluster_proportions <- function(ref_obj){
  proportions <- prop.table(table(ref_obj$orig.ident,ref_obj$seurat_clusters),margin=1)
  return(proportions)
}

calculate_mean_absolute_error <- function(actual_prop,estimated_proportion){
  avp <- as.vector(actual_prop)
  evp <- as.vector(as.matrix(estimated_proportion))
  return(mean(abs(avp - evp)))
}
