# Deconvolution Strategy Pattern
# Ingests a request to run deconvolution, returns a normalized result regardless of strategy
# @author jbard


#' run_deconvolution
#' This function takes in a reference single-cell object, bulk or pseudobulk query, and
#' requires a chosen deconvolution strategy. It seeks to call the external tool, and return
#' a standardized data table of sample-to-cell type predictions, regardless of implementation
#' @param strategy
#' @param ref_obj
#' @param query_obj
#' @return
#' @export
#' @examples
run_deconvolution <- function(ref_obj, query_obj, strategy){
  if(strategy="bayesprism"){
    run_bayesprism()
  } else if (strategy=="music") {
    run_music()
  } else if (strategy=="gedit"){
    run_gedit(ref_obj,query_obj)
  } else {
    warning("Invalid deconvolution strategy specified")
  }
}


