# Generate Cluster Statistics
# @author jbard
#

#' calc_cluster_proportions
#'
#' @return
#' @export
#'
#' @examples
calc_cluster_proportions <- function() {

}

### https://rdrr.io/github/statlab/permuter/src/R/utils.R
permute_within_rows <- function(x) {
  for (row in seq_len(nrow(x))) {
    x[row, ] <- sample(x[row, ])
  }
  return(x)
}

get_optimal_resolution <- function(drrsd_object){
  return(drrsd_object[head(n=1,order((drrsd_object$ACE_Random-drrsd_object$ACE),decreasing = T)),]$Resolution)
}

get_least_optimal_resolution <- function(drrsd_object){
  return(drrsd_object[head(n=1,order((drrsd_object$ACE_Random-drrsd_object$ACE),decreasing = F)),]$Resolution)
}
