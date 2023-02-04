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

