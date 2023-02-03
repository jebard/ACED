### Metrics
calculate_linear_regression <- function(actual_prop,estimated_proportion){
  avp <- as.vector(actual_prop)
  print(avp)
  evp <- as.vector(as.matrix(estimated_proportion))
  print(evp)
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

calculate_absolute_cell_error <- function(ref_obj=ref_obj,
                                          actual_prop=actual_prop,
                                          estimated_proportion=estimated_proportion){
  cells_per_sample <- as.vector(table(ref_obj$orig.ident)) ## calculate the number of cells-per-patient
  a.mat <- unclass(actual_prop) ### convert actual porpotion out of table object type
  b.mat <- as.matrix(estimated_proportion) ## gather estimated cells per cluster up
  a <- as.vector(a.mat * cells_per_sample) # multiply the actual proportion table against the total cells to get cells-per-cluster
  b <- as.vector(b.mat * cells_per_sample) # multiply the estimated proportion table against the total cells to get cells-per-cluster
  return(sum(abs(a-b))) ### take the absolute error rate of actual cells per cluster - estimated cells per cluster, and sum it up.
}

calculate_proportinal_correlation <- function(ref_obj=ref_obj,
                                          actual_prop=actual_prop,
                                          estimated_proportion=estimated_proportion){
  cells_per_sample <- as.vector(table(ref_obj$orig.ident)) ## calculate the number of cells-per-patient
  a.mat <- unclass(actual_prop) ### convert actual porpotion out of table object type
  b.mat <- as.matrix(estimated_proportion) ## gather estimated cells per cluster up
  a <- as.vector(a.mat * cells_per_sample) # multiply the actual proportion table against the total cells to get cells-per-cluster
  b <- as.vector(b.mat * cells_per_sample) # multiply the estimated proportion table against the total cells to get cells-per-cluster
  message(paste0("AMAT: ",a))
  message(paste0("BMAT: ",b))
  return(cor(a,b)) ### take the absolute error rate of actual cells per cluster - estimated cells per cluster, and sum it up.
}


