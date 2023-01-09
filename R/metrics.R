### Metrics
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
