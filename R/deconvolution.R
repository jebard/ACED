# Deconvolution Strategy Pattern
# Ingests a request to run cellular deconvolution, returns a normalized result regardless of strategy
# @author jbard

DRRSD <- function(ref_obj=ref_obj,query_obj=query_obj,start=0.01,stop=1,step=.05,algorithm="louvain",method="matrix"){
  values_mae = c();values_rse = c();values_smape = c();values_rmse = c()
  values_actual_zero = c();values_predicted_zero = c();clusters = c();resolution = c();
  values_ae = c();values_ae_cc = c();values_lm_res = c()
  values_ACE = c();values_ACE_random = c();values_MAE_random = c();values_PC = c()

  for (res in c(seq(from=start,to=stop,by=step))){
  #for (res in c(0.008,0.01,0.025,0.03,0.036,0.04,0.077,0.08,0.1,0.15,0.2,0.25,0.3,0.4,0.5)){#},0.6,0.8,1,1.2,1.5,2,2.5,3)){
    resolution_string <- paste0("integrated_snn_res.",res)
    #if (resolution_string %in% colnames(ref_obj@meta.data)) {
    #  message("Resolution already calculated for reference, skipping")
    #  ref_obj$seurat_clusters <- ref_obj[[resolution_string]]
    #} else {
      message(paste0("Calulating for reference, running ",algorithm," clustering"))
      if(algorithm=="leiden" || algorithm == 4){
        ref_obj <- FindClusters(ref_obj,resolution = res,algorithm=algorithm,verbose=T,method=method)
      } else {
      ref_obj <- FindClusters(ref_obj,resolution = res,algorithm=algorithm,verbose=T)
      }
    #}
    gedit_results <- evaluate_deconvolution(ref_obj,query_obj,"gedit")

    print(paste0("Res:",res,",",length(levels(ref_obj$seurat_clusters)),",",gedit_results))

#    return(c(MAE,RSE,SMAPE,RMSE,AVP_Z,
#             EVP_Z,AE,AE_CC,LM,ACE,
#             ACE_Random,MAE_RANDOM,ACE_Boot,PC))

    values_mae = c(values_mae,gedit_results[1])
    values_rse = c(values_rse,gedit_results[2])
    values_smape = c(values_smape,gedit_results[3])
    values_rmse = c(values_rmse,gedit_results[4])
    values_actual_zero = c(values_actual_zero,gedit_results[5])
    values_predicted_zero = c(values_predicted_zero,gedit_results[6])
    values_ae = c(values_ae,gedit_results[7])
    values_ae_cc = c(values_ae_cc,gedit_results[8])
    values_lm_res = c(values_lm_res,gedit_results[9])
    values_ACE = c(values_ACE,gedit_results[10])
    values_ACE_random = c(values_ACE_random,gedit_results[11])
    values_MAE_random = c(values_MAE_random,gedit_results[12])
    values_PC = c(values_PC,gedit_results[13])
    clusters = c(clusters,length(levels(ref_obj$seurat_clusters)))
    resolution = c(resolution,res)
    plot(values_ACE_random~resolution,col="red",ylim=c(0,max(values_ACE_random)))
    points(values_ACE_random-values_ACE~resolution,col="green")
    points(values_ACE~resolution,col="blue")
  }
  df <- data.frame("MAE"= values_mae, "RSE" = values_rse,
                   "SMAPE" = values_smape, "RMSE" = values_rmse,
                   "AZERO" = values_actual_zero, "PZERO"= values_predicted_zero,
                   "Clusters" = clusters,"AE"=values_ae,"AECC"=values_ae_cc,
                   "LMRES"=values_lm_res,"ACE"=values_ACE,"ACE_Random"=values_ACE_random,
                   "MAE_Random"=values_MAE_random,"Resolution"=resolution,
                   "PropCorr"=values_PC)
  return(df)
}

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
  message("Processing predictions")

  ### next we gather the various metrics and report back
  actual_proportion <- get_cluster_proportions(ref_obj)
  cluster_cell_counts <- get_cluster_cell_counts(ref_obj)
  random_proportions <- get_random_proportions(ref_obj)
  ## verify the row orders are equivalent
  print(estimated_proportions)
  if(rownames(estimated_proportions) == "all"){
    rownames(estimated_proportions) <- rownames(actual_proportion)
  }

  estimated_proportions <- estimated_proportions[rownames(actual_proportion),]
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
  #LM <- calculate_linear_regression(actual_proportion,estimated_proportions)
  LM <- calculate_cell_linear_regression(ref_obj,actual_proportion,estimated_proportions)
  ACE <- calculate_absolute_cell_error(ref_obj,actual_proportion,estimated_proportions)
  PC <- calculate_proportinal_correlation(ref_obj,actual_proportion,estimated_proportions)

   message("Bootstrapping 500 times a random background calculation")
  ACE_Boot <- c()
  for (boot in seq(1,500)){
  ACE_Boot <- rbind(ACE_Boot,
                    calculate_absolute_cell_error(ref_obj,
                                                  actual_proportion,
                                                  get_random_proportions(ref_obj)))
  }
  ACE_Random <- median(ACE_Boot)
  message("Bootstrapping the random ACE background calculation finished")

  message("Deconvolution results in: ",MAE,",",RSE,",",SMAPE,",",
          RMSE,",",AVP_Z,",",EVP_Z,",",AE,",",AE_CC,",",LM,",",
          ACE,",",ACE_Random,",",MAE_RANDOM,", PC:",PC)

  return(c(MAE,RSE,SMAPE,RMSE,AVP_Z,
           EVP_Z,AE,AE_CC,LM,ACE,
           ACE_Random,MAE_RANDOM,PC))
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

get_cluster_proportions <- function(ref_obj){
  proportions <- prop.table(table(ref_obj$orig.ident,ref_obj$seurat_clusters),margin=1)
  return(proportions)
}

get_cluster_cell_counts <- function(ref_obj){
  return(as.vector(table(ref_obj$orig.ident,ref_obj$seurat_clusters)))
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

PlotDRRSD <- function(df,xaxis="cluster"){
  if (xaxis == "cluster"){
  optimal_cluster = order((df$ACE_Random-df$ACE),decreasing = T)[1]
  optimal_cluster = df$Clusters[optimal_cluster]
  plot(df$ACE_Random~df$Clusters,col="red",
       ylim=c(min(c(df$ACE_Random,df$ACE_Random-df$ACE)),max(df$ACE_Random) + (max(df$ACE_Random) * .5)))
  points((df$ACE_Random-df$ACE)~df$Clusters,col="darkgreen")
  points(df$ACE~df$Clusters,col="blue")
  abline(v=optimal_cluster,h=max((df$ACE_Random-df$ACE)),lty=3,col="orange")
  legend("top",inset=c(-0.2,0),legend=c("Background ACE", "GEDIT3 ACE","DRRSD Score"),fill = c("red","blue","darkgreen"),xpd=TRUE)
  } else {
    optimal_cluster = order((df$ACE_Random-df$ACE),decreasing = T)[1]
    optimal_cluster = df$Resolution[optimal_cluster]
    plot(df$ACE_Random~df$Resolution,col="red",
         ylim=c(min(c(df$ACE_Random,df$ACE_Random-df$ACE)),max(df$ACE_Random) + (max(df$ACE_Random) * .5)))
    points((df$ACE_Random-df$ACE)~df$Resolution,col="darkgreen")
    points(df$ACE~df$Resolution,col="blue")
    abline(v=optimal_cluster,h=max((df$ACE_Random-df$ACE)),lty=3,col="orange")
    legend("top",inset=c(-0.2,0),legend=c("Background ACE", "GEDIT3 ACE","DRRSD Score"),fill = c("red","blue","darkgreen"),xpd=TRUE)
  }

}

