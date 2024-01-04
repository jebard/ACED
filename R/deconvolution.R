# Deconvolution Strategy Pattern
# Ingests a request to run cellular deconvolution, returns a normalized result regardless of strategy
# @author jbard

#' ACED Define a reference resolution for single-cell deconvolution
#'
#' @param ref_obj Seurat reference scRNA object
#' @param strategy Deconvolution strategy, default="gedit" or "music"
#' @param start minimum resolution to pass to Seurats FindClusters
#' @param stop maximum resolution to pass to Seurats FindClusters
#' @param step increment resolution to step from start to stop resolutions
#' @param algorithm Clustering strategy, "louvain" or "leiden"
#' @param method matrix or igraph strategys if using leiden clustering
#'
#' @return dataFrame of resulting metrics that ACED generates
#' @export ACED
#'
#' @examples ACED(Seurat.object, start = 0.05,stop=.7,step=0.01,algorithm="louvain",strategy="gedit")
ACED <- function(ref_obj=ref_obj,strategy="gedit",start=0.01,stop=1,
                  step=.05,algorithm="louvain",method="matrix"){

  #### instantiate lists to return
  query_obj <- ref_obj
  values_mae = c();values_rse = c();values_smape = c();values_rmse = c()
  values_actual_zero = c();values_predicted_zero = c();clusters = c();resolution = c();
  values_ae = c();values_ae_cc = c();values_lm_res = c()
  values_ACE = c();values_ACE_random = c();values_MAE_random = c();values_PC = c();
  values_background_stdev = c(); values_background_mean = c();

  #### for each of the resolutions we are going to be testing
  for (res in c(seq(from=start,to=stop,by=step))){
    #resolution_string <- paste0("integrated_snn_res.",res)
      message(paste0("Calulating resolution ",res," for reference. Clustering using the ",algorithm," clustering"))

        if(algorithm=="leiden" || algorithm == 4){
        ref_obj <- FindClusters(ref_obj,resolution = res,algorithm=algorithm,verbose=T,method=method)
      } else {
      ref_obj <- FindClusters(ref_obj,resolution = res,algorithm=algorithm,verbose=T)
      }
    #}
    gedit_results <- evaluate_deconvolution(ref_obj,query_obj,strategy,res)

    ##print(paste0("Res:",res,",",length(levels(ref_obj$seurat_clusters)),",",gedit_results))
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
    values_background_mean = c(values_background_mean,gedit_results[14])
    values_background_stdev = c(values_background_stdev,gedit_results[15])
    clusters = c(clusters,length(levels(ref_obj$seurat_clusters)))
    resolution = c(resolution,res)

    ## generate some on the fly plots
    plot(values_ACE_random~resolution,col="red",ylim=c(0,max(values_ACE_random + (values_ACE_random * .5))))
    arrows(x0=resolution, y0=values_background_mean-values_background_stdev,
           x1=resolution, y1=values_background_mean+values_background_stdev,
           code=3, angle=90, length=0.05,col="red",lty=2)
    points(values_ACE_random-values_ACE~resolution,col="darkgreen")
    points(values_ACE~resolution,col="blue")


  }
  df <- data.frame("MAE"= values_mae, "RSE" = values_rse,
                   "SMAPE" = values_smape, "RMSE" = values_rmse,
                   "AZERO" = values_actual_zero, "PZERO"= values_predicted_zero,
                   "Clusters" = clusters,"AE"=values_ae,"AECC"=values_ae_cc,
                   "LMRES"=values_lm_res,"ACE"=values_ACE,"ACE_Random"=values_ACE_random,"ACE_SCORE"=values_ACE_random-values_ACE,
                   "MAE_Random"=values_MAE_random,"Resolution"=resolution,
                   "PropCorr"=values_PC,"BGM"=values_background_mean,"BGSD"=values_background_stdev)
  PlotACED(df)
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
evaluate_deconvolution <- function(ref_obj, query_obj, strategy,res){
  ### first, launch the deconvolution analysis and get the results
  estimated_proportions <- NULL
  if(strategy=="bayesprism"){
    estimated_proportions <- run_bayesprism()
  } else if (strategy=="music") {
    estimated_proportions <- run_music(ref_obj)
  } else if (strategy=="gedit"){
    estimated_proportions <- run_gedit(ref_obj,query_obj,res)
  } else {
    estimated_proportions <- aced_glm(ref_obj)
    warning("Defaulting to ACeD pseudoinverse strategy to perform linear regression-based deconvolution")
  }
  message("Processing predictions")

  ### next we gather the various metrics and report back
  actual_proportion <- get_cluster_proportions(ref_obj)
  write.csv(file=paste0("ACED_ActProp",res,".csv"),actual_proportion)
  cluster_cell_counts <- get_cluster_cell_counts(ref_obj)
  random_proportions <- get_random_proportions(ref_obj)
  ## verify the row orders are equivalent
  print("Actual:")
  print(actual_proportion)

  if(length(rownames(estimated_proportions)) == 1 && strategy == "gedit"){
    print("Adjusting rownames because it is a single-sample in GEDIT")
    rownames(estimated_proportions) <- rownames(actual_proportion)
  }

  ### enforce the same column order
  estimated_proportions <- estimated_proportions[rownames(actual_proportion),]

  if(strategy == "music"){
    ### for some bizarre reason, MUSIC will return clusters out of order
    ### If a given cluster isn't in the subset, set estimated proportion to zeros
    #print(levels(ref_obj$seurat_clusters)[!(levels(ref_obj$seurat_clusters) %in% colnames(estimated_proportions))])
    for (clust in levels(ref_obj$seurat_clusters)[!(levels(ref_obj$seurat_clusters) %in% colnames(estimated_proportions))]) {
      print(paste0("Adding missing cluster: ",clust))
      estimated_proportions[,as.character(clust)] <- 0.0
    }
    estimated_proportions <- estimated_proportions[,colnames(actual_proportion)]
  }

  print("Estimated:")
  print(estimated_proportions)

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

  message("Bootstrapping a random background 500 times")
  ACE_Boot <- c()

  for (boot in seq(1,500)){
  ACE_Boot <- rbind(ACE_Boot,
                    calculate_absolute_cell_error(ref_obj,
                                                  actual_proportion,
                                                  #permute_within_rows(actual_proportion))) ### this line alternatively permutes the matrix
                                                  get_random_proportions(ref_obj)))
  }

  #print("Running up to 1000 permutations to compute backgrounds")
  #permuts <- permutations(n = length(as.numeric(get_cluster_proportions(ref_obj))),
  #                        r = length(as.numeric(get_cluster_proportions(ref_obj))),
  #                        v = as.numeric(get_cluster_proportions(ref_obj)))

  #if (nrow(permuts) > 1000){
  #permuts <- permuts[1:1000,]
  #}

  #for (boot in seq(1,nrow(permuts))){
  #  ACE_Boot <- rbind(ACE_Boot,
  #                    calculate_absolute_cell_error(ref_obj,
  #                                                  actual_proportion,
  #                                                  permuts[boot,]))
  #  }

  background_mean <- mean(ACE_Boot)
  background_stdev <- sd(ACE_Boot)
  ACE_Random <- mean(ACE_Boot)

  message("Deconvolution results in: ",MAE,",",RSE,",",SMAPE,",",
          RMSE,",",AVP_Z,",",EVP_Z,",",AE,",",AE_CC,",",LM,",",
          ACE,",",ACE_Random,",",MAE_RANDOM,", PC:",PC)

  return(c(MAE,RSE,SMAPE,RMSE,AVP_Z,
           EVP_Z,AE,AE_CC,LM,ACE,
           ACE_Random,MAE_RANDOM,PC,background_mean,background_stdev))
}

get_random_proportions <- function(ref_obj){
  actual_proportion <- get_cluster_proportions(ref_obj)
  # first get the number of samples in the object
  samples <- length(unique(ref_obj$orig.ident))
  # next get the number of clusters in the object
  clusts <- length(levels(ref_obj$seurat_clusters))
  m <- matrix(rnorm(samples * clusts,mean(actual_proportion),sd = 1), nrow=samples)
  #m <- matrix(rnorm(samples * clusts,mean(actual_proportion),sd = sd(actual_proportion)), nrow=samples)
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

clear_resolutions <- function(ref_obj){
  for (res in seq(from=0.001, to=3.00,by=.001)){
    resolution_string <- paste0("integrated_snn_res.",res)
    tryCatch(
      expr ={
      ref_obj[[resolution_string]] <- NULL
      }, error = function(e){
      })
  }
  return(ref_obj)
}

calculate_cluster_tree <- function(ref_obj){
  ref_obj <- BuildClusterTree(ref_obj)
  data.tree <- Tool(object = ref_obj, slot = "BuildClusterTree")
  return(ape::cophenetic.phylo(data.tree))
}

