library(Seurat)

### gedit Template Strategy Pattern
### GEDIT is a python based tool that solves deconvolution using  ...
###
### python2.7 GEDIT2.py -mix MixtureFullRefOrig.csv  -ref FullRefOrig.csv


run_gedit <- function(ref_obj,query_obj){
  message("Preparing the reference object for GEDIT")
  gedit_prep_reference(ref_obj)
  message("Preparing the query object for GEDIT")
  gedit_prep_query(query_obj)
  message("Running GEDIT deconvolution...")
  execute_gedit()
  message("Gathering GEDIT results")
  predictions <- gedit_gather_results()
  return(predictions)
}

gedit_prep_reference <- function(ref_obj){
  ### save out the reference GEDIT object
  Seurat::Idents(ref_obj) <- "seurat_clusters"
  write.csv(file="RefObj.csv",Seurat::AverageExpression(combined.seurat.sct,assays = "RNA",slot = "counts")$RNA)
}

gedit_prep_query <- function(query_obj){
  ### save out the mixture sets for each original sample
  Seurat::Idents(query_obj) <- "orig.ident"
  write.csv(file="MixtureQuery.csv",Seurat::AverageExpression(combined.seurat.sct,assays = "RNA",slot = "counts")$RNA)
}

execute_gedit <- function(){
  #my_env <- BasiliskEnvironment(envname="drrsd_gedit",
  #                             pkgname="DRRSD",
  #                              packages=c("random", "numpy")
  #)
  #proc <- basiliskStart(my_env)
  #on.exit(basiliskStop(proc))
  #basiliskRun(proc) {
  system("python $PWD/src/GEDIT/GEDITv3.0/GEDIT2.py -mix MixtureQuery.csv -ref RefObj.csv", TRUE)
  #}
}

gedit_gather_results <- function(){
  predictions = read.table(file="MixtureQuery.csv_RefObj.csv_50_Entropy_0.0_CTPredictions.tsv",header = TRUE, row.names = 1, sep = "\t")
  #truth <- prop.table(table(combined.seurat.sct$orig.ident,combined.seurat.sct$seurat_clusters),margin=1)
  #AE <- truth - predictions[rownames(truth),]
  return(predictions)
}

