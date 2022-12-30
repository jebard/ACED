### gedit Template Strategy Pattern
### GEDIT is a python based tool that solves deconvolution using  ...
###
### python2.7 GEDIT2.py -mix MixtureFullRefOrig.csv  -ref FullRefOrig.csv
library(basilisk)
library(BiocStyle)

gedit <- function(ref_obj,query_obj){
  message("Preparing the reference object for GEDIT")
  gedit_prep_reference()
  message("Preparing the query object for GEDIT")
  gedit_prep_query()
  message("Running GEDIT deconvolution...")
}

gedit_prep_reference <- function(ref_obj){
  ### save out the reference GEDIT object
  Idents(ref_obj) <- "seurat_clusters"
  write.csv(file="RefObj.csv",AverageExpression(combined.seurat.sct,assays = "RNA",slot = "counts")$RNA)
}

gedit_prep_query <- function(query_obj){
  ### save out the mixture sets for each original sample
  Idents(combined.seurat.sct) <- "orig.ident"
  write.csv(file="MixtureQuery.csv",AverageExpression(combined.seurat.sct,assays = "RNA",slot = "counts")$RNA)
}

execute_gedit <- function(){
  ?BasiliskEnvironment()
  my_env <- BasiliskEnvironment(envname="drrsd_gedit",
                                pkgname="DRRSD",
                                packages=c("random", "numpy")
  )
  proc <- basiliskStart(my_env)
  on.exit(basiliskStop(proc))
  basiliskRun(proc) {
    system("python GEDIT2.py -mix MixtureFullRefOrig.csv -ref FullRefOrig.csv", TRUE)
  }
}

gedit_gather_results(){
  predictions = read.table(file="MixtureQuery.csv_RefObj.csv_50_Entropy_0.0_CTPredictions.tsv",header = TRUE, row.names = 1, sep = "\t")
  truth <- prop.table(table(combined.seurat.sct$orig.ident,combined.seurat.sct$seurat_clusters),margin=1)
  AE <- truth - predictions[rownames(truth),]
}

