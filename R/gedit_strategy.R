library(Seurat)

### gedit Template Strategy Pattern
### GEDIT is a python based tool that solves deconvolution using  ...
###
### python2.7 GEDIT2.py -mix MixtureFullRefOrig.csv  -ref FullRefOrig.csv

run_gedit <- function(ref_obj,query_obj,res){
  message("Please cite GEDIT3 Nadel et al. https://doi.org/10.1093/gigascience/giab002")
  message("Preparing the reference object for GEDIT3")
  gedit_prep_reference(ref_obj,res)
  message("Preparing the query object for GEDIT3")
  gedit_prep_query(query_obj)
  message("Running GEDIT3 deconvolution...")
  execute_gedit(res)
  message("Gathering GEDIT3 results")
  predictions <- gedit_gather_results()
  return(predictions)
}

gedit_prep_reference <- function(ref_obj,res){
  ### save out the reference GEDIT object
  Seurat::Idents(ref_obj) <- "seurat_clusters"
  write.csv(file=paste0("RefObj.",res,".csv"),Seurat::AverageExpression(ref_obj,assays = "RNA",slot = "counts")$RNA)
}

gedit_prep_query <- function(query_obj){
  ### save out the mixture sets for each original sample
  Seurat::Idents(query_obj) <- "orig.ident"
  write.csv(file="MixtureQuery.csv",Seurat::AverageExpression(query_obj,assays = "RNA",slot = "counts")$RNA)
}

execute_gedit <- function(res){
  message("Running GEDIT3 using the following settings:")
  message(paste0(py_config()$python," ",package_info("DRRSD")$path,"/GEDIT3.py -mix $PWD/MixtureQuery.csv -ref $PWD/RefObj.",res,"csv -outFile $PWD/GEDIT_Deconv"))
  system(paste0(py_config()$python," ",package_info("DRRSD")$path,"/GEDIT3.py -mix $PWD/MixtureQuery.csv -ref $PWD/RefObj.",res,"csv -outFile $PWD/GEDIT_Deconv"),TRUE)
}

gedit_gather_results <- function(){
  ### TODO: Need to update this to a path within the DRRSD package where the output of deconv will go.
  predictions = read.table(file="GEDIT_Deconv_CTPredictions.tsv",header = TRUE, row.names = 1, sep = "\t")
  return(predictions)
}

