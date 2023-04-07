#' Title : run_music
#' Wrapped script to execute MuSiC deconvolution of a given dataset
#' @param music.input Seurat scRNA object containing the reference set
#' @param bulk.eset Bulk expressionSet object for deconvolution
#'
#' @return
#' @export
#'
#' @examples
run_music <- function(ref_object){
  message("Please cite MUSIC Wang et al. https://doi.org/10.1038/s41467-018-08023-x")
  message("Preparing the reference object for MUSIC")
  SC.eset <- music_prep_reference(music.input = ref_object)
  message("Preparing the query object for MUSIC")
  bulk.eset <- music_prep_query(ref_obj = ref_object)
  # Estimate cell type proportions
  print("Running MUSIC ...")
  estimated.prop = music_prop(bulk.eset = exprs(bulk.eset),
                              sc.sce = SC.eset,
                              clusters = 'cellType',
                              samples = 'SubjectName', verbose = T)
  return(as.data.frame(estimated.prop$Est.prop.weighted))
}

music_prep_reference <- function(music.input=NULL){
  scCounts = music.input@assays$RNA@counts
  individual.labels = music.input@assays$RNA@data@Dimnames[[2]]
  #individual.labels = colnames(MuWT.Pseudotime@assays$integrated@data)
  cell.type.labels = music.input$seurat_clusters
  music.input$names <- music.input$orig.ident

  sc.pheno <- data.frame(check.names=F, check.rows=F,
                         stringsAsFactors=F,
                         SubjectID = individual.labels,
                         SubjectName= music.input$orig.ident,
                         cellType=cell.type.labels)

  cell.type.char = as.character(sc.pheno[,1])

  metadata = data.frame(labelDescription = c('SubjectID','SubjectName','CellType'),
                        row.names = c('SubjectID','SubjectName','CellType'))
  print("Creating single-cell expression set")
  SC.eset =  ExpressionSet(assayData = data.matrix(scCounts),
                           phenoData =  new("AnnotatedDataFrame",
                                            data = sc.pheno,
                                            varMetadata = metadata))
  return(SC.eset)
}

music_prep_query <- function(ref_obj=NULL){
  Idents(ref_obj) <- "orig.ident"
  scrna.exp <- AverageExpression(ref_obj,assays = "RNA",slot = "counts")
  # row names identify features and column names identify samples.
  scrna.exp <- scrna.exp$RNA
  pseudobulk.eset = ExpressionSet(assayData = scrna.exp)
  return(pseudobulk.eset)
}
