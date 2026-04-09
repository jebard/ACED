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
  bulk.eset.mat <- Biobase::exprs(bulk.eset)
  # Estimate cell type proportions
  print("Running MUSIC ...")
  estimated.prop = music_prop(bulk.mtx = bulk.eset.mat,
                              sc.sce = SC.eset,
                              clusters = 'seurat_clusters',
                              samples = 'orig.ident', verbose = T)
  return(as.data.frame(estimated.prop$Est.prop.weighted))
}


music_prep_reference <- function(music.input=NULL){
  SC.sce = as.SingleCellExperiment(music.input,assay="RNA")
  return(SC.sce)
}

music_prep_reference_old <- function(music.input=NULL){
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
  # Detect Seurat version for slot/layer compatibility
  ver <- tryCatch(utils::packageVersion("Seurat"), error = function(e) "4.0.0")

  n_samples <- length(unique(ref_obj$orig.ident))

  if (n_samples == 1) {
    # Single-sample: Seurat v5 ignores group.by when all cells share one identity
    # and returns a matrix with an unpredictable column name. Bypass AverageExpression
    # entirely and compute rowMeans directly from the count matrix.
    message("music_prep_query: single-sample object detected -- ",
            "computing average expression directly from count matrix.")
    sample_name <- as.character(unique(ref_obj$orig.ident))
    if (ver >= "5.0.0") {
      counts <- as.matrix(Seurat::GetAssayData(ref_obj, assay = "RNA", layer = "counts"))
    } else {
      counts <- as.matrix(Seurat::GetAssayData(ref_obj, assay = "RNA", slot  = "counts"))
    }
    # Build G x 1 matrix with the correct sample name as column
    scrna.mat <- matrix(rowMeans(counts, na.rm = TRUE), ncol = 1,
                        dimnames = list(rownames(counts), sample_name))
  } else {
    # Multi-sample: use AverageExpression with correct slot/layer for Seurat version.
    # group.by="orig.ident" is used explicitly (not Idents) to avoid any active-ident
    # conflicts with single-sample objects.
    if (ver >= "5.0.0") {
      scrna.exp <- Seurat::AverageExpression(ref_obj, assays = "RNA",
                                             layer = "counts",
                                             group.by = "orig.ident")
    } else {
      scrna.exp <- Seurat::AverageExpression(ref_obj, assays = "RNA",
                                             slot = "counts",
                                             group.by = "orig.ident")
    }
    # Coerce to plain dense matrix -- Seurat v5 may return a sparse dgeMatrix
    # which Biobase::ExpressionSet() cannot accept (requires plain numeric matrix).
    scrna.mat <- as.matrix(scrna.exp$RNA)  # G x M
  }

  # Guard: replace any NA values before passing to ExpressionSet
  if (any(is.na(scrna.mat))) {
    scrna.mat[is.na(scrna.mat)] <- 0
    message("music_prep_query: NA values in query matrix replaced with 0.")
  }

  pseudobulk.eset <- Biobase::ExpressionSet(assayData = scrna.mat)
  return(pseudobulk.eset)
}
