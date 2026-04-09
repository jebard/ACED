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
  bulk.eset     <- music_prep_query(ref_obj = ref_object)
  bulk.eset.mat <- Biobase::exprs(bulk.eset)  # G x M matrix

  # Single-sample workaround: MuSiC uses colMeans() and colSums() internally,
  # which require at least 2 columns. When bulk.mtx has only 1 column, R may
  # drop the matrix class after division, producing a vector that fails colMeans.
  # Fix: duplicate the single sample so MuSiC operates on a G x 2 matrix.
  # Both columns are identical, so the returned proportions are identical --
  # we take the first row and return it as a single-sample result.
  single_sample <- ncol(bulk.eset.mat) == 1
  if (single_sample) {
    message("run_music: single-sample bulk matrix detected -- duplicating sample ",
            "to satisfy MuSiC multi-sample requirement. ",
            "Note: MuSiC's variance-weighting scheme is not meaningful for single samples; ",
            "consider strategy='lasso' or strategy='glm' for single-sample objects.")
    original_colname  <- colnames(bulk.eset.mat)
    bulk.eset.mat     <- cbind(bulk.eset.mat, bulk.eset.mat)
    colnames(bulk.eset.mat) <- c(original_colname,
                                 paste0(original_colname, "_dup"))
  }

  print("Running MUSIC ...")
  estimated.prop <- music_prop(bulk.mtx = bulk.eset.mat,
                               sc.sce   = SC.eset,
                               clusters = 'seurat_clusters',
                               samples  = 'orig.ident', verbose = T)

  result <- as.data.frame(estimated.prop$Est.prop.weighted)

  # If we duplicated the sample, drop the duplicate row before returning
  if (single_sample) {
    result <- result[1, , drop = FALSE]
    rownames(result) <- colnames(Biobase::exprs(bulk.eset))  # restore original name
  }

  return(result)
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
