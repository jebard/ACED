# DRRSD
Advances in cellular deconvolution algorithms have artfully leveraged single-cell resolution RNA sequencing datasets to predict the cellular composition of bulk tumors. Gene expression estimates at the cell-type level enable dissection of the tumor-microenvironment (TME) from large datasets like that of The Cancer Genome Atlas (TCGA), enabling analysis of the cellular-crosstalk within tumors at scale. However, our interpretation of the TME is contingent upon the accuracy and resolution of the reference single-cell dataset, and thus subsequently the deconvoluted gene expression profiles. In the case of head and neck squamous cell carcinomas (HNSCC), which are a diverse set of tumors with a range of differing etiological origins, the choice of reference scRNA is imperative when attempting to study the nuances within the TME. To address this challenge, we have developed am algorithm to establish an optimal, and label-free, scRNA reference resolution to improve deconvolutional robustness.

## Dependencies
1. GEDIT or MUSIC deconvolution algorithm and dependencies
2. Seurat R library and Seurat reference object to test
3. Metrics R library
4. Retriculate R library
5. Biobase R library
6. Devtools R library 

## INSTALL
install_github("jebard/DRRSD",force = T)

library(DRRSD)


## To Run DRRSD on a Seurat v3 or greater data object
drrsd.res <- DRRSD(Seurat.Object,start = 0.01,stop=.5,step=0.025)

#### Plotting the default resolution view
PlotDDRSD(drrsd.res,xaxis="resolution")

#### Plotting the optional cluster based view. Multiple resolutions hit the same cluster
PlotDDRSD(drrsd.res,xaxis="cluster")

#### Calculate and plot the most optimal UMAP clustering resolution
Seurat.Object <- FindClusters(Seurat.Object,resolution = get_optimal_resolution(drrsd.res))

DimPlot(Seurat.Object)

#### Calculate and plot the least optimal UMAP clustering resolution
Seurat.Object <- FindClusters(Seurat.Object,resolution = get_least_optimal_resolution(drrsd.res))
DimPlot(Seurat.Object)

