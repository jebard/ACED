# DRRSD

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

