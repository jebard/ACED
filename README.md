# ACED
Advances in cellular deconvolution algorithms have artfully leveraged single-cell resolution RNA sequencing datasets to predict the cellular composition of bulk tumors. Gene expression estimates at the cell-type level enable dissection of the tumor-microenvironment (TME) from large datasets like that of The Cancer Genome Atlas (TCGA), enabling analysis of the cellular-crosstalk within tumors at scale. However, our interpretation of the TME is contingent upon the accuracy and resolution of the reference single-cell dataset, and thus subsequently the deconvoluted gene expression profiles. In the case of head and neck squamous cell carcinomas (HNSCC), which are a diverse set of tumors with a range of differing etiological origins, the choice of reference scRNA is imperative when attempting to study the nuances within the TME. To address this challenge, we have developed an algorithm to identify the optimal resolution for a single-cell RNA reference set to enhance downstream deconvolutional robustness.

## Dependencies
1. GEDIT3 or MUSIC deconvolution algorithm and dependencies
2. Seurat R library and an scRNA Seurat reference object
3. Metrics R library
4. Reticulate R library
5. Biobase R library
6. Devtools R library 

## For GEDIT3 Deconvolution
GEDIT3 comes prepackages within ACED, so there is no need to download the actual algorithm itself. However GEDIT3 requires the python packages:

```
  random
  numpy
```

In addition, ACED assumes that your python installation is seen by your R. In R Studio this can be set by going to "Tools > Global Options > Python > Setting interpreter" You can test that this is working using the R pacakage "reticulate" function py_config()$python


## ACED INSTALL
```install_github("jebard/ACED",force = T)
library(ACED)
```


## To Run ACED on a Seurat v3 or greater data object

ACED currently assumes that your Seurat object sample names are stored in Seurat.Object$orig.ident. This field is set when creating the initial object. Please verify that your sample names are properly set in this field.
```
aced.res <- ACED(Seurat.Object,start = 0.01,stop=.5,step=0.1)
```
In testing, we have found that starting with a relative broad range between start and stop, with fairly big steps help identify general trends in the data, and then subsequent tests can fine tune the exact range and step size. 

#### Plotting the default resolution view
```
PlotACED(aced.res,title="ACED")
```
#### Calculate and plot the most optimal UMAP clustering resolution
```
Seurat.Object <- FindClusters(Seurat.Object,resolution = get_optimal_resolution(aced.res))
DimPlot(Seurat.Object)
```

#### Calculate and plot the least optimal UMAP clustering resolution
```
Seurat.Object <- FindClusters(Seurat.Object,resolution = get_least_optimal_resolution(aced.res))
DimPlot(Seurat.Object)
```

