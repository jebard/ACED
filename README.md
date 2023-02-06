# DRRSD

## INSTALL
install_github("jebard/DRRSD",force = T)
library(DRRSD)

## Dependencies
1. GEDIT or MUSIC deconvolution algorithms
2. Seurat and Seurat reference object
3. Metrics
4. Retriculate
5. Biobase

## To Run DRRSD on a Seurat v3 or greater data object
drrsd.res <- DRRSD(combined.seurat.sct,combined.seurat.sct,start = 0.01,stop=.5,step=0.025)
PlotDDRSD(drrsd.res,xaxis="cluster")
PlotDDRSD(drrsd.res,xaxis="resolution")
