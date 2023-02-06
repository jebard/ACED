# DRRSD

## INSTALL
install_github("jebard/DRRSD",force = T)

library(DRRSD)

## Dependencies
1. GEDIT or MUSIC deconvolution algorithm and dependencies
2. Seurat R library and Seurat reference object to test
3. Metrics R library
4. Retriculate R library
5. Biobase R library

## To Run DRRSD on a Seurat v3 or greater data object
drrsd.res <- DRRSD(Seurat.Object,start = 0.01,stop=.5,step=0.025)

PlotDDRSD(drrsd.res,xaxis="cluster")

PlotDDRSD(drrsd.res,xaxis="resolution")
