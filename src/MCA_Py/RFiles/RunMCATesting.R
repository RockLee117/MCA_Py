
if(!require("tidyverse")) install.packages("tidyverse")
if(!require("ggpubr")) install.packages("ggpubr")
devtools::install_github("RausellLab/CelliD", ref = "legacy")
library(CellID) # Note that the legacy version is called CellID and not CelliD

## Warning: replacing previous import 'data.table::shift' by 'tictoc::shift' when
## loading 'CelliD'
library(tidyverse)
library(ggpubr) #library for plotting



BaronMatrix   <- readRDS(url("https://storage.googleapis.com/cellid-cbl/BaronMatrix.rds"))
BaronMetaData <- readRDS(url("https://storage.googleapis.com/cellid-cbl/BaronMetaData.rds"))


#print(BaronMatrix@x)

data("HgProteinCodingGenes")
BaronMatrixProt <- BaronMatrix[rownames(BaronMatrix) %in% HgProteinCodingGenes,]
# print(BaronMatrixProt)
Baron <- CreateSeuratObject(counts = BaronMatrixProt, project = "Baron", min.cells = 5, meta.data = BaronMetaData)
# print(Baron)
Baron <- NormalizeData(Baron)
Baron <- ScaleData(Baron, features = rownames(Baron))



# this is to turn seurat object into readable object in python
# install.packages("remotes")
# remotes::install_github("mojaveazure/seurat-disk")
# library(Seurat)
# library(SeuratDisk)
# SaveH5Seurat(Baron, filename = "Baron.h5Seurat")

# Convert("Baron.h5Seurat", dest="h5ad")
######

#converting data that is actually passed into runmca()
# Baron_data is a dgcMatrix or a SPARSE MATRIX
Baron_data_assay = Baron@assays[["RNA"]]@data
dim(Baron_data_assay)
Baron_data_assay[6,6]
write.csv(Baron_data_assay, file = "Baron_assay_data.csv", row.names=FALSE)
#export rownames of the assay
Baron_assay_rownames = Baron_data_assay@Dimnames[[1]]
write.csv(Baron_assay_rownames, file = "Baron_assay_rownames.csv", row.names=FALSE)


install.packages("Rcpp")
install.packages("RcppArmadillo")
library(Rcpp)
library(RcppArmadillo)
# evalCpp("1+1")

BaronSmall <- Baron[1:400, 1:200]
#see input data
print(BaronSmall@assays[["RNA"]]@scale.data) # -> this is what is converted into the .h5ad and stored as X in the AnnData object in python 
print(BaronSmall@assays[["RNA"]]@data)  # -> this is what is passed into RunMCA() when calling the RunMCA.Seurat (assay)
# above is the assay for the dataset:
#Assay refers to the primary data matrix that contains the quantified gene expression levels for each cell in the experiment.

# tic toc used to see computation time
install.packages("tictoc")
# Load the tictoc package
library(tictoc)

source("R/mca.R") # if change code in mca.R must re-run this line to update changes here when calling RunMCA()
MCA <- RunMCA(BaronSmall)
MCA <- RunMCA(Baron)
MCA <- RunMCA(Baron, features = Baron_assay_rownames[1:60], nmcs=30)   # testing feature filter
DimPlotMC(MCA, reduction = "mca")

# IMPORTANT NOTE: if specifying a value for nmcs in RunMCA(), must provide same value for dims when calling GetCellGeneDistance()
#distance calculations
source("R/cell.R")
res = GetCellGeneDistance(MCA)
res = GetCellGeneDistance(MCA, dims=seq(30))
print(dim(res))
print(res)
write.csv(res, file = "dummy.csv")


