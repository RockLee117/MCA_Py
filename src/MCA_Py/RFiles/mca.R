#' Run Multiple Correspondence Analysis
#' 
#' @description RunMCA allows to compute the Multiple Corespondence Analysis on the single cell data contained in Seurat or SingleCellExperiment. 
#' MCA is a statistical technique close to PCA that provides a simultaneous 
#' representation of observations (e.g. cells) and variables (e.g. genes) in low-dimensional space.
#' The barycentric relation among cells and genes is a distinctive feature of MCA biplots 
#' and represents a major advantage as compared to other types of biplots such as those produced by Principal Component Analysis 
#' as well as over alternative low-dimensional transformations providing only cell projections. 
#' Thus, in the MCA biplot, analytical distances can be calculated not only between cells and between genes, 
#' but also between each cell and each gene in order to estimate its association. 
#' Thus, the closer a gene g is to a cell c, the more specific to such a cell it can be considered. 
#' Gene-to-cell distances can then be ranked for each individual cell, 
#' and the top-ranked genes may be regarded as a unique gene signature representing the identity card of the cell.
#' 
#'
#' @param X Seurat, SingleCellExperiment or matrix object
#' @param nmcs number of components to compute and store, default set to 30
#' @param features character vector of feature names. If not specified all features will be taken.
#' @param reduction.name name of the reduction default set to 'MCA' for SingleCellExperiment and mca
#' @param slot Which slot to pull expression data from? Default to logcounts for SingleCellExperiment and data for Seurat.
#' @param ... other aruments passed to methods
#'
#' @return Seurat or SCE object with MCA calculation stored in the reductions slot.
#' @importClassesFrom Seurat Seurat
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' @export
#'
#' @examples
#' seuratPbmc <- RunMCA(seuratPbmc, nmcs = 5)
RunMCA <- function(X, nmcs, features, reduction.name, slot, ...) {
    UseMethod("RunMCA", X)
}

#' @rdname RunMCA
#' @export
RunMCA.matrix <- function(X, nmcs = 50, features = NULL, reduction.name = "MCA", ...) {
    library(Rcpp)
    library(RcppArmadillo)
    
    sourceCpp("src/mca.cpp")

    # preprocessing matrix
    # ----------------------------------------------------
    if (!is.null(features)) {
        X <- X[features, ]
    }
    X <- as.matrix(X)
    X <- X[rowVars(X) != 0, ]
    X <- X[str_length(rownames(X)) > 0, ]
    X <- X[!duplicated(rownames(X)), ]
    cellsN <- colnames(X)
    featuresN <- rownames(X)
    tic()
    # print(cellsN)
    # print(featuresN)
    print(dim(X))
    message("Computing Fuzzy Matrix")
    MCAPrepRes <- MCAStep1(X)
    
    print("Z")
    print(dim(MCAPrepRes$Z))
    # print(MCAPrepRes$Z)
    print("Dc")
    print(dim(MCAPrepRes$Dc))
    # print(MCAPrepRes$Dc)
    
    # # Check if the file already exists, and remove it if it does
    # if (file.exists("Results_csv/R_MCAStep1_Z.csv")) {
    #   file.remove("Results_csv/R_MCAStep1_Z.csv")
    # }
    # if (file.exists("Results_csv/R_MCAStep1_Dc.csv")) {
    #   file.remove("Results_csv/R_MCAStep1_Dc.csv")
    # }
    # 
    # write.csv(MCAPrepRes$Z, file = "Results_csv/R_MCAStep1_Z.csv", row.names = FALSE)
    # write.csv(MCAPrepRes$Dc, file = "Results_csv/R_MCAStep1_Dc.csv", row.names = FALSE)
    
    toc()
    message("Computing SVD")
    tic()
    SVD <- irlba::irlba(A = MCAPrepRes$Z, nv = nmcs + 1, nu = 1)[seq(3)]
    
    # print("SVD results: ")
    # print("U")
    # print(dim(SVD$u))
    # # print(SVD$u)
    # print("S")
    # print(length(SVD$d))
    # # print(SVD$d)
    print("Vh")
    print(dim(SVD$v))
    # print(SVD$v)
    
    # if (file.exists("Results_csv/R_SVD_U.csv")) {
    #   file.remove("Results_csv/R_SVD_U.csv")
    # }
    # if (file.exists("Results_csv/R_SVD_S.csv")) {
    #   file.remove("Results_csv/R_SVD_S.csv")
    # }
    # if (file.exists("Results_csv/R_SVD_Vh.csv")) {
    #   file.remove("Results_csv/R_SVD_Vh.csv")
    # }
    # 
    # write.csv(SVD$u, file = "Results_csv/R_SVD_U.csv", row.names = FALSE)
    # write.csv(SVD$d, file = "Results_csv/R_SVD_S.csv", row.names = FALSE)
    # write.csv(SVD$v, file = "Results_csv/R_SVD_Vh.csv", row.names = FALSE)
    
    toc()
    message("Computing Coordinates")
    tic()
    MCA <- MCAStep2(Z = MCAPrepRes$Z, V = SVD$v[, -1], Dc = MCAPrepRes$Dc)
    component <- paste0(reduction.name, "_", seq(ncol(MCA$cellsCoordinates)))
    colnames(MCA$cellsCoordinates) <- component
    colnames(MCA$featuresCoordinates) <- component
    rownames(MCA$cellsCoordinates) <- cellsN
    rownames(MCA$featuresCoordinates) <- featuresN
    
    # print(length(rownames(MCA$cellsCoordinates)))
    # print(rownames(MCA$cellsCoordinates))
    # print(length(rownames(MCA$featuresCoordinates)))
    # print(rownames(MCA$featuresCoordinates))
    
    #print coordinates
    # print("Features Coordinates:")
    # print(dim(MCA$featuresCoordinates))
    # print(MCA$featuresCoordinates[199,])
    # # print("Cells Coordinates:")
    # # print(dim(MCA$cellsCoordinates))
    # print(MCA$cellsCoordinates[43,])
    print("Features Coordinates:")
    print(dim(MCA$featuresCoordinates))
    print(MCA$featuresCoordinates)
    print("Cells Coordinates:")
    print(dim(MCA$cellsCoordinates))
    print(MCA$cellsCoordinates)
    # 
    # if (file.exists("Results_csv/R_featuresCoordinates.csv")) {
    #   file.remove("Results_csv/R_featuresCoordinates.csv")
    # }
    # if (file.exists("Results_csv/R_cellsCoordinates.csv")) {
    #   file.remove("Results_csv/R_cellsCoordinates.csv")
    # }
    # 
    # write.csv(MCA$featuresCoordinates, file = "Results_csv/R_featuresCoordinates.csv", row.names = FALSE)
    # write.csv(MCA$cellsCoordinates, file = "Results_csv/R_cellsCoordinates.csv", row.names = FALSE)
    
    MCA$stdev <- SVD$d[-1]
    class(MCA) <- "MCA"
    toc()
    return(MCA)
}

#' @rdname RunMCA
#' @param assay Name of Assay MCA is being run on
#' @export
RunMCA.Seurat <- function(X, nmcs = 50, features = NULL, reduction.name = "mca", slot = "data", assay = DefaultAssay(X), ...) {
    InitAssay <- DefaultAssay(X)
    DefaultAssay(X) <- assay
    data_matrix <- as.matrix(GetAssayData(X, slot))
    #print(data_matrix)
    MCA <- RunMCA(X = data_matrix, nmcs = nmcs, features = features)
    geneEmb <- MCA$featuresCoordinates
    cellEmb <- MCA$cellsCoordinates
    # print("stdev")
    # print(length(MCA$stdev))
    # print(MCA$stdev)
    stdev <- MCA$stdev
    X <- setDimMCSlot(X = X, cellEmb = cellEmb, geneEmb = geneEmb, stdev = stdev, reduction.name = reduction.name)
    DefaultAssay(X) <- InitAssay
    return(X)
}

#' @rdname RunMCA
#' @export
RunMCA.SingleCellExperiment <- function(X, nmcs = 50, features = NULL, reduction.name = "MCA", slot = "logcounts", ...) {
    data_matrix <- as.matrix(SummarizedExperiment::assay(X, slot))
    MCA <- RunMCA(X = data_matrix, nmcs = nmcs, features = features, reduction.name = reduction.name)
    geneEmb <- MCA$featuresCoordinates
    cellEmb <- MCA$cellsCoordinates
    stdev <- MCA$stdev
    X <- setDimMCSlot(X, cellEmb = cellEmb, geneEmb = geneEmb, stdev = stdev, reduction.name = reduction.name)
    return(X)
}

# SetDimSlot
# --------------------------------------------------------------

#' SetDimSlot
#'
#' @description Integrate MCA in Seurat and SingleCellExperiment Dimensionlity reduction Slot.
#' It will set also a small parameter inside the dimensionality reduction object to signal if it is a MCA or not.
#'
#' @param X Seurat or SingleCellExperiment object
#' @param cellEmb cell coordinates returned by MCA
#' @param geneEmb feature coordinates returned by MCA
#' @param stdev eigen value returned by MCA
#' @param reduction.name name of the created dimensionlaity reduction, default set to 'mca' for Seurat and 'MCA' for SCE.
#' @param ... other arguments passed to methods
#'
#' @return Seurat or SingleCellExperiment object with MC stored in the reduction slot
setDimMCSlot <- function(X, cellEmb, geneEmb, stdev, reduction.name, ...) {
    UseMethod("setDimMCSlot", X)
}

#' @rdname setDimMCSlot
#' @param assay Seurat assay slot
setDimMCSlot.Seurat <- function(X, cellEmb, geneEmb, stdev = NULL, reduction.name = "mca", assay = DefaultAssay(X), ...) {
    colnames(cellEmb) <- paste0(reduction.name, "_", seq(ncol(cellEmb)))
    colnames(geneEmb) <- paste0(reduction.name, "_", seq(ncol(geneEmb)))
    DimReducObject <- CreateDimReducObject(embeddings = cellEmb, loadings = geneEmb, key = paste0(reduction.name, "_"), assay = assay)
    X@reductions[[reduction.name]] <- DimReducObject
    if (!is.null(stdev)) {
        X@reductions[[reduction.name]]@stdev <- sqrt(stdev)
    }
    X@reductions[[reduction.name]]@misc[["mca.flag"]] <- TRUE
    return(X)
}

#' @rdname setDimMCSlot
setDimMCSlot.SingleCellExperiment <- function(X, cellEmb, geneEmb, stdev = NULL, reduction.name = "MCA", ...) {
    reducedDim(X, reduction.name) <- cellEmb
    attr(reducedDim(X, reduction.name), "genesCoordinates") <- geneEmb
    attr(reducedDim(X, reduction.name), "mcaFlag") <- TRUE
    if (!is.null(stdev)) {
        attr(reducedDim(X, reduction.name), "stdev") <- stdev
    }
    return(X)
}
