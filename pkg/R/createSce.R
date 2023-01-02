#' Create SingleCellExperiment Object
#'
#' This function converts single cell expression data matrix and metadata to a
#' SingleCellExperiment object (sce).
#'
#' @param counts Optional. A gene by cell count matrix.
#' @param data Optional. A gene by cell expression matrix (can be normalized count).
#' @param colData Optional. A data frame storing sample metadata. The rownames (cell ids) of colData should match the colnames of counts/data matrix.
#' @param as.sparse Default: TRUE. Whether to convert counts/data matrix to sparse matrix format. Should be set to FALSE if matrices are non-sparse (not zero-inflated).
#' @return SingleCellExperiment Object.
#' @export
createSce <- function(counts = NULL, data = NULL, colData = NULL, as.sparse = TRUE){
  data.list <- list()
  if(!is.null(counts)){
    if(as.sparse) data.list[['counts']] <- Matrix::Matrix(as.matrix(counts),sparse = T) else data.list[['counts']] <- as.matrix(counts)
  }
  if(!is.null(data)){
    if(as.sparse) data.list[['data']] <- Matrix::Matrix(as.matrix(data),sparse = T) else data.list[['data']] <- as.matrix(data)
  }
  if(!is.null(colData)) sce <- SingleCellExperiment(data.list, colData = colData) else sce <- SingleCellExperiment(data.list)

  sce
}
