#' Aggregate
#'
#' Aggregate single cells to create pseudo-bulk samples.
#'
#' @param query Required. Query sce.
#' @param clusterid Required. The query Coldata attribute name. Cells will be aggregated within the clusters of 'clusterid'.
#' @param assay Default: 'counts'. On which assay of the query we perform aggregation.
#' @param size.factor Default: 1. Determine how many pseudo-bulk samples to be generated from each cluster. The number equals to "size.factor times cluster size".
#' @param pool.factor Default: 1  Determine how many cells to be pooled for aggregation in each cluster. The number equals to "pool.factor times cluster size"
#' @param nPool Optional. Pooling size of each pseudo-bulk sample. Will suppress pool.factor if specified.
#' @return A confusion matrix if trueAnno has been provided.
#' @export
aggregate <- function(query, clusterid, assay = 'counts', size.factor =1, pool.factor = 1, nPool = NULL){

  cluster <- colData(query)[,clusterid]
  cluname <- sort(unique(cluster))
  Ncluster <- length(cluname)
  G <- nrow(query)

  outList <- list()
  outCluster <- c()

  #Aggregate cells
  for(i in cluname){
    message(paste("\r Now aggregate", i), appendLF = F)
    subQuery <- Matrix::Matrix(assay(query, assay)[, cluster == i], sparse = T)
    nCell <- ncol(subQuery)
    nOut <- ceiling(nCell * size.factor)
    if(is.null(nPool)) nPool <- ceiling(nCell * pool.factor)

    subOut <- matrix(0, nrow = G, ncol = nOut)
    colnames(subOut) <- paste(i, '.', 1:nOut, sep = '')

    sampleIndex <- sample(1:nCell, nOut * nPool, replace = T) %>%
      matrix(ncol = nPool, nrow = nOut)

    if(nPool > 1){
      for(j in 1:nOut) subOut[,j] <- rowSums(subQuery[,sampleIndex[j,]])/nPool
    }else{
      for(j in 1:nOut) subOut[,j] <- subQuery[,sampleIndex[j,]]/nPool
    }

    outList[[i]] <- subOut
    outCluster <- c(outCluster, rep(i, nOut))
  }

  #Collapse the list
  outList <- unname(outList)
  outData <- do.call(cbind,outList)
  rownames(outData) <- rownames(query)
  sce <- SingleCellExperiment(list(data = outData))
  colData(sce)[,clusterid] <- outCluster

  #remove duplicated pseudo-bulk samples.
  dup <- duplicated(t(outData)) | duplicated(t(outData)[ncol(outData):1,])[ncol(outData):1]
  outData <- outData[,!dup]
  outCluster <- outCluster[!dup]

  message(paste('\n Sparsity after imputation is', round(mean(outData == 0),3)))
  sce

}
