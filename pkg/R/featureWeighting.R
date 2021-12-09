#' Feature Weighting
#'
#' Calculate Hellinger Distance (HD) between K-Mean discretized gene expression and reference
#' sample label. Genes with high HD are discriminant against sample labels. We can control whether
#' to filter the reference data by HD scores using the 'cut' and the 'n.cut' option.
#'
#' @param reference Required. Reference sce
#' @param clusterid Required. Variable name in the metadata that stores sample labels.
#' @param cut Default: TRUE Whether to filter the reference data to the gene sets with the 'n.cut' highest HD scores.
#' @param n.cut Default: 2000. How many genes to keep in the reference data.
#' @param assay2rank Default: 'data'. On the data of which assay we perform rank transformation and then HD calculation?
#' @param rowDataPrefix Default: 'HD_'. Prefix added to the column names of HD scores, which will be stored in the RowData of the reference sce.
#' @return SingleCellExperiment Object containing HD scores in the metadata and the RowData.
#' @export
featureWeighting <- function(reference, clusterid, cut = T, n.cut = 2000, assay2rank = 'data', rowDataPrefix = 'HD_'){

  rownames(reference) <-gsub("\\.",'-',rownames(reference))
  cluster <- colData(reference)[,clusterid]
  cluname <- sort(unique(cluster))
  n.cluster <- length(cluname)

  #remove genes with less than n.clusters expressed
  reference <- reference[apply(assay(reference, assay2rank),1,function(x) length(unique(x)))>=n.cluster,]

  n.sample <- ncol(reference)
  n.gene <- nrow(reference)

  #create rank assay
  reference <- rankTrans(reference, assay2rank)


  #discreteize reference data
  discRankDat <- apply(assay(reference, 'rank'), 1,
                       function(x) as.factor(kmeans(x, centers = n.cluster, iter.max = 100)$cluster))
  #calculate HD
  HD <- data.frame(t(apply(discRankDat,2,function(y) HellingerDist(cluster,y))))
  #generate metadata
  HD <- dplyr::mutate(HD, gene = rownames(HD),celltype = cluname[apply(HD, 1, which.max)],
                      mean = rowMeans(HD))
  #order the HD
  HD <- dplyr::arrange(HD, desc(mean))

  #Store HD in the metadata
  metadata(reference)[['HD']] <- HD

  colnames(HD) <- paste(rowDataPrefix,colnames(HD), sep = '')

  #filter the reference
  if(cut){
    if(n.cut > n.gene){
      warning(paste('Exceeding number of features were selected. Select a number smaller than', n.gene))
      n.cut <- n.gene
    }
    reference <- reference[HD[1:n.cut,n.cluster+1],]
    rowData(reference) <- cbind(rowData(reference), HD[1:n.cut,])

  }else{
    reference <- reference[HD[,n.cluster+1],]
    rowData(reference) <- cbind(rowData(reference),HD)

  }

  reference <- rankTrans(reference, 'rank')
  reference

}
