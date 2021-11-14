#' Sincast Capybara cell score
#'
#' Using restricted linear regression to predict cell identity continum.
#'
#' @param reference Required. Reference sce.
#' @param query Required. Query sce.
#' @param clusterid Required. Predict query cell identity based on the reference metadata attribute specified by 'clusterid'.
#' @param k Default: 5. Partition each reference cluster into k sub-clusters.
#' @param Npc Default: all. How many principle components to use for searching the nearest reference cluster centroid to each query cell.
#' @param w Optional. Assign weights to genes in restricted linear regression. Can be either a numerical vector of the same length as the number of genes or be the names of a reference metadata attribute.
#' @param colDataPrefix Default: 'Cb_'. Prefix added to the column names of Capybara scores, which will be stored in the ColData of the query sce.
#' @return Query sce object.
#' @export
SincastCapybara <- function(reference, query, clusterid, k = 5, Npc = NULL, w = NULL, colDataPrefix = 'Cb_'){
  cluster <- colData(reference)[,clusterid]
  cluname <- sort(unique(cluster))
  Ncluster <- length(cluname)
  G <- nrow(query)
  Nquery <- ncol(query)
  referencePC <- reducedDim(reference)
  queryPC <- reducedDim(query, 'referencePCA')
  if(is.null(Npc)) Npc <- ncol(reducedDim(reference))


  if(is.null(w)){
    w <- rep(1,G)

  }else if(length(w)==1){
    w <- rowData(reference)[,w]

  }else{
    if(length(w) != G){
      stop('Length of weights vector must be the same as the number of features')
    }

  }

  xlist <- list()
  pclist <- list()

  x <- matrix(ncol = Ncluster, nrow = G)
  neighbouPerCluster <- matrix(ncol = Ncluster, nrow = Nquery)
  cellScore <- matrix(ncol = Ncluster, nrow = Nquery)

  colnames(x) <- colnames(neighbouPerCluster) <- colnames(cellScore) <- cluname
  rownames(neighbouPerCluster) <- rownames(cellScore) <-  colnames(query)


  message("Now perform nearest sub-cluster centroid searching.")
  for(i in cluname){
    message(paste("\r \t Searching for Cluster", i), appendLF = F)
    refDataCluI <- matrix(t(assay(reference,'rank'))[cluster == i,], ncol = G)
    refPcCluI <- matrix(reducedDim(reference)[cluster == i,], ncol = Npc)
    NrefDataCluI<- nrow(refDataCluI)

    if(k > NrefDataCluI-1){
      pclist[[i]] <- refPcCluI
      xlist[[i]] <- refDataCluI

    }else{
      km <- cluster::pam(refPcCluI,k)
      pclist[[i]] <- km$medoids
      xlist[[i]] <- matrix(refDataCluI[km$id.med,], ncol = G)

    }

    if(k == 1){
      neighbouPerCluster[,i] <- 1

    }else{
      neighbouPerCluster[,i] <- RANN::nn2(pclist[[i]], reducedDim(query, 'referencePCA'), 1)$nn.idx

    }
  }
  message("Done")

  message("Now perform Restricted linear regression.")
  varExp <- c()
  pb = txtProgressBar(min = 0, max = Nquery, initial = 0, style = 3, width = 40)
  for(i in 1:Nquery){
    y <- assay(query, 'rank')[,i]
    for(j in cluname) x[,j] <- xlist[[j]][neighbouPerCluster[i,j],]

    Dmat <- crossprod(x, x * w)
    A <- rbind(-1,diag(Ncluster))
    b <- c(-1,rep(0,Ncluster))
    d <- crossprod(x, y * w)

    cellScore[i,] <- quadprog::solve.QP(Dmat = Dmat, dvec = d, Amat = t(A), bvec = b,meq = 0)$solution
    varExp[i] <- 1-sum(w*(y-x%*%cellScore[i,])^2)/sum((w*y)^2)

    setTxtProgressBar(pb,i)
  }
  close(pb)
  message("Done")

  cellScore[cellScore < 0] <- 0
  metadata(query)[['Capybara']] <- cellScore
  colnames(cellScore) <- paste(colDataPrefix,colnames(cellScore), sep = '')
  colData(query) <- cbind(query@colData,cellScore)

  query

}

