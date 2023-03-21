findk <- function(dist, k){

  if(k=='method1'){
    message('now searching for k via binary search')
    distRankMat <- apply(dist, 1, rank)

    upper.k <- N
    lower.k <- 0
    cur.k <- round(upper.k * 0.05)

    while(upper.k - lower.k > 1){
      adjTmp <- Matrix::Matrix(as.numeric(distRankMat <= cur.k), nrow = N, ncol = N, sparse = T)
      adjTmp <- adjTmp + t(adjTmp)
      isConnected <- igraph::is.connected(
        igraph::graph_from_adjacency_matrix(adjTmp, mode = 'undirected'))

      if(isConnected){
        upper.k <- cur.k
        cur.k <- round((lower.k + upper.k)/2)

      }else{
        lower.k <- cur.k
        cur.k <- round((lower.k + upper.k)/2)

      }
    }
    k <- upper.k

    message(paste('finish searching for k via binary search, k is', k))

  }else if(k == 'method2'){
    message('now searching for k via grid search')
    distRankMat <- apply(dist, 1, rank)

    cur.k <- 0.01*N

    for(i in 1:100){
      adjTmp <- Matrix::Matrix(as.numeric(distRankMat<= i*0.01*N), nrow = N, ncol = N, sparse = T)
      adjTmp <- adjTmp + t(adjTmp)
      if(mean((x %e% adjTmp) == 0) < 0.25) break

    }
    k <- round(i*0.01*N)

    message(paste('finish searching for k via binary search, k is', k))

  }

  k
}

symmetrization <- function(aff, norm){

  if(norm == 'probabilistic'){
    aff <- aff + t(aff) - aff * t(aff)

  }else if(norm == 'fuji'){
    tnorm <- aff %*% t(aff)
    x <- rowSums(aff)
    tconorm <- sapply(x, function(z) z+x) - tnorm
    aff <-0.5*(aff+t(aff))*(tnorm/tconorm)

  }else{
    stop('Wrong matrix symmetrization method specified')

  }

}

laplacian <- function(aff){

  Z <- rowSums(aff)
  Z.mat <- Z %*% t(Z)
  aff <- aff/Z.mat

}

medianScale <- function(X,Y){
  scale.factor <- apply(X, 1, function(x) median(x[x!=0]))/
    apply(replace(Y, X == 0, NA), 1, function(y) median(y, na.rm = T))
  Y <- Y * scale.factor
  Y[is.na(Y)] <- 0
  Y
}

find_sigma <- function(dk, a, k){
  lower <- 0
  upper <- Inf
  cur <- dk[k]
  while(T){
    psum <- sum(exp(-0.5*(dk/cur)^2))
    if(psum > a){
      upper <- cur
      cur <- (lower + cur)/2
    }else if(psum < a){
      if(is.infinite(upper)){
        lower <- cur
        cur <- 2 * cur
      }else{
        lower <- cur
        cur <- (upper + cur)/2
      }
    }
    if(abs(psum-a) < 1e-5) break
  }
  cur
}


#' Sincast imputation
#'
#' Perform Sincast imputation.
#'
#' @param query Required. Query sce.
#' @param dologScale Default: TRUE. Whether to log scale the query data.
#' @param doRpca Default: FALSE Whether to perform randomized PCA.
#' @param assay Default: 'data'. On which assay we perform imputation.
#' @param npc Default: 50. How many Principal components to compute on the query data.
#' @param scale Default: FALSE. Whether to Scale the query data for PCA.
#' @param k Default: 30. k neighborhood to infer adaptive Gaussian kernel.
#' @param umapdist Default: TRUE. Whether to scale Euclidean distances to distances beyond nearest neighbors as in the UMAP algorithm.
#' @param a Default: log2(k) and log(k/(k-1)) when setting umapdist to TRUE and FALSE respectively.
#'        a < 1 represents the probability of a cell communicating with its kth nearest neighbor.
#'        a > 1 represents the sum of probabilities of a cell communicating with its k nearest neighbors.
#' @param knn Default: k. A cell can only be connected to knn neighbors when t is set 1.
#' @param doLaplacian Default: TRUE. Whether to perform Laplacian normalization on the affinity matrix.
#' @param norm Default: 'probabilistic'. How to symmetrize the affinity matrix. Default is Probabilistic t-norm. The other option is 'fuji' (Fuzzy Jaccard index).
#' @param t Default: 3. Diffusion time. The power of Markov transition matrix.
#' @param col.by Optional. Color query cells by which metadata attribute.
#' @param colors Optional. Color Scheme of col.by. Should be a named vector using color codes as values and labels of 'col.by' as names.
#' @param text Optional. Additional annotation of query cells.
#' @param vis.dc Default: TRUE. Whether to visualize the diffusion component extracted from the Markov transition matrix.
#' @param saveGrpah Default: FALSE. Whether to save the igraph model
#' @return Query sce object.
#' @export
FastSincastImp <- function(query, dologScale = T, doRpca = T, assay = 'data', npc = 50, scale = F, k = 30, umapdist = TRUE,
                       a = NULL, knn = NULL, doLaplacian = TRUE, norm = 'probabilistic', t = 3, col.by = NULL,
                       colors = NULL, text = NULL, vis.dc = FALSE, saveGrpah = FALSE){

  if(dologScale) x <- Matrix::Matrix(log(assay(query, assay)+1),sparse = T) else x <- assay(query, assay)
  N <- ncol(x); G <- nrow(x)
  out <- list()

  message('Now perform PCA')
  r <- Sincast::pca(t(x), scale = scale, rank = npc, doRpca = doRpca)$x
  message('Finish PCA')

  message('Now construct affinity matrix')
  message('\t Calcualting distance')

  k <- findk(out$dist, k) #kappa
  if(is.null(a)){
    if(umapdist) a <- log2(k) else a <- log(k/(k-1)) # probability
  }
  if(is.null(knn)) knn <- k #maximum neighborhood size

  out$nn_res <- RANN::nn2(r, k = max(knn,k) + 1)

  message('\t Scaling distance')
  if(umapdist){
    out$scaled_dist <- out$nn_res$nn.dists - out$nn_res$nn.dists[,2]
    out$scaled_dist[,1] <- 0
  }else{
    out$scaled_dist <- NULL
  }

  message('\t Calculating band width')
  out$sigma <- c()
  for(i in 1:N){
    tmp <- if(umapdist) out$scaled_dist[i,2:(k+1)] else out$nn_res$nn.dists[i,2:(k+1)]
    if(a < 1){
      out$sigma[i] <- tmp[k]/sqrt(-2*log(a))
    }else if(a > 1){
      out$sigma[i] <- find_sigma(tmp, a, k)
    }else{
      out$sigma[i] <- Inf
    }
  }
  names(out$sigma) <- colnames(query)

  if(umapdist){
    tmp <- exp(-0.5 * (out$scaled_dist/out$sigma)^2)
  }else{
    tmp <- exp(-0.5 * (out$nn_res$nn.dists/out$sigma)^2)
  }

  out$aff <- matrix(0, N, N)
  for(i in 1:N) out$aff[i, out$nn_res$nn.idx[i, 2:(knn+1)]] <- tmp[i, 2:(knn+1)]
  out$aff <- Matrix::Matrix(out$aff)
  out$aff <- symmetrization(out$aff, norm = norm)

  if(doLaplacian){message('\t Laplacian normalization'); out$aff <- laplacian(out$aff)}
  message('Finish construct affinity matrix')

  message('Now impute')
  out$p <- out$aff/rowSums(out$aff) #diffusion operator
  assay(query, 'SincastImpData') <- as.matrix(x)

  message('Diffusing')
  for(i in 1:t) assay(query, 'SincastImpData', withDimnames = F) <- assay(query, 'SincastImpData', withDimnames = F) %*% out$p

  message('Scaling')
  assay(query, 'SincastImpData') <- medianScale(as.matrix(x), assay(query, 'SincastImpData'))
  message('Finish impute')

  sparsity <- mean(assay(query, 'SincastImpData')==0)
  message(paste('Sparsity after imputation is', round(sparsity,3)))

  #make plots
  if(vis.dc){

    #plotly 3d map
    s <- RSpectra::eigs(out$p,k =20, method = 'LR')
    s$values <- Re(s$values); s$vectors <- Re(s$vectors)
    s$dc <-  sweep(s$vectors, 2, s$values^t, '*')
    rownames(s$dc) <- colnames(query); colnames(s$dc) <- paste('DC',1:20,sep='')

    if(!is.null(col.by)) color <- colData(query)[,col.by] else  color <- NULL
    fig <- plotly::plot_ly() %>% plotly::add_trace(data = data.frame(s$dc),
                                                   x = ~DC2, y = ~DC3, z = ~DC4,
                                                   color = color,
                                                   colors = colors,
                                                   text = text,
                                                   marker = list(size = 5)) %>% plotly::layout(legend = list(orientation = 'h'))

    suppressMessages(suppressWarnings(print(fig)))


    #igraph network
    g <- igraph::graph_from_adjacency_matrix(out$aff, mod = 'undirected', weighted = T,diag = FALSE)
    c_scale <- colorRamp(RColorBrewer::brewer.pal(9,'Blues'))
    igraph::E(g)$color = apply(c_scale(igraph::E(g)$weight/max(igraph::E(g)$weight)),
                               1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255))

    if(!is.null(col.by)){
      clustMed <- plyr::ddply(data.frame(group = color,id = colnames(query),s$dc[,2:3]),'group',
                              function(x) x[cluster::pam(x[,3:4],k = 1)$id.med,])$id

      plot(layout = s$dc[,2:3],
           g, vertex.size= 0.01,
           vertex.label.color = if(is.null(colors)) as.numeric(as.factor(color)) else colors[color],
           vertex.frame.color = adjustcolor("white", alpha.f = 0),
           vertex.color = if(is.null(colors)) as.numeric(as.factor(color)) else colors[color],
           vertex.label = ifelse(colnames(query)%in%clustMed, color,''),
           edge.width= 1)

    }else{
      plot(layout = s$dc[,2:3],
           g, vertex.size= 0.01,vertex.label.cex = 0.01,
           vertex.frame.color = adjustcolor("white", alpha.f = 0),
           vertex.color = 'steelblue',
           vertex.label = rep('',N),
           edge.width= 1)

    }

    out$SincastDMapMod <- s
    if(saveGrpah) out$graphMod <- g
    out$fig <- fig
    reducedDim(query, 'SincastDMap') <- s$dc
  }

  metadata(query)[['SincastImp']] <- out

  query

}
