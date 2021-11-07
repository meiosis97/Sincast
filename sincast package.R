library(SingleCellExperiment)
library(zeallot)
library(plotly)
library(Rcpp)
sourceCpp("H:/sincast paper/mat.cpp")


load('H:/sincast query dataset/Villani etal.Rdata')
load('H:/sincast reference dataset/Rdata/Myeloid atlas-DC subset.Rdata')

query <-  SingleCellExperiment(list(data = Matrix::Matrix(as.matrix(query.data),sparse = T)), colData = query.annotation)
reference <- SingleCellExperiment(list(data = Matrix::Matrix(as.matrix(atlas.data),sparse = T)), colData = atlas.annotation)

################# cpp multiply ##############
`%e*%` <- function(x, y) eigenMapMatMult(x, y)

assayNames(reference)
################# color hue #################
rank.trans <- function(sce, assay = 'data'){
  if(!assay%in%assayNames(sce)){
    warnings('Assay not found, use the first assay to rank transform.')
    assay(sce, 'rank') <- (apply(assay(sce),2,rank, ties.method = 'min')-1)/(nrow(sce) - 1)
  }else{
    assay(sce, 'rank') <- (apply(assay(sce, assay),2,rank, ties.method = 'min')-1)/(nrow(sce) - 1)
    
  }
  sce  

}


################# color hue #################
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


################# HD #################
HellingerDist <- function(x, y){
  
  t1 <- table(x, y)
  tot <- colSums(t1)
  t2 <- -sweep(t1,2,tot)
  t1 <- t1/rowSums(t1)
  t2 <- t2/rowSums(t2)
  sqrt(rowSums((sqrt(t1)-sqrt(t2))^2))
  
}


################# feature weighting #################
feature.weighting <- function(reference, clusterid, cut = T, n.cut = 2000, assay2rank = 'data', assay2weight = 'rank'){
  
  rownames(reference) <-gsub("\\.",'-',rownames(reference))
  cluster <- colData(reference)[,clusterid]
  cluname <- sort(unique(cluster))
  n.cluster <- length(cluname)
  
  #create rank assay
  reference <- rank.trans(reference, assay2rank)
  
  #remove genes with less than n.clusters expressed
  reference <- reference[rowSums(assay(reference, assay2weight)!=0)>n.cluster,]
  
  n.sample <- ncol(reference)
  n.gene <- nrow(reference)
  
  #rerank
  reference <- rank.trans(reference, assay2rank)
  
  
  #discreteize reference data
  discRankDat <- apply(assay(reference, assay2weight), 1, 
                      function(x) as.factor(kmeans(x, centers = n.cluster, iter.max = 30)$cluster))
  #calculate HD
  HD <- data.frame(t(apply(discRankDat,2,function(y) HellingerDist(cluster,y))))
  #generate metadata
  HD <- dplyr::mutate(HD, gene = rownames(HD),
                      mean = rowMeans(HD), celltype = cluname[apply(HD, 1, which.max)])
  #order the HD
  HD <- arrange(HD, desc(mean))
  
  
  #filter the reference
  if(cut){
    if(n.cut > n.gene){
      warning(paste('Exceeding number of features were selected. Select a number smaller than', n.gene))
      n.cut <- n.gene
    }
    reference <- reference[HD$gene[1:n.cut],]
    rowData(reference) <- HD[1:n.cut,]
    
  }else{
    reference <- reference[HD$gene,]
    rowData(reference) <- HD
    
  }
  
  reference <- rank.trans(reference, 'rank')
  metadata(reference)[['HD']] <- HD
  reference
  
}


###########################filter data###################################
filter.data <- function(reference, query, min.cell = 10, max.sparsity = 0.98, referenceAssay = 'rank', queryAssay = 'data'){
  
  rownames(query) <-gsub("\\.",'-',rownames(query))
  
  #filter cells with sparsity rate larger than max.sparsity pct
  keep <- colSums(assay(query,queryAssay) == 0, na.rm =T)/nrow(query) < max.sparsity
  query <- query[,keep]
  
  #filter genes that expressed in at least min.cell
  query <- query[rowSums(assay(query,queryAssay) !=0, na.rm =T)> min.cell,]
  
  intersect <- intersect(rownames(reference),rownames(query))
  reference <- reference[intersect,]
  query <- query[intersect,]
  
  #re-rank the reference
  reference <- rank.trans(reference, referenceAssay) 
  
  return(list(reference = reference, query = query))
  
}


###########################pca##################################
pca <- function(x, retx = TRUE, center = TRUE, scale. = FALSE, rank = 100){
    x <- as.matrix(x)
    x <- scale(x, center = center, scale = scale.)
    cen <- attr(x, "scaled:center")
    sc <- attr(x, "scaled:scale")
    if(any(sc == 0))
      stop("cannot rescale a constant/zero column to unit variance")
    s <- RSpectra::svds(x, k = 100, nu = 0)
    s$d <- s$d / sqrt(max(1, nrow(x) - 1))
    dimnames(s$v) <-
      list(colnames(x), paste0("PC", seq_len(ncol(s$v))))
    r <- list(sdev = s$d, rotation = s$v,
              center = if(is.null(cen)) FALSE else cen,
              scale = if(is.null(sc)) FALSE else sc)
    if (retx) r$x <- x %e*% s$v
    rownames(r$x) <- rownames(x); colnames(r$x) <- colnames(s$v)
    totVar <- sum(x^2)/max(1, nrow(x)-1)
    r$Explainedvar <- r$sdev^2/totVar
    r$cumExplainedvar <- cumsum(r$sdev^2/totVar)
    
    r
}

###########################make atlas##################################
make.atlas <- function(reference, npc = 100, scale = FALSE,
                       vis.atlas = T, col.by = NULL, colors = NULL, text = NULL, rank = TRUE, assay = 'data'){
  
  #rank transformation and PCA
  if(rank) reference <- rank.trans(reference, assay)
  r <- pca(t(assay(reference, 'rank')),
           scale = scale, rank = npc)
  
  #make plots
  if(!is.null(col.by)) color <- colData(reference)[,col.by] else  color <- NULL

  fig <- plot_ly() %>% add_trace(data = data.frame(r$x),
                                 x = ~PC1, y = ~PC2, z = ~PC3,
                                 color = color,
                                 colors = colors,
                                 text = text,
                                 marker = list(size = 5)) %>% layout(legend = list(orientation = 'h'))
  
  if(vis.atlas) suppressMessages(suppressWarnings(print(fig)))
  
  r$fig <- fig
  metadata(reference)[['PCA']] <- r
  reducedDim(reference, 'PCA') <- r$x
  
  reference
  
}




################# pdist #################
pdist <- function(tmat){
  # @param tmat A non-negative matrix with samples by features
  # @reference http://r.789695.n4.nabble.com/dist-function-in-R-is-very-slow-td4738317.html
  mtm <- Matrix::tcrossprod(tmat)
  sq <- rowSums(tmat^2)
  out0 <- outer(sq, sq, "+") - 2 * mtm
  out0[out0 < 0] <- 0
  
  sqrt(out0)
}


################# sincast #################
findk <- function(dist, k){
  distRankMat <- apply(dist, 1, rank)
  
  if(k=='method1'){
    message('now searching for k via binary search')
    
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
    
    cur.k <- 0.01*N
    
    for(i in 1:100){
      adjTmp <- Matrix::Matrix(as.numeric(distRankMat<= i*0.01*N), nrow = N, ncol = N, sparse = T)
      adjTmp <- adjTmp + t(adjTmp)
      if(mean((x %e*% adjTmp) == 0) < 0.25) break
      
    }
    k <- round(i*0.01*N)
    
    message(paste('finish searching for k via binary search, k is', k))
  
  }
    
  k
}


scaleDistance <- function(dist){
  
  for(i in 1:nrow(dist)){
    r <-  sort(dist[i,],partial=2)[2]
    dist[i,] <- dist[i,]-r
    dist[i,][dist[i,] <= 0] <- 0
  }
  dist
  
}


symmetrization <- function(aff, norm){
  
  if(norm == 'probablistic'){
    aff <- aff + t(aff) - aff * t(aff)
    
  }else if(norm == 'fuji'){
    tnorm <- aff %e*% t(aff)
    x <- rowSums(aff)
    tconorm <- sapply(x, function(z) z+x) - tnorm
    aff <-0.5*(aff+t(aff))*(tnorm/tconorm)
    
  }else{
    stop('Wrong matrix symmetrization method specified')
    
  }
  
}


laplacian <- function(aff){

  Z <- rowSums(aff)
  Z.mat <- Z %e*% t(Z)
  aff <- aff/Z.mat

}


medianScale <- function(x,y){
  for(i in 1:nrow(x)){
    cellExpressed <- x[i,]!=0
    scale.factor <- median(x[i,cellExpressed])/median(y[i,cellExpressed])
    y[i,] <- y[i,] *scale.factor
  }
  y
  
}

sincastImp <- function(query, logScale = T, assay = 'data', npc = 50, scale = F, k = 30, umapdist = TRUE,
                       a = NULL, knn = NULL, doLaplacian = TRUE, norm = 'fuji', t = 3, col.by = NULL,
                       colors = NULL, text = NULL, vis.dc = TRUE, saveGrpah = FALSE){

  if(logScale) x <- Matrix::Matrix(log(assay(query, assay)+1),sparse = T) else x <- assay(query, assay)
  N <- ncol(x); G <- nrow(x)
  out <- list()
  
  message('Now perform PCA')
    r <- pca(t(x), scale = scale, rank = npc)
  message('Finish PCA')
  
  message('Now construct affinity matrix')
    message('\t Calcualting distance')
    out$dist <- pdist(r$x)
    k <- findk(out$dist, k) #kappa
    if(is.null(a)) a <- log(k/(k-1)) #probability
    if(is.null(knn)) knn <- k #maximum neighborhood size
    
    message('\t Scaling distance')
    if(umapdist) out$dist <- scaleDistance(out$dist)
  
    message('\t Calculating band width')
    out$sigma <- c()
    for(i in 1:N) out$sigma[i] <- sort(out$dist[i,], partial = k)[k]/ sqrt(-2*log(a))
    
    out$aff <- exp(-0.5*(out$dist/out$sigma)^2)
    for(i in 1:N) out$aff[i,][rank(-out$aff[i,]) > knn] <- 0
    
    out$aff <- symmetrization(out$aff,norm = norm)
  
    if(doLaplacian){message('\t Laplacian normalization'); out$aff <- laplacian(out$aff)}
  message('Finish construct affinity matrix')
  
  message('Now impute')
    out$p <- out$aff/rowSums(out$aff) #diffusion operator
    assay(query, 'SincastImpData') <- as.matrix(x)
    message('Diffusing')
    for(i in 1:t) assay(query, 'SincastImpData') <- assay(query, 'SincastImpData') %e*% out$p
    
    message('Scaling')
    assay(query, 'SincastImpData') <- medianScale(x, assay(query, 'SincastImpData'))
  message('Finish impute')
  
  sparsity <- mean(assay(query, 'SincastImpData')==0)
  message(paste('Sparsity after imputation is', round(sparsity,3)))
  
  #make plots
  s <- RSpectra::eigs(out$p,k =20, method = 'LR')
  s$values <- Re(s$values); s$vectors <- Re(s$vectors)
  s$dc <-  sweep(s$vectors, 2, s$values^t, '*')
  rownames(s$dc) <- colnames(query); colnames(s$dc) <- paste('DC',1:20,sep='')
  
  if(!is.null(col.by)) color <- colData(query)[,col.by] else  color <- NULL
  fig <- plot_ly() %>% add_trace(data = data.frame(s$dc),
                                 x = ~DC2, y = ~DC3, z = ~DC4,
                                 color = color,
                                 colors = colors,
                                 text = text,
                                 marker = list(size = 5)) %>% layout(legend = list(orientation = 'h'))
  
  if(vis.dc) suppressMessages(suppressWarnings(print(fig)))
  
  #make network
  g <- igraph::graph_from_adjacency_matrix(out$aff, mod = 'undirected', weighted = T,diag = FALSE)
  c_scale <- colorRamp(RColorBrewer::brewer.pal(9,'Blues'))
  igraph::E(g)$color = apply(c_scale(igraph::E(g)$weight/max(igraph::E(g)$weight)),
                             1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255))
  if(vis.dc){
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
  }
  
  out$SincastDMapMod <- s
  if(saveGrpah) out$graphMod <- g
  out$fig <- fig
  
  reducedDim(query, 'SincastDMap') <- s$dc
  metadata(query)[['SincastImp']] <- out

  query

}


################# post scaling #################
postScale <- function(query, preImpAssay = 'data', postImpAssay = 'SincastImpData', plotTrend = T, dologScale = T, k = NULL,
                          returnObsVarEst = FALSE){
  
  if(dologScale) y <- Matrix::Matrix(log(assay(query, preImpAssay)+1),sparse = T) else y <- assay(query, preImpAssay)
  N <- ncol(query)
  G <- nrow(query)
  r <- ppoints(N)
  q <- qnorm(r)
  w <- rowMeans(assay(query, preImpAssay)!=0)
  
  #Scaling by Median
  x <- medianScale(y, assay(query, postImpAssay))
  
  #gene wise mean and variance estimate on expressed genes
  message('Genewise Mean and Variance estimation on imputed data')
  mu <-c()
  var <- c()
  for(i in 1:G){
    z <- sort(x[i,])
    lmod <- lm(z[z>0]~q[z>0])
    mu[i] <- lmod$coefficients[1]
    var[i] <- lmod$coefficients[2]^2
  }
  names(mu) <-  names(var) <- names(w) <- rownames(query)
  message('Done')
  
  #estimate mean and variance trend
  suppressWarnings(log.mu <- sort(log(mu)))
  log.var <- log(var)[names(log.mu)]
  w <- w[names(log.mu)]
  
  message("Now perform GAM fit.")
  if(is.null(k)){
    k <- 3
    GAMmod <- mgcv::gam(log.var~s(log.mu, k = k, bs = 'cr'))
    while(mgcv::k.check(GAMmod)[,4] < 0.05){
      k <- k + 2
      GAMmod <- mgcv::gam(log.var~s(log.mu, k = k, bs = 'cr'))
      
    } 
    
  }else{
    GAMmod <- mgcv::gam(log.var~s(log.mu, k = k, bs = 'cr'))
  }
  
  if(plotTrend){
    p <- ggplot() + geom_point(aes(log.mu,log.var, col = w))+
      scale_color_gradientn('Zero proportion',colours = rainbow(100)[1:70]) + geom_path(aes(log.mu, predict(GAMmod,data.frame(log.mu))),size = 1.5, linetype = 'dashed') + 
      theme_bw() + xlab('Log-Mean') + ylab('Log-Variance') +   theme(text = element_text(size=15))
    print(p)
    
  }
  message(paste("Finish regress. k is", k))
  
  message(paste("Perform observation wise variance estimation"))
  x[x == 0] <- NA
  sigma <- log(x)
  for(i in 1:N){
    sigma[,i] <- exp(predict(GAMmod, data.frame(log.mu = sigma[,i])))
  }
  
  x[is.na(x)] <- 0
  e <- (y-x)^2
  message(paste("Done"))
  
  lambda<- colMeans(e/(e+sigma),na.rm = T)
  assay(query, 'SincastScaledData') <-sweep(y,2,lambda,"*") +  sweep(x,2,1-lambda,"*") 
  metadata(query)[['SincastScale']] <- list(GAMmod = GAMmod, lambda = lambda, mu = mu, var = var)
  if(returnObsVarEst){metadata(query)[['SincastScale']]$e <- e; metadata(query)[['SincastScale']]$sigma <- sigma}
  
  query
  
}


################# project #################
projectDefault <- function(x, load, center, scale = F){
  pcs <- scale(t(x), center = center, scale= scale) %e*% load
  colnames(pcs) <- paste('PC', 1:ncol(load), sep = '')
  pcs
  
}


project <- function(reference, query, assay = 'SincastScaledData'){
  query <- rank.trans(query, assay)
  reducedDim(query, paste('referencePCA')) <- projectDefault(assay(query, 'rank'), load = reference@metadata$PCA$rotation,
                 center = reference@metadata$PCA$center, scale = reference@metadata$PCA$scale)
  
  query
}

















colors <- atlas.color$color
names(colors) <- atlas.color$celltype

reference <- feature.weighting(reference, 'celltype')
c(reference, query) %<-% filter.data(reference, query)
reference <- make.atlas(reference = reference, col.by = 'celltype', colors = colors, vis.atlas = T)
query <- sincastImp(query, col.by = 'celltype', t = 3)
query <- postScale(query)
query <- project(reference, query)


plot_ly() %>% add_trace(data = data.frame(reducedDim(query,'referencePCA')),
                               x = ~PC1, y = ~PC2, z = ~PC3,
                               color = query$celltype,
                               marker = list(size = 5)) %>% layout(legend = list(orientation = 'h'))
