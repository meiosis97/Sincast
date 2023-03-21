find.elbow <- function(x, y){
  n <- length(x)
  scaled.y <- y * max(x)/max(y)
  r <- -(scaled.y[1] - scaled.y[n])/(x[1]-x[n])
  c <- - scaled.y[1] - r * x[1]
  best.index <- which.max(abs(r * x +  scaled.y + c))
  best.index
}


#' PCA
#'
#' Standard Principal component analysis
#'
#' @param x Required. data.
#' @param retx Default: TRUE.
#' @param center Default: TRUE.
#' @param scale. Default: FALSE.
#' @param rank Default: Find by the elbow method.
#' @return pca object.
#' @export
pca<- function(x, retx = TRUE, center = TRUE, scale. = FALSE, rank = NULL,
               doRpca = FALSE){
  x <- scale(x, center = center, scale = scale.)
  cen <- attr(x, "scaled:center")
  sc <- attr(x, "scaled:scale")
  x <- Matrix::Matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  min_dim <- min(n,p)
  if(is.null(rank)) rank <- min_dim
  if(n*p > 500^2 & rank < 0.8*min_dim & is.null(doRpca)) doRpca <- T
  if(any(sc == 0))
    stop("cannot rescale a constant/zero column to unit variance")
  if(rank == min_dim){
    message("Find the number of PC to use via the elbow method.")
    s <- svd(x)
    rank <- find.elbow(1:length(s$d),s$d)
    message(paste("The best PC found is", rank))
  }else{
    s <- if(doRpca) rsvd::rsvd(x, k = rank) else RSpectra::svds(x, k = rank)
  }
  dimnames(s$v) <-
    list(colnames(x), paste0("PC", seq_len(ncol(s$v))))
  r <- list(sdev = s$d[1:rank]/sqrt(max(1, nrow(x) - 1)), rotation = s$v[,1:rank],
            center = if(is.null(cen)) FALSE else cen,
            scale = if(is.null(sc)) FALSE else sc)

  if (retx) r$x <- sweep(s$u[,1:rank],2, s$d[1:rank], '*')
  dimnames(r$x) <- list(rownames(x), colnames(s$v)[1:rank])
  totVar <- sum(x^2)/max(1, nrow(x)-1)
  r$Explainedvar <- r$sdev^2/totVar
  r$cumExplainedvar <- cumsum(r$Explainedvar)

  r
}
