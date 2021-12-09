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
pca <- function(x, retx = TRUE, center = TRUE, scale. = FALSE, rank = NULL){
  x <- as.matrix(x)
  x <- scale(x, center = center, scale = scale.)
  cen <- attr(x, "scaled:center")
  sc <- attr(x, "scaled:scale")
  if(any(sc == 0))
    stop("cannot rescale a constant/zero column to unit variance")
  if(is.null(rank)){
    message("Find the number of PC to use via the elbow method.")
    s <- svd(x, nu = 0)
    rank <- find.elbow(1:length(s$d),s$d)
    message(paste("The best PC found is", rank))
  }else{
    s <- RSpectra::svds(x, k = rank, nu = 0)
  }
  s$d <- s$d / sqrt(max(1, nrow(x) - 1))
  dimnames(s$v) <-
    list(colnames(x), paste0("PC", seq_len(ncol(s$v))))
  r <- list(sdev = s$d, rotation = s$v[,1:rank],
            center = if(is.null(cen)) FALSE else cen,
            scale = if(is.null(sc)) FALSE else sc)
  if (retx) r$x <- x %e% s$v[,1:rank]
  dimnames(r$x) <- list(rownames(x), colnames(s$v)[1:rank])
  totVar <- sum(x^2)/max(1, nrow(x)-1)
  r$Explainedvar <- r$sdev^2/totVar
  r$cumExplainedvar <- cumsum(r$sdev^2/totVar)

  r
}
