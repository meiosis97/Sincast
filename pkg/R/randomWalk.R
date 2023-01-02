medianScale <- function(X,Y){
  scale.factor <- apply(X, 1, function(x) median(x[x!=0]))/
    apply(replace(Y, X == 0, NA), 1, function(y) median(y, na.rm = T))
  Y <- Y * scale.factor
  Y[is.na(Y)] <- 0
  Y
}


#' Diffusion random walk
#'
#' Applying the diffusion operator computed by sincastImp to impute data. Can be used
#' for imputing full data with diffusion opperators computed on filtered features.
#'
#' @param sce Required. A sce object, on which sincastImp has been applied. The
#'  diffusion operator compted by sincastImp will be applied to facilitate diffusion random
#'  walk.
#' @param data Required. A gene by cell expression matrix to be diffused. Expect matched columns
#' between the sce object and the data.
#' @param logScale Default: TRUE. Whether to logscale the data.
#' @param t Default: 3. Diffusion time. The power of Markov transition matrix.
#' @return A imputed data matrix
#' @export
randomWalk <- function(sce, data, logScale = T, t = 3){
  if(logScale) data <- log(data+1)
  impdata <- as.matrix(data)

  message('Diffusing')
  for(i in 1:t) impdata <- impdata %e% t(sce@metadata$SincastImp$p)

  message('Scaling')
  impdata <- medianScale(as.matrix(data), impdata)
  message('Finish impute')

  sparsity <- mean(impdata==0)
  message(paste('Sparsity after imputation is', round(sparsity,3)))

  dimnames(impdata) <- dimnames(data)
  impdata
}
