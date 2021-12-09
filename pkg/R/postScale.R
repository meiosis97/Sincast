#' Post imputation data scaling.
#'
#' Perform Sincast post imputation data scaling.
#'
#' @param query Required. Query sce.
#' @param preImpAssay Default: 'data'. Which assay in the query sce corresponds to the original data matrix.
#' @param postImpAssay Default: 'SincastImpData'. Which assay in the query sce corresponds to the imputed data matrix.
#' @param plotTrend Default: TRUE. Whether to plot the genewise mean versus variance trend estimation.
#' @param dologScale Default: TRUE. Whether to perform log scaling on the original data matrix (Which should be done if the imputation also performed logscaling)..
#' @param k Optional. Basis dimension for GAM.
#' @param returnObsVarEst Default: FALSE. Whether to return observational variance estimation.
#' @return Query sce object.
#' @export
postScale <- function(query, preImpAssay = 'data', postImpAssay = 'SincastImpData', plotTrend = TRUE, dologScale = TRUE, k = NULL,
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
    p <- ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(log.mu,log.var, col = w))+
      ggplot2::scale_color_gradientn('Zero proportion',colours = rainbow(100)[1:70]) +
      ggplot2::geom_path(ggplot2::aes(log.mu, predict(GAMmod,data.frame(log.mu))),size = 1.5, linetype = 'dashed') +
      ggplot2::theme_bw() + ggplot2::xlab('Log-Mean') + ggplot2::ylab('Log-Variance') + ggplot2::theme(text = ggplot2::element_text(size=15))
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
