#' Hillinger Distance
#'
#' Calculate Hellinger Distance between two categorical variables
#'
#'
#' @param x Required. a categorical variable.
#' @param y Required. a categorical variable.
#' @return Hillinger Distance.
#' @export
HellingerDist <- function(x, y){

  t1 <- table(x, y)
  tot <- colSums(t1)
  t2 <- -sweep(t1,2,tot)
  t1 <- t1/rowSums(t1)
  t2 <- t2/rowSums(t2)
  sqrt(rowSums((sqrt(t1)-sqrt(t2))^2))

}
