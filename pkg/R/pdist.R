#' Distance
#'
#' Fast calculation of Euclidean distance matrix.
#'
#' @param tmat Required. Data matrix.
#' @return distance matrix.
#' @export
pdist <- function(tmat){
  # @param tmat A non-negative matrix with samples by features
  # @reference http://r.789695.n4.nabble.com/dist-function-in-R-is-very-slow-td4738317.html
  mtm <- Matrix::tcrossprod(tmat)
  sq <- rowSums(tmat^2)
  out0 <- outer(sq, sq, "+") - 2 * mtm
  out0[out0 < 0] <- 0

  sqrt(out0)
}
