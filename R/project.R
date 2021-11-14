#' Default PCA projection
#'
#' Perform PCA projection.
#'
#' @param x Required. A feature by observation data matrix.
#' @param load Required. Loading matrix.
#' @param center Default: FALSE Centering factors of the reference.
#' @param scale Default: FALSE. Scaling factors of the reference.
#' @return PC scores
#' @export
projectDefault <- function(x, load, center = FALSE, scale = FALSE){
  pcs <- scale(t(x), center = center, scale= scale) %e*% load
  colnames(pcs) <- paste('PC', 1:ncol(load), sep = '')
  pcs

}




#' PCA projection
#'
#' Perform PCA projection on sce object.
#'
#' @param reference Required. Reference sce with PCA performed.
#' @param query Required. Query sce.
#' @param assay Default: 'SincastScaledData'. The query assay to be projected onto the reference atlas.
#' @return PC scores
#' @export
project <- function(reference, query, assay = 'SincastScaledData'){
  query <- rankTrans(query, assay)
  reducedDim(query, paste('referencePCA')) <- projectDefault(assay(query, 'rank'), load = reference@metadata$PCA$rotation,
                                                             center = reference@metadata$PCA$center, scale = reference@metadata$PCA$scale)

  query
}
