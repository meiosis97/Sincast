#' Relative count normalization
#'
#' This function normalizes count matrix in the sce object, generates a separate
#' relative count matrix stored in the data slot of the sce.
#'
#' @param sce Required. A sce object containing count matrix.
#' @param scale.factor Default: 1e6. Will produce Count Per Million transformed data if set default.
#' @return SingleCellExperiment Object.
#' @export
rcTransform <- function(sce, scale.factor = 1e6){
  f <- scale.factor/colSums(assay(sce, 'counts'))
  assay(sce, 'data') <- sweep(assay(sce, 'counts'), 2, f, '*')

  sce
}
