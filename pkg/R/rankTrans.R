#' Rank transformation
#'
#' This function rank transforms a specified data slot in sce, generates a separate
#' rank expression matrix stored in the rank slot of the sce.
#'
#'
#' @param sce Required. A sce object.
#' @param assay Default: 'data'. On which data assay we perform rank transformation.
#' @return SingleCellExperiment Object.
#' @export
rankTrans <- function(sce, assay = 'data'){
  if(!assay%in%assayNames(sce)){
    warnings('Assay not found, use the first assay to rank transform.')
    assay(sce, 'rank') <- (apply(assay(sce),2,rank, ties.method = 'min')-1)/(nrow(sce) - 1)
  }else{
    assay(sce, 'rank') <- (apply(assay(sce, assay),2,rank, ties.method = 'min')-1)/(nrow(sce) - 1)

  }
  sce

}

