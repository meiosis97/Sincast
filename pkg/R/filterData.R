#' Filter data
#'
#' Quality control the query data and filter the reference and query to their
#' common gene set.
#'
#' @param reference Required. Reference sce.
#' @param query Required. Query sce.
#' @param min.cell Default: 10. Query genes expressed in less then 'min.cell' cells are filtered.
#' @param max.sparsity Default: 0.98. Cells with proportion of zeros larger than 'max.sparsity' are filtered.
#' @param referenceAssay Default: 'rank'. On which reference assay we perform rank transformation.
#' @param queryAssay Default: 'data'. On which query assay we perform quality controls.
#' @return reference and query sce objects.
#' @export
filterData <- function(reference, query, min.cell = 10, max.sparsity = 0.98, referenceAssay = 'rank', queryAssay = 'data'){

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
  reference <- rankTrans(reference, referenceAssay)

  return(list(reference = reference, query = query))

}
