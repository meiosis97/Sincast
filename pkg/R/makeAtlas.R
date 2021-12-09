#' Make Atlas
#'
#' Perform PCA on the ranked reference data to build atlas.
#'
#' @param reference Required. Reference sce.
#' @param npc Default: 100. Number of Principal Components to build for the atlas.
#' @param scale Default: FALSE. Whether to scale the reference features.
#' @param vis.atlas Default: TRUE. Whether to visualize the atlas.
#' @param col.by Optional. Color reference samples by which metadata attribute.
#' @param colors Optional. Color Scheme of col.by. Should be a named vector using color codes as values and labels of 'col.by' as names.
#' @param text Optional. Additional annotation of reference samples.
#' @param rank Default: TRUE. Whether to perform rank transformation on the reference data.
#' @param assay Default: 'data'. On which reference assay we perform rank transformation.
#' @return Reference sce.
#' @export
makeAtlas <- function(reference, npc = NULL, scale = FALSE,
                      vis.atlas = T, col.by = NULL, colors = NULL, text = NULL, rank = TRUE, assay = 'data'){

  #rank transformation and PCA
  if(rank) reference <- rankTrans(reference, assay)
  r <- pca(t(assay(reference, 'rank')),
           scale = scale, rank = npc)

  #make plots
  if(!is.null(col.by)) color <- colData(reference)[,col.by] else  color <- NULL

  fig <- plotly::plot_ly() %>% plotly::add_markers(data = data.frame(r$x),
                                                 x = ~PC1, y = ~PC2, z = ~PC3,
                                                 color = color,
                                                 colors = colors,
                                                 text = text,
                                                 marker = list(size = 5)) %>% plotly::layout(legend = list(orientation = 'h'))

  if(vis.atlas) suppressMessages(suppressWarnings(print(fig)))

  r$fig <- fig
  metadata(reference)[['PCA']] <- r
  reducedDim(reference, 'PCA') <- r$x

  reference

}



