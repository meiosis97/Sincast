#' Visualize projection
#'
#' Using plotly 3D to visualize query projection on the reference atlas.
#'
#' @param reference Required. Reference sce.
#' @param query Required. Query sce.
#' @param colReference.by Optional. Color reference samples by which reference metadata attribute.
#' @param colQuery.by Optional. Color quey samples by which query metadata attribute.
#' @param referenceColors Optional. Color Scheme of colReference.by Should be a named vector using color codes as values and labels of 'colReference.by' as names.
#' @param queryColors Optional. Color Scheme of colQuery.by Should be a named vector using color codes as values and labels of 'colQuery.by' as names.
#' @param referenceText Optional. Additional annotation of reference samples
#' @param queryText Optional. Additional annotation of query cells.
#' @return A plotly object.
#' @export
visProjection <- function(reference, query, colReference.by = NULL, colQuery.by = NULL, referenceColors = NULL, queryColors = NULL,
                          referenceText = NULL, queryText = NULL){

  colors <- c(referenceColors, queryColors)


  if(!is.null(colReference.by)){
    referenceColor <- colData(reference)[,colReference.by]
  }else{
    referenceColor <- 'reference';colors['reference'] <- 'steelblue'
  }
  if(is.null(referenceColors)&!is.null(colReference.by)){
    colors[sort(unique(referenceColor))] <- suppressWarnings(colorRampPalette(RColorBrewer::brewer.pal(12,'Paired'))(length(unique(referenceColor))))
  }


  if(length(colQuery.by) == 1 | is.null(colQuery.by)){
    if(!is.null(colQuery.by)){
      queryColor <- colData(query)[,colQuery.by]
    }else{
      queryColor <- 'query'; colors['query'] <- 'lightgrey'
    }

    if(is.null(queryColors)&!is.null(colQuery.by)&!is.numeric(queryColor)){
      colors[sort(unique(queryColor))] <- suppressWarnings(colorRampPalette(RColorBrewer::brewer.pal(12,'Set2'))(length(unique(queryColor))))
    }

  }else if(is.numeric(colQuery.by)){
    queryColor <- colQuery.by
  }else{
    stop('Wrong colQuery.by specification')

  }



  fig <- plotly::plot_ly() %>% plotly::add_markers(data = data.frame(reducedDim(reference, 'PCA')),
                                                   x = ~PC1, y = ~PC2, z = ~PC3,
                                                   text = referenceText,
                                                   color = referenceColor,
                                                   colors = colors,
                                                   marker = list(size = 5))

  if(!is.numeric(queryColor)){
    fig <- fig %>% plotly::add_markers(data = data.frame(reducedDim(query, 'referencePCA')),
                                       x = ~PC1, y = ~PC2, z = ~PC3,
                                       color = queryColor,
                                       colors = colors,
                                       text = queryText,
                                       marker = list(size = 3,symbol =  ~'x')) %>% plotly::layout(legend = list(orientation = 'h'))

  }else{
    fig <- fig %>% plotly::add_markers(x = ~PC1, y = ~PC2, z = ~PC3,
                                       data = data.frame(reducedDim(query, 'referencePCA')),
                                       name = 'query',
                                       text = paste('Value: ', queryColor,'\n',queryText, sep=''),
                                       marker = list(symbol =  ~'x', size = 3, colorscale='Viridis',reversescale = F,
                                                     color = ~queryColor,
                                                     colorbar=list(
                                                       title='', x = -0.15, y= 5
                                                     )))

  }


  suppressWarnings(suppressMessages(fig))

}


