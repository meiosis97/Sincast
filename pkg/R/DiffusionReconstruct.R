#' Diffusion Reconstruction
#'
#' Diffusion reconstruction of the PCA projection landscape.
#'
#' @param reference Required. Reference sce.
#' @param query Required. Query sce.
#' @param Npc Default: all. How many principle components to use for diffusion map.
#' @param a Default: 2. The diffusion bandwidth will be set to 'a' times the maximum distance between reference sample paires.
#' @param returnSce Default: TRUE. Whether to return sce objects.
#' @param returnDiffuseObj Default: FALSE. Whether to return the diffusion map model.
#' @param nPool Optional. Pooling size of each pseudo-bulk sample. Will suppress pool.factor if specified.
#' @param colReference.by Optional. Color reference samples by which reference metadata attribute.
#' @param colQuery.by Optional. Color quey samples by which query metadata attribute.
#' @param referenceColors Optional. Color Scheme of colReference.by Should be a named vector using color codes as values and labels of 'colReference.by' as names.
#' @param queryColors Optional. Color Scheme of colQuery.by Should be a named vector using color codes as values and labels of 'colQuery.by' as names.
#' @param referenceText Optional. Additional annotation of reference samples
#' @param queryText Optional. Additional annotation of query cells.
#' @return A plotly object.
#' @export
DiffusionReconstruct <- function(reference, query, Npc, a = 2, returnSce = TRUE, returnDiffuseObj = FALSE,
                                 colReference.by = NULL, colQuery.by = NULL, referenceColors = NULL, queryColors = NULL,
                                 referenceText = NULL, queryText = NULL, ...){

  PCAcombined <- rbind(reducedDim(reference, 'PCA'), reducedDim(query, 'referencePCA'))
  #Compute bandwidth and perform diffusion embedding.
  D <- pdist(PCAcombined)
  m <- max(D[colnames(reference),colnames(reference)])
  DiffuseObj <- diffusionMap::diffuse(D = D, eps.val = a * m, ...)

  #Extract diffusion map
  DM <- DiffuseObj$X
  rownames(DM) <- rownames(PCAcombined)
  colnames(DM) <- paste('DC',1:ncol(DM),sep='')

  #Put DMs to sces
  reducedDim(reference, 'DM') <-  DM[colnames(reference),]
  reducedDim(query, 'DM') <-  DM[colnames(query),]

  #Ploter
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


  fig <- plotly::plot_ly() %>% plotly::add_markers(data = data.frame(reducedDim(reference, 'DM')),
                                                   x = ~DC1, y = ~DC2, z = ~DC3,
                                                   text = referenceText,
                                                   color = referenceColor,
                                                   colors = colors,
                                                   marker = list(size = 5))

  if(!is.numeric(queryColor)){
    fig <- fig %>% plotly::add_markers(data = data.frame(reducedDim(query, 'DM')),
                                       x = ~DC1, y = ~DC2, z = ~DC3,
                                       color = queryColor,
                                       colors = colors,
                                       text = queryText,
                                       marker = list(size = 3,symbol =  ~'x')) %>% plotly::layout(legend = list(orientation = 'h'))

  }else{
    fig <- fig %>% plotly::add_markers(x = ~DC1, y = ~DC2, z = ~DC3,
                                       data = data.frame(reducedDim(query, 'DM')),
                                       name = 'query',
                                       text = paste('Value: ', queryColor,'\n',queryText, sep=''),
                                       marker = list(symbol =  ~'x', size = 3, colorscale='Viridis',reversescale = F,
                                                     color = ~queryColor,
                                                     colorbar=list(
                                                       title='', x = -0.15, y= 5
                                                     )))

  }

  print(fig)

  if(returnSce & returnDiffuseObj) return(invisible(list(reference = reference, query = query, DiffuseObj = DiffuseObj, fig = fig)))

  if(returnSce) return(invisible(list(reference = reference, query = query, fig =fig))) else invisible(fig)

}


