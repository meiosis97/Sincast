#' Plot Capybara heatmap
#'
#' A function to plot Capybara scores predicted on the query cells.
#'
#' @param query Required. Query sce.
#' @param trueAnno Optional. The known annotation of the query cells.
#' @param cluCols Default: FALSE. Whether to perform hierachical clustering	by columns (query cells).
#' @param Npc Default: FALSE. Whether to perform hierachical clustering	by rows (reference labels).
#' @param gap Default: TRUE. Whether to add gaps between columns to separate query clusters.
#' @return A confusion matrix if trueAnno has been provided.
#' @export
CapybaraHeatmap <- function(query, trueAnno = NULL, cluCols = FALSE, cluRows = FALSE, gap = T, ...){

  score <- metadata(query)[['Capybara']][colnames(query),]
  N <- nrow(score)
  if(is.null(trueAnno)){
    predAnno <- colnames(score)[apply(score,1,which.max)]
    names(predAnno) <- rownames(score)
    #create column annotation
    colAnno<- data.frame(Predicted = predAnno, row.names = rownames(score))
    predAnno <- sort(predAnno)
    score <- score[names(predAnno),]

    if(!cluCols){
      #Find position of first cell in each cluster
      colLab <- data.frame(index = 1:N, anno =predAnno)
      colLabFirst <- plyr::ddply(colLab, 'anno', function(x) x[1,1])
      #First cells are labeled, otherwise put blank.
      colLab <- rep('',N)
      colLab[colLabFirst$V1] <- as.character(colLabFirst$anno)
      showColnames <- T

    }else{
      colLab <- NULL
      showColnames <- F

    }

    if(gap) gapsCol <-colLabFirst$V1[-1]-1 else gapsCol <- NULL

    pheatmap::pheatmap(t(score), color = viridis::viridis(20), annotation_col = colAnno, cluster_cols = cluCols,
                       cluster_rows = cluRows, labels_col = colLab, show_colnames = showColnames, gaps_col = gapsCol, ...)

  }else{
    anno <- colData(query)[,trueAnno]
    predAnno <- colnames(score)[apply(score,1,which.max)]
    confusion <- table(predAnno, anno)
    #create column annotation
    colAnno <- data.frame(Predicted = predAnno, True = anno)
    rownames(colAnno) <- rownames(score)
    names(anno) <- rownames(score)
    anno <- sort(anno)
    score <- score[names(anno),]

    if(!cluCols){
      #Find position of first cell in each cluster
      colLab <- data.frame(index = 1:N, anno =anno)
      colLabFirst <- plyr::ddply(colLab, 'anno', function(x) x[1,1])
      #First cells are labeled, otherwise put blank.
      colLab <- rep('',N)
      colLab[colLabFirst$V1] <- as.character(colLabFirst$anno)
      showColnames <- T

    }else{
      colLab <- NULL
      showColnames <- F

    }

    if(gap) gapsCol <-colLabFirst$V1[-1]-1 else gapsCol <- NULL

    pheatmap::pheatmap(t(score), color = viridis::viridis(20), annotation_col = colAnno, cluster_cols = cluCols,
                       cluster_rows = cluRows, labels_col = colLab, show_colnames = showColnames, gaps_col = gapsCol, ...)

    confusion

  }

}
