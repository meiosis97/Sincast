% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CapybaraHeatmap.R
\name{CapybaraHeatmap}
\alias{CapybaraHeatmap}
\title{Plot Capybara heatmap}
\usage{
CapybaraHeatmap(
  query,
  trueAnno = NULL,
  cluCols = FALSE,
  cluRows = FALSE,
  gap = T,
  ...
)
}
\arguments{
\item{query}{Required. Query sce.}

\item{trueAnno}{Optional. The known annotation of the query cells.}

\item{cluCols}{Default: FALSE. Whether to perform hierachical clustering    by columns (query cells).}

\item{gap}{Default: TRUE. Whether to add gaps between columns to separate query clusters.}

\item{Npc}{Default: FALSE. Whether to perform hierachical clustering    by rows (reference labels).}
}
\value{
A confusion matrix if trueAnno has been provided.
}
\description{
A function to plot Capybara scores predicted on the query cells.
}
