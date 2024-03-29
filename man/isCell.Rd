% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quality_control.R
\name{isCell}
\alias{isCell}
\title{isCell}
\usage{
isCell(
  obj,
  num.test = 20000,
  num.tn5 = NULL,
  num.ref = 1000,
  background.cutoff = 100,
  min.pTSS = 0.2,
  min.FRiP = 0.2,
  min.pTSS.z = -2,
  min.FRiP.z = -2,
  verbose = F
)
}
\arguments{
\item{obj}{Object output from generateMatrix. Requires counts and meta data (buildMetaData) slots populated.}

\item{num.test}{Number of barcodes to query (ranked by total # of unique Tn5 insertions). Set to NULL to select barcodes
to test using obj minimum total number Tn5 insertions. Defaults to 20,000}

\item{num.tn5}{Set the minimum number of Tn5 insertions to select test cells. Overridden by num.test. Defaults to NULL.}

\item{num.ref}{Number of cells to use as the cell bulk reference (top X cells based on # unique Tn5 insertions)}

\item{background.cutoff}{Maximum unique Tn5 insertions to use for selecting barcodes for the background reference set}

\item{min.pTSS}{Minimum per cent TSS for good cells. Defaults to 0.2.}

\item{min.FRiP}{Minimum fraction reads in peaks for good cells. Defaults to 0.2.}

\item{min.pTSS.z}{Minimum z-score from per cent TSS for good cells. Defaults to -2.}

\item{min.FRiP.z}{Minimum z-score from fraction reads in peaks for good cells. Defaults to -2.}

\item{verbose}{Default to False. Set to TRUE to progress messages.}
}
\description{
This function compares each barcode to obj background and cell bulk reference to identify
barcodes representing ambient DNA or broken nuclei. Three columns are added to the metadata:
background, cellbulk, is_cell, which reflect the correlation of the cell's chromatin profile with
the background reference, the correlation with the predicted cellbulk reference, and whether a
barcode is predicted to be a cell (1 = cell, 0 = background noise).
}
