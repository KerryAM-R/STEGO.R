% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mergin_harmony.R
\name{merging_ClusTCR2}
\alias{merging_ClusTCR2}
\title{Merging the ClusTCR2 files}
\usage{
merging_ClusTCR2(
  directory = "1_ClusTCR",
  chain = "AG",
  output_dir = "2_ClusTCR2"
)
}
\arguments{
\item{directory}{location of the multiple clusTCR2 files from step 1. default is "1_ClusTCR"}

\item{chain}{Either select AG or BD. The AG are in the same column so that the A will pair with B and G with a D.}

\item{output_dir}{location of where the file will be saved. default is "2_ClusTCR2"}
}
\value{
A merged file for clustering step based in ClusTCR2
}
\description{
This function is to aid im merging multiple Seurat object, which will then need to undergo harmony merging after this step. We have included a loop that merges two Seurat objects at a time to ensure time efficiency  (exponentially gets smaller from the previous sequential merging process)
}
