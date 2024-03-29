% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Update_auto_prior.R
\name{selected_clonotypes}
\alias{selected_clonotypes}
\title{Selected clones function for the automated process}
\usage{
selected_clonotypes(sc = sc, restrict_to_expanded = F)
}
\arguments{
\item{sc}{Add the merged and annotated file}

\item{Samp_col}{Sample column name with the default being "Sample_Name"}

\item{V_gene_sc}{column name vdj_gene_cdr3_AG_BD}

\item{colourtype}{colour types include: default, hcl.colors, topo.colors, heat.colors, terrain.colors, rainbow, random}

\item{NA_col_analysis}{If NA's are present it will make them coloured as grey90}

\item{font_type}{from the extrafont package; default: Times New Roman}

\item{title_size}{Size of the graph plot in 20}

\item{text_size}{Size of the text 12}

\item{Legend_size}{Size of the text 12}

\item{height_px}{Height of plot 1200}

\item{resolution_px}{Width of plot 144}

\item{save_file}{save file to director or you can store as an object}
}
\description{
This section is to add the factor order for the graphs
}
