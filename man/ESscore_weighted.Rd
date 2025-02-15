% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GSEA.R
\name{ESscore_weighted}
\alias{ESscore_weighted}
\title{Weighted Enrichment Score Calculation}
\usage{
ESscore_weighted(labels.list, correl.vector = NULL)
}
\arguments{
\item{labels.list}{A binary vector indicating membership in a gene set (1 for inclusion, 0 for exclusion).}

\item{correl.vector}{A vector of correlation values (e.g., gene expression correlations).}
}
\value{
The weighted enrichment score (ES) for the given labels and correlation vector.
}
\description{
Calculates the weighted enrichment score (ES) for a given set of labels and correlation vector.
}
