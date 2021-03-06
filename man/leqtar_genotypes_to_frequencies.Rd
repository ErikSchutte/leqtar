% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/leqtar_genotypes_to_frequencies.R
\name{leqtar_genotypes_to_frequencies}
\alias{leqtar_genotypes_to_frequencies}
\title{leqtar_genotypes_to_frequencies}
\usage{
leqtar_genotypes_to_frequencies(snps, VERBOSE = F)
}
\arguments{
\item{snps}{the genotype file content (genotypeData)}

\item{VERBOSE}{verbosity level, defaults to false}
}
\value{
genotypeData as numeric matrix.
}
\description{
Turns genotypes into frequencies to Matrix eQTL is able to perform linear regression analysis.
}
