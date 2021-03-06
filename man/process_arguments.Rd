% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/leqtar.R
\name{process_arguments}
\alias{process_arguments}
\title{Processes the user input arguments}
\usage{
process_arguments(run_name, genotypeFile, genotypePositionFile, phenotypeFile,
  phenotypePositionFile, covariateFile, geneNames, useModel, output_dir,
  genoToFreq, forceRun)
}
\arguments{
\item{run_name}{Paramter that defines the name of the current run.}

\item{genotypeFile}{A file/object containing genotypes.}

\item{genotypePositionFile}{A file/object containing the positions for each SNP.}

\item{phenotypeFile}{A file/object containing phenotypes.}

\item{phenotypePositionFile}{A file/object containing the positions for each Gene.}

\item{covariateFile}{A file/object contaning covariates for each sample in the genotype- and phenotype-File.}

\item{geneNames}{A file containing the gene names that corrospond to Ensemble ID's or any other ID.}

\item{useModel}{A string representing the model that should be used for QTL mapping.}

\item{output_dir}{A path were the output from Leqtar is stored.}

\item{genoToFreq}{A boolean flag, when set to 'TRUE' genotypes are converted to frequencies.}

\item{forceRun}{A boolean flag, when set to 'TRUE' an already existing run can be overwritten.}
}
\value{
Arguments List containing all processed arguments.
}
\description{
Processes the user input arguments
}
