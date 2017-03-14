#' @title leqtar
#' @author Erik Schutte
#' @seealso https://github.com/ErikSchutte/leqtar/issues
#' @
# Function    : leqtar
# Description : Performs linear regression analysis to test the association between
#               genotype and expression data.
# Input       : A genotype file/data.frame, expression file/data.frame.
# Optional    : A covariate file/data.frame.
# Output      : Outputs the results to /results/data/

# Build time and version ------------------------------------------------------------
build_time <- as.character(Sys.time())
build_version <- "0.1.0"
cat("Building package", build_version, "on", build_time, "\n")

# leqtar main function --------------------------------------------------------------
#' Main leqtar functiono
#'
#' @param genotypeFile the genotype file for all samples.
#' @param expressionFile the expresion file for all samples.
#' @param output_dir the output directory where results are stored, defaults to your home folder.
#' @param covariateFile the covariates you want to use in your association analysis.
#' @note genotypeFile and expressionFile are both required, the output_dir is set automatically and the covariateFile is optional.
leqtar <- function(genotypeFile = NULL, expressionFile = NULL, output_dir = NULL, covariateFile = NULL) {

  packageStartupMessage("
  [INFO] leqtar stands for Linear eQTL analysis in R
  [INFO] Thanks for using this package, if you find any bugs please report them on https://github.com/ErikSchutte/leqtar/issues
  [INFO] Package version ", build_version, "
  [INFO] This package was build on ", build_time,
                        "\n")

  # Processes arguments, returns a list with all arguments.
  arguments <- process_arguments(genotypeFile, expressionFile, output_dir, covariateFile)

}
#' Processes the user input arguments
#'
#' @param genotypeFile the genotype file.
#' @param expressionFile the expression file.
#' @param output_dir the output dir, either default or user specified.
#' @param covariateFile the covariates, either zero or defined.
#' @return arguments object, containing all processed arguments.
#' @importFrom "utils" "modifyList"
process_arguments <- function(genotypeFile, expressionFile, output_dir, covariateFile) {

  # Create an arguments object
  arguments <- list(genotype=NULL, expression=NULL, covariates=NULL,
                    output=NULL, valid=FALSE)

  # Check if the genotype file is provided/exists.
  if ( is.null(genotypeFile) ) {
    stop("[STOP] No genotype file provided, provide a genotype file..\n")
  } else {
    if ( file.exists(genotypeFile) ) {
      arguments <- modifyList(arguments, list(genotype = genotypeFile) )
    } else {
      stop("[STOP] The genotype file you provided does not exist..\n")
    }
  }

  # Check if the expression file is provided/exists.
  if ( is.null(expressionFile) ) {
    stop("[STOP] No expression file provided, provide a genotype file..\n")
  } else {
    if ( file.exists(expressionFile) ) {
      arguments <- modifyList(arguments, list(expression = expressionFile) )
    } else {
      stop("[STOP] The expression file you provided does not exist..\n")
    }
  }

  if ( is.null(covariateFile) ) {
    message("No covariateFile specified, moving on..\n")
  } else {
    if ( file.exists(covariateFile) ) {
      arguments <- modifyList(arguments, list(covariates = covariateFile) )
    } else {
      stop("[STOP] The covariate file you provided does not exist..\n")
    }
  }

  # Check if the output directory is specified, if not create the output directory.
  if ( is.null(output_dir) ) {
    message("No output_dir specified, using default '~/leqtar/' as output directory..")
    home = file.path("~/")
    leqtar_out = file.path("leqtar/")

    if ( !dir.exists( file.path(home, leqtar_out) ) ) {
      message("\\___ The output direcotry does not yet exist, creating ", as.character( paste(home, leqtar_out, sep="") ), "..\n" )
      dir.create( file.path(home, leqtar_out) )
      dir.create( file.path(home, leqtar_out, "data/") )
      dir.create( file.path(home, leqtar_out, "images/") )
    } else {
      message("\\___ The output directory ", as.character( paste(home, leqtar_out, sep="") ), " already exists, moving on..\n")
    }

    # Path to output directory.
    output_dir = file.path(home, leqtar_out)
    arguments <- modifyList(arguments, list(output = output_dir) )
  }

  arguments <- modifyList(arguments, list(valid = TRUE) )
  return(arguments)
}

