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

# On load functionality -------------------------------------------------------------
.onLoad <- function(libname, pkgname) {
  op <- options()
  op.devtools <- list(
    devtools.path = "~/git/R-dev",
    devtools.install.args = "",
    devtools.name = "Erik Schutte",
    devtools.desc.author = '"Erik Schutte <schutte.erik@hotmail.com>"',
    devtools.desc.license = "What license is it under?",
    devtools.desc.suggests = NULL,
    devtools.desc = list()
  )
  toset <- !(names(op.devtools) %in% names(op))
  if(any(toset)) options(op.devtools[toset])

  invisible()
}

# leqtar main function --------------------------------------------------------------
leqtar <- function(genotypeFile = NULL, expressionFile = NULL, output_dir = NULL, covariateFile = NULL) {

  packageStartupMessage("
  [INFO] leqtar stands for Linear eQTL analysis in R
  [INFO] Thanks for using this package, if you find any bugs, please report them to schutte.erik@hotmail.com
  [INFO] Package version ", build_version, "
  [INFO] This package was build on ", build_time,
                        "\n")

  arguments <- check_arguments(genotypeFile, expressionFile, output_dir, covariateFile)
  print(arguments)
}

check_arguments <- function(genotypeFile, expressionFile, output_dir, covariateFile) {

  # Create an arguments object
  arguments <- list()

  # Check if the genotype file is provided/exists.
  if ( is.null(genotypeFile) ) {
    stop("[STOP] No genotype file provided, provide a genotype file..\n")
  } else {
    if ( file.exists(genotypeFile) ) {
      append(arguments, genotype=genotypeFile)
    } else {
      stop("[STOP] The genotype file you provided does not exist..\n")
    }
  }

  # Check if the expression file is provided/exists.
  if ( is.null(expressionFile) ) {
    stop("[STOP] No expression file provided, provide a genotype file..\n")
  } else {
    if ( file.exists(expressionFile) ) {
      append(arguments, expression=expressionFile)
    } else {
      stop("[STOP] The expression file you provided does not exist..\n")
    }
  }

  if ( is.null(covariateFile) ) {
    message("No covariateFile specified, moving on..\n")
  } else {
    if ( file.exists(covariateFile) ) {
      append(arguments, covariates=covariateFile)
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
    append(arguments, output_dir)
  }

  append(arguments, valid=TRUE)
  return(TRUE)
}

