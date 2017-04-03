#' @include zzz.R

# leqtar
# Erik Schutte
# https://github.com/ErikSchutte/leqtar/issues
# Function    : leqtar
# Description : Performs linear regression analysis to test the association between
#               genotype and expression data
# Input       : A genotype file/data.frame, expression file/data.frame
# Optional    : A covariate file/data.frame
# Output      : Outputs the results to /results/data/

# Build time and version ------------------------------------------------------------
build_time <- as.character(Sys.time())
build_version <- "0.1.0"
cat("Building package", build_version, "on", build_time, "\n")

# leqtar main function --------------------------------------------------------------
#' Main leqtar function
#'
#' @export
#' @param genotypeFile the genotype file for all samples
#' @param expressionFile the expresion file for all samples
#' @param output_dir the output directory where results are stored, defaults to your home folder
#' @param covariateFile the covariates you want to use in your association analysis
#' @param genoToFreq a variable that turns on conversion from genotypes i.e. 'AC' to a frequency for linear regression analysis.
#' @note genotypeFile and expressionFile are both required, the output_dir is set automatically and the covariateFile is optional
leqtar <- function(genotypeFile = NULL, expressionFile = NULL, covariateFile = NULL, output_dir = NULL, run_name = NULL, genoToFreq = F, forceRun = F) {

  message("[INFO] leqtar stands for Linear eQTL analysis in R",
                        "\n[INFO] Thanks for using this package, if you find any bugs please report them on https://github.com/ErikSchutte/leqtar/issues",
                        "\n[INFO] Package version ", build_version,
                        "\n[INFO] This package was build on ", build_time)

  # Processes arguments, returns a list with all arguments.
  arguments <- process_arguments(genotypeFile, expressionFile, covariateFile, output_dir, run_name, genoToFreq, forceRun)

  # Parse arguments to leqtar_analysis.
  dataFiles <- leqtar_process_files(arguments)

  # Parse data from the input files for analysis.
  leqtar_analysis(dataFiles, arguments)

  # Process results.
  leqtar_results(arguments)
}
# process_arguments function -----------------------------------------------------
#' Processes the user input arguments
#'
#' @param genotypeFile the genotype file
#' @param expressionFile the expression file
#' @param output_dir the output dir, either default or user specified
#' @param covariateFile the covariates, either zero or defined
#' @param genoToFreq a flag that converts genotypes AA AT TT to frequencies 0 1 2
#' @return arguments object, containing all processed arguments
#' @importFrom "utils" "modifyList"
process_arguments <- function(genotypeFile, expressionFile, covariateFile, output_dir, run_name, genoToFreq, forceRun) {

  message("[INFO] ----------#----------")
  # Bind the arguments variable.
  arguments <- list(genotype=NULL, expression=NULL, covariates=NULL,
                    output=NULL, run_name=NULL, genoToFreq=FALSE, forceRun = F, valid=FALSE)

  # Check if the genotype file is provided/exists.
  if ( is.null(genotypeFile) ) {
    stop("[STOP] No genotype file provided, provide a genotype file..")
  } else if ( is.matrix(genotypeFile) | is.data.frame(genotypeFile) ) {
    arguments <- modifyList(arguments, list(genotype = genotypeFile) )
  } else {
    if ( file.exists(genotypeFile) ) {
      arguments <- modifyList(arguments, list(genotype = genotypeFile) )
    } else {
      stop("[STOP] The genotype file you provided: ", as.character(genotypeFile), " does not exist..")
    }
  }

  # Check if the expression file is provided/exists.
  if ( is.null(expressionFile) ) {
    stop("[STOP] No expression file provided, provide a genotype file..")
  } else if ( is.matrix(expressionFile) || is.data.frame(expressionFile) ) {
    arguments <- modifyList(arguments, list(expression = expressionFile) )
  } else {
    if ( file.exists(expressionFile) ) {
      arguments <- modifyList(arguments, list(expression = expressionFile) )
    } else {
      stop("[STOP] The expression file you provided: ", as.character(expressionFile), " does not exist..")
    }
  }

  if ( is.null(covariateFile) ) {
    message("[INFO] No covariateFile specified, moving on..")
  } else if ( is.matrix(covariateFile) || is.data.frame(covariateFile) ) {
    arguments <- modifyList(arguments, list(covariates = covariateFile) )
  } else {
    if ( file.exists(covariateFile) ) {
      arguments <- modifyList(arguments, list(covariates = covariateFile) )
    } else {
      stop("[STOP] The covariate file you provided: ", as.character(covariateFile), " does not exist..")
    }
  }

  # Check if the output directory is specified, if not create the output directory.
  if ( is.null(output_dir) ) {
    home <- getwd()
    leqtar_out = file.path("leqtar_out", fsep=.Platform$file.sep)
    message("[INFO] No output_dir specified, using default ", as.character( file.path(home, leqtar_out, fsep=.Platform$file.sep) ), " as output directory..")

    if ( !dir.exists( file.path(home, leqtar_out, fsep=.Platform$file.sep) ) ) {
      message("\\___   The output direcotry does not yet exist, creating ", as.character( file.path(home, leqtar_out, fsep=.Platform$file.sep) ), ".." )
      dir.create( file.path(home, leqtar_out, fsep=.Platform$file.sep) )

      if ( !dir.exists( file.path( home, leqtar_out, run_name) ) ) {
        message("\\___   Creating directory for current run: ", run_name, "..")
        dir.create( file.path(home, leqtar_out, run_name, fsep=.Platform$file.sep) )
      } else if ( forceRun == T && dir.exists( file.path( home, leqtar_out, run_name) ) ) {
        message("\\___   Overwrite directory for current run: ", run_name, "..")
        unlist( file.path( home, leqtar_out, run_name), recursive = T )
        dir.create( file.path(home, leqtar_out, run_name, fsep=.Platform$file.sep) )
      } else {
        stop("[STOP] The directory with run name: ", run_name, " already exists, please choose a different name or use forceRun = TRUE to overwrite data..")
      }

      # Merge paths.
      leqtar_out <- file.path(home, leqtar_out, run_name, fsep=.Platform$file.sep)

      if ( !dir.exists( file.path(leqtar_out, "data", fsep=.Platform$file.sep) ) ) {
        message("\\___   Creating subdirectory 'data' for run: ", run_name, "..")
        dir.create( file.path(leqtar_out, "data", fsep=.Platform$file.sep) )
      }
      if ( !dir.exists( file.path(leqtar_out, "data", fsep=.Platform$file.sep) ) ) {
        message("\\___   Creating subdirectory 'data' for run: ", run_name, "..")
        dir.create( file.path(leqtar_out, "images", fsep=.Platform$file.sep) )

        if ( !dir.exists( file.path(leqtar_out, "data", "manhattan", fsep=.Platform$file.sep) ) ) {
          message("\\___   Creating subdirectory 'manhattan' for run: ", run_name, "..")
          dir.create( file.path(leqtar_out, "images", "manhattan", fsep=.Platform$file.sep) )
        }
      }
      if ( !dir.exists( file.path(leqtar_out, "data", fsep=.Platform$file.sep) ) ) {
        message("\\___   Creating subdirectory 'data' for run: ", run_name, "..")
        dir.create( file.path(leqtar_out, "info", fsep=.Platform$file.sep) )
      }



    } else {
      message("\\___   The output directory ", as.character( file.path(home, leqtar_out, fsep=.Platform$file.sep) ), " already exists, moving on..")
    }

    # Path to output directory.
    output_dir = file.path(home, leqtar_out, fsep=.Platform$file.sep)
    arguments <- modifyList(arguments, list(output = output_dir) )
  } else {
    leqtar_out = file.path("leqtar_out", fsep=.Platform$file.sep)
    message("[INFO] output_dir specified, using ", as.character( file.path( output_dir, leqtar_out, fsep=.Platform$file.sep) ), " as output directory..")

    if ( !dir.exists( file.path(output_dir, leqtar_out, fsep=.Platform$file.sep) ) ) {
      message("\\___   The output direcotry does not yet exist, creating ", as.character( file.path(output_dir, leqtar_out, fsep=.Platform$file.sep) ), ".." )
      dir.create( file.path(output_dir, leqtar_out, fsep=.Platform$file.sep) )
      dir.create( file.path(output_dir, leqtar_out, "data", fsep=.Platform$file.sep) )
      dir.create( file.path(output_dir, leqtar_out, "images", fsep=.Platform$file.sep) )
    } else {
      message("\\___   The output directory ", as.character( file.path(output_dir, leqtar_out, fsep=.Platform$file.sep) ), " already exists, moving on..")
    }

    # Path to output directory
    output_dir = file.path(output_dir, leqtar_out, fsep=.Platform$file.sep)
    arguments <- modifyList( arguments, list(output = output_dir) )
  }

  if ( genoToFreq == F ) {
    message("[INFO] Expecting genotypes that are already converted to numbers..")
    arguments <- modifyList( arguments, list(genoToFreq = genoToFreq) )
  } else {
    message("[INFO] Leqtar will try to turn genotypes (i.e. 'AA', 'AT', 'TT') into frequencies")
    arguments <- modifyList( arguments, list(genoToFreq = genoToFreq) )
  }

  arguments <- modifyList(arguments, list(valid = TRUE) )
  message("[INFO] All arguments are processed..")
  if (arguments$valid == T) {
    message("\\___   Status: OK, moving on..")
  } else {
    stop("\\___   Status: INVALID..\n Check your arguments, if you are convinced these are correct contact me on github.\n Shoot in an issue and don't forget your stacktrace() output..")
  }
  return(arguments)
}

# get_os function -----------------------------------------------------
#' Determine the current operating system.
#'
#' @return os as string, identifies the current os
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  return(c(.Platform$file.sep, tolower(os)))
}
