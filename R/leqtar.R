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
#' Leqtar
#'
#' Leqtar uses a linear regression models for eQTL analysis. It has a multitude of options, play around with the test data that is availble to get a good idea of your optimal settings.
#' Leqtar tries to do much of the normal pre-processing work for you. Leqtar corrects sample differences, dimensions problems, 'cell' type differences and other things.
#' Please note to always double check your results and if there are any issues, please refer to https://github.com/ErikSchutte/leqtar/issues. This is a work in progress, I expect there to be bugs.
#'
#' @export
#' @param genotypeFile [REQUIRED] A Genotype file. SNP names are expected as colnames and sample names as rownames.
#' @param phenotypeFile [REQUIRED] A Phenotype file. Sample names are expected as colnames and stimulations/genes as rownames.
#' @param run_name [REQUIRED] A name that will be used for the current 'run'.
#' @param genotypePositionFile [OPTIONAL] A file containing the positions for each SNP. The colnames are 'snp', 'chr' and	'pos'.
#' @param phenotypePositionFile [OPTIONAl] A file containing the positions for each gene.
#' @param covariateFile [OPTIONAL] A Covariate file containing covariates. Sample names are expected as colnames and covariates as rownames.
#' @param output_dir [OPTIONAL] A relative path from your current working directory or an absolute path to a specific directory. The results will be stored in this location.
#' @param genoToFreq [OPTIONAL] Turns on conversion from genotypes i.e. 'AC' to a frequency for linear regression analysis.
#' @param forceRun [OPTIONAL] Normally Leqtar perserves data, by turning this to `TRUE` runs that already exist will be overwritten.
#' @note For a complete view of how to run and use Leqtar, please visit https://github.com/ErikSchutte/leqtar.
leqtar <- function(genotypeFile = NULL, phenotypeFile = NULL, covariateFile = NULL, output_dir = NULL, run_name = NULL, genoToFreq = F, forceRun = F) {

  message("[INFO] leqtar stands for Linear eQTL analysis in R",
          "\n[INFO] Thanks for using this package, if you find any bugs please report them on https://github.com/ErikSchutte/leqtar/issues",
          "\n[INFO] Package version ", build_version,
          "\n[INFO] This package was build on ", build_time)

  # Processes arguments, returns a list with all arguments.
  arguments <- process_arguments(genotypeFile, phenotypeFile, covariateFile, output_dir, run_name, genoToFreq, forceRun)

  # Parse arguments to leqtar_analysis.
  dataFiles <- leqtar_process_files(arguments)

  # Parse data from the input files for analysis.
  leqtar_analysis(dataFiles, arguments)

  # Process results.
  leqtar_results(arguments)

  message("[INFO] --------DONE!--------")
}
# process_arguments function -----------------------------------------------------
#' Processes the user input arguments
#'
#' @param genotypeFile the genotype file
#' @param phenotypeFile the expression file
#' @param output_dir the output dir, either default or user specified
#' @param covariateFile the covariates, either zero or defined
#' @param genoToFreq a flag that converts genotypes AA AT TT to frequencies 0 1 2
#' @return arguments object, containing all processed arguments
#' @importFrom "utils" "modifyList"
process_arguments <- function(genotypeFile, phenotypeFile, covariateFile, output_dir, run_name, genoToFreq, forceRun) {

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
  if ( is.null(phenotypeFile) ) {
    stop("[STOP] No expression file provided, provide a genotype file..")
  } else if ( is.matrix(phenotypeFile) || is.data.frame(phenotypeFile) ) {
    arguments <- modifyList(arguments, list(expression = phenotypeFile) )
  } else {
    if ( file.exists(phenotypeFile) ) {
      arguments <- modifyList(arguments, list(expression = phenotypeFile) )
    } else {
      stop("[STOP] The expression file you provided: ", as.character(phenotypeFile), " does not exist..")
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

  if ( is.null(run_name) ) {
    stop("[STOP] A name for the current run is required, preferably a unique one. See '?leqtar' for 'forceRun'..")
  } else {
    arguments <- modifyList(arguments, list(run_name = run_name) )
  }

  # Check if the output directory is specified, if not create the output directory.
  sub_arguments <- init_basic_directories(output_dir, run_name, arguments, forceRun)
  output_dir <- sub_arguments$output_dir
  arguments <- sub_arguments$arguments

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
  message("[INFO] ----------#----------")
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

# create intial directories ----------------------------------------------------
#' Creates initial file structure for leqtar output.
#' @param output_dir the output_dir variable set by the user.
#' @param run_name variable name for current run.
#' @param arguments the empty/partially empty arguments variable used for processing.
#' @return output_dir contains the output path for the current run
#' @return arguments modified arguments list.
init_basic_directories <- function(output_dir, run_name, arguments, forceRun) {

  leqtar_out = file.path("leqtar_out", fsep=.Platform$file.sep)

  if ( is.null(output_dir) ) {
    output_dir <- getwd()
    message("[INFO] No output_dir specified, using default ", as.character( file.path(output_dir, leqtar_out, run_name, fsep=.Platform$file.sep) ), " as output directory..")
  } else {
    message("[INFO] output_dir specified, using ", as.character( file.path( output_dir, leqtar_out, run_name, fsep=.Platform$file.sep) ), " as output directory..")
  }

  if ( !dir.exists( file.path(output_dir, leqtar_out, fsep=.Platform$file.sep) ) ) {
    message("\\___   The output directory 'leqtar_out' does not yet exists, creating leqtar_out..")
    dir.create( file.path(output_dir, leqtar_out, fsep=.Platform$file.sep) )
  }

  if ( !dir.exists( file.path(output_dir, leqtar_out, run_name, fsep=.Platform$file.sep) ) ) {
    message("\\___   The output direcotry does not yet exist, creating ", as.character( file.path(output_dir, leqtar_out, run_name, fsep=.Platform$file.sep) ), ".." )
    dir.create( file.path(output_dir, leqtar_out, run_name, fsep=.Platform$file.sep) )

  } else if ( forceRun == T && dir.exists( file.path( output_dir, leqtar_out, run_name) ) ) {
    message("\\___   Overwrite directory for current run: ", as.character( file.path(output_dir, leqtar_out, run_name, fsep=.Platform$file.sep) ), "..")
    unlink( file.path(output_dir, leqtar_out, run_name, fsep=.Platform$file.sep), recursive = T )
    dir.create( file.path(output_dir, leqtar_out, run_name, fsep=.Platform$file.sep) )

  } else {
    stop("[STOP] The directory with run name: ", as.character( file.path(output_dir, leqtar_out, run_name, fsep=.Platform$file.sep) ), " already exists..\n\\___   Please choose a different name or use 'forceRun = TRUE' to overwrite data..")
  }

  # Merge paths.
  leqtar_out <- file.path(output_dir, leqtar_out, run_name, fsep=.Platform$file.sep)

  if ( !dir.exists( file.path(leqtar_out, "data", fsep=.Platform$file.sep) ) ) {
    message("\\___   Creating subdirectory 'data' for run: ", run_name, "..")
    dir.create( file.path(leqtar_out, "data", fsep=.Platform$file.sep) )
  }
  if ( !dir.exists( file.path(leqtar_out, "data", "images", fsep=.Platform$file.sep) ) ) {
    message("\\___   Creating subdirectory 'images' for run: ", run_name, "..")
    dir.create( file.path(leqtar_out, "images", fsep=.Platform$file.sep) )

    if ( !dir.exists( file.path(leqtar_out, "data", "images", "manhattan", fsep=.Platform$file.sep) ) ) {
      message("\\___   Creating subdirectory 'manhattan' for run: ", run_name, "..")
      dir.create( file.path(leqtar_out, "images", "manhattan", fsep=.Platform$file.sep) )
    }
    if ( !dir.exists( file.path(leqtar_out, "data", "images", "genotype", fsep=.Platform$file.sep) ) ) {
      message("\\___   Creating subdirectory 'genotype' for run: ", run_name, "..")
      dir.create( file.path(leqtar_out, "images", "genotype", fsep=.Platform$file.sep) )
    }
  }
  if ( !dir.exists( file.path(leqtar_out, "tables", fsep=.Platform$file.sep) ) ) {
    message("\\___   Creating subdirectory 'tables' for run: ", run_name, "..")
    dir.create( file.path(leqtar_out, "tables", fsep=.Platform$file.sep) )
  }

  # Path to output directory.
  output_dir = file.path(leqtar_out, fsep=.Platform$file.sep)
  arguments <- modifyList(arguments, list(output = output_dir) )
  return( list( output_dir=output_dir, arguments=arguments) )
}
