#' @include zzz.R

# leqtar
# Erik Schutte
# https://github.com/ErikSchutte/leqtar/issues
# Function    : leqtar
# Description : Performs linear regression analysis to test the association between
#               genotype and phenotype data
# Input       : A genotype file/data.frame, phenotype file/data.frame
# Optional    : A covariate file/data.frame
# Output      : Outputs the results to /results/data/

# Build time and version ------------------------------------------------------------
build_time <- as.character(Sys.time())
build_version <- "0.1.1"
cat("Building package", build_version, "on", build_time, "\n")

# leqtar main function --------------------------------------------------------------
#' Leqtar
#'
#' Leqtar uses a linear regression models for eQTL analysis. It has a multitude of options, play around with the test data that is availble to get a good idea of your optimal settings.
#' Leqtar tries to do much of the normal pre-processing work for you. Leqtar corrects sample differences, dimensions problems, 'cell' type differences and other things.
#' Please note to always double check your results and if there are any issues, please refer to https://github.com/ErikSchutte/leqtar/issues. This is a work in progress, I expect there to be bugs.
#'
#' @export
#' @param run_name [REQUIRED] - A name that will be used for the current 'run'. This name has to be unique to perserve data, if you want to overwrite a run see 'forceRun'.
#' @param genotypeFile [REQUIRED] - A Genotype file. SNP names are expected as colnames and sample names as rownames.
#' @param phenotypeFile [REQUIRED] - A Phenotype file. Sample names are expected as colnames and stimulations/genes as rownames.
#' @param useModel [REQUIRED] - A model used for QTL analysis. This can be 'linear', 'anova' or 'linear_cross'. This defaults to the 'linear' model.
#' @param genotypePositionFile [OPTIONAL] - A file containing the positions for each SNP. The colnames are 'snp', 'chr' and	'pos'.
#' @param phenotypePositionFile [OPTIONAL] - A file containing the positions for each gene.
#' @param covariateFile [OPTIONAL] - A Covariate file containing covariates. Sample names are expected as colnames and covariates as rownames.
#' @param geneNames [OPTIONAL] - A file containing the gene names that corrospond to Ensemble ID's or any other ID.
#' @param output_dir [OPTIONAL] - A relative path from your current working directory or an absolute path to a specific directory. The results will be stored in this location. This defaults to your current working directory.
#' @param genoToFreq [OPTIONAL] - Turns on conversion from genotypes i.e. 'AC' to a frequency for linear regression analysis. The default is set to 'FALSE'.
#' @param forceRun [OPTIONAL] - Normally Leqtar perserves data, by turning this to `TRUE` runs that already exist will be overwritten. The default is set to 'FALSE'.
#' @note For a complete view of how to run and use Leqtar, please visit https://github.com/ErikSchutte/leqtar.
leqtar <- function(run_name = NULL, genotypeFile = NULL, phenotypeFile = NULL, useModel = 'linear',
                    genotypePositionFile = NULL, phenotypePositionFile = NULL, covariateFile = NULL,
                   geneNames = NULL, output_dir = NULL,  genoToFreq = FALSE, forceRun = FALSE) {

  message("[INFO] leqtar stands for Linear eQTL analysis in R",
          "\n[INFO] Thanks for using this package, if you find any bugs please report them on https://github.com/ErikSchutte/leqtar/issues",
          "\n[INFO] Package version ", build_version,
          "\n[INFO] This package was build on ", build_time)

  # Processes arguments, returns a list with all arguments.
  arguments <- process_arguments(run_name, genotypeFile, genotypePositionFile, phenotypeFile, phenotypePositionFile,
                                 covariateFile, geneNames, useModel, output_dir, genoToFreq, forceRun)

  # Parse arguments to leqtar_analysis.
  arguments <- leqtar_process_files(arguments)

  # Parse data from the input files for analysis.
  leqtar_analysis(arguments)

  # Return used settings.
  message("[INFO] --------DONE!--------")
  return(arguments)

}
# process_arguments function -----------------------------------------------------
#' Processes the user input arguments
#'
#' @param run_name Paramter that defines the name of the current run.
#' @param genotypeFile A file/object containing genotypes.
#' @param genotypePositionFile A file/object containing the positions for each SNP.
#' @param phenotypeFile A file/object containing phenotypes.
#' @param phenotypePositionFile A file/object containing the positions for each Gene.
#' @param covariateFile A file/object contaning covariates for each sample in the genotype- and phenotype-File.
#' @param geneNames A file containing the gene names that corrospond to Ensemble ID's or any other ID.
#' @param useModel A string representing the model that should be used for QTL mapping.
#' @param output_dir A path were the output from Leqtar is stored.
#' @param genoToFreq A boolean flag, when set to 'TRUE' genotypes are converted to frequencies.
#' @param forceRun A boolean flag, when set to 'TRUE' an already existing run can be overwritten.
#' @return Arguments List containing all processed arguments.
#' @importFrom "utils" "modifyList"
#' @importFrom "stringr" "str_match"
process_arguments <- function(run_name, genotypeFile, genotypePositionFile, phenotypeFile, phenotypePositionFile,
                              covariateFile, geneNames, useModel, output_dir, genoToFreq, forceRun) {

  message("[INFO] ----------#----------")
  message("[INFO] Processing arguments..")
  message("[INFO] ----------#----------")
  # Bind the arguments variable.
  arguments <- list(run_name = NULL, genotype = NULL, genotypePosition = NULL, genotypeData = NULL, genotypePositionData = NULL, genotypeUnconvertedData = NULL,
                    phenotype = NULL, phenotypePosition = NULL, phenotypeData = NULL, phenotypePositionData = NULL, covariates=NULL,
                    covariatesData = NULL, geneNames = NULL, useModel = NULL, output = NULL, genoToFreq = FALSE,
                    forceRun = FALSE)

  if ( is.null(run_name) ) {
    stop("[STOP] A name for the current run is required, preferably a unique one. See '?leqtar' for 'run_name'..")
  } else {
    arguments <- modifyList(arguments, list(run_name = run_name) )
  }

  # Check if the genotype file is provided/exists.
  if ( is.null(genotypeFile) ) {
    stop("[STOP] No genotype file provided, provide a genotype file..")
  } else if ( is.matrix(genotypeFile) | is.data.frame(genotypeFile) ) {
    arguments <- modifyList(arguments, list(genotypeData = genotypeFile) )
  } else {
    if ( file.exists(genotypeFile) ) {
      arguments <- modifyList(arguments, list(genotype = genotypeFile) )
    } else {
      stop("[STOP] The genotype file you provided: ", as.character(genotypeFile), " does not exist..")
    }
  }

  # Check if the genotype position file is provided/exists.
  if ( is.null(genotypePositionFile) ) {
    stop("[STOP] No genotype position file provided, provide a genotype position file..")
  } else if ( is.matrix(genotypePositionFile) | is.data.frame(genotypePositionFile) ) {
    arguments <- modifyList(arguments, list(genotypePositionData = genotypePositionFile) )
  } else {
    if ( file.exists(genotypePositionFile) ) {
      arguments <- modifyList(arguments, list(genotypePosition = genotypePositionFile) )
    } else {
      stop("[STOP] The genotype position file you provided: ", as.character(genotypePositionFile), " does not exist..")
    }
  }

  # Check if the phenotype file is provided/exists.
  if ( is.null(phenotypeFile) ) {
    stop("[STOP] No phenotype file provided, provide a phenotype file..")
  } else if ( is.matrix(phenotypeFile) || is.data.frame(phenotypeFile) ) {
    arguments <- modifyList(arguments, list(phenotypeData = phenotypeFile) )
  } else {
    if ( file.exists(phenotypeFile) ) {
      arguments <- modifyList(arguments, list(phenotype = phenotypeFile) )
    } else {
      stop("[STOP] The phenotype file you provided: ", as.character(phenotypeFile), " does not exist..")
    }
  }

  # Check if the phenotype file is provided/exists.
  if ( is.null(phenotypePositionFile) ) {
    stop("[STOP] No phenotype position file provided, provide a phenotype position file..")
  } else if ( is.matrix(phenotypePositionFile) || is.data.frame(phenotypePositionFile) ) {
    arguments <- modifyList(arguments, list(phenotypePositionData = phenotypePositionFile) )
  } else {
    if ( file.exists(phenotypePositionFile) ) {
      arguments <- modifyList(arguments, list(phenotypePosition = phenotypePositionFile) )
    } else {
      stop("[STOP] The phenotype position file you provided: ", as.character(phenotypePositionFile), " does not exist..")
    }
  }

  if ( is.null(covariateFile) ) {
    message("[INFO] No covariateFile specified, moving on..")
  } else if ( is.matrix(covariateFile) || is.data.frame(covariateFile) ) {
    arguments <- modifyList(arguments, list(covariatesData = covariateFile) )
  } else {
    if ( file.exists(covariateFile) ) {
      arguments <- modifyList(arguments, list(covariates = covariateFile) )
    } else {
      stop("[STOP] The covariate file you provided: ", as.character(covariateFile), " does not exist..")
    }
  }

  if ( is.null(useModel) ) {
    stop("[STOP] You've overwritten the default model and set it to NULL. Please use a model or leave it at it's default value 'linear'..")

  } else if ( tolower(useModel) == "linear" ) {
    message("[INFO] Model: 'linear' (default)..")
    arguments <- modifyList( arguments, list(useModel = useModel) )

  } else if ( tolower(useModel) == "anova" ) {
    message("[INFO] Model: 'anova'..")
    arguments <- modifyList( arguments, list(useModel = useModel) )

  } else if ( tolower(useModel) == "linear_cross" ) {
    message("[INFO] Model: 'linear_cross'..")
    arguments <- modifyList( arguments, list(useModel = useModel) )

  } else {
    stop("[STOP] Model: ", as.character(useModel), " is not recognized by Leqtar..")
  }

  # Check for gene names files -----
  if (is.null(geneNames) ) {
    message("[INFO] geneName: missing, using package default. This might result in missing values.")
    gencode_names <- leqtar::gencode_names
    arguments <- modifyList(arguments, list(geneNames = gencode_names) )
  } else {
    message("[INFO] geneName: ", as.character(geneNames), "..")
    arguments <- modifyList(arguments, list(geneNames = geneNames) )
  }

  # Check how Force Run is set -----
  if ( forceRun == F ) {
    message("[INFO] Forced: ", as.character(forceRun), "..")
  } else if ( forceRun == T ) {
    message("[INFO] Forced: ", as.character(forceRun), "..")
  } else {
    stop("[STOP] The flag 'forceRun' has to be a boolean..")
  }
  arguments <- modifyList( arguments, list(forceRun = forceRun) )

  # Check if the output directory is specified, if not create the output directory -----
  if ( is.null(output_dir) ) {
    arguments$output <- getwd()
  } else {
    arguments$output <- output_dir
  }

  arguments <- init_basic_directories(arguments)

  # Check how genoToFreq is set -----
  if ( genoToFreq == F ) {
    message("[INFO] Genotype conversion: ", as.character(genoToFreq), " ..")
  } else if ( genoToFreq == T) {
    message("[INFO] Genotype conversion: ", as.character(genoToFreq), " ..")
  } else {
    stop("[STOP] The flag 'genoToFreq' has to be a boolean..")
  }
  arguments <- modifyList( arguments, list(genoToFreq = genoToFreq) )
  message("[INFO] ----------#----------")
  message("[INFO] Processing arguments OK..")
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
#' @param arguments the empty/partially empty arguments variable used for processing.
#' @return arguments the partially empty arguments variables used for processing.
init_basic_directories <- function(arguments) {

  leqtar_out = file.path("leqtar_out", fsep=.Platform$file.sep)
  message("[INFO] Output directory: ", as.character( file.path(arguments$output, leqtar_out, arguments$run_name, fsep=.Platform$file.sep) ), "..")

  if ( !dir.exists( file.path(arguments$output, leqtar_out, fsep=.Platform$file.sep) ) ) {
    message("\\___   The output directory: ", as.character(leqtar_out), ", does not yet exists, creating directory..")
    dir.create( file.path(arguments$output, leqtar_out, fsep=.Platform$file.sep) )
  }

  if ( !dir.exists( file.path(arguments$output, leqtar_out, arguments$run_name, fsep=.Platform$file.sep) ) ) {
    message("\\___   The output directory: ", as.character(arguments$run_name), ", does not yet exist, creating directory.." )
    dir.create( file.path(arguments$output, leqtar_out, arguments$run_name, fsep=.Platform$file.sep) )

  } else if ( arguments$forceRun == T && dir.exists( file.path( arguments$output, leqtar_out, arguments$run_name) ) ) {
    message("\\___   Overwrite directory for current run: ", as.character( file.path(arguments$output, leqtar_out, arguments$run_name, fsep=.Platform$file.sep) ), "..")
    unlink( file.path(arguments$output, leqtar_out, arguments$run_name, fsep=.Platform$file.sep), recursive = T )
    dir.create( file.path(arguments$output, leqtar_out, arguments$run_name, fsep=.Platform$file.sep) )

  } else {
    stop("[STOP] The directory with run name: ", as.character( file.path(arguments$output, leqtar_out, arguments$run_name, fsep=.Platform$file.sep) ), " already exists..\n  \\___   Please choose a different name or use 'forceRun = TRUE' to overwrite data..")
  }

  # Merge paths.
  leqtar_out <- file.path(arguments$output, leqtar_out, arguments$run_name, fsep=.Platform$file.sep)

  if ( !dir.exists( file.path(leqtar_out, "data", fsep=.Platform$file.sep) ) ) {
    message("\\___   Creating subdirectory 'data' for run: ", arguments$run_name, "..")
    dir.create( file.path(leqtar_out, "data", fsep=.Platform$file.sep) )
  }
  if ( !dir.exists( file.path(leqtar_out, "data", "images", fsep=.Platform$file.sep) ) ) {
    message("\\___   Creating subdirectory 'images' for run: ", arguments$run_name, "..")
    dir.create( file.path(leqtar_out, "images", fsep=.Platform$file.sep) )

    if ( !dir.exists( file.path(leqtar_out, "data", "images", "manhattan", fsep=.Platform$file.sep) ) ) {
      message("\\___   Creating subdirectory 'manhattan' for run: ", arguments$run_name, "..")
      dir.create( file.path(leqtar_out, "images", "manhattan", fsep=.Platform$file.sep) )
    }
    if ( !dir.exists( file.path(leqtar_out, "data", "images", "genotype", fsep=.Platform$file.sep) ) ) {
      message("\\___   Creating subdirectory 'genotype' for run: ", arguments$run_name, "..")
      dir.create( file.path(leqtar_out, "images", "genotype", fsep=.Platform$file.sep) )
    }
  }
  if ( !dir.exists( file.path(leqtar_out, "tables", fsep=.Platform$file.sep) ) ) {
    message("\\___   Creating subdirectory 'tables' for run: ", arguments$run_name, "..")
    dir.create( file.path(leqtar_out, "tables", fsep=.Platform$file.sep) )
  }

  # Path to output directory.
  arguments <- modifyList(arguments, list(output = leqtar_out) )
  return( arguments )
}
