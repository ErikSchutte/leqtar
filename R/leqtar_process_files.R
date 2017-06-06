# Function    : leqtar_process_files
# Input       : processed arguments
# Output      : Data objects
# Note to self: 1 genotype, 2 phenotype, 3 covariates, 4 output, 5 valid

# leqtar_process_files -----------------------------------------------------------
#' leqtar_process_files function
#'
#' Processes files based on their extensions.
#' Checks the order of column- and row- names and dimensions of the data sets.
#' If all checks pass, continue with genotype conversion.
#'
#' @param arguments path to file supplied by user.
#' @return content of data files.
#' @importFrom "gdata" "read.xls"
#' @importFrom "stringr" "str_split"
#' @importFrom "gtools" "mixedorder"
#' @importFrom "gtools" "mixedsort"
#' @importFrom "utils" "read.table"
#' @importFrom "utils" "read.csv"
#' @note Hard requirement, no dots should be present in the file name, except for the extension.
leqtar_process_files <- function(arguments) {

  # Check for data or file path --------------
  # Set flags
  message("[INFO] Processing Files..")
  message("[INFO] ----------#----------")
  message("[INFO] Reading files..")
  # Check paramters for file's or objects --------------------
  phenotype_file_content <- check_object_or_file(arguments$phenotype, arguments$phenotypeData, "Phenotype")
  phenotype_position_content <- check_object_or_file(arguments$phenotypePosition, arguments$phenotypePositionData, "Phenotype positions")

  genotype_file_content <- check_object_or_file(arguments$genotype, arguments$genotypeData, "Genotype")
  genotype_position_content <- check_object_or_file(arguments$genotypePosition, arguments$genotypePositionData, "Genotype positions")

  covariates_file_content <- check_object_or_file(arguments$covariates, arguments$covariatesData, "Covariates")
  message("[INFO] ----------#----------")
  message("[INFO] Processing Files OK..")

  # Check dimensions -------------------
  message("[INFO] ----------#----------")
  message("[INFO] Checking dimensions..")
  message("[INFO] ----------#----------")
  dim_genotype <- dim(genotype_file_content)
  dim_phenotype <- dim(phenotype_file_content)
  warnings <- 0
  if ( !is.null(covariates_file_content) ) {
    dim_covariates <- dim(covariates_file_content)

    if (dim_genotype[2] != dim_covariates[2]) {
      message("[WARN] The number of samples in your genotype and covariate files do not match..",
              "\n\\___   leqtar will try to correct for this sample indifference..")
      warnings <- warnings + 1
    } else if (dim_genotype[2] != dim_phenotype[2]) {
      message("[WARN] The number of samples in your genotype and phenotype files do not match..",
              "\n\\___   leqtar will try to correct for this sample indifference..")
      warnings <- warnings + 1
    } else if (dim_covariates[2] != dim_phenotype[2]) {
      message("[WARN] The number of samples in your covariates and phenotype files do not match..",
              "\n\\___   leqtar will try to correct for this sample indifference..")
      warnings <- warnings + 1
    }
  } else {
    if (dim_genotype[2] != dim_phenotype[2]) {
      message("[WARN] The number of samples in your genotype and phenotype files do not match..",
              "\n\\___   leqtar will try to correct for this sample indifference..")
      warnings <- warnings + 1
    }
  }
  message("[INFO] ----------#----------")
  if ( warnings > 0 ) {
    message("[INFO] Checking dimensions OK with: ", as.character(warnings), " warning(s) ..")
  } else {
    message("[INFO] Checking dimensions OK..")
  }


  # Changing genotypes to frequencies -------------------------
  message("[INFO] ----------#----------")
  message("[INFO] Checking genotype data..")
  message("[INFO] ----------#----------")
  if ( class( as.vector(genotype_file_content[1,1]) ) == "character" && arguments$genoToFreq == F ) {
    stop("[STOP] Detected characters in genotype data. If you want leqtar to change them to\n  \\___   frequencies, set argument 'genoToFreq=T'..")

  } else if ( class( as.vector(genotype_file_content[1,1]) ) == "character" && arguments$genoToFreq == T ) {
    message("[INFO] Genotype conversion: ", as.character(arguments$genoToFreq), ", conversing genotypes..")
    sub_arguments <- leqtar_genotypes_to_frequencies(genotype_file_content)
    genotype_file_content <- sub_arguments$genotypeConverted
    genotype_file_content_unconverted <- sub_arguments$genotypeNotConverted

  } else if ( class( as.vector(genotype_file_content[1,1]) )  == "integer" && arguments$genoToFreq == F ||
              class( as.vector(genotype_file_content[1,1]) )  == "numeric" && arguments$genoToFreq == F ) {
    message("[INFO] Genotype conversion: ", as.character(arguments$genoToFreq), ", but genotypes already conversed..")

  } else if ( class( as.vector(genotype_file_content[1,1]) )  == "integer" && arguments$genoToFreq == T ||
              class( as.vector(genotype_file_content[1,1]) )  == "numeric" && arguments$genoToFreq == T ) {
    message("[INFO] Genotype conversion: ", as.character(arguments$genoToFreq), ", but genotypes already conversed..")

  } else if ( class( as.vector(genotype_file_content[1,1]) ) == "factor" ) {
    stop("[STOP] Factor variables are not yet supported..")

  } else {
    stop("[STOP] Unexpected error, your genotype file is probably incorrect. If this is a persistent error,
         report the issue in the github issue tracker..")
  }
  message("[INFO] ----------#----------")
  message("[INFO] Checking genotype data OK..")


  # Check column name order ---------------------------------
  message("[INFO] ----------#----------")
  message("[INFO] Checking sample names..")
  message("[INFO] ----------#----------")
  # Bind 'covariates_samples'
  covariates_samples <- NULL
  genotype_samples <- colnames(genotype_file_content)
  phenotype_samples <- colnames(phenotype_file_content)
  if ( !is.null(covariates_file_content) ) {
    covariates_samples <- colnames(covariates_file_content)
  }

  if ( !is.null(covariates_file_content) ) {

    if ( length(genotype_samples) == length(phenotype_samples) &&
         length(genotype_samples) == length(covariates_samples) ) {
      if ( all(genotype_samples == phenotype_samples) && all(genotype_samples == covariates_samples) ) {
        message("[INFO] Checking sample names OK..")
      }
    }
    if ( length( intersect(genotype_samples, phenotype_samples) ) == length(genotype_samples) &&
         length( intersect(genotype_samples, covariates_samples) ) == length(genotype_samples) ) {
      message("[WARN] Sample names OK, but are in the wrong order..\n\\___   Re-ordering samples..")
      stop("[STOP] This is not yet implemented!")
    } else {
      message("[WARN] Some samples do not co-exists in all files..\n\\___   Re-ordering samples and trying to run anyway..")

      # Index the different samples for the covariates
      covariates_in_phenotype <- covariates_samples[which(covariates_samples %in% phenotype_samples)]
      covariates_in_genotype <- covariates_samples[which(covariates_samples %in% genotype_samples)]

      # Index the different samples for the covariates
      genotype_in_covariates <- genotype_samples[which(genotype_samples %in% covariates_samples)]
      genotype_in_phenotype <- genotype_samples[which(genotype_samples %in% phenotype_samples)]

      # Index the different samples for the covariates
      phenotype_in_genotype <- phenotype_samples[which(phenotype_samples %in% genotype_samples)]
      phenotype_in_covariates <- phenotype_samples[which(phenotype_samples %in% covariates_samples)]

      # Detect differences between samples.
      if ( length(covariates_samples) != length(covariates_in_phenotype) ) {
        message("[WARN] Detected different amount of samples between covariates file and the phenotype file..")
      } else if ( length(covariates_samples) != length(covariates_in_genotype) ) {
        message("[WARN] Detected different amount of samples between covariates file and the genotype file..")
      }
      if ( length(genotype_samples) != length(genotype_in_covariates) ) {
        message("[WARN] Detected different amount of samples between genotype file and the covariate file..")
      } else if ( length(genotype_samples) != length(genotype_in_phenotype) ) {
        message("[WARN] Detected different amount of samples between genotype file and the phenotype file..")
      }
      if ( length(phenotype_samples) != length(phenotype_in_genotype) ) {
        message("[WARN] Detected different amount of samples between phenotype file and the genotype file..")
      } else if ( length(phenotype_samples) != length(phenotype_in_covariates) ) {
        message("[WARN] Detected different amount of samples between phenotype file and the covariate file..")
      }

      # Define coexisting samples.
      coexistingSamples <- intersect( intersect(covariates_samples, genotype_samples), phenotype_samples)

      # Create subsets.
      genotype_file_content <- genotype_file_content[,coexistingSamples, drop=F]
      phenotype_file_content <- phenotype_file_content[,coexistingSamples, drop=F]
      covariates_file_content <- covariates_file_content[,coexistingSamples, drop=F]

      # Re-order data.
      genotype_file_content <- genotype_file_content[,mixedsort( colnames(genotype_file_content) ), drop=F]
      phenotype_file_content <- phenotype_file_content[,mixedsort( colnames(phenotype_file_content) ), drop=F]
      covariates_file_content <- covariates_file_content[,mixedsort( colnames(covariates_file_content) ), drop=F]
      if (arguments$genoToFreq == T) {
        genotype_file_content_unconverted <- genotype_file_content_unconverted[,coexistingSamples, drop=F]
        genotype_file_content_unconverted <- genotype_file_content_unconverted[, mixedsort( colnames(genotype_file_content_unconverted) ), drop=F]
      }

      # Output changes.
      message("[WARN] Initial number of phenotype samples: ", length(phenotype_samples),
              "\n\\___   Initial number of genotype samples: ", length(genotype_samples),
              "\n\\___   Initial number of covariate samples: ", length(covariates_samples),
              "\n\\___   Excluded: ", as.character( length(phenotype_samples) - length( colnames(phenotype_file_content) ) ),
              " samples from the phenotype data.",
              "\n\\___   Exlcuded: ", as.character( length(genotype_samples) - length( colnames(genotype_file_content) ) ),
              " samples from the genotype data.",
              "\n\\___   Excluded: ", as.character( length(covariates_samples) - length( colnames(covariates_file_content) ) ),
              " samples from the covariates data.")
    }
  } else {
    if ( length(genotype_samples) == length(phenotype_samples) ) {
      if ( all(genotype_samples == phenotype_samples) ) {
        message("[INFO] Checking sample names OK..")
      } else if ( length( intersect(genotype_samples, phenotype_samples) ) == length( unique(genotype_samples) ) ) {
        message("[WARN] Sample names OK, but are in the wrong order..\n\\___   Re-ordering samples..")
        stop("[STOP] This is not yet implemented!")
      }
    } else {
      message("[WARN] Some samples do not co-exists in both files..\n\\___   Re-ordering samples and trying to run anyway..")

      # Define coexisting samples.
      coexistingSamples <- intersect( genotype_samples, phenotype_samples)

      # Set subsets of data.
      phenotype_file_content <- phenotype_file_content[,coexistingSamples, drop=F]
      genotype_file_content <- genotype_file_content[,coexistingSamples, drop=F]

      # Re-order data.
      phenotype_file_content <- phenotype_file_content[,mixedsort( colnames(phenotype_file_content) ), drop=F]
      genotype_file_content <- genotype_file_content[,mixedsort( colnames(genotype_file_content) ), drop=F]
      if ( arguments$genoToFreq == T ) {
        genotype_file_content_unconverted <- genotype_file_content_unconverted[,coexistingSamples, drop=F]
        genotype_file_content_unconverted <- genotype_file_content_unconverted[, mixedsort( colnames(genotype_file_content_unconverted) ), drop=F]
      }
      # Output changes.
      message("[WARN] Initial number of phenotype samples: ", length(phenotype_samples),
              "\n\\___   Initial number of genotype samples: ", length(genotype_samples),
              "\n\\___   Excluded: ", as.character( length(phenotype_samples) - length( colnames(phenotype_file_content) ) ),
              " samples from the phenotype data.",
              "\n\\___   Exlcuded: ", as.character( length(genotype_samples) - length( colnames(genotype_file_content) ) ),
              " samples from the genotype data.",
              "\n\\___   New number of phenotype samples: ", length( colnames(phenotype_file_content) ),
              "\n\\___   New number of genotype samples: ", length ( colnames(genotype_file_content) ) )
    }
  }
  message("[INFO] ----------#----------")
  message("[INFO] Checking sample names OK..")
  message("[INFO] ----------#----------")

  # Checking phenotype data -----------------
  message("[INFO] Checking phenotype data..")
  message("[INFO] ----------#----------")
  if ( class( as.vector(phenotype_file_content[1,1]) ) == "character" ) {
    message("[WARN] Detected characters in phenotype data.\n\\___   Conversing to integers/numeric values..")

    # Save the number of NA's before conversion.
    numberOfNABefore <- sum( is.na(phenotype_file_content) )

    # Change the type of values in the data.frame
    phenotype_file_content <- as.matrix(phenotype_file_content)
    suppressWarnings( class(phenotype_file_content) <- "double" )

    # Check the number of Na's after conversion
    numberOfNAAfter <- sum( is.na(phenotype_file_content) )

    message("\\___   Number of NA's in phenotype data before conversion: ", numberOfNABefore,
            "\n\\___   Number of NA's in phenotype data after conversion: ", numberOfNAAfter)

  }
  message("[INFO] ----------#----------")
  message("[INFO] Checking phenotype data OK..")
  message("[INFO] ----------#----------")
  message("[INFO] Checking additional files..")
  message("[INFO] ----------#----------")

  # Genotype Position Data -----
  if ( !is.null(genotype_position_content) ) {
    expected_genotypePosDataCols <- c("snps", "chr", "pos")
    genotypePosDataCols <- colnames( arguments$genotypePositionData )
    if ( length( genotypePosDataCols ) == length( expected_genotypePosDataCols ) ) {
      if ( all( genotypePosDataCols == expected_genotypePosDataCols ) ) {
        message("[INFO] Columns genotype position: OK..")
      } else {
        stop("[STOP] Please check the leqtar help, the column names of the genotype position file have to match 'snps chr pos'.." )
      }
    } else {
      stop("[STOP] Please check the leqtar help, there are either too many or too few columns in the genotype position file..",
           "\n  \\___   Expected number of columns in the genotype position file: ", length( expected_genotypePosDataCols ),
           "\n  \\___   Observed number of columns in the genotype position file: ", length( genotypePosDataCols ) )

    }
  }


  # Phenotype Position Data -----
  if ( !is.null(phenotype_position_content) ) {
    expected_phenotypePosDataCols <- c("geneid", "chr", "s1", "s2")
    phenotypePosDataCols <- colnames(arguments$phenotypePositionData)
    if ( length( phenotypePosDataCols ) == length( expected_phenotypePosDataCols ) ) {
      if ( all( phenotypePosDataCols == expected_phenotypePosDataCols ) ) {
        message("[INFO] Columns phenotype position: OK..")
      } else {
        stop("[STOP] Please check the leqtar help, the column names of the phenotype position file have to match 'geneid chr s1 s2'.." )
      }
    } else {
      stop("[STOP] Please check the leqtar help, there are either too many or too few columns in the phenotype position file..",
           "\n  \\___   Expected number of columns in the phenotype position file: ", length( expected_phenotypePosDataCols ),
           "\n  \\___   Observed number of columns in the phenotype position file: ", length( phenotypePosDataCols ) )
    }
  }


  message("[INFO] ----------#----------")
  message("[INFO] Checking additional files OK..")
  message("[INFO] ----------#----------")
  message("[INFO] Processing files OK..")
  message("[INFO] ----------#----------")

  # Store content in global argument object.
  arguments$genotypeData <- genotype_file_content
  arguments$genotypePositionData <- genotype_position_content
  arguments$phenotypeData <- phenotype_file_content
  arguments$phenotypePositionData <- phenotype_position_content
  arguments$covariatesData <- covariates_file_content

  if ( arguments$genoToFreq == T ) {
    arguments$genotypeUnconvertedData <- genotype_file_content_unconverted
  }
  return( arguments )
}
# check_object_or_file -------------------
#' check_object_or_file
#'
#' Checks wether the given parameter is a file path or an R object. Thus determines wether to read it or not.
#'
#' @param path_argument the argument that should contain the path.
#' @param data_argument the argument that contains the object.
#' @param name_argument variable name to specify genotype phenotype covariaties or w/e in the messages.
#' @return If the path_argument is indeed given return the file content, else return the object.
check_object_or_file <- function(path_argument, data_argument, name_argument) {
  # If not required
  if ( is.null(path_argument) & is.null(data_argument) ) {
    message("[INFO] Not using ", name_argument, " argument..")
    return( NULL )
  }
  # Check if Object
  else if ( is.null(path_argument) ) {
    if ( is.matrix(data_argument) | is.data.frame(data_argument) ) {
      message("[INFO] ", name_argument, " argument: Object..")
      file_content <- data_argument
      return( file_content )
    } else {
      message("[WARN] ", name_argument, " argument is not a matrix or a data.frame.. ")
    }
  # Check if File
  } else {

    if ( is.null(data_argument) ) {
      if ( is.character(path_argument) ) {
        message("[INFO] ", name_argument, " argument: File path..")
        extension <- unlist(str_split(path_argument, "\\."))[2]
        file_content <- read_files(extension, path_argument)
        return( file_content )
      } else {
        stop("[STOP] ", name_argument, " argument has to be a file path.")
      }
    } else {
      stop("[STOP] Dev note: ", name_argument, " data and Genotype cannot both be assigned.")
    }
  }


}

# read_files function ---------------------------------------------
#' read_files function
#'
#' Tries to read in files determined by their file extensions.
#' Always returns the content of the given file.
#'
#' @param file_extension the file extension
#' @param file_path the file path
#' @return the file content
read_files <- function(file_extension, file_path) {
  file_content <- tryCatch({
    message("[INFO] File path: ", file_path, "\n\\___   File extension: ", file_extension)
    if (file_extension == "txt") {
      read.table(file_path, stringsAsFactors = F, header = T, sep="\t", row.names=NULL)
    } else if (file_extension == "Rdata" | file_extension == "RData") {
      tmp_env <- new.env()
      load( file_path, tmp_env )
      file_content = get( ls( tmp_env )[1], envir=tmp_env )
      rm(tmp_env)
      file_content
    } else if (file_extension == "xlsx" | file_extension == "xls") {
      read.xls(file_path)
    } else if (file_extension == "csv") {
      read.csv(file_path, header = T, sep=",")
    } else if (file_extension == "tsv") {
      read.table(file_path, stringsAsFactors = F, header = T, sep="\t")
    } else {
      stop("[STOP] File extension: ", file_extension, " not supported!\n\\___   Currently supports RData, txt, xlsx, csv and tsv files.")
    }
  },
  error=function(condition) {
    message("[ERR] Could not read in file, critical error..\n\\___   Reason:\n\\___   ", condition)
    return(NULL)
  },
  warning=function(condition) {
    message("[WARN] Could not read in file..\n\\___   Reason:\n\\___   ", condition)
    return(NULL)
  },
  finally={
    message("[INFO] File: ", file_path, " read correctly..")
  })
  return(file_content)
}
