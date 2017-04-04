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
  genotypeIsFilePath <- TRUE
  phenotypeIsFilePath <- TRUE
  if ( !is.null(arguments$covariates) ) {
    covariateIsFilePath <- TRUE
  }

  if ( is.matrix(arguments$genotype) | is.data.frame(arguments$genotype) ) {
    message("[INFO] Genotype file is an object, not file path.")
    genotypeIsFilePath <- FALSE
  }

  if ( is.matrix(arguments$phenotype) | is.data.frame(arguments$phenotype) ) {
    message("[INFO] phenotype file is an object, not file path.")
    phenotypeIsFilePath <- FALSE
  }

  if ( !is.null(arguments$covariates) ) {
    if ( is.matrix(arguments$covariates) | is.data.frame(arguments$covariates) ) {
      message("[INFO] Covariates file is an object, not file path.")
      covariateIsFilePath <- FALSE
    }
  }

  # Depending on path or object read or set content -------
  if ( all(genotypeIsFilePath, phenotypeIsFilePath) ) {

    # Extract file extensions --------------------
    message("[INFO] ----------#----------")
    message("[INFO] Checking extensions..")
    genotype_file_extension <- unlist(str_split(arguments$genotype, "\\."))[2]
    phenotype_file_extension <- unlist(str_split(arguments$phenotype, "\\."))[2]

    #Visibly bind 'covariates_file_extension'
    covariates_file_extension <- NULL
    if ( !is.null(arguments$covariates) ) {
      if ( covariateIsFilePath ) {
        covariates_file_extension <- unlist(str_split(arguments$covariates, "\\."))[2]
      }
    }
    message("[INFO] Checking extensions OK..")

    # Determine read methods.
    message("[INFO] ----------#----------")
    message("[INFO] Reading files..")
    # genotype --------------------
    genotype_file_content <- read_files(genotype_file_extension, arguments$genotype)

    # phenotype -------------------
    phenotype_file_content <- read_files(phenotype_file_extension, arguments$phenotype)

    # covariates -------------------
    #Visibly bind 'covariate_file_content'
    covariates_file_content <- NULL
    if ( !is.null(arguments$covariates) ) {
      if ( covariateIsFilePath ) {
        covariates_file_content <- read_files(covariates_file_extension, arguments$covariates)
      }

    }
  } else {
    genotype_file_content <- arguments$genotype
    phenotype_file_content <- arguments$phenotype
    if ( !is.null(arguments$covariates) ) {
      covariates_file_content <- arguments$covariates
    }
  }
  message("[INFO] Reading files OK..")

  # Check dimensions -------------------
  message("[INFO] ----------#----------")
  message("[INFO] Checking dimensions..")
  dim_genotype <- dim(genotype_file_content)
  dim_phenotype <- dim(phenotype_file_content)
  warnings <- 0
  if ( !is.null(arguments$covariates) ) {
    dim_covariates <- dim(covariates_file_content)

    if (dim_genotype[1] != dim_covariates[2]) {
      message("[WARN] The number of samples in your genotype and covariate files do not match..",
              "\n\\___   leqtar will try to correct for this sample indifference..")
      warnings <- warnings + 1
    } else if (dim_genotype[1] != dim_phenotype[2]) {
      message("[WARN] The number of samples in your genotype and phenotype files do not match..",
              "\n\\___   leqtar will try to correct for this sample indifference..")
      warnings <- warnings + 1
    } else if (dim_covariates[2] != dim_phenotype[2]) {
      message("[WARN] The number of samples in your covariates and phenotype files do not match..",
              "\n\\___   leqtar will try to correct for this sample indifference..")
      warnings <- warnings + 1
    }
  } else {
    if (dim_genotype[1] != dim_phenotype[2]) {
      message("[WARN] The number of samples in your genotype and phenotype files do not match..",
              "\n\\___   leqtar will try to correct for this sample indifference..")
      warnings <- warnings + 1
    }
  }
  if ( warnings > 0 ) {
    message("[INFO] Checking dimensions OK with: ", as.character(warnings), " warning(s) ..")
  } else {
    message("[INFO] Checking dimensions OK..")
  }


  # Changing genotypes to frequencies -------------------------
  message("[INFO] ----------#----------")
  message("[INFO] Checking genotype data..")

  if ( class( as.vector(genotype_file_content[1,1]) ) == "character" && arguments$genoToFreq == F ) {
    stop("[STOP] Detected characters in genotype data. If you want leqtar to change them to\n\\___   frequencies, set argument 'genoToFreq=T'..")

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
  message("[INFO] Checking genotype data OK..")


  # Check column name order ---------------------------------
  message("[INFO] ----------#----------")
  message("[INFO] Checking sample names..")

  # Bind 'covariates_samples'
  covariates_samples <- NULL
  genotype_samples <- colnames(genotype_file_content)
  phenotype_samples <- colnames(phenotype_file_content)
  if ( !is.null(arguments$covariates) ) {
    covariates_samples <- colnames(covariates_file_content)
  }

  if ( !is.null(arguments$covariates) ) {

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
      genotype_file_content <- genotype_file_content[,mixedsort (colnames(genotype_file_content) ), drop=F]
      phenotype_file_content <- phenotype_file_content[,mixedsort( colnames(phenotype_file_content) ), drop=F]
      covariates_file_content <- covariates_file_content[,mixedsort( colnames(covariates_file_content) ), drop=F]

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

  # Checking phenotype data -----------------
  message("[INFO] Checking phenotype data..")
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
  message("[INFO] Checking phenotype data OK..")
  message("[INFO] Processing files OK..")
  message("[INFO] ----------#----------")

  # Store content in global argument object.
  arguments$genotypeData <- genotype_file_content
  arguments$phenotypeData <- phenotype_file_content

  if ( !is.null(arguments$covariates) ) {
    arguments$covariatesData <- covariates_file_content
  }

  return( arguments )
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
      stop("[STOP] File extension: ", file_extension, " not supported!\nCurrently supports RData, txt, xlsx, csv and tsv files.")
    }
  },
  error=function(condition) {
    message("[ERR] Could not read in file, critical error..\n[ERR] Reason:\n[ERR] ", condition)
    return(NULL)
  },
  warning=function(condition) {
    message("[WARNING] Could not read in file..\n[WARNING] Reason:\n[WARNING] ", condition)
    return(NULL)
  },
  finally={
    message("[INFO] File: ", file_path, " read correctly..")
  })
  return(file_content)
}
