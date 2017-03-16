# Function    : leqtar_process_files
# Input       : processed arguments
# Output      : Data objects
# Note to self: 1 genotype, 2 expression, 3 covariates, 4 output, 5 valid

# leqtar_process_files -----------------------------------------------------------
#' leqtar_process_files function
#'
#' Processes files based on their extensions.
#' Checks the order of column- and row- names and dimensions of the data sets.
#' If all checks pass, continue with genotype conversion.
#'
#' @param arguments path to file supplied by user.
#' @return content of data files.
#' @import "gdata"
#' @import "stringr"
#' @importFrom "utils" "read.table"
#' @importFrom "utils" "read.csv"
#' @note Hard requirement, no dots should be present in the file name, except for the extension.
leqtar_process_files <- function(arguments) {

  # Extract file extensions --------------------
  message("[INFO] Determine extensions..")
  genotype_file_extension <- unlist(str_split(arguments$genotype, "\\."))[2]
  expression_file_extension <- unlist(str_split(arguments$expression, "\\."))[2]
  if ( !is.null(arguments$covariates) ) {
    covariates_file_extension <- unlist(str_split(arguments$covariates, "\\."))[2]
  }
  message("[INFO] File extensions determined..")

  # Determine read methods.
  message("[INFO] Reading files..")
  # genotype --------------------
  genotype_file_content <- read_files(genotype_file_extension, arguments$genotype)

  # expression -------------------
  expression_file_content <- read_files(expression_file_extension, arguments$expression)

  # covariates -------------------
  if ( !is.null(arguments$covariates) ) {
    covariate_file_content <- read_files(covariate_file_extension, arguments$covariates)
  }
  message("[INFO] Files read..")

  # Check dimensions -------------------
  message("\n[INFO] Checking dimensions..")
  dim_genotype <- dim(genotype_file_content)
  dim_expression <- dim(expression_file_content)
  if ( !is.null(arguments$covariates) ) {
    dim_covariates <- dim(covariates_file_content)

    if (dim_genotype != dim_covariates) {
      stop("[STOP] Dimensions of your genotype and covariate files do not match!")
    } else if (dim_genotype != dim_expression) {
      stop("[STOP] Dimensions of your genotype and expression files do not match!")
    } else if (dim_covariates != dim_expression) {
      stop("[STOP] Dimensions of your covariates and expression files do not match!")
    }
  } else {
    if (dim_genotype != dim_expression) {
      stop("[STOP] Dimension of your genotype and expression files do not match!")
    }
  }
  message("[INFO] Dimensions match..")

  # Check column name order ---------------------------------
  message("[INFO] Checking column names..")
  genotype_colnames <- colnames(genotype_file_content)
  expression_colnames <- colnames(expression_file_content)
  if ( !is.null(arguments$covariates) ) {
    covariates_colnames <- colnames(covariates_file_content)
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
    message("[INFO] File extension: ", file_extension, "\n[INFO] File path: ", file_path)
    if (file_extension == "txt") {
      read.table(file_path, stringsAsFactors = F, header = T, sep="\t")
    } else if (file_extension == "Rdata" | file_extension == "RData") {
      tmp_env <- new.env()
      load( file_path, tmp_env )
      file_content = get( ls( tmp_env )[1], envir=tmp_env )
      rm(tmp_env)
    } else if (file_extension == "xlsx" | file_extension == "xls") {
      read.xls(file_path)
    } else if (file_extension == "csv") {
      read.csv(file_path, header = T, sep=",")
    } else if (file_extension == "tsv") {
      read.table(file_path, stringsAsFactors = F, header = T, sep="\t")
    } else {
      stop("[STOP] File extension: ", file_extension, " not supported!\nCurrently supports txt, xlsx, csv and tsv files.")
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
