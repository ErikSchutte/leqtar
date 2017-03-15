# Function    : leqtar_process_files
# Input       : processed arguments
# Output      : Data objects
# Note to self: 1 genotype, 2 expression, 3 covariates, 4 output, 5 valid

# leqtar_process_files -----------------------------------------------------------
#' leqtar_process_files function
#'
#' @param arguments path to file supplied by user.
#' @return content of data files.
#' @import "gdata"
#' @import "stringr"
#' @importFrom "utils" "read.table"
#' @importFrom "utils" "read.csv"
leqtar_process_files <- function(arguments) {

  # Extract file extensions --------------------
  message("\n[INFO] Determine extensions..")
  genotype_file_extension <- unlist(str_split(arguments$genotype, "\\."))[2]
  expression_file_extension <- unlist(str_split(arguments$expression, "\\."))[2]
  if ( !is.null(arguments$covariates) ) {
    covariates_file_extension <- unlist(str_split(arguments$covarites, "\\."))[2]
  }
  message("\n[INFO] File extensions determined..")
  # Determine read methods.
  message("\n[INFO] Reading files..")
  # genotype --------------------
  if (genotype_file_extension == "txt") {
    genotype_file_content <- read.table(arguments$genotype, stringsAsFactors = F, header = T, sep="\t")
  } else if (genotype_file_extension == "Rdata" | genotype_file_extension == "RData") {
    tmp_env <- new.env()
    load( arguments$genotype, tmp_env )
    genotype_file_extension = get( ls( tmp_env )[1], envir=tmp_env )
    rm(tmp_env)
  }else if (genotype_file_extension == "xlsx") {
    genotype_file_content <- read.xls(arguments$genotype)
  } else if (genotype_file_extension == "csv") {
    genotype_file_content <- read.csv(arguments$genotype, header = T, sep=",")
  } else if (genotype_file_extension == "tsv") {
    genotype_file_content <- read.table(arguments$genotype, stringsAsFactors = F, header = T, sep="\t")
  } else {
    stop("[STOP] File extension: ", genotype_file_extension, " not supported!\nCurrently supports txt, xlsx, csv and tsv files.")
  }
  # expression -------------------
  if (expression_file_extension == "txt") {
    expression_file_content <- read.table(arguments$expression, stringsAsFactors = F, header = T, sep="\t")
  } else if (expression_file_content == "Rdata" | expression_file_content == "RData") {
    tmp_env <- new.env()
    load( arguments$expression, tmp_env )
    expression_file_content = get( ls( tmp_env )[1], envir=tmp_env )
    rm(tmp_env)
  } else if (expression_file_extension == "xlsx") {
    expression_file_content <- read.xls(arguments$expression)
  } else if (expression_file_extension == "csv") {
    expression_file_content <- read.csv(arguments$expression, header = T, sep=",")
  } else if (expression_file_extension == "tsv") {
    expression_file_content <- read.table(arguments$expression, stringsAsFactors = F, header = T, sep="\t")
  } else {
    stop("[STOP] File extension: ", expression_file_extension, " not supported!\nCurrently supports txt, xlsx, csv and tsv files.")
  }
  # covariates -------------------
  if ( !is.null(arguments$covariates) ) {
    if (covariates_file_extension == "txt") {
      covariates_file_content <- read.table(arguments$covariates, stringsAsFactors = F, header = T, sep="\t")
    } else if (covariates_file_extension == "Rdata" | covariates_file_extension == "RData") {
      tmp_env <- new.env()
      load( arguments$covariates, tmp_env )
      covariates_file_content = get( ls( tmp_env )[1], envir=tmp_env )
      rm(tmp_env)
    } else if (covariates_file_extension == "xlsx") {
      covariates_file_content <- read.xls(arguments$covariates)
    } else if (covariates_file_extension == "csv") {
      covariates_file_content <- read.csv(arguments$covariates, header = T, sep=",")
    } else if (covariates_file_extension == "tsv") {
      covariates_file_content <- read.table(arguments$covariates, stringsAsFactors = F, header = T, sep="\t")
    } else {
      stop("[STOP] File extension: ", covariates_file_extension, " not supported!\nCurrently supports txt, xlsx, csv and tsv files.")
    }
  }
  message("\n[INFO] Files read..")

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
  message("\n[INFO] Dimensions match..")

  # Check column name order ---------------------------------
  message("\n[INFO] Checking column names..")
  genotype_colnames <- colnames(genotype_file_content)
  expression_colnames <- colnames(expression_file_content)
  if ( !is.null(arguments$covariates) ) {
    covariates_colnames <- colnames(covariates_file_content)
  }


}
