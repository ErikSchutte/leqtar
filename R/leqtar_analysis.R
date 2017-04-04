# Function    : leqtar_analysis
# Description : Acts as a wrapper for Matrix eQTL.
# Input       : Processed arguments
# Output      : Outputs the results to /results/data/

# leqtar_analysis ------------------------------------------------------------
#' leqtar_analysis function
#'
#' Uses linear regression analysis to identify the influence of a specific genotype on expression.
#'
#' @param arguments parsed and processed arguments
#' @import "MatrixEQTL"
leqtar_analysis <- function(arguments) {
  message("[INFO] Initializing linear regression analysis..")
  # Settings ------
  # Cis window
  cisDist <- 5e5

  # Statistical model
  useModel <- modelLINEAR

  # Covariates
  if ( !is.null(arguments$covariates) ) {
    covariates <- arguments$covariatesData
  } else {
    covariates <- character()
  }

  # P-value threshold
  pvOutputThreshold <- 0.5

  # Error covariance
  errorCovariance <- numeric()

  # Output file
  output_file_name <- tempfile()

  # Set genotype variables for analysis ----------------
  snps = SlicedData$new()
  snps$CreateFromMatrix( as.matrix(arguments$genotypeData) )
  snps$fileDelimiter <- "\t";      # the TAB character
  snps$fileOmitCharacters <- "NA"; # denote missing values;
  snps$fileSkipRows <- 1;          # one row of column labels
  snps$fileSkipColumns <- 1;       # one column of row labels
  snps$fileSliceSize <- 2000;      # read file in pieces of 2,000 rows

  # Set expression variables for analysis ----------------
  expr = SlicedData$new()
  expr$CreateFromMatrix( as.matrix(arguments$phenotypeData) )
  expr$fileDelimiter <- "\t";      # the TAB character
  expr$fileOmitCharacters <- "NA"; # denote missing values;
  expr$fileSkipRows <- 1;          # one row of column labels
  expr$fileSkipColumns <- 1;       # one column of row labels
  expr$fileSliceSize <- 2000;      # read file in pieces of 2,000 rows

  # Set covariates variables for analysis ----------------
  if ( !is.null(arguments$covariates) ) {
    covs = SlicedData$new()
    covs$CreateFromMatrix( as.matrix(arguments$covariatesData) )
    covs$fileDelimiter <- "\t";      # the TAB character
    covs$fileOmitCharacters <- "NA"; # denote missing values;
    covs$fileSkipRows <- 1;          # one row of column labels
    covs$fileSkipColumns <- 1;       # one column of row labels
    covs$fileSliceSize <- 2000;      # read file in pieces of 2,000 rows
  } else {
    covs = SlicedData$new()
    covs$fileDelimiter <- "\t";      # the TAB character
    covs$fileOmitCharacters <- "NA"; # denote missing values;
    covs$fileSkipRows <- 1;          # one row of column labels
    covs$fileSkipColumns <- 1;       # one column of row labels
    covs$fileSliceSize <- 2000;      # read file in pieces of 2,000 rows
  }


  # Analysis -----------------------
  message("[INFO] Running Linear Analysis..")

  # Supress matrix eqtl output
  sink("/dev/null")
  me = Matrix_eQTL_engine(
    snps = snps,
    gene = expr,
    cvrt = covs,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)
  sink()
  # Remove temp file.
  unlink(output_file_name)

  message("[INFO] Detected eQTLs: ", me$all$neqtls, " in ", unlist(me$time.in.sec[[1]]), " seconds..")

  # Save output.
  run_name <- paste(arguments$run_name, ".Rdata", sep="")
  save(me, file= file.path( arguments$output, "data", run_name, fsep=.Platform$file.sep) )

}
