library(leqtar)
context("Testing input arguments are valid arguments per leqtar's input standard..")

# Set files
genotypeFile = file.path("..", "..", "data", "genotype_test_data.RData", fsep=.Platform$file.sep)
expressionFile = file.path("..", "..", "data", "expression_test_data.RData", fsep=.Platform$file.sep)
covariateFile = file.path("..", "..", "data", "covariate_test_data.RData", fsep=.Platform$file.sep)
output_dir = "~"

# Test case 1 -----------------------------
test_that("arguments are correctly processed when parsing a genotypeFile and an expressionFile..", {
  arguments <- leqtar:::process_arguments(genotypeFile = genotypeFile, expressionFile = expressionFile,
                                          covariateFile = NULL, output_dir = NULL)

  # Test if genotype data is loaded correctly.
  expect_true( file.exists( file.path( arguments$genotype) ), TRUE )
  # Test if expression data is loaded correctly.
  expect_true( file.exists( file.path( arguments$expression) ), TRUE )
  # Test if the output folder is created.
  expect_true( dir.exists( file.path( getwd(), "leqtar", fsep=.Platform$file.sep) ), TRUE )
  # Test if the process_arguments function validated the arguments.
  expect_true(arguments$valid, TRUE)

  # clean-up directories.
  unlink( file.path(getwd(), "leqtar", fsep=.Platform$file.sep), recursive = T)
})

# Test case 2 -----------------------------
test_that("arguments are correctly processed when parsing a genotypeFile, expressionFile and a covariatesFile", {
  arguments <- leqtar:::process_arguments(genotypeFile = genotypeFile, expressionFile = expressionFile,
                                          covariateFile = covariateFile, output_dir = NULL)

  # Test if genotype data is loaded correctly.
  expect_true( file.exists( file.path( arguments$genotype) ), TRUE )
  # Test if expression data is loaded correctly.
  expect_true( file.exists( file.path( arguments$expression) ), TRUE )
  # Test if covariate data is loaded correctly.
  expect_true( file.exists( file.path( arguments$covariate) ), TRUE )
  # Test if the output folder is created.
  expect_true( dir.exists( file.path( getwd(), "leqtar", fsep=.Platform$file.sep) ), TRUE )
  # Test if the process_arguments function validated the arguments.
  expect_true(arguments$valid, TRUE)

  # clean-up directories.
  unlink( file.path(getwd(), "leqtar", fsep=.Platform$file.sep), recursive = T)
})

# Test case 3 ----------------------------
test_that("arguments are correctly processed when parsing a genotypeFile, expressionFile and a covariatesFile", {
  arguments <- leqtar:::process_arguments(genotypeFile = genotypeFile, expressionFile = expressionFile,
                                          covariateFile = covariateFile, output_dir = output_dir)

  # Test if genotype data is loaded correctly.
  expect_true( file.exists( file.path( arguments$genotype) ), TRUE )
  # Test if expression data is loaded correctly.
  expect_true( file.exists( file.path( arguments$expression) ), TRUE )
  # Test if covariate data is loaded correctly.
  expect_true( file.exists( file.path( arguments$covariate) ), TRUE )
  # Test if the output folder is created.
  expect_true( dir.exists( file.path( getwd(), "leqtar", fsep=.Platform$file.sep) ), TRUE )
  # Test if the output variable is set & created.
  expect_true( dir.exists( file.path( arguments$output) ), TRUE )
  # Test if the subdirectories for the output directory is created
  expect_true( dir.exists( file.path( arguments$output, "images", fsep=.Platform$file.sep) ), TRUE )
  expect_true( dir.exists( file.path( arguments$output, "data", fsep=.Platform$file.sep) ), TRUE )
  # Test if the process_arguments function validated the arguments.
  expect_true(arguments$valid, TRUE)

  # clean-up directories.
  unlink( file.path(getwd(), "leqtar", fsep=.Platform$file.sep), recursive = T)
})
