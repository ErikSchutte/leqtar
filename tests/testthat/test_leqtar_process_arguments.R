library(leqtar)
library(testthat)
library(rprojroot)
context("Testing input arguments are valid arguments per leqtar's input standard..")

# Set package root
root <- is_testthat
root_file <- root$make_fix_file()

# Set arguments.
genotypeFile <- root_file("test_data", "genotype_test_data.RData")
expressionFile <- root_file("test_data", "expression_test_data.RData")
covariateFile <- root_file("test_data", "covariate_test_data.RData")
output_dir <- file.path( getwd() )
genoToFreq = F

# Test case 1 -----------------------------
test_that("arguments are correctly processed when parsing a genotypeFile and an expressionFile..", {
  arguments <- leqtar:::process_arguments(genotypeFile = genotypeFile, expressionFile = expressionFile,
                                          covariateFile = NULL, output_dir = NULL, genoToFreq)

  # Test if genotype data is loaded correctly.
  expect_true( file.exists( file.path(arguments$genotype) ), info = "Note: genotype file should exist in this scenario." )
  # Test if expression data is loaded correctly.
  expect_true( file.exists( file.path(arguments$expression) ), info = "Note: expression file should exist in this scenario." )
  # Test if the output folder is created.
  expect_true( dir.exists( file.path(arguments$output, fsep=.Platform$file.sep) ), info = "Note: leqtar_out should exist in this scenario" )
  # Test if the genoToFreq flag works.
  expect_false(arguments$genoToFreq, info = "Note: genoToFreq should be FALSE in this scenario")
  # Test if the process_arguments function validated the arguments.
  expect_true(arguments$valid, info = "Note: valid should be TRUE in this scenario")

  # clean-up directories.
  unlink( file.path(arguments$output, fsep=.Platform$file.sep), recursive = T)
})

# Test case 2 -----------------------------
test_that("arguments are correctly processed when parsing a genotypeFile, expressionFile and a covariatesFile", {
  arguments <- leqtar:::process_arguments(genotypeFile = genotypeFile, expressionFile = expressionFile,
                                          covariateFile = covariateFile, output_dir = NULL, genoToFreq)

  # Test if genotype data is loaded correctly.
  expect_true( file.exists( file.path(arguments$genotype) ), info = "Note: genotype file should exist in this scenario." )
  # Test if expression data is loaded correctly.
  expect_true( file.exists( file.path(arguments$expression) ), info = "Note: expression file should exist in this scenario." )
  # Test if covariate data is loaded correctly.
  expect_true( file.exists( file.path(arguments$covariate) ), info = "Note: covariate file should exist in this scenario." )
  # Test if the output folder is created.
  expect_true( dir.exists( file.path(arguments$output, fsep=.Platform$file.sep) ), info = "Note: leqtar_out should exist in this scenario" )
  # Test if the genoToFreq flag works.
  expect_false(arguments$genoToFreq, info = "Note: genoToFreq should be FALSE in this scenario")
  # Test if the process_arguments function validated the arguments.
  expect_true(arguments$valid, info = "Note: valid should be TRUE in this scenario")

  # clean-up directories.
  unlink( file.path(arguments$output, fsep=.Platform$file.sep), recursive = T)
})

# Test case 3 ----------------------------
test_that("arguments are correctly processed when parsing a genotypeFile, an expressionFile, a covariatesFile and a specific output directory", {
  arguments <- leqtar:::process_arguments(genotypeFile = genotypeFile, expressionFile = expressionFile,
                                          covariateFile = covariateFile, output_dir = output_dir, genoToFreq)

  # Test if genotype data is loaded correctly.
  expect_true( file.exists( file.path(arguments$genotype) ), info = "Note: genotype file should exist in this scenario." )
  # Test if expression data is loaded correctly.
  expect_true( file.exists( file.path(arguments$expression) ), info = "Note: expression file should exist in this scenario." )
  # Test if covariate data is loaded correctly.
  expect_true( file.exists( file.path(arguments$covariate) ), info = "Note: covariate file should exist in this scenario." )
  # Test if the output folder is created.
  expect_true( dir.exists( file.path(arguments$output, fsep=.Platform$file.sep) ), info = "Note: leqtar should exist in this scenario" )
  # Test if the output variable is set & created.
  expect_true( dir.exists( file.path(arguments$output) ), info = "Note: leqtar should exist in this scenario" )
  # Test if the subdirectories for the output directory is created
  expect_true( dir.exists( file.path(arguments$output, "images", fsep=.Platform$file.sep) ), info = "Note: leqtar_out/images should exist in this scenario" )
  expect_true( dir.exists( file.path(arguments$output, "data", fsep=.Platform$file.sep) ), info = "Note: leqtar_out/data should exist in this scenario" )
  # Test if the genoToFreq flag works.
  expect_false(arguments$genoToFreq, info = "Note: genoToFreq should be FALSE in this scenario")
  # Test if the process_arguments function validated the arguments.
  expect_true(arguments$valid, info = "Note: valid should be TRUE in this scenario")

  # clean-up directories.
  unlink( file.path(arguments$output, fsep=.Platform$file.sep), recursive = T )
})

# Test case 4 ----------------------------
test_that("arguments are correctly processed when parsing a genotypeFile, an expressionFile, a covariatesFile, a specific output directory and the genoToFreq flag", {
  arguments <- leqtar:::process_arguments(genotypeFile = genotypeFile, expressionFile = expressionFile,
                                          covariateFile = covariateFile, output_dir = output_dir, genoToFreq)

  # Test if genotype data is loaded correctly.
  expect_true( file.exists( file.path(arguments$genotype) ), info = "Note: genotype file should exist in this scenario." )
  # Test if expression data is loaded correctly.
  expect_true( file.exists( file.path(arguments$expression) ), info = "Note: expression file should exist in this scenario." )
  # Test if covariate data is loaded correctly.
  expect_true( file.exists( file.path(arguments$covariate) ), info = "Note: covariate file should exist in this scenario." )
  # Test if the output folder is created.
  expect_true( dir.exists( file.path(arguments$output, fsep=.Platform$file.sep) ), info = "Note: leqtar should exist in this scenario" )
  # Test if the output variable is set & created.
  expect_true( dir.exists( file.path(arguments$output) ), info = "Note: leqtar should exist in this scenario" )
  # Test if the subdirectories for the output directory is created
  expect_true( dir.exists( file.path(arguments$output, "images", fsep=.Platform$file.sep) ), info = "Note: leqtar_out/images should exist in this scneario" )
  expect_true( dir.exists( file.path(arguments$output, "data", fsep=.Platform$file.sep) ), info = "Note: leqtar_out/data should exist in this scenario" )
  # Test if the genoToFreq flag works.
  expect_false(arguments$genoToFreq, info = "Note: genoToFreq should be FALSE in this scenario")
  # Test if the process_arguments function validated the arguments.
  expect_true(arguments$valid, info = "Note: valid should be TRUE in this scenario")

  # clean-up directories.
  unlink( file.path(arguments$output, fsep=.Platform$file.sep), recursive = T)
})

