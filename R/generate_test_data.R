#' Generate test data
#'
#' Generates test data files, mainly for debugging and testing.
generate_test_data <- function() {

  test_sample_names <- lapply(seq_along(1:22), function(i) {
    paste("sample", i, sep="")
  })
  test_snp_names <- lapply(seq_along(1:50), function(i) {
    paste("rs00000", i, sep="")
  })
  test_gene_names <- lapply(seq_along(1:1000), function(i) {
    paste("ENSG0000", i, sep="")
  })

  return( list( test_sample_names=unlist(test_sample_names),
             test_snp_names=unlist(test_snp_names),
             test_gene_names=unlist(test_gene_names) ) )
}

test_cases <- generate_test_data()
test_cases$test_sample_names
