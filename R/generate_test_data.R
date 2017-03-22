#' Generate test data
#'
#' Generates test data files, mainly for debugging and testing.
generate_test_data <- function() {

  test_sample_names <- lapply(seq(101, 300, 1), function(i) {
    paste("sample", i, sep="")
  })
  test_snp_names <- lapply(seq(101, 600, 1), function(i) {
    paste("rs00000", i, sep="")
  })
  test_gene_names <- lapply(seq(1001, 6000, 1), function(i) {
    paste("ENSG0000", i, sep="")
  })

  return( list( test_sample_names=unlist(test_sample_names),
                test_snp_names=unlist(test_snp_names),
                test_gene_names=unlist(test_gene_names) ) )
}

# Generate samples, snps and genes.
test_cases <- generate_test_data()


# Genrate genotype data.
snps_geno <- sapply(test_cases$test_snp_names, function(snp) {
  pattern = stringi::stri_rand_strings(1, 2, pattern = "[ACTG]")
  data.frame(snp=stringi::stri_rand_strings(200, 2, pattern = paste("[",pattern,"]", sep="") ) )
})

genotype_test_data <- t(data.frame(snps_geno))
colnames(genotype_test_data) <- test_cases$test_sample_names
# dim(genotype_test_data)
# length(test_cases$test_snp_names)
# length(test_cases$test_sample_names)
# head(genotype_test_data)

# Generate expression data.
expression_test_data <- matrix(rbinom(10*1000, 3000, .5), ncol=200, nrow=5000)
rownames(expression_test_data) <- test_cases$test_gene_names
colnames(expression_test_data) <- test_cases$test_sample_names
# head(expression_test_data)
# tail(expression_test_data)
# plot(density(expression_test_data))

# Generate covariates.
covariate_test_data <- t(data.frame( replicate(200, sample(18:75, 1, rep=T) ), replicate(200, sample(0:1, 1, rep=T) ) ))
colnames(covariate_test_data) <- test_cases$test_sample_names
rownames(covariate_test_data) <- c("age", "gender")
# head(covariate_test_data)

# Save data.
# save(genotype_test_data, file="~/git/leqtar/data/genotype_test_data.RData")
# save(expression_test_data, file="~/git/leqtar/data/expression_test_data.RData")
# save(covariate_test_data, file="~/git/leqtar/data/covariate_test_data.RData")

