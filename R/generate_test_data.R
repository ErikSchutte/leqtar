#' Generate test data
#'
#' Generates test data files, mainly for debugging and testing.
generate_test_data <- function() {

  # Create sample names
  test_sample_names <- lapply(seq(1, 30, 1), function(i) {
    paste("sample", i, sep="")
  })

  # Load default data.
  phpo <- leqtar::gencode_position
  gepo <- leqtar::gepo

  #Genrate genotype data.
  snps <- sample(gepo$snp, 30)
  snps_geno <- sapply(snps, function(snp) {
    pattern = stringi::stri_rand_strings(1, 2, pattern = "[ACTG]")
    data.frame(snp=stringi::stri_rand_strings(30, 2, pattern = paste("[",pattern,"]", sep="") ) )
  })

  genotype_test_data <- t(data.frame(snps_geno))
  colnames(genotype_test_data) <- test_sample_names
  rownames(genotype_test_data) <- snps
  genotype_test_locations <- gepo[which(gepo$snp %in% snps),]
  colnames(genotype_test_locations) <- c("snps", "chr", "pos")

  # Generate phenotype data.
  phenotype_genes <- sample(phpo$geneid, 50)
  phenotype_test_data <- matrix(rbinom(10*300, 10, .5), ncol=30, nrow=50)
  rownames(phenotype_test_data) <- phenotype_genes
  colnames(phenotype_test_data) <- test_sample_names
  phenotype_test_locations <- phpo[which(phpo$geneid %in% phenotype_genes),]

  # Generate covariates.
  covariate_test_data <- t(data.frame( replicate(30, sample(18:75, 1, rep=T) ), replicate(30, sample(0:1, 1, rep=T) ) ))
  colnames(covariate_test_data) <- test_sample_names
  rownames(covariate_test_data) <- c("age", "gender")

  # Save data.
  devtools::use_data(genotype_test_data, overwrite = T)
  devtools::use_data(genotype_test_locations, overwrite = T)
  devtools::use_data(phenotype_test_data, overwrite = T)
  devtools::use_data(phenotype_test_locations, overwrite = T)
  devtools::use_data(covariate_test_data, overwrite = T)
  # save(genotype_test_data, file="~/git/leqtar/data/genotype_test_data.RData", compress = "xz")
  # save(phenotype_test_data, file="~/git/leqtar/data/phenotype_test_data.RData", compress = "xz")
  # save(covariate_test_data, file="~/git/leqtar/data/covariate_test_data.RData", compress = "xz")
  # save(genotype_test_locations, file="~/git/leqtar/data/genotype_test_locations.RData", compress = "xz")
  # save(phenotype_test_locations, file="~/git/leqtar/data/phenotype_test_locations.RData", compress = "xz")


}
