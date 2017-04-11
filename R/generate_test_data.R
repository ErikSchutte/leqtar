#' Generate test data
#'
#' Generates test data files, mainly for debugging and testing.
generate_test_data <- function() {

  test_sample_names <- lapply(seq(1, 30, 1), function(i) {
    paste("sample", i, sep="")
  })
  test_snp_names <- lapply(seq(1, 30, 1), function(i) {
    paste("rs00000", i, sep="")
  })
  test_gene_names <- lapply(seq(11, 6000, 1), function(i) {
    paste("ENSG0000", i, sep="")
  })

  return( list( test_sample_names=unlist(test_sample_names),
                test_snp_names=unlist(test_snp_names),
                test_gene_names=unlist(test_gene_names) ) )
}

# Generate samples, snps and genes.
test_cases <- generate_test_data()
#
# gepo <- read.table("~/git/leqtar/data-raw/modified/mergedSNPPositions.txt", header=F, stringsAsFactors = F)
# phpo <- leqtar::gencode_position
#
# colnames(gepo) <- c("chr", "pos", "snp")
# gepo <- gepo[, c("snp", "chr", "pos")]
# gepo <- gepo[1:50,]
# devtools::use_data(gepo, overwrite = T)
#
# #Genrate genotype data.
# gepo <- leqtar::gepo
# snps <- sample(gepo$snp, 30)
# snps_geno <- sapply(snps, function(snp) {
#   pattern = stringi::stri_rand_strings(1, 2, pattern = "[ACTG]")
#   data.frame(snp=stringi::stri_rand_strings(30, 2, pattern = paste("[",pattern,"]", sep="") ) )
# })
#
# genotype_test_data <- t(data.frame(snps_geno))
# colnames(genotype_test_data) <- test_cases$test_sample_names
# rownames(genotype_test_data) <- snps
# dim(genotype_test_data)
# length(test_cases$test_snp_names)
# length(test_cases$test_sample_names)
# head(genotype_test_data)
#
# # Generate phenotype data.
# gn <- leqtar::gencode_names
# phenotype_test_data <- matrix(rbinom(10*300, 10, .5), ncol=30, nrow=50)
# rownames(phenotype_test_data) <- sample(gn$gene_id, 50)
# colnames(phenotype_test_data) <- test_cases$test_sample_names
# # head(phenotype_test_data)
# # tail(phenotype_test_data)
# # plot(density(phenotype_test_data))
#
# # Generate covariates.
# covariate_test_data <- t(data.frame( replicate(30, sample(18:75, 1, rep=T) ), replicate(30, sample(0:1, 1, rep=T) ) ))
# colnames(covariate_test_data) <- test_cases$test_sample_names
# rownames(covariate_test_data) <- c("age", "gender")
# # head(covariate_test_data)
#
# phenotype_test_locations <- phpo[which(sample(phpo$geneid, 50) %in% phpo$geneid),]
# genotype_test_locations <- gepo[which(snps %in% gepo$snp),]

# # Generate genotype position file.
# pos <- round( seq( from = 1, to = 100000, length.out = 30 ) )
# chr <- sort( rep_len( seq( from = 1, to = 22 ) , 30) )
# genotype_test_locations <- data.frame( snp = test_cases$test_snp_names,
#                                        chr = chr,
#                                        pos = pos)
# # Generate phenotype position file.
# s1 <- round( seq( from = 1, to = 4600000000, length.out = 50 ) )
# s2 <- round( seq( from = 400, to = 4600000000, length.out = 50 ) )
# chr <- sort( rep_len( seq( from = 1, to = 22 ) , 50 ) )
# phenotype_test_locations <- data.frame( geneid = test_cases$test_gene_names,
#                                         chr = chr,
#                                         s1 = s1,
#                                         s2 = s2)



# # Save data.
# devtools::use_data(genotype_test_data, overwrite = T)
# devtools::use_data(genotype_test_locations, overwrite = T)
# devtools::use_data(phenotype_test_data, overwrite = T)
# devtools::use_data(phenotype_test_locations, overwrite = T)
# devtools::use_data(covariate_test_data, overwrite = T)
# save(genotype_test_data, file="~/git/leqtar/data/genotype_test_data.RData", compress = "xz")
# save(phenotype_test_data, file="~/git/leqtar/data/phenotype_test_data.RData", compress = "xz")
# save(covariate_test_data, file="~/git/leqtar/data/covariate_test_data.RData", compress = "xz")
# save(genotype_test_locations, file="~/git/leqtar/data/genotype_test_locations.RData", compress = "xz")
# save(phenotype_test_locations, file="~/git/leqtar/data/phenotype_test_locations.RData", compress = "xz")
