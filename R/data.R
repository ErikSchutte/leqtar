#' Example of a genotype data set.
#'
#' This data set contains some genotypes.
#'
#' @format A data frame with 200 rows and 500 variables:
#' \describe{
#'   \item{rs000001}{Example SNP with rs number 1.}
#'   \item{rs000002}{Example SNP with rs number 2.}
#'   \item{rs000003}{Example SNP with rs number 3.}
#'   ...
#' }
"genotype_test_data"

#' Example of a genotype position data set.
#'
#' This data set contains some genotypes.
#'
#' @format A data frame with 200 rows and 3 variables:
#' \describe{
#'   \item{snp}{SNP IDs.}
#'   \item{chr}{The chr location for each SNP.}
#'   \item{pos}{The position for each SNP.}
#'   ...
#' }
"genotype_test_locations"

#' Example of a phenotype position data set.
#'
#' This data set contains some expression values. These values represent raw read counts.
#'
#' @format A data frame with 5000 rows and 200 variables:
#' \describe{
#'   \item{sample1}{Sample 1 with expression in raw read counts per gene.}
#'   \item{sample2}{Sample 2 with expression in raw read counts per gene.}
#'   \item{sample3}{Sample 3 with expression in raw read counts per gene.}
#'   ...
#' }
"phenotype_test_data"

#' Example of a phenotype data set.
#'
#' This data set contains some expression values. These values represent raw read counts.
#'
#' @format A data frame with 5000 rows and 4 variables:
#' \describe{
#'   \item{geneid}{Gene IDs.}
#'   \item{chr}{The chr location for each gene.}
#'   \item{s1}{The start position for each gene.}
#'   \item{s2}{The stop position for each gene.}
#' }
"phenotype_test_locations"

#' Example of a covariate data set.
#'
#' This data set contains two covariates, age and gender. Gender is denoted as 0 for female and 1 for male.
#'
#' @format A data frame with 2 rows and 200 variables:
#' \describe{
#'   \item{sample1}{Sample 1 with covariates such as gender and time.}
#'   \item{sample2}{Sample 2 with covariates such as gender and time.}
#'   \item{sample3}{Sample 3 with covariates such as gender and time.}
#'   ...
#' }
"covariate_test_data"

#' Gencode gene id's and gene names.
#'
#' This data set contains a column with gencode gene ids and gencode gene names.
#'
#' @format A data frame with 60252 rows and 2 variables
#' \describe{
#'   \item{gene_id}{Gencode gene id's.}
#'   \item{gene_name}{Corrosponding gene name.}
#' }
"gencode_names"

#' Gene positions.
#'
#' This data set contain the positions for each gene id on the genome.
#'
#' @format A data frame with 60252 rows ans 4 variables
#' \describe{
#'   \item{geneid}{The Id's.}
#'   \item{chr}{located on chromosome .. .}
#'   \item{s1}{start position of the gene.}
#'   \item{s2}{stop position of the gene.}
#' }
"gencode_position"

#' Genotype positions.
#'
#' This data set contains the positions for each SNP.
#'
#' @format A data frame with 5493617 rows and 3 variables.
#' \describe{
#'   \item{snp}{The rs id for SNPs.}
#'   \item{chr}{The chromosome location.}
#'   \item{pos}{The position on the genomte.}
#' }
"gepo"
