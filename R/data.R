#' Example of a genotype data set.
#'
#' This data set contains some genotypes.
#'
#' @format A data frame with 200 rows and 500 variables:
#' \describe{
#'   \item{rs000001}{Example SNP with rs number 1}
#'   \item{rs000002}{Example SNP with rs number 2}
#'   \item{rs000003}{Example SNP with rs number 3}
#'   ...
#' }
"genotype_test_data"

#' Example of an expression data set.
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
"expression_test_data"

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
