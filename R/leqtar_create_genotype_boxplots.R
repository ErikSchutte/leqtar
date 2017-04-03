# leqtar_create_genotype_boxplots ----------------------------
#' leqtar_create_genotype_boxplots
#'
#' Using the output from leqtar_analysis this function will try to create genotype boxplots for the top x associations found by Matrix eQTL.
#' @param output_data the folder containing the output from leqtar_analysis
#' @param output_img the folder whre genotype boxplots should go. Note it contains subfolders for manhattan plots and for genotype boxplots.
leqtar_create_genotype_boxplots <- function(output_data, output_img) {

  message("[INFO] ----------#----------")
  message("[INFO] Creating genotype boxplots..")
  list.files( file.path( output_data ) )
}
