# leqtar_create_manhattan_plots -------------------------
#' leqtar_create_manhattan_plots
#'
#' Creates manhattan plots from the data file produced by leqtar_analysis.
#'
#' @param arguments the processed input arguments.
#' @param output_data the folder containing the output from leqtar_analysis.
#' @param output_img the folder whre genotype boxplots should go. Note it contains subfolders for manhattan plots and for genotype boxplots.
#' @importFrom "gtools" "mixedorder"
#' @importFrom "stringr" "str_split"
#' @importFrom "qqman" "manhattan"
#' @export
leqtar_create_manhattan_plots <- function(arguments) {
  message("[INFO] ----------#----------")
  message("[INFO] Creating manhattan plots..")

  # Set paths.
  output_dir <- arguments$output
  output_data <- file.path( output_dir, "data", fsep = .Platform$file.sep)
  output_img <- file.path( output_dir, "images", fsep = .Platform$file.sep)
  output_tbl <- file.path( output_dir, "tables", fsep = .Platform$file.sep)

  # Get data
  load( file.path("..", "data", "modified_data_files", "uniqueShuffeledMergedSNPPositions.RData", fsep = .Platform$file.sep) )
  stimulations <- list.files( file.path( "..", "results", "cytokineQTLs", "log2", fsep = .Platform$file.sep), pattern = ".Rdata", full.names = T)

  lapply(stimulations, function(stimulation) {
    tmp_env <- new.env()
    load(stimulation, envir = tmp_env)
    stimulation_data <- get( ls( tmp_env)[1], envir = tmp_env )
    rm(tmp_env)

    # Subset results
    cytokines <- stimulation_data$all$eqtls[which(stimulation_data$all$eqtls$pvalue < 0.05),]

    # Should be the same length as the cytokine file.
    length(which(stimulation_data$all$eqtls$snps %in% uniqueShuffeledSNPPositions[,1]))

    # Order on snps.
    cytokines <- cytokines[mixedorder(cytokines$snps),]
    uniqueShuffeledSNPPositions <- uniqueShuffeledSNPPositions[mixedorder(uniqueShuffeledSNPPositions[,1]),]

    # Get snp info from positions file.
    index <- which(uniqueShuffeledSNPPositions[,1] %in% cytokines$snps)

    print( head( cytokines ) )
    print( head( uniqueShuffeledSNPPositions ) )
    # re-aragne the info for plot.
    newdf <- cbind.data.frame(uniqueShuffeledSNPPositions[index,], cytokines$pvalue)

    # Set name of stimulation
    path_blocks <- unlist(str_split(stimulation, .Platform$file.sep))

    fileName <- path_blocks[length(path_blocks)]

    fileName <- unlist(str_split(fileName, ":[0-9]{2}_"))[2]

    fileName <- unlist(str_split(fileName, "Rdata"))[1]
    plotTitle <- fileName
    tableName <- paste(fileName, "tsv", sep="")
    fileName <- paste(fileName, "png", sep="")


    colnames(newdf) <- c("SNP", "CHR", "BP", "P")

    # genome wide significance
    genwide_sig_snps <- newdf[which(newdf$P < 5e-08 ),]
    write.table(genwide_sig_snps, file= file.path("..", "results", "cytokineQTLs", "log2", "manhattan_plots", tableName, fsep=.Platform$file.sep),
                quote = F, sep="\t" )

    # Open device.
    png( file.path("..", "results", "cytokineQTLs", "log2", "manhattan_plots", fileName, fsep=.Platform$file.sep) )

    # Generate plots
    manhattan(newdf, chr = "CHR", bp = "BP", p = "P", snp = "SNP", main = plotTitle,
              cex = 0.5, cex.axis = 0.8, col = c("blue4", "orange3") )

    # Close device
    dev.off()
    cat("Processed: ", as.character(fileName), "\n" )

  })

}
