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
#' @import "qqman"
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
  #load( file.path("..", "data", "modified_data_files", "uniqueShuffeledMergedSNPPositions.RData", fsep = .Platform$file.sep) )
  #stimulations <- list.files( file.path( "..", "results", "cytokineQTLs", "log2", fsep = .Platform$file.sep), pattern = ".Rdata", full.names = T)
  # Get files.
  output <- list.files( file.path( output_data ), pattern = paste0(arguments$run_name, ".R[Dd]{1}ata") )

  lapply(output, function(stimulation) {
    tmp_env <- new.env()
    load(stimulation, envir = tmp_env)
    stimulation_data <- get( ls( tmp_env)[1], envir = tmp_env )
    rm(tmp_env)

    # Subset results
    cytokines <- stimulation_data$all$eqtls[which(stimulation_data$all$eqtls$pvalue < 0.05),]
    print("Subset significant cQTLs..")

    # Sort cytokines on snps
    cytokines <- cytokines[mixedorder(as.character(cytokines$snps)),]

    # Sort genomewide snps
    genomewideSortedSNPPositions <- genomewideSortedSNPPositions[which(as.character(genomewideSortedSNPPositions$snps) %in% as.character(cytokines$snps)),]
    genomewideSortedSNPPositions <- genomewideSortedSNPPositions[mixedorder(as.character(genomewideSortedSNPPositions$snps)),]

    # Create data.tables
    cytokines <- as.data.table(cytokines)
    genomewideSortedSNPPositions <- as.data.table(genomewideSortedSNPPositions)
    print(head(cytokines))
    print("*****")
    print(head(genomewideSortedSNPPositions))
    # Merge data.tables
    setkey(cytokines, snps)
    setkey(genomewideSortedSNPPositions, snps)
    newdf <- merge(cytokines, genomewideSortedSNPPositions, by.x = "snps", by.y = "snps", all.x = T)

    newdf <- newdf[, c("snps", "chr", "pos", "pvalue")]
    print(colnames(newdf))
    # Set name of stimulation
    print("Fabricating new names..")
    path_blocks <- unlist(str_split(stimulation, .Platform$file.sep))

    fileName <- path_blocks[length(path_blocks)]

    fileName <- unlist(str_split(fileName, ":[0-9]{2}_"))[2]

    fileName <- unlist(str_split(fileName, "Rdata"))[1]
    plotTitle <- fileName
    print(fileName)
    tableNameGenomewideSignificant <- paste(fileName, "genome_wide.tsv", sep="")
    tableNameHighlySignificant <- paste(fileName, "highly_signif.tsv", sep = "")
    fileName <- paste(fileName, "png", sep="")


    colnames(newdf) <- c("SNP", "CHR", "BP", "P")
    print(colnames(newdf))
    print(newdf)
    # genome wide significance
    print("Writing tables..")
    genwide_sig_snps <- newdf[which(newdf$P < 5e-08 ),]
    highly_significant_snps <- newdf[which(newdf$P < 1e-05 & newdf$P > 5e-08),]
    write.table(genwide_sig_snps, file= file.path("..", "results", type, "figures", "manhattan", tableNameGenomewideSignificant, fsep=.Platform$file.sep),
                quote = F, sep="\t" )
    write.table(highly_significant_snps, file= file.path("..", "results", type, "figures", "manhattan", tableNameHighlySignificant, fsep=.Platform$file.sep),
                quote = F, sep = "\t")

    # Open device.
    print("Creating plots..")
    png( file.path("..", "results", type, "figures", "manhattan", fileName, fsep=.Platform$file.sep) )

    # Generate plots
    manhattan(newdf, chr = "CHR", bp = "BP", p = "P", snp = "SNP", main = plotTitle,
              cex = 0.5, cex.axis = 0.8, col = c("blue4", "orange3") )

    # Close device
    dev.off()
    cat("Processed: ", as.character(fileName), "\n" )

  })

}
