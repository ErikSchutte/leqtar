# leqtar_create_genotype_boxplots ----------------------------
#' leqtar_create_genotype_boxplots
#'
#' Using the output from leqtar_analysis this function will try to create genotype boxplots for the top x associations found by Matrix eQTL.
#'
#' @param arguments the processed input arguments.
#' @importFrom "gtools" "mixedorder"
#' @importFrom "gtools" "mixedsort"
#' @import "ggplot2"
#' @importFrom "reshape2" "melt"
#' @importFrom "stringr" "str_split"
#' @export
leqtar_create_genotype_boxplots <- function(arguments) {

  message("[INFO] ----------#----------")
  message("[INFO] Creating genotype boxplots..")
  message("[INFO] ----------#----------")

  # Set paths.
  output_dir <- arguments$output
  output_data <- file.path( output_dir, "data", fsep = .Platform$file.sep)
  output_img <- file.path( output_dir, "images", fsep = .Platform$file.sep)
  output_tbl <- file.path( output_dir, "tables", fsep = .Platform$file.sep)
  # Get files.
  output <- list.files( file.path( output_data ), pattern = paste0(arguments$run_name, ".R[Dd]{1}ata") )

  # Get result files.
  result_set <- lapply(output, function(result_file) {

    # Load data.
    tmp_env <- new.env()
    load( file.path( output_data, output, fsep = .Platform$file.sep ), envir = tmp_env)
    f.data <- get( ls(tmp_env)[1], envir = tmp_env)
    rm(tmp_env)

    # Extract all qtls.
    qtls <- f.data$all$eqtls

    # Get significant qtls.
    treshold <- 5e-8
    qtls.subset <- qtls[which(qtls$pvalue < treshold),, drop =F]
    message("[INFO] Total QTLs: ", dim(qtls)[1], "\n\\___   Treshold Significant QTLs (< ", as.character(treshold), "): ", dim(qtls.subset)[1], "..")
    message("[INFO] Preparing QTL tables..")
    if (dim(qtls.subset)[1] > 0 ) {
      qtls.subset <- prepare_df_qtls(qtls.subset, arguments)
    } else {
      message("[INFO] Not able to write and visualize QTLs using this treshold..\n\\___   Showing only top one..")
      qtls.subset <- qtls[which(qtls$pvalue == range(qtls$pvalue)[1]),, drop = F]
      qtls.subset <- prepare_df_qtls(qtls.subset, arguments)
      message("[INFO] Writing all qtls to file..")
    }


      # Iterate over df.
      ll <- apply(qtls.subset, 1, function(qtl) {
        # Set each row as a dataframe instead of a vector
        qtl <- data.frame(t(qtl))

        # Retrieve all genotypes from all samples for the specified snps.
        if ( arguments$genoToFreq == T ) {
          genotypes <- arguments$genotypeUnconvertedData
          nCols <- ncol(genotypes) - 1
          genotypes <- genotypes[which( rownames(genotypes) == as.character(qtl$snps) ), 1:nCols, drop = F]

        } else {
          genotypes <- arguments$genotypeData
          genotypes <- genotypes[which( rownames(genotypes) == as.character(qtl$snps) ),]
        }

        # Transpose and to data frame.
        genotypes <- t(data.frame(genotypes))

        # Filter NA's.
        if (length( which(genotypes == "00") ) > 0) {
          genotypes[which(genotypes=="00")] <- NA
          samples.with.na <- which(is.na(genotypes)==T)
          genotypes = t(data.frame(genotypes[,-samples.with.na]))
        } else if ( length( which(is.na(genotypes) == T) ) > 0 ) {
          samples.with.na <- which(is.na(genotypes)==T)
          genotypes = t(data.frame(genotypes[,-samples.with.na]))
        }
        # colnames(genotypes) = gsub(colnames(genotypes), pattern="\\.[0-9]+", replacement="")

        # Melt the transposed genotypes.
        genotypes.melt <- suppressMessages( melt( t(genotypes) ) )

        # Position current gene in the expression file.
        gene.expression <- arguments$phenotypeData
        gene.expression <- gene.expression[ which( rownames( gene.expression ) == qtl$gene ), , drop =F]


        # Add a column to the gene expression dataframe and fill it with time points for later use.
        df.ge <- data.frame(gene.expression)

        colnames(df.ge) <- gsub(colnames(df.ge), pattern = "X", replacement = "")
        df.ge.melt <- suppressMessages(melt(df.ge))

        # Set sample names.
        if (length(names(gene.expression) > 0) ) {
          df.ge.melt[,1] <- names(gene.expression)
        }

        # Bind the dataframes together.
        df.melt <- cbind.data.frame(expression=df.ge.melt[,2],
                                    sample=df.ge.melt[,1],
                                    genotypes=genotypes.melt)

        # Create an order for the genotypes.
        if ( length( levels( df.melt$genotypes.value ) ) > 3 ) {
          df.melt[ which( df.melt$genotypes.value == levels(df.melt$genotypes.value)[[3]] ), "genotypes.value"] <- paste( rev( unlist( stringr::str_split( levels(df.melt$genotypes.value)[[3]], "") ) ), collapse = "" )
          df.melt$genotypes.value <- factor(df.melt$genotypes.value)
          genotype.levels <- levels(df.melt$genotypes.value)
        } else {
          genotype.levels <- levels(df.melt$genotypes.value)
        }
        if ( !is.null(arguments$genotypeUnconvertedData) ) {
          major.allele = strsplit(as.character(qtl$alleles),"_")[[1]][2]
          genotype.order <- c()

          if ( length(genotype.levels) == 2 ) { # Only 2 levels for genotypes.

            one <- strsplit(genotype.levels,"")[[1]] # first level
            two <- strsplit(genotype.levels,"")[[2]] # second level

            if (major.allele == one[1] && major.allele == one[2]) { # first level is homozygous major allele.
              # print("one is major")
              major = genotype.levels[1]
            } else if (major.allele == one[1] && major.allele != one[2] | major.allele != one[1] && major.allele == one[2]) { # first level is heteryzygous
              # print("one is hetero")
              heter = genotype.levels[1]
              genotype.order <- c(genotype.order, one)
            } else {
              print("WARNING: Third genotype should not exists here.")
            }

            if (major.allele == two[1] && major.allele == two[2]) { # first level is homozygous major allele.
              # print("two is major")
              major = genotype.levels[2]
            } else if (major.allele == two[1] && major.allele != two[2] | major.allele != two[1] && major.allele == two[2]) { # first level is heteryzygous
              # print("two is heter")
              heter = genotype.levels[2]
            } else {
              print("WARNING: Third genotype should not exists here.")
            }

            genotype.order <- c(major, heter)

          } else { # 3 or 4 levels for genotypes.
            # print(genotype.levels)

            one <- strsplit(genotype.levels,"")[[1]] # first level
            two <- strsplit(genotype.levels,"")[[2]] # second level
            three <- strsplit(genotype.levels,"")[[3]] #third level
            if (major.allele == one[1] && major.allele == one[2]) { # first level is major allele.
              # print("one is major")
              major = genotype.levels[1]
            } else if (major.allele == one[1] && major.allele != one[2] | major.allele != one[1] && major.allele == one[2]) { # first level is heteryzygous
              # print("one is hetero")
              heter = genotype.levels[1]
              genotype.order <- c(genotype.order, one)
            } else {
              # print("one is minor")
              minor = genotype.levels[1]
            }

            if (major.allele == two[1] && major.allele == two[2]) { # first level is homozygous major allele.
              # print("two is major")
              major = genotype.levels[2]
            } else if (major.allele == two[1] && major.allele != two[2] | major.allele != two[1] && major.allele == two[2]) { # first level is heteryzygous
              # print("two is heter")
              heter = genotype.levels[2]
            } else {
              minor = genotype.levels[2]
            }

            if (major.allele == three[1] && major.allele == three[2]) { # first level is homozygous major allele.
              # print("three is major")
              major = genotype.levels[3]
            } else if (major.allele == three[1] && major.allele != three[2] | major.allele != three[1] && major.allele == three[2]) { # first level is heteryzygous
              # print("three is heter")
              heter = genotype.levels[3]
            } else {
              minor = genotype.levels[3]
            }

            genotype.order <- c(major, heter, minor)
          }

          # Order the factor levels for the genotypes.
          df.melt$genotypes.value <- factor(df.melt$genotypes.value, levels=genotype.order)

          # Order data on factor levels.
          df.melt <- df.melt[order(df.melt$sample),]

          # Save ggplot in variable and plot.
          mi <- min( df.melt$expression ) - 1
          ma <- max(df.melt$expression ) + 1
          p <- ggplot(data=df.melt, aes(x=genotypes.value, y=expression, group=genotypes.value) ) +
            geom_boxplot(aes( fill=genotypes.value), outlier.shape=NA ) +
            geom_point( position=position_jitter(width=0.15),colour = "darkgrey") +
            coord_cartesian( ylim = c( mi,ma ) ) +
            ggtitle( paste(qtl$gene.name, " - ", qtl$snps, sep = ""),
                     subtitle = paste( "P-value: ", qtl$pvalue ) ) +
            theme( plot.title = element_text( size = rel(1.6), hjust = 0.5 ),
                   plot.subtitle = element_text(size = rel(1), hjust = 0.5 ) ) +
            xlab(paste("Genotypes",sep="")) + ylab("Norm. read count")

          p + scale_fill_discrete( name="Genotypes",
                                   labels=paste( names( table( df.melt$genotypes.value ) ),"(", table( df.melt$genotypes.value ), ")", sep ="") )


          suppressMessages(ggsave( filename=paste( output_img, "/genotype/", qtl$gene.name, "_", qtl$snps,".pdf", sep=""), plot=last_plot(), device = "pdf"))
        } else {
          df.melt <- cbind.data.frame(df.melt, rounded=round(df.melt$genotypes.value))
          df.melt <- df.melt[which(!is.na(df.melt$expression)),, drop = F]
          p <- ggplot(data=df.melt, aes(x=rounded, y=expression, group = rounded ) ) +
            geom_boxplot(aes( fill=rounded), outlier.shape=NA) +
            geom_point( position=position_jitter(width=0.15), colour = "darkgrey") +
            coord_cartesian( ylim = c( (min(df.melt$expression) -1 ), (max(df.melt$expression) +1) ) ) +
            ggtitle( paste( qtl$gene, " - ", qtl$snps, sep = ""),
                     subtitle = paste( "P-value: ", qtl$pvalue ) ) +
            theme( plot.title = element_text( size = rel(1.6), hjust = 0.5),
                   plot.subtitle = element_text(size= rel(1), hjust = 0.5) ) +
            xlab(paste("Rounded genotype dosages", sep = "") ) + ylab("Log2 cytokine levels")
          #breaks was as.vector( as.numeric( names(table( round( df.melt$genotypes.value ) ) ) ) )

          p + scale_fill_continuous(name="Rounded genotypes", breaks = as.vector( as.numeric( names(table( round( df.melt$genotypes.value ) ) ) ) ),
                                   labels=paste( "Group: ", names( table( round(df.melt$genotypes.value) ) )," (#", table( round(df.melt$genotypes.value) ), ")", sep ="") ) +
            scale_x_discrete( labels = paste( names( table( round(df.melt$genotypes.value) ) ),"(", table( round(df.melt$genotypes.value) ), ")", sep =""),
                               breaks =  as.vector( as.numeric( names(table( round( df.melt$genotypes.value ) ) ) ) ) )
          # scale_x_discrete( labels = paste( names( table( round(df.melt$genotypes.value) ) ),"(", table( round(df.melt$genotypes.value) ), ")", sep =""),
          #                   breaks =  as.vector( as.numeric( names(table( round( df.melt$genotypes.value ) ) ) ) ) )
          suppressMessages(ggsave( filename=paste( output_img, "/genotype/", qtl$gene, "_", qtl$snps,".pdf", sep=""), plot=last_plot(), device = "pdf"))
        }


      })


      # Save result tables.
      message("[INFO] Saving result tables..")
      write.table(qtls.subset, file = file.path(output_tbl, paste0(arguments$run_name, "_all_eqtls_", as.character(treshold), ".tsv"), fsep = .Platform$file.sep),
                  sep="\t", quote = F, row.names = F)

      message("[INFO] ----------#----------")
      message("[INFO] Creating genotype boxplots.. OK")
      message("[INFO] ----------#----------")

    message("[INFO] --------DONE!--------")
  })

}

# prepare_df_qtls ----------------
#' prepare_df_qtls
#'
#' prepares the matrix qtl output data.frame and transforms it slightly. Adding information
#' to easily create genotype plots.
#'
#' @param qtl a single qtl from matrix qtl.
#' @param arguments the processed input arguments.
#' @importFrom "data.table" "as.data.table" "setkey"
prepare_df_qtls <- function(qtl, arguments) {
  #Set temp df.
  tmp.df <- NULL

  # Col order boolean
  colOrder <- c(0, 0, 0)
  # Data tables ------------
  # Convert to data.tables
  qtl <-  as.data.table(qtl)
  if ( !is.null(arguments$geneNames ) ) {
    geneNames <- as.data.table( arguments$geneNames )
  }

  # Snp positions -------------
  if ( !is.null(arguments$genotypePositionData) ) {
    # If genotype.location file is given, add it to the table.
    genotype.loc <- as.data.table( arguments$genotypePositionData )

    # Create keys + merge
    setkey(qtl, snps)
    setkey(genotype.loc, snps)
    colNames <- colnames(qtl)
    qtl <- merge(qtl, genotype.loc, all.x=TRUE, by.x = "snps", by.y = "snps")
    colnames(qtl) <- c(colNames, "snps.chr", "snps.pos")
    # print("******* snps ********")
    # print(qtl)
    colOrder[1] <- 1
  }

  # Gene IDs -------------
  if ( !is.null(arguments$phenotypePositionData) ) {
    # If phenotype.location file is given, add it to the table.
    phenotype.loc <- as.data.table( arguments$phenotypePositionData )

    # Create keys + merge
    setkey(qtl, gene)
    setkey(phenotype.loc, geneid)
    colNames <- colnames(qtl)
    qtl <- merge(qtl, phenotype.loc, all.x = TRUE, by.x = "gene", by.y = "geneid")
    colnames(qtl) <- c(colNames, "gene.chr", "gene.start", "gene.end")
    colnames(qtl)[1] <- "gene"
    colnames(qtl)[2] <- "snps"
    # print("************** genes IDs *************")
    # print(qtl)

    # Gene Names -------------
    setkey(qtl, gene)
    setkey(geneNames, gene_id)

    qtl <- merge(qtl, geneNames, all.x = TRUE, by.x = "gene", by.y = "gene_id")
    # print("************** genes names *************")
    # print(qtl)
    colnames(qtl)[length(colnames(qtl))] <- "gene.name"
    colOrder[2] <- 3
  }

  # Minor alleles ------------
  if ( !is.null( arguments$genotypeUnconvertedData) ) {
    snps <- cbind.data.frame(snps=rownames(arguments$genotypeUnconvertedData), arguments$genotypeUnconvertedData)
    snps <- as.data.table( snps )
    cols <- c("snps", "MAF")
    snps <- snps[, cols, with=F]
    colnames(snps) <- c("snps", "alleles")

    # Set keys for qtls and for snps.
    setkey(qtl, snps)
    setkey(snps, snps)

    # Subset the snps data.frame before the join.
    qtl <- merge(qtl, snps, all.x=TRUE)
    # print("******* minor_alleles ********")
    # print(qtl)
    colOrder[3] <- 5
  }


  if ( sum(colOrder) == 0) {
    stop("[STOP] Dev error, create an issue in github..")
  } else if ( sum(colOrder) == 1 ) {
    # Re order tmp df cols.
    cols <- c("snps", "snps.chr", "snps.pos", "gene", "statistic", "pvalue", "FDR", "beta")
  } else if ( sum(colOrder) == 3 ) {
    cols <- c("snps", "gene", "gene.name", "gene.chr", "gene.start", "gene.end", "statistic", "pvalue", "FDR", "beta")
  } else if ( sum(colOrder) == 5 ) {
    cols <- c("snps", "alleles", "gene", "statistic", "pvalue", "FDR", "beta")
  } else if ( sum(colOrder) == 4 ) {
    cols <- c("snps", "snps.chr", "snps.pos", "gene", "gene.name", "gene.chr", "gene.start", "gene.end", "statistic", "pvalue", "FDR", "beta")
  } else if ( sum(colOrder) == 6 ) {
    cols <- c("snps", "snps.chr", "snps.pos", "alleles", "gene", "statistic", "pvalue", "FDR", "beta")
  } else if ( sum(colOrder) == 8 ) {
    cols <- c("snps", "alleles", "gene", "gene.name", "gene.chr", "gene.start", "gene.end",  "statistic", "pvalue", "FDR", "beta")
  } else if ( sum(colOrder) == 9 ) {
    cols <- c("snps", "snps.chr", "snps.pos", "alleles", "gene", "gene.name", "gene.chr", "gene.start", "gene.end", "statistic", "pvalue", "FDR", "beta")
  }
  tmp.df <- qtl[,cols, with=F]
  # print(head(tmp.df))
  # Return tmp df -------------
  return(tmp.df)
}
