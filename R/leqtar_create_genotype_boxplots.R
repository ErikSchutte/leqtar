# leqtar_create_genotype_boxplots ----------------------------
#' leqtar_create_genotype_boxplots
#'
#' Using the output from leqtar_analysis this function will try to create genotype boxplots for the top x associations found by Matrix eQTL.
#' @param output_data the folder containing the output from leqtar_analysis
#' @param output_img the folder whre genotype boxplots should go. Note it contains subfolders for manhattan plots and for genotype boxplots.
#' @import gtools
#' @import ggplot2
#' @import reshape2
leqtar_create_genotype_boxplots <- function(output_data, output_img) {

  message("[INFO] ----------#----------")
  message("[INFO] Creating genotype boxplots..")

  # Get files.
  output <- list.files( file.path( output_data ) )

  # Load data.
  tmp_env <- new.env()
  load( file.path( output_data, output, fsep = .Platform$file.sep ), envir = tmp_env)
  output.data <- get( ls(tmp_env)[1], envir = tmp_env)
  rm(tmp_env)

  # Extract all eqtls.
  eqtls <- output.data$all$eqtls

  # Get significant eqtls.
  eqtls.05 <- eqtls[which(eqtls$pvalue < 0.05),]

  message("[INFO] FURTHER IMPLEMENTATION REQUIRED")


  # ## Set path to sqtl result files.
  # data.folder = file.path("~/Dropbox/Erik Schutte Internship 2016/Results/eQTLs/")
  # eqtl.files = list.files(data.folder, pattern = "*.Rdata", recursive=T)
  # eqtl.files.cis = eqtl.files[grep(eqtl.files, pattern="cis")]
  # eqtl.files.cis.no.individual <- eqtl.files.cis[c(1,2,5)]
  #
  # ## Load expression data
  # load("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/gsTcell_gene_expr_vst_88samples_gencode.Rdata")
  #
  # ## Genotype file + Location file
  # load("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/gsTcell_genotype_frequencies_transformed.Rdata")
  # load("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/gsTcell_genotypes_transformed.Rdata")
  # genotype.data <- snps.genotype.t #snps.t.gencode OR snps.t
  # genotype.loc <- read.table("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/snp.txt",
  #                            header=T, stringsAsFactors = F)
  #
  # ## Loop over Rdata files.
  # lapply(eqtl.files.cis.no.individual, function(f) {
  #
  #   ## Create new environmet for data.
  #   tmp_env <- new.env()
  #
  #   ## Load Rdata file.
  #   load( file.path( paste(data.folder, f, sep="") ), tmp_env )
  #
  #   ## Test, remove when done.
  #   # load( file.path( paste(data.folder, eqtl.files[1], sep="") ), tmp_env )
  #   # f = eqtl.files[1]
  #
  #   ## Get data loaded, destory temp environment.
  #   f.data = get( ls( tmp_env )[1], envir=tmp_env )
  #   rm(tmp_env)
  #
  #   ## Set name of file.
  #   f.name = sub(f, pattern=".*/", replacement="")
  #   f.name = sub(f.name, pattern="\\.Rdata", replacement="")
  #
  #   ## Determine wether it's over all, interaction, or individual time points.
  #   if ( length( grep(f, pattern = "interaction" ) ) > 0 ) {
  #     f.type = 'interaction'
  #   } else if ( length( grep(f, pattern = "all" ) ) > 0 ) {
  #     f.type = "all"
  #   } else if ( length( grep(f, pattern = "individual" ) ) > 0 ) {
  #     f.type = "individual"
  #   } else {
  #     print("Unexpected: something is wrong!\nCheck your input files.")
  #   }
  #
  #   if ( length( grep(f, pattern = "cis" ) ) > 0 )  {
  #     f.acting = "cis"
  #   } else {
  #     f.acting = "trans"
  #   }
  #
  #   ## Determine wether the current file is specified with or without covariates.
  #   if ( length( grep(f, pattern="no_covs" ) ) > 0 ) {
  #     f.cov = "no_covs"
  #   } else {
  #     f.cov = "covs"
  #   }
  #
  #   ## If sub directories do not yet exist, create them for the plots.
  #   image.base = file.path("~/Dropbox/Erik Schutte Internship 2016/Results/eQTLs/")
  #   if ( !dir.exists( file.path( paste( image.base, f.type, "/", f.acting, "/", f.cov, "/", sep="") ) ) ) {
  #     dir.create( file.path( paste( image.base, f.type, "/", f.acting, "/", f.cov, "/", sep="") ), recursive = TRUE )
  #   }
  #
  #   ## Create table of written eqtls.
  #   subset <- 1:10
  #   subsetName <- paste( "top_10_", f.name, sep="")
  #   write.table(f.data, file=paste( image.base, f.type, "/", f.acting, "/", subsetName, ".tsv", sep=""), quote=F, sep="\t")
  #
  #   ## f.type determines sample size, 'all' and 'interaction' are based on 88 samples and 'individual' on 22 samples.
  #   ll <- apply(f.data, 1, function(eqtl) {
  #     # Set each row as a dataframe instead of a vector
  #     eqtl <- data.frame(t(eqtl))
  #
  #     ## Test, remove after
  #     # eqtl <- f.data[1,]
  #
  #     # Retrieve all genotypes from all samples for the specified snps.
  #     genotypes <- snps.genotype.t[as.character(eqtl$snps),1:22]
  #
  #     # Transpose and to data frame.
  #     genotypes <- t(data.frame(genotypes))
  #
  #     if ( f.type != "individual" ) {
  #       # Multiply the genotypes 4 times.
  #       genotypes <- cbind.data.frame(genotypes, genotypes, row.names = NULL) # 44x / 42x
  #       genotypes <- cbind.data.frame(genotypes, genotypes, row.names = NULL) # 88x / 84x
  #       genotypes <- genotypes[mixedsort(names(genotypes))]
  #       colnames(genotypes) <- gsub(colnames(genotypes), pattern="\\.[0-9]+", replacement="")
  #
  #       # Filter NA's.
  #       if (length( which(genotypes == "00") ) > 0) {
  #         genotypes[which(genotypes=="00")] <- NA
  #         samples.with.na <- which(is.na(genotypes)==T)
  #         genotypes = genotypes[,-samples.with.na]
  #       }
  #     } else {
  #       # Filter NA's.
  #       if (length( which(genotypes == "00") ) > 0) {
  #         genotypes[which(genotypes=="00")] <- NA
  #         samples.with.na <- which(is.na(genotypes)==T)
  #         genotypes = t(data.frame(genotypes[,-samples.with.na]))
  #       }
  #     }
  #
  #     colnames(genotypes) = gsub(colnames(genotypes), pattern="\\.[0-9]+", replacement="")
  #
  #
  #     # Melt the transposed genotypes.
  #     genotypes.melt <- melt(t(genotypes))
  #
  #     # Position current gene in the expression file.
  #     if ( f.type != "individual" ) {
  #       gene.expression <- ge.gencode[ which( eqtl$gene == rownames( ge.gencode ) ), ]
  #       if ( length( names(genotypes) ) < 88 ) {
  #         gene.expression <- gene.expression[-samples.with.na]
  #       }
  #     } else {
  #       times = list( t0=seq( 1, length( colnames(ge.gencode) ), 4 ),
  #                     t10=seq( 2, length( colnames(ge.gencode) ), 4 ),
  #                     t30=seq( 3, length( colnames(ge.gencode) ), 4 ),
  #                     t180=seq( 4, length( colnames(ge.gencode) ), 4 ) )
  #       gene.expression <- ge.gencode[ which( eqtl$gene == rownames( ge.gencode ) ), times[eqtl$t.interval][[1]] ]
  #       if ( length( colnames(genotypes) ) < 22 ) {
  #         gene.expression <- gene.expression[-samples.with.na]
  #       }
  #     }
  #
  #     # Add a column to the gene expression dataframe and fill it with time points for later use.
  #     df.ge <- data.frame(gene.expression)
  #
  #     if ( f.type != "individual" ) {
  #       df.ge <- cbind.data.frame(df.ge, time=c(0))
  #       df.ge[seq(1,length(gene.expression),4),2] <- "t0"
  #       df.ge[seq(2,length(gene.expression),4),2] <- "t10"
  #       df.ge[seq(3,length(gene.expression),4),2] <- "t30"
  #       df.ge[seq(4,length(gene.expression),4),2] <- "t180"
  #     } else {
  #       df.ge <- cbind.data.frame(df.ge, time=eqtl$t.interval)
  #     }
  #
  #     df.ge.melt <- melt(df.ge)
  #
  #     # Set sample names.
  #     df.ge.melt[,2] <- names(gene.expression)
  #
  #     # Bind the dataframes together.
  #     df.melt <- cbind.data.frame(timepoints=df.ge.melt[,1],
  #                                 expression=df.ge.melt[,3],
  #                                 sample=df.ge.melt[,2],
  #                                 genotypes=genotypes.melt)
  #
  #     # Create an order for the timepoints.
  #     timepoints.order <- c("t0","t10","t30","t180")
  #
  #     # Order the factor levels for the timepoints.
  #     df.melt$timepoints <- factor(df.melt$timepoints, levels=timepoints.order)
  #
  #     # Create an order for the genotypes.
  #     genotype.levels <- levels(df.melt$genotypes.value)
  #     major.allele = strsplit(as.character(eqtl$minor_major),"_")[[1]][2]
  #     genotype.order <- c()
  #     if (length(genotype.levels) == 2) { # Only 2 levels for genotypes.
  #
  #       one <- strsplit(genotype.levels,"")[[1]] # first level
  #       two <- strsplit(genotype.levels,"")[[2]] # second level
  #
  #       if (major.allele == one[1] && major.allele == one[2]) { # first level is homozygous major allele.
  #         # print("one is major")
  #         major = genotype.levels[1]
  #       } else if (major.allele == one[1] && major.allele != one[2] | major.allele != one[1] && major.allele == one[2]) { # first level is heteryzygous
  #         # print("one is hetero")
  #         heter = genotype.levels[1]
  #         genotype.order <- c(genotype.order, one)
  #       } else {
  #         print("WARNING: Third genotype should not exists here.")
  #       }
  #
  #       if (major.allele == two[1] && major.allele == two[2]) { # first level is homozygous major allele.
  #         # print("two is major")
  #         major = genotype.levels[2]
  #       } else if (major.allele == two[1] && major.allele != two[2] | major.allele != two[1] && major.allele == two[2]) { # first level is heteryzygous
  #         # print("two is heter")
  #         heter = genotype.levels[2]
  #       } else {
  #         print("WARNING: Third genotype should not exists here.")
  #       }
  #
  #       genotype.order <- c(major, heter)
  #
  #     } else { # 3 levels for genotypes.
  #       one <- strsplit(genotype.levels,"")[[1]] # first level
  #       two <- strsplit(genotype.levels,"")[[2]] # second level
  #       three <- strsplit(genotype.levels,"")[[3]] #third level
  #       if (major.allele == one[1] && major.allele == one[2]) { # first level is major allele.
  #         # print("one is major")
  #         major = genotype.levels[1]
  #       } else if (major.allele == one[1] && major.allele != one[2] | major.allele != one[1] && major.allele == one[2]) { # first level is heteryzygous
  #         # print("one is hetero")
  #         heter = genotype.levels[1]
  #         genotype.order <- c(genotype.order, one)
  #       } else {
  #         # print("one is minor")
  #         minor = genotype.levels[1]
  #       }
  #
  #       if (major.allele == two[1] && major.allele == two[2]) { # first level is homozygous major allele.
  #         # print("two is major")
  #         major = genotype.levels[2]
  #       } else if (major.allele == two[1] && major.allele != two[2] | major.allele != two[1] && major.allele == two[2]) { # first level is heteryzygous
  #         # print("two is heter")
  #         heter = genotype.levels[2]
  #       } else {
  #         minor = genotype.levels[2]
  #       }
  #
  #       if (major.allele == three[1] && major.allele == three[2]) { # first level is homozygous major allele.
  #         # print("three is major")
  #         major = genotype.levels[3]
  #       } else if (major.allele == three[1] && major.allele != three[2] | major.allele != three[1] && major.allele == three[2]) { # first level is heteryzygous
  #         # print("three is heter")
  #         heter = genotype.levels[3]
  #       } else {
  #         minor = genotype.levels[3]
  #       }
  #
  #       genotype.order <- c(major, heter, minor)
  #     }
  #
  #     # Order the factor levels for the genotypes.
  #     df.melt$genotypes.value <- factor(df.melt$genotypes.value, levels=genotype.order)
  #
  #     # Order data on factor levels.
  #     df.melt <- df.melt[order(df.melt$timepoints),]
  #     df.melt <- df.melt[order(df.melt$sample),]
  #
  #     # Save ggplot in variable and plot.
  #     mi <- min( df.melt$expression ) - 1
  #     ma <- max(df.melt$expression ) + 1
  #     p <- ggplot(data=df.melt, aes(x=genotypes.value, y=expression, group=genotypes.value) ) +
  #       geom_boxplot(aes( fill=genotypes.value), outlier.shape=NA ) +
  #       geom_point( position=position_jitter(width=0.15),colour = "darkgrey") +
  #       coord_cartesian( ylim = c( mi,ma ) ) +
  #       ggtitle( paste(eqtl$genenames, " - ", eqtl$snps, sep = ""),
  #                subtitle = paste( "P-value: ", eqtl$pvalue ) ) +
  #       theme( plot.title = element_text( size = rel(1.6), hjust = 0.5 ),
  #              plot.subtitle = element_text(size = rel(1), hjust = 0.5 ) ) +
  #       xlab(paste("Genotypes",sep="")) + ylab("Norm. read count")
  #     if ( f.type != "individual" ) {
  #       p + scale_fill_discrete( name="Genotypes",
  #                                labels=paste( names( table( df.melt$genotypes.value ) ),"(", table( df.melt$genotypes.value )/4, ")", sep ="") ) +
  #         facet_wrap( ~ timepoints, scales="free")
  #     } else {
  #       p + scale_fill_discrete( name="Genotypes",
  #                                labels=paste( names( table( df.melt$genotypes.value ) ),"(", table( df.melt$genotypes.value ), ")", sep ="") ) +
  #         facet_wrap( ~ timepoints, scales="free")
  #     }
  #
  #     ggsave( filename=paste( image.base, f.type, "/", f.acting, "/", f.cov, "/", eqtl$origin, "_eQTL_", eqtl$t.interval, "_", eqtl$genenames, "_", eqtl$snps,
  #                             "_", eqtl$minor_major, "_", f.cov, ".pdf", sep=""), plot=last_plot(), device = "pdf")
  #   })
  # })

}