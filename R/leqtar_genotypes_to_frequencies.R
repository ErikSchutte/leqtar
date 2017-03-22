#' leqtar_genotypes_to_frequencies
#'
#' Turns genotypes into frequencies to Matrix eQTL is able to perform linear regression analysis.
#'
#' @param genotypeData the genotype file content
#' @return genotypeData as numeric matrix.
leqtar_genotypes_to_frequencies <- function(snps, VERBOSE = F) {

  # Transpose data.frame
  snps <- t(snps)

  # Add a row to the data frame with MAF where the alleles will be stored.
  tmp.df <- t(data.frame(MAF=rep("",length(colnames(snps)))))
  snps.genotype <- snps # Copy of snps that won't be changed to frequencies.
  snps.genotype <- rbind(snps.genotype, tmp.df) # after re-arangement bind alleles to genotypes.

  ## Convert genotypes to 0, 1 and 2. 0 ref/ref, 1 ref/mut and 2 mut/mut.
  # Keeps track of the snps that are causing noise in the data and are not representable.
  snps.discarded <- c()
  snps.discarded.pos <- c()
  snps.discarded.counter <- 0
  # For every column in the genotype data.
  for (i in 1:ncol(snps)) {

    # Verbose
    if ( VERBOSE == TRUE ) {
      # Verbose - Print current column and column name.
      cat("#",i,"- current column: ",colnames(snps)[i],"\n")
    }
    # Save genotype and genotype count for each column.
    alleles.all <- table(snps[,i])

    # A minimum of 3 samples for a genotype is required to effectively calculate correlation.
    # Therefore we are removing the rows where this is not the case.
    if ( length(alleles.all) < 3 & min(alleles.all) < 2) {

      # Set an index for the SNP, so that it can be removed later.
      snps.discarded.pos <- c(snps.discarded.pos,i)

      # Store the discarded SNP for later review.
      snps.discarded <- c(snps.discarded, colnames(snps)[i])
      snps.discarded.counter <- snps.discarded.counter + 1
    }

    # List with all the major, minor, and heterozygote genotypes and their respecitve counts.
    allele.list <- list(major <- list(genotype = NULL, count = 0),
                        minor <- list(genotype = NULL, count = 0),
                        hetero <- list(genotype =c(), count = c()),
                        other <- list(genotype = NULL, count = 0))

    # Temp value for checking which genotype occurs the most.
    highest_count_homozygote <- 0

    # For every allele that the genotype is counted for.
    for (a in 1:length(alleles.all)) {

      # Extract the genotype.
      geno <- names(alleles.all[a])

      # Verbose
      if ( VERBOSE == TRUE ) {
        # Verbose - Print current allele.
        cat("Current allele = ", geno,'\n')
      }

      # Extract the number of occurrences.
      count <- alleles.all[[a]]

      # Split the genotype
      alleles <- strsplit(geno, "")

      #Check if the first allele is the same as the second
      # thus checking for heterozygote or homozygote genotype.
      if (identical(alleles[[1]][1],alleles[[1]][2]) && !is.na(geno)) {

        # Verbose
        if ( VERBOSE == TRUE ) {
          # Verbose - Print if alleles are identical, and thus are homozygous.
          cat('Alleles are identical, sorted as homozygote\n')
        }

        # If the highest_count_homozygote == 0 (the first iteration).
        if (highest_count_homozygote == 0) {

          # Verbose
          if ( VERBOSE == TRUE ) {
            # Verbose - Print which genotype and count is added as major allele.
            cat("Adding ",geno," with count ", count, "as major allele\n")
          }

          # For the first iteration, the first genotype's occurences are the most occured
          # genotype and will be added on the first position.
          highest_count_homozygote <- count
          allele.list[[1]]$genotype <- geno
          allele.list[[1]]$count <- count

          # If the current genotype has more occurences than the lastly noted genotype,
          # it is automatically the major allele on the first position.
        } else if (highest_count_homozygote < count) {

          # Verbose
          if ( VERBOSE == TRUE ) {
            # Verbose - Print which genotype and count is added as new major allele,
            # and which genotype and count previously occupied that spot.
            cat("Adding ",geno," with count ", count, "as major allele\nPreviously ",
                allele.list[[1]]$genotype," with count ", allele.list[[1]]$count,'\n')
          }

          # Altering the position from the previous current major allele to minor allele
          # in the allele.list.
          allele.list[[2]]$genotype <- allele.list[[1]]$genotype
          allele.list[[2]]$count <- allele.list[[1]]$count

          # Setting the new major allele on the first position.
          highest_count_homozygote <- count
          allele.list[[1]]$genotype <- geno
          allele.list[[1]]$count <- count

        } else {
          # If the highest occuring genotype already found has more occurences,
          # the resulting genotype will be the minor allele which is the 2nd position
          # in the list.
          allele.list[[2]]$genotype <- geno
          allele.list[[2]]$count <- count

          # Verbose
          if ( VERBOSE == TRUE ) {
            # Verbose - Print which genotype and count is added as minor allele.
            cat("Adding ",geno," with count ", count, "as minor allele\n")
          }
        }

      } else {
        if (identical(geno,NA)) {
          # Separating the '00'/NA's from the data, these will be represented as 9.
          allele.list[[4]]$genotype <- geno
          allele.list[[4]]$count <- count

        } else {
          # If the genotype is not '00' it is automatically a heterozygote.
          allele.list[[3]]$genotype <- c(allele.list[[3]]$genotype, geno)
          allele.list[[3]]$count <- c(allele.list[[3]]$count, count)

          # Verbose
          if ( VERBOSE == TRUE ) {
            # Verbose - Print if genes are not identical, and thus are heterozygote.
            cat('alleles are not identical, heterozygote\n')
            # Verbose - Print which genotype and count is added as heterozygous.
            cat("Adding ",geno," with count ", count, "as heterozygote allele\n")
          }
        }
      }
    }

    # For every row in the genotype data.
    for (j in 1:nrow(snps.genotype) ) {

      # First we want to create a dict that takes the inverse of the MAJOR alleles we put there,
      # since we do not save all the minor alleles.
      #list()

      # If J is smaller than the maximum lengt of rows, e.g. all the samples minus the last row
      # where the MAF will be stored.
      if (j < (nrow(snps.genotype) ) ) {
        if (VERBOSE == T) {
          cat("Current row is: ",j,"\n")
        }

        # Check the current snp agains the determined major, minor and heterozygote genotypes,
        # and determine the allele frequency.
        # print(allele.list)
        snp.freq <- change_genotype_to_frequency(snps[j,i],allele.list,VERBOSE)

        # Verbose.
        if ( VERBOSE == TRUE ) {
          # Verbose - Print the current snip' genotype and the corrosponding frequency.
          cat("Changing current snp: ",snps[j,i]," to: ", snp.freq,"\n")
        }

        # Save the SNP frequency at previous occupied genotype location.
        snps[j,i] <- as.numeric(snp.freq)
      }
      # J is equal to the length of rows for the data frame snps. meaning it is hte last row that will contain
      # the MAF.
      else {
        if (VERBOSE == T) {
          cat("Current row is: ",j,"\n")
        }
        # We don't always have the minor allele, so instead we look at the major and heterozygous genotypes.
        # E.g. all 22 samples have the following genotyeps for snp X, TT (19x) CT (3x) CC(0x).
        # Because we can't count CC we have no genotype CC and thus have to 'create' it.
        major.allele = strsplit(allele.list[[1]]$genotype,"")[[1]] # Contains both alleles, e.g. T and T.

        if ( length(allele.list[[3]]$genotype) < 1 && length(allele.list[[2]]$genotype) < 1 ) {

          # In the off case that there are TT(22x) CT(0x) CC(0x)
          minor.allele <- "?"
        } else if ( length(allele.list[[3]]$genotype) < 1 && length(allele.list[[2]]$genotype) > 0 ) {

          # In the off case that there are TT(17x) CT(0x) CC(5x)
          minor <- allele.list[[2]]$genotype
          minor.alleles <- strsplit(minor,"")[[1]]
          minor.allele <- minor.alleles[1]

        } else {

          heter.allele = strsplit(allele.list[[3]]$genotype,"")[[1]] # Contains both alleles, e.g. C and T.

          if ( heter.allele[1] == major.allele[1] ) {

            # If the first allele, in our example C matches the major allele T we know that he minor allele must be
            # the second allele of the heterozygous alleles.
            minor.allele <- heter.allele[2]

          } else if ( heter.allele[2] == major.allele[1] ) {

            # If the second allele, in our exmaple T matches the major allele T we know that the minior allele must be
            # the first allele of the heterozgyous alleles.
            minor.allele <- heter.allele[1]
          }
        }

        # Save the genotypes in the data frame as MINOR / MAJOR.
        alleles <- paste(minor.allele, "_", major.allele[1], sep="")
        # print(alleles)
        snps.genotype[j,i] <- alleles
      }
    }
  }

  # Remove discarded snps from data set.
  snps <- snps[-snps.discarded.pos,]
  snps.freq <- snps.freq[,-snps.discarded.pos]
  # Convert genotype matrix to numeric values for Matrix eQTL analaysis.
  suppressMessages(
    class(snps) <- "numeric"
  )

  # Transfer the genotype matrix so that the columns 'align' with the gene expression matrix.
  snps.t <- t(snps)
  snps.genotype.t <- t(snps.genotype)
}

# Change_allele_to_frequencies --------------
#' change_allele_to_requencies
#'
#' Function to identify the genotype and determine wheter it is major, minor or heterozygote.
#' @param x current cell.
#' @param allele.list list containg the current alleles per snp.
#' @param VERBOSE verbosity true/false.
change_genotype_to_frequency <- function (x, allele.list, VERBOSE) {
  if ( identical(x, allele.list[[1]]$genotype) ) {
    # The genotype is major, return 0.
    if (VERBOSE == T) {
      cat(x, "is identical to: ", allele.list[[1]]$genotype,"\n")
    }
    return(0)

  } else if ( identical(allele.list[[4]]$genotype, x) ) {
    # The genotype is NA or 00, return NA.
    if (VERBOSE == T) {
      cat(x, "is identical to: ", allele.list[[4]]$genotype,"\n")
    }
    return('NA')

  } else if ( identical(x, allele.list[[3]]$genotype[1]) | identical(x,allele.list[[3]]$genotype[2]) ) {
    # If the genotype is heterozygote, return 1.
    if (VERBOSE == T) {
      cat(x, "is identical to:", allele.list[[3]]$genotype[1], " or ", allele.list[[3]]$genotype[2],"\n")
    }
    return(1)

  } else if ( identical(x, allele.list[[2]]$genotype) ) {
    # If the genotype is one of the minor alleles, return 2.
    if (VERBOSE == T) {
      cat(x, "is identical to: ", allele.list[[2]]$genotype,"\n")
    }
    return(2)

  }
  else {
    # If the genotype is not major, hetero or minor it is NA.
    if (VERBOSE == T) {
      cat(x, "is identical to: ", NA,"\n")
    }
    return(NA)
  }
}
