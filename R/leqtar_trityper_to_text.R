# Leqtar
# Erik Schutte
# https://github.com/ErikSchutte/leqtar/issues
# Function    : leqtar_trityper_to_text
# Description : Converts Trityper data to dosage text fiels.

#' leqtar_trityper_to_text
#'
#' Conversst Trityper output to a human readable dosage matrix.
#'
#' @note Only works on Linux/MacOSX
#' @param xmx Java Xmx size in gb.
#' @param xms Java Xms size in gb.
#' @param cohort The name of the cohort for conversion.
#' @param snp_location A path to file that contains SNPs that will be included in the conversion. Usually this is in path/to/trityper/chr_i/SNPs.txt.
#' @param trityper_location A path to the trityper files.
#' @param sample_location A path to a file containing included samples. Useally this is in path/to/trityper/chr_i/Individuals.txt.
#' @param output_location A path to an output directory.
leqtar_trityper_to_text <- function(xmx = 30 , xms = 30,
                                    cohort, snp_location,
                                    trityper_location, sample_location,
                                    output_location) {

  # Parameters rJava.
  parameters <- paste("-Xmx", xmx, "g -Xms", xms, "g", sep = "")

  # Init rJava.
  rJava::.jinit()
  .jcall( file.path("java", "lang", "System", fsep = .Platform$file.sep ),
          "S", "getProperty", "java.runtime.version" )
  dir_Genotype_java <- paste0( "Genotype-IO-1.0.6-SNAPSHOT-jar-with-dependencies.jar" )
  .jinit( classpath = dir_Genotype_java, parameters = parameters )

  # Does output location exist? If not than create all files in requested output path.
  ifelse( !dir.exists( output_location ), dir.create( output_location, recursive = T ),
          message( "[INFO] ", output_location, " exists.." ) )

  # Does trityper_location exists?
  ifelse( !dir.exists( trityper_location ), { message( "[ERR] ", trityper_location, " does not exist!" ); q( save='no' ) },
          message( "[INFO] ", trityper_location, " exists.." ) )

  # Does sample file exist?
  ifelse( !file.exists( sample_location ), { message( "[ERR] ", sample_location, " does not exist!" ); q( save='no' ) },
          message( "[INFO] ", sample_location, " exists.." ) )

  # Does included_snps exist?
  ifelse( !file.exists( snp_location ), { message( "[ERR] ", snp_location, " does not exist!" ); q( save='no' ) },
          message( "[INFO] ", snp_location, " exists.." ) )

  # Create sampleIndex
  sampleIndex <- try( scan( sample_location, "character" ) )

  # Get samples
  samples <- scan( file = sample_location, what = "character" )

  # included snps
  genotype <- read.csv( snp_location, header = T, stringsAsFactors = F,
                        check.names = F )
  genotype <- genotype[,1]
  genotype_total <- length( genotype )
  genotype_unique <- unique( genotype )
  genotype_unique_total <- length( genotype_unique )

  # For the cohort.
  base <- trityper_location



}
#' loadGenotypeData
#'
#' Loads genotype data using GenotypeHarmonizer.
#'
#' @param basePath basePath, i.e. often working directory of script or executable.
#' @param dataType TRITYPER, VCFFOLDER. Conversion for GenotypeHarmonizer | Genotype-IO
#' @param cacheSize Cache size.
#' @return returns loaded genotype data.
loadGenotypeData <- function( basePath, dataType, cacheSize=1000,
                              variantFilter = .jnull( class = file.path( "org", "molgenis", "genotype", "variantFilter", "VariantFilter", fsep = .Platform$file.sep ) ),
                              sampleFilter = j.null( class = file.path( "org", "molgenis", "genotype", "sampleFilter", "SampleFilter", fsep = .Platform$file.sep ) ) ) {
  dataType <- toupper(dataType)
  genotypeDataFormat <- .jcall ( file.path( "org", "molgenis", "genotype", "RandomAccessGenotypeDataReaderFormats", fsep = .Platform$file.sep),
                                 file.path( "Lorg", "molgenis", "genotype", "RandomAccessGenotypeDataReaderFormats;", fsep = .Platform$file.sep),
                                 "valueOf", dataType)
  return( .jcall( genotypeDataFormat, file.path( "Lorg", "molgenis", "genotype", "RandomAccessGenotypeData;", fsep = .Platform$file.sep),
                  "createFilteredGenotypeData", basePath, as.integer(cacheSize), variantFilter, sampleFilter) )


}

#' complementSNP
#'
#' Calculates the complement for input allele.
#'
#' @param allele input allele.
#' @return returns complement from allele.
complementSNP <- function(allele){
  if(nchar(allele) == 1){
    if(allele == "A"){
      return("T")
    }else if(allele == "T"){
      return("A")
    } else if(allele == "G"){
      return("C")
    } else if(allele == "C"){
      return("G")
    } else{
      return("not Avail")
    }
  }
}

#' inverseDosage
#'
#' inverts dosage according to vector.
#'
#' @param dosageVector vector of dosage values per rsID.
#' @return inverse dosage vector.
inverseDosage <- function(dosageVector){
  return(((dosageVector -2) * -1))
}
