print(getwd())
## Libs
library(rJava)
library(data.table)

## Functions
# Load Genotype Data
loadGenotypeData <- function(
                        basePath,
                        dataType,
                        cacheSize=1000,
                        variantFilter = .jnull(class = "org/molgenis/genotype/variantFilter/VariantFilter"),
                        sampleFilter = .jnull(class= "org/molgenis/genotype/sampleFilter/SampleFilter")
){
  dataType <- toupper(dataType)
  genotypeDataFormat <- .jcall("org/molgenis/genotype/RandomAccessGenotypeDataReaderFormats", "Lorg/molgenis/genotype/RandomAccessGenotypeDataReaderFormats;","valueOf", dataType)
  return(.jcall(genotypeDataFormat, "Lorg/molgenis/genotype/RandomAccessGenotypeData;", "createFilteredGenotypeData", basePath, as.integer(cacheSize), variantFilter, sampleFilter))
}

# Complement SNP
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

# Inverse Dosage
inverseDosage <- function(dosageVector){
        return(((dosageVector -2) * -1))
}

## Code
# Init
.jinit()
.jcall("java/lang/System", "S", "getProperty", "java.runtime.version")
configurationTable <- read.csv( "source/config.csv", header = T, stringsAsFactors = F )
.jinit(classpath=as.character(configurationTable[3,2]), parameters=paste("-Xmx", configurationTable[2,2], "g -Xms", configurationTable[1,2] ,"g", sep = "") )

# Get arguments
args <- commandArgs(TRUE)

if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      The R Script
      Arguments:
      --arg1= Name of Cohort 
      --arg2= path to a csv file which contains a column with SNPS names (colname should be SNP)
      --arg3= path to TryTyper directory
      --arg4= path to csv file which contains samples to be included
      --arg4= path to results folder
      --help= print this text
 
      Example:
      Rscript source/getDosages.R COHORT path/to/variants path/to/trityperdir path/to/sample_file path/to/output_directory \n\n")
  q(save="no")
} 

# Set variables
cohortName <- args[1]
variantsFile <- args[2]
trityperDir <- args[3]
samplesFile <- args[4]
outputDir <- args[5]

completeSNPs <- list()
minorAlleleList <- list()
completeSNPInfo <- list()

variants <- read.csv(variantsFile, header = T, stringsAsFactors = F, check.names = F)
variants <- variants[,1]
numberOfVariants <- length(variants)

sampleIndex <- try(scan(samplesFile, "character"))
samples <- scan(file=file.path( trityperDir, "Individuals.txt"), what = "character")

# Loading the filter parameters, this includes the list of SNPs from the previews step. 
variantFilter <-  .jcast(.jnew("org/molgenis/genotype/variantFilter/VariantIdIncludeFilter",
									variants), ### list of SNPs
									"org/molgenis/genotype/variantFilter/VariantFilter")
genotypeData <- loadGenotypeData(basePath= trityperDir, dataType= "TriTyper", variantFilter= variantFilter)

# Get correct positions per SNP
mappings <- fread( file.path( trityperDir, "SNPMappings.txt" ), data.table=FALSE)
mappedSNPs <- mappings[match(variants, mappings[,3]),]
if(sum(is.na(mappedSNPs[,1])) >=1){
	mappedSNPs <- mappedSNPs[-which(is.na(mappedSNPs[,1])),]
}
acceptedNumberOfVariants <- nrow(mappedSNPs)
completeSNPs <- mappedSNPs[,3]
rownames(mappedSNPs) <- mappedSNPs[,3]
completeSNPInfo <- mappedSNPs[,c(3,1,2)]

# Results 
start_dsg <- Sys.time()
sampleDosages <- apply(mappedSNPs[1:1000,], 1, function(snp) {
	
	snp <- t(as.data.frame(snp))
	snp <- .jcall(genotypeData, "Lorg/molgenis/genotype/variant/GeneticVariant;","getSnpVariantByPos",
			as.character(snp[,1]), ## snp chromosome needs to be as character
			as.integer(snp[,2]) ## snp position as integer
		)
	if(.jcall(snp, "Z", "isSnp")) {
		sampleDosage <- round((.jcall(snp, "[F", "getSampleDosages")), digits= 5)
		return(sampleDosage)
	}
})

rownames(sampleDosages) <- samples
sampleDosages <- t(sampleDosages)

# Save dosages.
save( sampleDosages, file = paste(outputDir, "dosages.RData", sep = "" ) )
write.table( sampleDosages, file = paste( outputDir, "dosages.txt", sep = "" ) )

end_dsg <- Sys.time()
start_vrt <- Sys.time()
variantInfo <- apply(mappedSNPs[1:1000,], 1, function(snp) {
	snp <- t(as.data.frame(snp))
	snp.id <- snp[,3]
	snp <- .jcall(genotypeData, "Lorg/molgenis/genotype/variant/GeneticVariant;","getSnpVariantByPos",
                        as.character(snp[,1]), ## snp chromosome needs to be as character
                        as.integer(snp[,2]) ## snp position as integer
                )
	if(.jcall(snp, "Z", "isSnp")) {
		tmpAlleles <- .jcall (snp, "Lorg/molgenis/genotype/Alleles;", "getVariantAlleles")
                minorAllele <- unlist(strsplit(tmpAlleles$toString(), "\\", fixed=TRUE)) [1]
                refAllele <- unlist(strsplit(tmpAlleles$toString(), "\\",fixed=TRUE))[2]
		MAF <- format(.jcall (snp, "D", "getMinorAlleleFrequency"),digits=5)
                HwePvalue <- format(.jcall (snp, "D", "getHwePvalue"), digits=5)
                callRate <- format(.jcall (snp, "D", "getCallRate"), digits=5)
                biAllelic <- .jcall (snp, "Z", "isBiallelic")
		
		
		variants <- data.frame(minorAllele = minorAllele,
				refAllele = refAllele,
				MAF = MAF,
				HwePvalue = HwePvalue,
				callRate = callRate,
				biAllelic = biAllelic)
		rownames(variants) <- snp.id
		colnames(variants) <- c("minorAllele", "refAllele", "MAF", "HwePvalue", "callRate", "biAllelic")
		return( variants )
        }
})
end_vrt <- Sys.time()

message("Run Times:",
	"\nDosages: ",as.character(end_dsg - start_dsg),
	"\nVariant: ", as.character(end_vrt - start_vrt))
	
allVariantInfo <- do.call(rbind, variantInfo)

save( allVariantInfo, file = paste(outputDir, "allVariantInfo.RData", sep = "" ) )
write.table( allVariantInfo, file = paste( outputDir, "allVariantInfo.txt", sep = "" ) )
