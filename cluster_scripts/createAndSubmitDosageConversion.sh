#!/bin/bash

## PATHS
JOB_DIR=/groups/umcg-wijmenga/tmp02/projects/umcg-eschutte/jobs/

for i in {1..22};
do
	(printf "#!/bin/bash
#SBATCH --job-name=Forester_TrityperToDosage_$i.run
#SBATCH --output="$JOB_DIR"Forester_TrityperToDosage_$i.out
#SBATCH --error="$JOB_DIR"Forester_TrityperToDosage_$i.err
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=80gb
#SBATCH --nodes=1
#SBATCH --qos=regular
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=30L

## Modules
module load R

## VARIABLES
COHORT=Foresters

## PATHS
TRITYPER_DIR=/groups/umcg-wijmenga/tmp02/projects/umcg-eschutte/Foresters_Tritryper_030417/
TRITYPER_OUT='/groups/umcg-wijmenga/tmp02/projects/umcg-eschutte/'\$COHORT'_Dosages/'

echo -e 'Started Trityper to dosage conversion'

echo -e 'Check output directories'
if [[ -d \$TRITYPER_OUT ]];
then
	echo 'Output directory: '\$TRITYPER_OUT' exists'
else
	mkdir -p $TRITYPER_OUT
        echo 'Output directory: '\$TRITYPER_OUT' created'
fi

# Define chromosome
CHR='chr_'$i
echo -e 'Processing '\$CHR

# Directories for current chr
INPUT_DIR=\$TRITYPER_DIR\$CHR'/'
OUTPUT_DIR=\$TRITYPER_OUT\$CHR'/'

# If output directory exists, remove and remake.
if [[ -d \$OUTPUT_DIR ]];
then
	rm -r \$OUTPUT_DIR
	mkdir \$OUTPUT_DIR
else
	mkdir \$OUTPUT_DIR
fi
	
# Pre-define variables
variants=\$INPUT_DIR'SNPs.txt'
uniqueVariants=\$INPUT_DIR'uniqueVariants.csv'
duplicateVariants=\$INPUT_DIR'duplicateVariants.csv'
samples=\$INPUT_DIR'Individuals.txt'
	
# Make unique variant file
echo -e 'Creating unique variants'
echo -e 'SNP' > \$uniqueVariants
sort \$variants | uniq -u >> \$uniqueVariants
sort \$variants | uniq -d >> \$duplicateVariants
	
# Call R script for Trityper-Dosage conversion
echo -e 'Convert Trityper-Dosage'
Rscript source/getDosages.R \$COHORT \$uniqueVariants \$INPUT_DIR \$samples \$OUTPUT_DIR
") > $JOB_DIR'chr_'$i.sh
sbatch $JOB_DIR'chr_'$i.sh
done
	
	
