# leqtar
Linear eQTL analysis in R

Please note that this is a work in progress. This readme is not yet fully complete but will be updated along the way. If there are any
bugs or problems that you might spot, please report them in the issue tracker.

## 1. Install

In order to install leqtar you need `devtools`. If not already installed go ahead and enter `install.packages("devtools")`
in your RStudio terminal or in R using the terminal, else load using `library(devtools)`.

After `devtools` is installed and loaded using the library function enter `install_github("ErikSchutte/leqtar")`. This defaults to the master branch. If you want the latest update from the dev branch you have to specify `install_github("ErikSchutte/leqtar", ref = "dev")`.
For more information about `install_github` please use `?install_github` in the terminal.

Once the package is installed you can load it using `library(leqtar)`.

## 2. Data

I expect the user to be ignorant of the required steps in acquiring genotype and phenotype data. Supposidly you know what you are doing,
the phenotype data can be derived from aligning Fastq files with HISAT2 or STAR. Resulting BAM files can be counted with HTSeq or Feature counts.
Raw phenotype data is now availble, but do note that you probably have to normalize this data first. You can either just do a simple log transformation or a Variance Stabilizing Transformation (VST).
I would refer to `?DESeq2` for VST and the log transformation is in the basic R packages.

Genotype data can be anything from dosage files to literal genotypes i.e. 'AA/AT/TT'.

## 3. Examples

Leqtar is designed to simplify association analysis, pre-process data, post-processing and visualization of data all in one go.
Let's define some usage scenario's before we get started, the Leqtar library provides a basic non-sense data set. Currently leqtar
has been tested using assocations between cytokines and genotypes and gene expression and genotypes, being cytokine QTLs and expression QTLs respectively.

When you load leqtar with `library(leqtar)` three data sets are immediatly availble for use.
These are `genotype_test_data.RData`, `expression_test_data.RData` and `covariate_test_data.RData`.
Because these are included in the package they are automatically loaded. But don't worry! You don't need to transform any of your data.
Leqtar accepts a multitude of input. Data objects, RData files, .csv files, .tsv files, .txt files and .xls files are all supported!

### 3.1 Examples - Basic analysis

We want to know the effect of a genotype on a specific phenotype. We would also like to see visualize the effect of each significant assocation.
To do this we simply just call the following little script.

```
library(leqtar)

leqtar(genotypeFile = genotype_test_data, expressionFile = expression_test_data, covariateFile = covariate_test_data,
       run_name = "test01")
```

It's as easy as that. Please do note that the genotypeFile, expressionFile and run_name are required parameters. For more info see `?leqtar`.

When the option `output_dir` is not given, leqtar will automatically store the output in your current working directory. To see what your current working directory is in R,
you can enter `getwd()` in the terminal or in RStudio's terminal.

Let's assume your current working directory is your home folder, i.e. `"~/"` for Linux and OSX users.
Leqtar should create a directory structure starting with `leqtar_out` as the root of all output. 
The full path would look something like `/Users/someone/leqtar_out`. Below `leqtar_out` a sub directory is created for your run, called `test01`.
In the `leqtar_out/test01/data/` folder an .RData file is stored with the output of linear regression analysis which you can ofcourse always review with
`load( file.path( "path/to/leqtar_out/test01/data/test01.Rdata" ) )`.

Currently there are two types of visualizations that are done. In `leqtar_out/test01/images/` are two folders, genotypes and manhattan.
`/genotypes/` will contain genotype boxplots for each significant association found using linear regression analysis. `/manhattan/` shows the significance of
your SNPs according to the associations found by linear regression analsyis.

Last but not least there is a folder called `tables` below your run folder. This folder contains a variety of tables containing important information
that is picked up during visualization. The manhattan plot for instance would show genomewide signficant SNPs, and they are reported in a table here. You get the idea.






