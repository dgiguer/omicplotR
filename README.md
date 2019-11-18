# omicplotR
An R package to visualize high-throughput sequencing data. Click on the wiki for more!  
[![DOI](https://zenodo.org/badge/101769044.svg)](https://zenodo.org/badge/latestdoi/101769044)

Read more about omicplotR here: 

> [__Giguere, DJ, Macklaim, JM, Lieng, BY, Gloor, GB.__ omicplotR: visualizing omic datasets as compositions. BMC Bioinformatics 20, 580 (2019)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3174-x)

# Table of content

* [Introduction](#introduction)
* [Installation](#installation)
* [Quick usage](#quick-usage)
* [Submit bug report](#submit-bug-report)

## Introduction

As input, omicplotR takes the following: 
* A data file and optionally a metadata file
* The data file is a count table of sequencing reads with samples as columns and features by rows
* An additional column of taxonomic identifiers may be included. Column name must be 'taxonomy', identifiers must contain at least 4 levels separated by a semi-colon (Bacteria;Firmicutes;Bacilli;Lactobacillales)

Here is an example of how your data file should look:
<p align="center"><img src="https://raw.githubusercontent.com/wiki/dgiguer/omicplotR/www/example_data.png" alt="Data" width="600"></p>

* The metadata file must contain samples by rows and identifiers by columns

Here is an example of how your metadata should look like:

<p align="center"><img src="https://raw.githubusercontent.com/wiki/dgiguer/omicplotR/www/example_metadata.png" alt="Metadataata" width="600"></p>

Reasons to use omicplotR: 
* It allows you to explore your raw sequencing count data without coding in R
* It performs interactive filtering
* It allows you to easily colour PCA biplots according to metadata
* It generates interactive effect plots
* It allows easy download of publically available data from the EBI MGnify portal

Reasons to **not** use omicplotR: 
* You have thousands of samples and tens of thousands of features (it will be quicker using the command line)

## Installation

The first thing to check is if you have the current version of [Bioconductor](http://bioconductor.org). Look at the `BiocManager` [vignette](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html) for instructions on installing or upgrading to 3.10.

You can check the version of Bioconductor that is install if you already have `BiocManager` installed.

```
BiocManager::version()
```

Once you have ensured that you have the correct version of Bioconductor, install `omicplotR`: 

```
BiocManager::install("omicplotR")
```

The development version can be installed directly from Github: 

```
install.packages("devtools")
devtools::install_github("dgiguer/omicplotR")
```

## Quick usage

The quickest way to run `omicplotR` is:

```
library(omicplotR)
omicplotr.run()
```

This will pop up a window in your default browser. 

To view the tutorials, visit [the wiki](https://github.com/dgiguer/omicplotR/wiki).

## Submit bug report

If omicplotR is not working as expected, please submit an issue that includes the following information. This will help me help you!

```
Please make sure to fill in this template when submitting an issue. Thank you for taking the time to submit an issue!
- [x] Please check off boxes below like this.

What is the error message that appears in the R console?

- [ ] Check that you are using the newest release of R.
- [ ] Check that you are using the newest release of omicplotR.

What OS are you on? 

- [ ] Windows
- [ ] macOS
- [ ] Linux

What version of your OS are you on?

What is the expected behaviour? 

What is the actual behaviour? 

What does the error look like? (submit screenshot of browser if possible)
```

