# AUCell
AUCell is an R-package to analyze the state of gene-sets in single-cell RNA-seq data (i.e. identify cells with active gene signatures).







*Note:* This is a version in development. The main version of the package will be submitted Bioconductor. 

To install this BETA version, you can run the following commands from R:
```
# devtools::install_github("aertslab/AUCell") # Does not build the tutorial (faster)
devtools::install_github("aertslab/AUCell", build_vignettes=TRUE)

# You might need to install these packages first:
install.packages(c("devtools", "data.table", "zoo")) # Recommended
source("https://bioconductor.org/biocLite.R"); biocLite("BiocInstaller")
install.packages(c("DT")); biocLite("GEOquery") # For tutorial
```

A tutorial (vignette) is included in the package.
An HTML version of the [tutorial](http://scenic.aertslab.org/tutorials/AUCell.html) is available at http://scenic.aertslab.org.

