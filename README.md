# AUCell
AUCell is an R-package to analyze the state of gene-sets in single-cell RNA-seq data (i.e. identify cells with active gene signatures).

** Availability**: 

- The newest stable version of AUCell is available in *[Bioconductor](https://bioconductor.org/packages/AUCell)*. 
Unless it is for using with SCENIC, we recommend to install AUCell from Bioconductor.
The package also contains a **tutorial** [(vignette)](https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html) with information on how to run AUCell and how to interprete its results.

- **SCENIC** (v 0.1.5) requires [AUCell 0.99.5](http://scenic.aertslab.org/downloads/Rpackages/AUCell_0.99.5.tar.gz). More details on how to install this specific version are provided in SCENIC tutorial. 
*(Using newer versions might result in errors related to missing functions/classes, as some of these have been renamed).*

- This Github repository is meant mostly for development. "Use at your own risk" :-). Developers are welcome to contribute. 

You might also need to install these depencencies first:
```
source("https://bioconductor.org/biocLite.R"); biocLite("BiocInstaller")
biocLite(c("data.table", "zoo", "doMC", "doRNG", "mixtools", "GEOquery", 
"SummarizedExperiment", "DT", "plotly", "NMF", "d3heatmap"))
```
