
# DeadOrAlive

Analysis of High-throughput Colony Forming Unit Assays

## Installation

``` r
install.packages("adimpro")
package <- "https://cran.r-project.org/src/contrib/Archive/PET/PET_0.5.1.tar.gz"
fileLocation <- tempfile()
download.file(package, fileLocation)
install.packages(fileLocation, type="source", repos=NULL)
install.packages("BiocManager")
BiocManager::install("EBImage")
install.packages(c("jpeg", "tiff", "logging", "ggplot2"))
package <- "https://cran.r-project.org/src/contrib/Archive/gitter/gitter_1.1.1.tar.gz"
fileLocation <- tempfile()
download.file(package, fileLocation)
install.packages(fileLocation, type="source", repos=NULL)
install.packages(c("zoo", "magick", "gplots", "RColorBrewer", "rmarkdown", "toOrdinal", "cobs", "rootSolve", "ggplot2"))
install.packages("devtools")
devtools::install_github("JohnTownsend92/DeadOrAlive", build_vignettes=TRUE)
```

## Tutorial

``` r
vignette("DeadOrAliveVignette", package="DeadOrAlive")
```
