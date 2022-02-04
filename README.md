# STIE
Spatial Transcriptome Image and Expression integrative deconvolution enables single cell level spatial analysis

## Description
STIE is a novel computational method tailored for spatial transcriptome data analysis, which integrated both pathology image and gene expression to perform cell type deconvolution, enabling single cell-level spatial gene expression anlayiss.


## Dependencies on packages
-  imagemagick - [imagemagick](http://www.imagemagick.org/script/install-source.php)
-  ImageJ/Fiji - [ImageJ/Fiji](https://imagej.net/software/fiji/downloads)
-  R package: [magick](https://cran.r-project.org/web/packages/magick/vignettes/intro.html); [EBImage](https://bioconductor.org/packages/release/bioc/html/EBImage.html); [Seurat](https://satijalab.org/seurat/articles/install.html); [foreach](https://cran.r-project.org/web/packages/foreach/); [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)

## Installation
### 1) Obtaining a recent version of R
Clear instructions for different version can be found here:
http://cran.fhcrc.org/

### 2) Install the dependent R packages
```
# install R package of Biostrings. 
# try http:// if https:// URLs are not supported
> source("https://bioconductor.org/biocLite.R")
> biocLite("Biostrings")

# install packages for parallel computating
> install.packages(c("foreach","doMC","magick"))

```

### 3) Install the STIE R package
```
git clone https://github.com/zhushijia/STIE.git
R CMD INSTALL -l userFolder STIE
```
Alternatively, use [devtools](https://github.com/hadley/devtools) package
```
> install.packages("devtools")
> library(devtools)
> install_github("zhushijia/STIE")
```

## Tutorial
   See our [wiki](https://github.com/zhushijia/STIE/wiki)

## Citation
 [(link)](asdfads)
