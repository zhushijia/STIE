# STIE
Spatial Transcriptome Image and Expression integration enables single-cell level spatial transcriptomics data analysis

## Description
STIE is a novel computational method tailored for spatial transcriptomics data analysis, which integrated spot level gene expression, nuclear segmentation, and nuclear morphology to perform cell type deconvolution/convolution and clustering, therefore enabling the single-cell level spatial transcriptomics anlayiss.
![Figure1](https://user-images.githubusercontent.com/5418417/182406531-3f623ed0-41ad-484c-9c77-f1707d2fc34c.jpg)


## Dependencies on packages
STIE has been tested on GNU/Linux but should run on all major operating systems. STIE depends on the following packages:
-  imagemagick - [imagemagick (>=7.1.0)](http://www.imagemagick.org/script/install-source.php)
-  ImageJ/Fiji - [ImageJ/Fiji](https://imagej.net/software/fiji/downloads)
-  R package: [magick](https://cran.r-project.org/web/packages/magick/vignettes/intro.html); [EBImage](https://bioconductor.org/packages/release/bioc/html/EBImage.html); [Seurat](https://satijalab.org/seurat/articles/install.html); [CellChat](https://github.com/sqjin/CellChat); [quadprog](https://cran.r-project.org/web/packages/quadprog/index.html); 

## Installation
### 1) Obtain R (>=3.6)
Clear instructions for different version can be found here:
http://cran.fhcrc.org/

### 2) Install the dependent R packages
```
# install R packages of computing
> install.packages(c("quadprog"))

# install magick
> install.packages("magick")

# install EBImage
> if (!require("BiocManager", quietly = TRUE))
>     install.packages("BiocManager")
> BiocManager::install("EBImage")

# install CellChat
> BiocManager::install("ComplexHeatmap")
> devtools::install_github("sqjin/CellChat")

# install Seurat
> install.packages('Seurat')
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
