# STIE
Spatial Transcriptome Image and Expression integration enables single-cell level spatial transcriptomics data analysis

## Description
STIE is a novel computational method tailored for spatial transcriptome data analysis, which integrated both pathology image and gene expression to perform cell type deconvolution, enabling single cell-level spatial gene expression anlayiss.
![Figure1_STIE_flowchart2](https://user-images.githubusercontent.com/5418417/168800169-94375fdf-9e42-40b5-976e-b707365cf4ca.jpg)


## Dependencies on packages
-  imagemagick - [imagemagick (>=7.1.0)](http://www.imagemagick.org/script/install-source.php)
-  ImageJ/Fiji - [ImageJ/Fiji](https://imagej.net/software/fiji/downloads)
-  R package: [magick](https://cran.r-project.org/web/packages/magick/vignettes/intro.html); [EBImage](https://bioconductor.org/packages/release/bioc/html/EBImage.html); [Seurat](https://satijalab.org/seurat/articles/install.html); [CellChat](https://github.com/sqjin/CellChat); [quadprog](https://cran.r-project.org/web/packages/quadprog/index.html); 

## Installation
### 1) Obtain R (>=3.6)
Clear instructions for different version can be found here:
http://cran.fhcrc.org/

### 2) Install the dependent R packages
```
# install R packages of computating
> install.packages(c("foreach","doParallel","quadprog"))

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
devtools::install_github("zhushijia/STIE""
                         ,ref="main"
                         ,auth_token = "My script"
                         )
```


## Tutorial
   See our [wiki](https://github.com/zhushijia/STIE/wiki)

## Citation
 [(link)](asdfads)
