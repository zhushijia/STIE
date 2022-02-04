# STIE
Spatial Transcriptome Image and Expression integrative deconvolution enables single cell level spatial analysis

## Description
STIE is a novel hybrid sequencing-based alternative splicing identification method, which is specially tailored for long read RNAseq correction. It was used in the paper [(E Flaherty*, S Zhu*, et. al., Neuronal impact of patient-specific aberrant NRXN1α splicing, Nature Genetics, 2019)](https://www.nature.com/articles/s41588-019-0539-z). 


## Dependencies on packages
-  STARlong & STARshort - [STAR](https://github.com/alexdobin/STAR)
-  samtools - [samtools](http://samtools.sourceforge.net/)
-  bedtools - [bedtools](http://bedtools.readthedocs.io/en/latest/)
-  R package: [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html); [foreach](https://cran.r-project.org/web/packages/foreach/); [doMC](https://cran.r-project.org/web/packages/doMC/); [reshape](https://cran.r-project.org/web/packages/reshape/index.html)


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
> install.packages(c("foreach","doMC","reshape"))

```

### 3) Install the STAR2bSMRT R package
```
git clone https://github.com/zhushijia/STAR2bSMRT.git
R CMD INSTALL -l userFolder STAR2bSMRT
```
Alternatively, use [devtools](https://github.com/hadley/devtools) package
```
> install.packages("devtools")
> library(devtools)
> install_github("zhushijia/STAR2bSMRT")
```

## Tutorial
   See our [wiki](https://github.com/zhushijia/STAR2bSMRT/wiki)

## Citation
Erin Flaherty*, Shijia Zhu*, et al., Neuronal impact of patient-specific aberrant NRXN1α splicing, Nature Genetics, 51:1679–1690, 2019. [(link)](https://www.nature.com/articles/s41588-019-0539-z)
