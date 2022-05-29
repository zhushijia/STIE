library(Matrix)
library(Seurat)

#####################################################################################
######## scRNAseq data
#####################################################################################

refdir <- "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq"
setwd(refdir)
counts <- readMM("count_matrix_sparse.mtx")
genes <- read.delim("count_matrix_genes.tsv",header=F)
cell_ids <- read.delim("count_matrix_barcodes.tsv",header=F)
rownames(counts) <- as.character(genes[,1])
colnames(counts) <- as.character(cell_ids[,1])

meta_data = read.csv("metadata.csv",header=T,row.names=1)
cell_types <- meta_data$celltype_major
names(cell_types) <- rownames(meta_data) # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- meta_data$nCount_RNA
names(nUMI) <- rownames(meta_data) # create nUMI named list

#####################################################################################
######### downsampling to 15%
#####################################################################################
uni_celltypes = unique(cell_types)
set.seed(1234567)
tmp = lapply(uni_celltypes, function(x) {
    index = which(cell_types==x)
    sample( index, ceiling(length(index)*0.05) )
} )
selected = c()
for(k in 1:length(tmp)) selected = c(selected, tmp[[k]])
counts = counts[,selected]
cell_types = cell_types[selected]
#####################################################################################
sc_cnt = data.frame( cell=colnames(counts), t(counts) )
sc_mta = data.frame( cell=names(cell_types), bio_celltypeXXX=cell_types )
setwd("/archive/SCCC/Hoshida_lab/s184554/Code/github/stereoscope/data/HumanBreastCancer")
write.table( sc_cnt, "sc_cnt.tsv", sep="\t", quote=F, col.names=T, row.names=F)
write.table( sc_mta, "sc_mta.tsv", sep="\t", quote=F, col.names=T, row.names=F)


#####################################################################################
######## spatial data
#####################################################################################
datadir <- '/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/Visium_FFPE_Human_Breast_Cancer'
matrix_paths <- paste( datadir, "outs/filtered_feature_bc_matrix.h5", sep="/")
tissue_paths <- paste( datadir, "outs/spatial/tissue_positions_list.csv", sep="/")
counts <- as.matrix(Read10X_h5(matrix_paths))
coords <- read.csv(tissue_paths, header=F, row.names=1)
coords <- coords[,c(4:5)]
colnames(coords) = c("xcoord", "ycoord")

st_cnt = data.frame( cell=colnames(counts), t(counts) )
setwd("/archive/SCCC/Hoshida_lab/s184554/Code/github/stereoscope/data/HumanBreastCancer")
write.table( st_cnt, "st_cnt.tsv", sep="\t", quote=F, col.names=T, row.names=F)




#####################################################################################
module load python/3.8.x-anaconda
source activate myVE.python3.8
module load gcc/6.3.0

cd /archive/SCCC/Hoshida_lab/s184554/Code/github/stereoscope/res

stereoscope run \
--sc_cnt ../data/HumanBreastCancer/sc_cnt.tsv \
--sc_labels ../data/HumanBreastCancer/sc_mta.tsv \
--st_cnt ../data/HumanBreastCancer/st_cnt.tsv \
-o HumanBreastCancer




