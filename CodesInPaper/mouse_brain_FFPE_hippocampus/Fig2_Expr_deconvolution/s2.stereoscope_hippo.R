library(Matrix)
library(Seurat)
library(data.table)
#####################################################################################
######## scRNAseq data
#####################################################################################

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/BroadInstitute_SingleCell")
count = data.frame(fread("DATA_MATRIX_LOG_TPM.txt"))
rownames(count) = as.character(count$GENE)
colnames(count) = gsub("[.]","_",colnames(count))
count = count[,-1]
#count = exp(count) - 1
annot = read.delim("CLUSTER_AND_SUBCLUSTER_INDEX.txt",sep='\t',skip=1,header=T)
annot$TYPE = gsub("-","_",as.character(annot$TYPE))
annot$TYPE = gsub("[.]","_",as.character(annot$TYPE))
colnames(annot)[2] = "subclass"
rownames(annot) = as.character(annot$TYPE)

all( colnames(count) %in% as.character(annot$TYPE) )
annot = annot[ match(colnames(count), as.character(annot$TYPE) ) , ]
all( colnames(count) == rownames(annot) )

index = which( !as.character(annot$subclass)%in%c("Ependymal","Non") )
counts = count[,index]
cell_types = annot[index,]
#####################################################################################
sc_cnt = data.frame( cell=colnames(counts), t(counts) )
sc_mta = data.frame( cell=as.character(cell_types$TYPE), bio_celltypeXXX=as.character(cell_types$subclass) )

setwd("/archive/SCCC/Hoshida_lab/s184554/Code/github/stereoscope/data/MouseBrainHippo")
write.table( sc_cnt, "sc_cnt.tsv", sep="\t", quote=F, col.names=T, row.names=F)
write.table( sc_mta, "sc_mta.tsv", sep="\t", quote=F, col.names=T, row.names=F)


#####################################################################################
######## spatial data
#####################################################################################
datadir <- '/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/Visium_FFPE_Mouse_Brain'
matrix_paths <- paste( datadir, "outs/filtered_feature_bc_matrix.h5", sep="/")
tissue_paths <- paste( datadir, "outs/spatial/tissue_positions_list.csv", sep="/")
counts <- as.matrix(Read10X_h5(matrix_paths))
coords <- read.csv(tissue_paths, header=F, row.names=1)
coords <- coords[,c(4:5)]
colnames(coords) = c("xcoord", "ycoord")

st_cnt = data.frame( cell=colnames(counts), t(counts) )
setwd("/archive/SCCC/Hoshida_lab/s184554/Code/github/stereoscope/data/MouseBrainHippo")
write.table( st_cnt, "st_cnt.tsv", sep="\t", quote=F, col.names=T, row.names=F)




#####################################################################################
module load python/3.8.x-anaconda
source activate myVE.python3.8
module load gcc/6.3.0

cd /archive/SCCC/Hoshida_lab/s184554/Code/github/stereoscope/res

stereoscope run \
--sc_cnt ../data/MouseBrainHippo/sc_cnt.tsv \
--sc_labels ../data/MouseBrainHippo/sc_mta.tsv \
--st_cnt ../data/MouseBrainHippo/st_cnt.tsv \
-o MouseBrainHippo




