deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")

########################################################################
######## STIE
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
load("BreastCancer_spot_BIC_new_get_cells_on_spot.RData")
stie = result[['2']]


setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq")
y = load("Wu_etal_2021_BRCA_scRNASeq_DWLS_Signature.RData")
ord = c("Plasmablasts", "Bcells", "Tcells", "Myeloid", "CAFs", "Endothelial", "PVL", "NormalEpithelial", "CancerEpithelial")
sig_corr = cor(Signature[,ord],method='pearson')


cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.8*spot_radius)
stie <- STIE(ST_expr, Signature, cells_on_spot, features, 
               lambda=0, steps=30, min_cells=-1, 
               known_signature=TRUE, known_cell_types=FALSE)

stie_prop = table( stie$cells_on_spot$spot, stie$cell_types )
stie_prop = data.frame( t(apply(stie_prop,1,function(x)x/sum(x))) )
stie_prop$NormalEpithelial = 0
stie_prop = stie_prop[,intersect( ord, colnames(stie_prop) )]
stie_corr = cor(stie_prop, method='pearson')

########################################################################
######## DWLS
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/Visium_FFPE_Human_Breast_Cancer/DWLS_on_Wu_etal_2021_BRCA_Signature")
x = load("Visium_FFPE_Human_Breast_Cancer_DWLS_on_Wu_etal_2021_BRCA_Signature.RData" )
#prop_mat = prop_mat[, !colnames(prop_mat)%in%c("Myeloid","PVL")]
dwls_prop = prop_mat
dwls_prop = dwls_prop[, intersect( ord, colnames(dwls_prop) )]
dwls_corr = cor(dwls_prop, method='pearson')

########################################################################
######## SPOTlight
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/SPOTlight")
spotlight_ls <- readRDS(file = "spotlight_ls.rds" )
nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]
rownames(decon_mtrx) = rownames(ST_expr)
spotlight_prop = decon_mtrx[,-ncol(decon_mtrx)]
spotlight_prop = spotlight_prop[,intersect( ord, colnames(spotlight_prop) )]
spotlight_corr = cor(spotlight_prop, method='pearson')


setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
corr.write = function(corr, filename) {
    diag(corr) = NA
    corr2 = data.frame(id=rownames(corr), signif(corr,2) )
    write.table(corr2, filename, sep="\t", col.names=T, row.names=F, quote=F)
}
corr.write(stie_corr,"breastcancer_cell_type_colocalizalition_stie.txt")
corr.write(dwls_corr,"breastcancer_cell_type_colocalizalition_dwls.txt")
corr.write(spotlight_corr,"breastcancer_cell_type_colocalizalition_spotlight.txt")
######################################################

pdf("breastcancer_cell_type_colocalizalition.pdf", w=6, h=7)

par(mfrow=c(2,2))
x = list(dwls = dwls_corr[upper.tri(dwls_corr)], 
         spotlight = spotlight_corr[upper.tri(spotlight_corr)], 
         stie = stie_corr[upper.tri(stie_corr)] )
boxplot(x,las=3, ylim=range(do.call(c,x),na.rm=T), outline=F )
stripchart(x,method='jitter',add=T,vertical=T)

######################################################
par(mfrow=c(2,2))
x = sig_corr[upper.tri(sig_corr)]
y = dwls_corr[upper.tri(dwls_corr)]
z = stie_corr[upper.tri(stie_corr)]
plot( x, y, pch=16 )
abline( lm( y~x ) )
plot( x, z, pch=16 )
abline( lm( z~x ) )

cor.test( stie_corr[upper.tri(stie_corr)], sig_corr[upper.tri(sig_corr)] )
cor.test( dwls_corr[upper.tri(dwls_corr)], sig_corr[upper.tri(sig_corr)] )

dev.off()





