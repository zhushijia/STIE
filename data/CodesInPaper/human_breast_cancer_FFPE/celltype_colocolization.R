################################################################################################
########## load data
################################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq")
y = load("Wu_etal_2021_BRCA_scRNASeq_DWLS_Signature.RData")
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/Visium_FFPE_Human_Breast_Cancer/DWLS_on_Wu_etal_2021_BRCA_Signature")
x = load("Visium_FFPE_Human_Breast_Cancer_DWLS_on_Wu_etal_2021_BRCA_Signature.RData" )

prop_mat = prop_mat[, !colnames(prop_mat)%in%c("Myeloid","PVL")]
Signature = Signature[, match(colnames(prop_mat),colnames(Signature))]

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
z = load("BreastCancer_spot_BIC_new_get_cells_on_spot.RData")

################################################################################################
######### cell type colocolization
################################################################################################
rmse = sapply(score,function(x) mean(sqrt(x$mse)))
f = function(x) (x-min(x))/(max(x)-min(x))
plot( ratio, f(mse), type="l")
lines( ratio, f(cell_count) )


method = 'pearson'
c1 = cor(prop_mat, method=method)
c2 = cor(score[['6']]$celltypes_on_spot, method=method)
cs = cor(Signature, method=method)

index = upper.tri(c1)
wilcox.test( c(c1[index])^2, c(c2[index])^2 , 'greater', paired=T)
t.test( c(c1[index])^2, c(c2[index])^2 , 'greater', paired=T)

c1^2 - c2^2

par(mfrow=c(2,2))
plot( c(c2[index])^2, c(c1[index])^2)
abline(b=1,a=0)

boxplot( c(c1[index])^2 , c(c2[index])^2 )

plot( c1[index], cs[index] )
plot( c2[index], cs[index] )

# to test linear relationship
# both are associated with signature correlation, but c1 is driven by it,  
# while c2 is not
cor.test( c1[index], cs[index], method='pearson' )
cor.test( c2[index], cs[index], method='pearson' )

