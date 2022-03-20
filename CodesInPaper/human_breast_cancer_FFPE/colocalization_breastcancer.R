deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")

########################################################################
######## STIE
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
load("BreastCancer_spot_BIC_new_get_cells_on_spot.RData")
stie = result[['2']]
stie_prop = table( stie$cells_on_spot$spot, stie$cell_types )
stie_prop = t(apply(stie_prop,1,function(x)x/sum(x)))
stie_prop = stie_prop[,order(colnames(stie_prop))]
stie_corr = cor(stie_prop)

########################################################################
######## DWLS
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/Visium_FFPE_Human_Breast_Cancer/DWLS_on_Wu_etal_2021_BRCA_Signature")
x = load("Visium_FFPE_Human_Breast_Cancer_DWLS_on_Wu_etal_2021_BRCA_Signature.RData" )
prop_mat = prop_mat[, !colnames(prop_mat)%in%c("Myeloid","PVL")]
dwls_prop = prop_mat
dwls_prop = dwls_prop[, order(colnames(dwls_prop))]
dwls_corr = cor(dwls_prop)

########################################################################
######## SPOTlight
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/SPOTlight")
spotlight_ls <- readRDS(file = "spotlight_ls.rds" )
nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]
rownames(decon_mtrx) = rownames(ST_expr)
spotlight_prop = decon_mtrx[,-ncol(decon_mtrx)]
spotlight_prop = spotlight_prop[,order(colnames(spotlight_prop))]
spotlight_corr = cor(spotlight_prop)


boxplot(abs(stie_corr)[upper.tri(stie_corr)], 
        abs(dwls_corr)[upper.tri(stie_corr)], 
        abs(spotlight_corr)[upper.tri(stie_corr)])


melt_corr = function(corr)
{
    corr2 = signif(corr,2)
    corr2[!lower.tri(corr2)] = 0
    corr2 = melt(corr2)
    colnames(corr2) = c("source", "target", "weight")
    corr2 = subset(corr2, weight!=0)
    corr2 = corr2[order(corr2$weight),]
    print(corr2)
    corr2
}

plot_interaction_circos( melt_corr(dwls_corr) )
plot_interaction_circos( melt_corr(spotlight_corr) )
plot_interaction_circos( melt_corr(stie_corr) )

