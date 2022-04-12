deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_hippocampus/parameters_hippo.R")

########################################################################
######## STIE
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/STIE")
load("MouseBrainHippo_spot_BIC_new_cells_on_spot.RData")
stie = result[['2']]
stie_prop = table( stie$cells_on_spot$spot, stie$cell_types )
stie_prop = t(apply(stie_prop,1,function(x)x/sum(x)))
stie_prop = stie_prop[,order(colnames(stie_prop))]
stie_corr = cor(stie_prop)

########################################################################
######## DWLS
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/BroadInstitute_SingleCell")
x = load("Visium_FFPE_Mouse_Brain_hippocampus_DWLS_on_BroadInstitute_Signature.RData" )
dwls_prop = prop_mat
dwls_prop = dwls_prop[, order(colnames(dwls_prop))]
dwls_corr = cor(dwls_prop)

########################################################################
######## SPOTlight
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/SPOTlight")
spotlight_ls <- readRDS(file = "spotlight_ls.rds" )
nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]
spotlight_prop = decon_mtrx[,1:6]
spotlight_prop = spotlight_prop[,order(colnames(spotlight_prop))]
spotlight_corr = cor(spotlight_prop)


boxplot(abs(stie_corr)[upper.tri(stie_corr)], 
        abs(dwls_corr)[upper.tri(stie_corr)], 
        abs(spotlight_corr)[upper.tri(stie_corr)])


boxplot((stie_corr)[upper.tri(stie_corr)], 
        (dwls_corr)[upper.tri(stie_corr)], 
        (spotlight_corr)[upper.tri(stie_corr)])

