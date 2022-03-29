
result = list()

source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_kidney_FFPE/parameter_MouseKidney.R")
unit = spot_radius*2 / 55 / args$x_scale
result$mKidney = list( morphology_fts=morphology_fts, unit=unit)

source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")
unit = spot_radius*2 / 55 / args$x_scale
result$hBreastCancer = list( morphology_fts=morphology_fts, unit=unit)

source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_cortex/parameters_cortex.R")
unit = spot_radius*2 / 55 / args$x_scale
result$mBrain = list( morphology_fts=morphology_fts, unit=unit)

result = result[c("mKidney", "hBreastCancer", "mBrain")]

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
pdf("morphology_3datasets.pdf",h=5,w=3)
par(mfrow=c(2,2))
cols = c('steelblue','darkgrey','darkorange')
area = lapply(result, function(x) (x$morphology_fts$Area)/(x$unit)/(x$unit) )
major = lapply(result, function(x) (x$morphology_fts$Major)/(x$unit) )
minor = lapply(result, function(x) (x$morphology_fts$Minor)/(x$unit) )
ecc = lapply(result, function(x) (x$morphology_fts$Eccentricity) )
boxplot(area, outline=F, las=3, col=cols, main="Area" )
boxplot(ecc, outline=F, las=3, col=cols, main="Eccentricity" )
boxplot(major, outline=F, las=3, col=cols, main="Major" )
boxplot(minor, outline=F, las=3, col=cols, main="Minor" )
dev.off()


