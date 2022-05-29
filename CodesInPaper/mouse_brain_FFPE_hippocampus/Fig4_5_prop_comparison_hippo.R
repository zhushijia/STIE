library(Seurat)
matrix_paths = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/Visium_FFPE_Mouse_Brain/outs/filtered_feature_bc_matrix.h5"
count <- as.data.frame(t(as.matrix(Read10X_h5(matrix_paths))))

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/BroadInstitute_SingleCell")
Signature = read.delim("Major_cell_types_marker_genes.txt", header=T, row.names=1)
Signature = Signature[,!colnames(Signature)%in%c("Ependymal.cells")]

ST_expr = count[,colnames(count)%in%rownames(Signature)]
Signature = as.matrix(Signature)[rownames(Signature)%in%colnames(ST_expr),]
colnames(Signature) = gsub("Pyramidal.neurons.|.interneurons|.like.cells|Granule.cells.", "", colnames(Signature))
Signature = Signature[,order(colnames(Signature))]
Signature = Signature[rownames(Signature) %in% colnames(count) ,  ]


ord = colnames(Signature)
sig_corr = cor(Signature[,ord],method='pearson')


###########################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/STIE")
load("MouseBrainHippo_spots.RData")
stie = results[['2.5']]
stie_prop = table( stie$cells_on_spot$spot, stie$cell_types )
stie_prop = data.frame( t(apply(stie_prop,1,function(x)x/sum(x))) )
stie_prop$GABAergic = 0
stie_prop = stie_prop[,ord]
stie_corr = cor(stie_prop)
stie_corr[is.na(stie_corr)] = 0

###########################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/DWLS")
x = load("Visium_FFPE_Mouse_Brain_hippocampus_DWLS_on_BroadInstitute_Signature.RData")
dwls_prop = prop_mat
dwls_prop = dwls_prop[, ord]
dwls_corr = cor(dwls_prop)

###########################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/SPOTlight")
spotlight_ls <- readRDS(file = "spotlight_ls_log.rds" )
nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]
spotlight_prop = decon_mtrx
spotlight_prop = spotlight_prop[, colnames(spotlight_prop)!="res_ss"]
spotlight_prop = spotlight_prop[, ord]
rownames(spotlight_prop) = rownames(dwls_prop)
spotlight_corr = cor(spotlight_prop)

###########################################################################################
setwd("/archive/SCCC/Hoshida_lab/s184554/Code/github/stereoscope/res/MouseBrainHippo/st_cnt")
stereoscope_prop = read.delim("W.2022-05-05000218.951870.tsv", sep="\t", header=T, row.names=1)
stereoscope_prop = stereoscope_prop[, order(colnames(stereoscope_prop))]
colnames(stereoscope_prop) = gsub("[.]","",colnames(stereoscope_prop))
stereoscope_prop = stereoscope_prop[,ord]
stereoscope_corr = cor(stereoscope_prop)

###########################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/RCTD")
x = load("myRCTD.RData" )
rctd_prop = as.matrix(norm_weights)
rctd_prop = rctd_prop[, order(colnames(rctd_prop))]
colnames(rctd_prop) = gsub("[ ]|-", "", colnames(rctd_prop))
rctd_prop = rctd_prop[,ord]
rctd_corr = cor(rctd_prop)

###########################################################################################
######## correlation with signature similarity
###########################################################################################

cor.test( sig_corr[upper.tri(sig_corr)], dwls_corr[upper.tri(dwls_corr)] )
cor.test( sig_corr[upper.tri(sig_corr)], stereoscope_corr[upper.tri(stereoscope_corr)] )
cor.test( sig_corr[upper.tri(sig_corr)], spotlight_corr[upper.tri(spotlight_corr)] )
cor.test( sig_corr[upper.tri(sig_corr)], rctd_corr[upper.tri(stie_corr)] )
cor.test( sig_corr[upper.tri(sig_corr)], stie_corr[upper.tri(stie_corr)] )


setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/STIE")

pdf("deconvolution_method_comparison_hippo.pdf", w=10, h=5)

f  = function(x,y,title) {
    plot(x,y,ylim=c(-1,1),pch=16,col='black',main=title)
    z = lm(y~x)
    abline(z, col='darkred')
}
par(mfrow=c(2,5))
f( sig_corr[upper.tri(sig_corr)], dwls_corr[upper.tri(dwls_corr)], "dwls_corr" )
f( sig_corr[upper.tri(sig_corr)], stereoscope_corr[upper.tri(stereoscope_corr)], "stereoscope_corr" )
f( sig_corr[upper.tri(sig_corr)], spotlight_corr[upper.tri(spotlight_corr)], "spotlight_corr" )
f( sig_corr[upper.tri(sig_corr)], rctd_corr[upper.tri(stie_corr)], "rctd_corr" )
f( sig_corr[upper.tri(sig_corr)], stie_corr[upper.tri(stie_corr)], "stie_corr" )

par(mfrow=c(1,4))
corrs = list(dwls = dwls_corr[upper.tri(dwls_corr)], 
             stereoscope = stereoscope_corr[upper.tri(stereoscope_corr)], 
             spotlight = spotlight_corr[upper.tri(spotlight_corr)], 
             rctd = rctd_corr[upper.tri(stie_corr)],
             stie = stie_corr[upper.tri(stie_corr)] )
boxplot(corrs,las=3, ylim=range(do.call(c,corrs),na.rm=T), outline=F )
stripchart(corrs,method='jitter',add=T,vertical=T)


###########################################################################################
######## RMSE
###########################################################################################

calculate_rmse = function( B, X=Signature, Y=t(ST_expr) ) {
    B = B[ rownames(B) %in% colnames(Y), ]
    Y = Y[, colnames(Y) %in% rownames(B)]
    Y = Y[, match(rownames(B), colnames(Y))]
    Y = Y[match(rownames(X), rownames(Y)), ]
    X = X[, colnames(X) %in% colnames(B)]
    B = B[, match(colnames(X), colnames(B))]
    mse = sapply(1:ncol(Y), function(i) {
        y = as.numeric(Y[, i])
        b = as.numeric(B[i, ])
        bx = X %*% b
        t = sum(y)/sum(bx)
        mean((y - t * bx)^2)
    })
    sqrt(mse)
}


stie_rmse = calculate_rmse( B=stie_prop, X=Signature, Y=t(ST_expr) )
dwls_rmse = calculate_rmse( B=dwls_prop, X=Signature, Y=t(ST_expr) )
spotlight_rmse = calculate_rmse( B=spotlight_prop, X=Signature, Y=t(ST_expr) )
stereoscope_rmse = calculate_rmse( B=stereoscope_prop, X=Signature, Y=t(ST_expr) )
rctd_rmse = calculate_rmse( B=rctd_prop, X=Signature, Y=t(ST_expr) )

rmse = list(NULL=NULL,
            dwls = dwls_rmse, 
            stereoscope = stereoscope_rmse, 
            spotlight = spotlight_rmse, 
            rctd = rctd_rmse,
            stie = stie_rmse,
            NULL=NULL)

par(mfrow=c(2,5))
errplot( rmse, ylab="RMSE", line=F, ylim=c(3.1,3.8) )

dev.off()


###########################################################################################
library(Seurat)
library(ggplot2)

data_dir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/Visium_FFPE_Mouse_Brain/outs"
hipo_st = Load10X_Spatial(data.dir = data_dir)
se_obj = hipo_st
slice <- names(se_obj@images)[1]
x_scale = se_obj@images[[slice]]@scale.factors$lowres
coordinates = data.frame(se_obj@images[[slice]]@coordinates)
img_path = paste0(data_dir,"/spatial/tissue_lowres_image.png")
colors = c("magenta", "blue", "green", "black", "cyan","orange")

prop_data = list(dwls = dwls_prop, 
                stereoscope = stereoscope_prop, 
                spotlight = spotlight_prop, 
                rctd = rctd_prop,
                stie = stie_prop)

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/STIE")
pdf("deconvolution_method_comparison_hippo_prop_on_image.pdf")

for(i in 1:length(prop_data))
{
    #i = i+1
    coordinates_i = coordinates[ match( rownames(prop_data[[i]]), rownames(coordinates) ),  ]
    ST_scatterpie_overlay( prop_data[[i]], coordinates_i, 
                           x_scale, cell_types_all=ord, 
                           img_path=img_path, scatterpie_alpha=1, 
                           pie_scale = 0.8, colors=colors )
    
}

dev.off()
