deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

load("parameter_BreastCancer.R")

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

results = result

############################################################
########## deep learning features
############################################################
library(data.table)
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/")
cells = data.frame(fread("Breast_Cancer_FFPE_information_model40_image20.csv"))
dl_cell_feature = with(cells, data.frame(
                    celltypes=nucleus_class,
                    probability=probability,
                    X=centroid_x0,
                    Y=centroid_y0,
                    Area=area,
                    Major=major_axis_length,
                    Minor=minor_axis_length,
                    Eccentricity=eccentricity,
                    Orientation=orientation,
                    Solidity=solidity) )
dl_cell_feature = subset(dl_cell_feature, celltypes!="blood")
dl_cell_feature = dl_cell_feature[order(dl_cell_feature$probability,decreasing=T),]
###########################################################
######### 
############################################################

result = results[['2.8']]

stie_cell_types = result$cell_types
stie_cell_types = stie_cell_types[ match(unique(names(stie_cell_types)), names(stie_cell_types)) ]

stie_cell_feature = morphology_fts[ match( names(stie_cell_types), rownames(morphology_fts) ), ]
stie_cell_contour = cell_info$cell_contour[ match( names(stie_cell_types), names(cell_info$cell_contour) ) ]

###########################################################
######### 
############################################################
dl_cell_feature2 = subset(dl_cell_feature,probability>0.5)
overlap_index = lapply( stie_cell_contour, function(x) {
    i = with( dl_cell_feature2, which( X>min(x[,1]) & X<max(x[,1]) & Y>min(x[,2]) & Y<max(x[,2]) ) )
    cat(i,"\n")
    i
} )

dl_cell_types = sapply(overlap_index, function(i) {
    cc = ""
    if(length(i)>0) {
        j = i[which.max(dl_cell_feature2$probability[i])]
        cc = as.character(dl_cell_feature2$celltypes[j])
    }
    cc
}) 

table(dl_cell_types, stie_cell_types)


dl_cell_types2 = dl_cell_types
dl_cell_types2[dl_cell_types%in%c("tumor", "necrosis", "ductal epithelium")] = "Epithelia"
dl_cell_types2[dl_cell_types%in%c("lymphocyte", "macrophage")] = "Immune"
dl_cell_types2[dl_cell_types%in%c("stroma")] = "Stroma"

stie_cell_types2 = stie_cell_types
stie_cell_types2[stie_cell_types%in%c("CancerEpithelial", "NormalEpithelial")] = "Epithelia"
stie_cell_types2[stie_cell_types%in%c("Tcells", "Bcells", "Plasmablasts")] = "Immune"
stie_cell_types2[stie_cell_types%in%c("CAFs", "Endothelial")] = "Stroma"

table(dl_cell_types2, stie_cell_types2)
x = table(dl_cell_types2, stie_cell_types2)
sum(diag(x[-1,]))/sum(x[-1,])
sum(diag(x[-1,]))/sum(x)


stie_cell_types3 = paste0("DeepLearning_",dl_cell_types2, "_vs_STIE_",stie_cell_types2)
names(stie_cell_types3) = names(stie_cell_types2)
stie_cell_types3 = stie_cell_types3[ stie_cell_types2!=dl_cell_types2 ]

###########################################################
######### 
############################################################
if(0)
{
    contour3 = contour[ match(names(stie_cell_types3),names(contour)) ]
    
    plot_sub_image(im=im, w=3000, h=3000, xoff=9000, yoff=9000, 
                   x_scale=1, 
                   contour=contour3, cell_types=stie_cell_types3, plot_spot=F, plot_cell=T  )
    
    plot_sub_image(im=im, 
                   x_scale=args$x_scale, 
                   contour=contour3, cell_types=stie_cell_types3, plot_spot=F, plot_cell=T )
    
}

############################################################
###########################################################
############################################################
cells_on_spot <- result$cells_on_spot
Signature = result$Signature

cells_on_spot3 = cells_on_spot[ as.character(cells_on_spot$cell_id) %in% names(stie_cell_types3), ]
cells_on_spot_celltypes3 = stie_cell_types3[ match( as.character(cells_on_spot3$cell_id), names(stie_cell_types3) ) ]

tt = table(cells_on_spot3$spot, cells_on_spot_celltypes3)

STdata = load_STdata( "Visium_FFPE_Human_Breast_Cancer", args$spaceranger_count_dir, umi_cutoff=0, is_normalized=TRUE ) 
normalized_count = STdata$normalized_matrix[[1]]
#count = STdata$matrix[[1]]
tt = tt[rownames(tt)%in%rownames(normalized_count), ]
ST_expr = normalized_count[match( rownames(tt), rownames(normalized_count) ),]

Signature = Signature[rownames(Signature) %in% colnames(ST_expr),]
ST_expr = ST_expr[ , match( rownames(Signature), colnames(ST_expr) ) ]



cor( tt, ST_expr[,"CD4"])
cor( tt, ST_expr[,"CD34"])

cor( tt, ST_expr[,"EPCAM"])

cor( tt, ST_expr[,"CD34"])

cor( tt, Signature)


(i = i+1)
index = order(tt[,i],decreasing=T)[1:10]
merged_expr = t(ST_expr[index,]) %*% apply(as.matrix.noquote(tt[index,]),2,as.numeric)


cor(Signature,as.numeric(merged_expr[,i]),method='spearman')
cat(colnames(merged_expr)[i],"\n")
x = solveOLS(S = Signature, B=merged_expr[,i], scaled = T)
sort(signif(x,3),decreasing=T)


X = data.frame(log2(Signature+1))
Y = as.numeric(log2(merged_expr[,i]+1))
z = lm(Y~Bcells+Tcells+Plasmablasts, data=X )
summary(z)$adj.r.squared # immune

z = lm(Y~Endothelial+CAFs, data=X )
summary(z)$adj.r.squared # stroma

z = lm(Y~CancerEpithelial+NormalEpithelial, data=X )
summary(z)$adj.r.squared # Epi

z = lm(Y~Bcells+Tcells+Plasmablasts + Endothelial+CAFs, data=X )
summary(z)$coef




