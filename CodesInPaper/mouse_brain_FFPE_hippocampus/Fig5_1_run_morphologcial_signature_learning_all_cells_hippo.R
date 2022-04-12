deconvolution = TRUE
source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_hippocampus/parameters_hippo3.R")

############################################################
## deconvolution on multiorgan
############################################################

cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)

result <- STIE(ST_expr, Signature, cells_on_spot, features, 
               lambda=0, steps=30, 
               known_signature=TRUE, known_cell_types=FALSE)

################################################################################

cells_on_spot$cell_types = result$cell_types
spot_id = as.character(cells_on_spot$spot)
cell_id = as.character(cells_on_spot$cell_id)
Signature = Signature[rownames(Signature)%in%colnames(ST_expr),]
ST_expr2 = t( ST_expr[ rownames(ST_expr)%in%spot_id, match(rownames(Signature),colnames(ST_expr)) ] )

coefs = table( as.character(cells_on_spot$spot), as.character(cells_on_spot$cell_types) )
coefs = coefs[match( colnames(ST_expr2), rownames(coefs) ), ]
Signature2 = t( apply( ST_expr2, 1, function(x) solveNNLS( coefs, as.matrix(x), scaled=F ) ) )

cc = data.frame( celltypes = c( "CA1", "CA2",  "CA3", "DG", "GABAerig", "Glia"), 
                      cols = c( "magenta", "blue", "green", "black", "cyan", "orange") )

Signature = Signature[, intersect(as.character(cc$celltypes),colnames(Signature)) ]
Signature2 = Signature2[, intersect(as.character(cc$celltypes),colnames(Signature2)) ]

prop2 = apply( Signature2, 2, function(x) solveNNLS(Signature,x) )  
prop_dwls2 = dwls(S=Signature, B=Signature2)
corr2 = cor(log2(Signature2+1),log2(Signature+1))


setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/STIE")
write.table( data.frame(celltype=rownames(corr2),corr2), "hippo_signature_learning_corr.txt", sep="\t", col.names=T, row.names=F, quote=F )

pdf("hippo_signature_learning.pdf")
par(mfrow=c(2,2))
cols = as.character(cc$cols)[ match( colnames(prop2), as.character(cc$celltypes) ) ]
barplot(prop2, col=cols, las=3, border='white')
barplot(prop_dwls2, col=cols, las=3, border='white')
heatmap(corr2,Colv=NA,Rowv=NA)
dev.off()


