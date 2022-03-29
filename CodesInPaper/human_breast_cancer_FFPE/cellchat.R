cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2


library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")
nPatterns = 6
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")


selectK(cellchat, pattern = "incoming")
nPatterns = 6
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")


cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)










features.info4 = NULL
X = STIE_result$Signature
for(i in 1:nrow(X)) 
{
	cat(i,"\n")
	x = X[i, ]
	features = rownames(X)[i]
	for(j in 1:ncol(X))
	{
		clusters = colnames(X)[j]
		pct.1  =  x[j]
		pct.2  =  mean(x[-j])
		logFC  =  log( pct.1/pct.2 )
		pvalues = t.test(log(x[-j]+1), mu=log(x[j]+1), alternative='less')$p.val
		res = data.frame( clusters, features, logFC, pct.1, pct.2, pvalues )
		features.info4 = rbind(features.info4, res)
	}
}



X = STIE_result$Signature
features.info4 = do.call(rbind, lapply( 1:nrow(X), function(i) {
	cat(i,"\n")
	x = X[i, ]
	features = rownames(X)[i]
	do.call(rbind, lapply( 1:length(x), function(j) {
		clusters = colnames(X)[j]
		pct.1  =  x[j]
		pct.2  =  mean(x[-j])
		diff  =  log( pct.1/pct.2 )
		pvalues = t.test(log(x[-j]+1), mu=log(x[j]+1), alternative='less')$p.val
		data.frame( clusters, features, diff, pct.1, pct.2, pvalues )
	}))
}))

logfc.threshold = 0.25
min.pct = 0.1

x = subset(features.info4, logFC>0.25 & pct.1>0.1 & pvalues<0.01 )
table(x$clusters)


thres=0.05/19743
x3 = subset(f3, diff>0 & pct.1>0.5        & pvalues<thres ); mean(x3$features%in%rownames(truth))
x4 = subset(f4, diff>0 & log(pct.1+1)>0.5 & pvalues<1e-2 ); mean(x4$features%in%rownames(truth))

x3 = do.call(rbind, lapply(1:6,function(k) {
    tmp = subset(f3, clusters==as.character(k))
    tmp = tmp[order(tmp$pvalues),]
    ind = with(tmp, which( diff>0 & pct.1>0.1 & pvalues<0.05/19743) )
    if(length(ind)>20) {
        tmp[ind,]
    } else {
        subset(tmp, pvalues<1e-2)
    }
} ))

x4 = do.call(rbind, lapply(1:6,function(k) {
    tmp = subset(f4, clusters==as.character(k))
    tmp = tmp[order(tmp$pvalues),]
    ind = with(tmp, which( diff>0 & log(pct.1+1)>0.5 & pvalues<1e-2 ) )
    if(length(ind)>20) {
        tmp[ind,]
    } else {
        tmp[1:20,]
    }
} ))



dim(x3)
table(x3$clusters)

dim(x4)
table(x4$clusters)


k = k+1
kc = as.character(k)
subset(x4,clusters==kc)$features %in% subset(x3,clusters==kc)$features

y3 = subset(x3,clusters==kc); y3 = y3[order(y3$pvalues), ]
y4 = subset(x4,clusters==kc); y4 = y4[order(y4$pvalues), ]

a = setdiff( y3$features, y4$features )
subset(f3, features%in%a & clusters==k)

intersect(a,rownames(truth))

b = y3$features[1:min(200,nrow(y3))]
c = y4$features[1:min(200,nrow(y4))]
c(sum(b%in%rownames(truth)), length(b) )
c( sum(c%in%rownames(truth)), length(c) )




