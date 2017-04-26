library(calibrate)
library(gplots)
library(genefilter)
library(dendextend)
library("beeswarm")


clustering_variety <- function(outdir,data,sample_data,method="") {

	#cluster on correlation distance

	diss <- 1 - cor(data,use="complete.obs",method=method)
	diss.2 <- as.dist(diss)
	dend <- as.dendrogram(hclust(diss.2,method="average"))

	#original labels
	png("./Cluster/cluster-all-samples-correlation-distance.png",width=1200,height=600)
	par(mar=c(11,7,7,7),cex=1.5)
	plot(dend,main="Hierarchical clustering on all genes",ylab="Correlation distance")
	dev.off()

	#clsuter on euclidean distance

	diss <- dist(t(data))
	dend <- as.dendrogram(hclust(diss,method="average"))

	png("./Cluster/cluster-all-samples-euclidean.png",width=1200,height=600)
	par(mar=c(11,7,7,7),cex=1.5)
	plot(dend,main="Hierarchical clustering on all genes",ylab="Euclidean distance")
	dev.off()

	#cluster on top 5000 more variable genes

	vars <- rowVars(data)
	t1 <- sort(vars,index.return=T,decreasing=T)
	data2 <- data[t1$ix,]
	data2 <- data2[1:5000,]

	diss <- dist(t(data2))
	dend <- as.dendrogram(hclust(diss,method="average"))

	png("./Cluster/cluster-all-samples-euclidean-top5000.png",width=1200,height=600)
	par(mar=c(11,7,7,7),cex=1.5)
	plot(dend,main="Hierarchical clustering on all genes (top 5000 genes by variance)",xlab="")
	dev.off()

	#cluster on top 1000 more variable gene

	vars <- rowVars(data)
	t1 <- sort(vars,index.return=T,decreasing=T)
	data2 <- data[t1$ix,]
	data2 <- data2[1:1000,]

	diss <- dist(t(data2))
	dend <- as.dendrogram(hclust(diss,method="average"))

	png("./Cluster/cluster-all-samples-euclidean-top1000.png",width=1200,height=600)
	par(mar=c(11,7,7,7),cex=1.5)
	plot(dend,main="Hierarchical clustering on all genes (top 1000 genes by variance)",xlab="")
	dev.off()

}

perform_pcas <- function(outdir,data,sample_data,sample_labels) {

	re <- prcomp(t(data))

	pcs=re$sdev**2
	pcs = pcs/sum(pcs) * 100

	pdf(paste(outdir,"/PC-percent-variance.pdf",sep=""))
	par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE,cex=1.0)
	barplot(
		pcs,
		names.arg=unlist(lapply(1:length(pcs),function(x) return(paste("PC",x,sep="")))),
		las=2,
		ylab="Percentage of variance (%)",
		ylim=c(0,100)
	)
	dev.off()

	#save loadings
	write.csv(re$rotation,paste(outdir,"/PCA-loadings.csv",sep=""))

	write.csv(re$x,paste(outdir,"/PCA-PCs.csv",sep=""))

	colors <- rainbow(max(sample_labels[[1]]),start=0.1,end=0.9)[sample_labels[[1]]]

	legend_colors = rainbow(max(sample_labels[[1]]),start=0.1,end=0.9)


	#pca, color samples by cohort
	for (i in 1:NUMBER_PCS) {
		for (j in 1: NUMBER_PCS) {
			if (i == j) {
				next
			}
			pdf(paste(outdir,"/PC",i,"PC",j,".pdf",sep=""))
			par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE,cex=1.5)
			plot(
				re$x[,i],
				re$x[,j],
				col=colors,
				pch=16,
				ylab=paste("PC",j,sep=""),
				xlab=paste("PC",i,sep=""),
				main="PCA on all samples",
				ylim=c(1.2*min(re$x[,j]),1.2*max(re$x[,j])),
				xlim=c(1.2*min(re$x[,i]),1.2*max(re$x[,i]))
			)
			textxy(re$x[,i],re$x[,j],names(re$x[,i]),cex=0.4)
			legend(x="topright",names(sample_labels[[2]]),cex=0.85,inset=c(-0.55,0),pch=16,col=legend_colors,title="Group")
			dev.off()
		}
	}
}
