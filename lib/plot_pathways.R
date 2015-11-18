library(gplots)
source("/Users/charlesmurphy/Desktop/Research/lib/gsea.R")

plot_pathway <- function(fpkm,fpkm_threshold=0.1,gmt,pathway,samples,prefix,main) {
	fpkm <- as.matrix(read.csv(fpkm,header=TRUE,row.names=1))

	gmt <- parse_gmt(gmt=gmt)

	genes = gmt[[pathway]]
	genes.2 <- c()
	for (i in 1:length(genes)) {
		if (genes[i] %in% row.names(fpkm)) {
			genes.2 <- c(genes.2,genes[i])
		}
	}
	fpkm <- fpkm[genes.2,samples]
	fpkm[fpkm<=fpkm_threshold] <- 0.0
	fpkm <- fpkm[rowSums(fpkm)>0,]

	png(paste(prefix,".png",sep=""),width=1000,height=1000)
	heatmap.2(log(fpkm+0.1),trace="none",keysize=1.0,scale="row",margins=c(8,12),cexRow=0.25,cexCol=1.2,col= redgreen(75),key.title="Gene expression",main=main)
	dev.off()
	write.csv(fpkm,paste(prefix,".csv",sep=""))
}
