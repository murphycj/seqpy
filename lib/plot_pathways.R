library(gplots)
source("/Users/charlesmurphy/Desktop/Research/lib/seqpy/lib/gsea.R")

plot_pathway <- function(data,threshold=0.1,gmt,pathway,samples,prefix,main,margins=c(8,12),cexRow=0.25,cexCol=1.2,log_data=T,scale="row") {
	data <- as.matrix(read.csv(data,header=TRUE,row.names=1,check.names=F))

	gmt <- parse_gmt(gmt=gmt)

	genes = gmt[[pathway]]
	genes.2 <- c()
	for (i in 1:length(genes)) {
		if (genes[i] %in% row.names(data)) {
			genes.2 <- c(genes.2,genes[i])
		}
	}

	if (length(samples)==0) {
		samples = colnames(data)
	}

	data <- data[genes.2,samples]
	data[data<=threshold] <- 0.0
	data <- data[rowSums(data)>0,]

	if (log_data) {
		data = log(data+0.1)
	}

	png(paste(prefix,".png",sep=""),width=1000,height=1000)
	heatmap.2(
		data,
		trace="none",
		keysize=1.0,
		scale=scale,
		margins=margins,
		cexRow=cexRow,
		cexCol=cexCol,
		col=greenred(75),
		key.title="Gene expression",
		main=main
	)
	dev.off()
	write.csv(data,paste(prefix,".csv",sep=""))
}

plot_multiple_pathways <- function(fpkm,fpkm_threshold=0.1,gmt,pathways,samples,prefix,main) {
	#
	# plot the genes involved in multiple pathways
	#

	fpkm <- as.matrix(read.csv(fpkm,header=TRUE,row.names=1))

	gmt <- parse_gmt(gmt=gmt)

	genes.2 <- c()
	for (j in pathways) {
		genes = gmt[[j]]
		for (i in 1:length(genes)) {
			if (genes[i] %in% row.names(fpkm)) {
				genes.2 <- c(genes.2,genes[i])
			}
		}
	}
	genes.2 <- unique(genes.2)
	fpkm <- fpkm[genes.2,samples]
	fpkm[fpkm<=fpkm_threshold] <- 0.0
	fpkm <- fpkm[rowSums(fpkm)>0,]

	png(paste(prefix,".png",sep=""),width=1000,height=1000)
	heatmap.2(log(fpkm+0.1),trace="none",keysize=1.0,scale="row",margins=c(8,12),cexRow=0.25,cexCol=1.2,col= redgreen(75),key.title="Gene expression",main=main)
	dev.off()
	write.csv(fpkm,paste(prefix,".csv",sep=""))
}


plot_pathway_leading_edge <- function(fpkm,fpkm_threshold=0.1,gmt,pathway,samples,prefix,main,rnk,gsea_results) {
	#
	# plot a heatmap of the genes in the leading edge of a pathway
	#

	rnk <- read.table(rnk,sep='\t',header=F)
	rnk <- rnk[order(rnk[,2],decreasing=T),]

	gsea_results <- read.table(gsea_results,header=T,row.names=1,sep="\t")

	rank_at_max <- gsea_results[pathway,"RANK.AT.MAX"]
	rnk <- rnk[1:rank_at_max,]

	fpkm <- as.matrix(read.csv(fpkm,header=TRUE,row.names=1))

	gmt <- parse_gmt(gmt=gmt)

	genes = gmt[[pathway]]
	genes.2 <- c()
	for (i in 1:length(genes)) {
		if ((genes[i] %in% row.names(fpkm)) & (genes[i] %in% rnk[,1])) {
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


plot_pathways_leading_edge <- function(fpkm,fpkm_threshold=0.1,gmt,pathways,samples,prefix,main,rnk,gsea_results) {
	#
	# plot a heatmap of the genes in the leading edge of a pathway
	#

	rnk <- read.table(rnk,sep='\t',header=F)
	rnk <- rnk[order(rnk[,2],decreasing=T),]

	gsea_results <- read.table(gsea_results,header=T,row.names=1,sep="\t")

	fpkm <- as.matrix(read.csv(fpkm,header=TRUE,row.names=1))

	gmt <- parse_gmt(gmt=gmt)

	genes.2 <- c()
	for (j in pathways) {
		genes = gmt[[j]]
		rank_at_max <- gsea_results[j,"RANK.AT.MAX"]
		rnk.2 <- rnk[1:rank_at_max,]
		for (i in 1:length(genes)) {
			if ((genes[i] %in% row.names(fpkm)) & (genes[i] %in% rnk.2[,1])) {
				genes.2 <- c(genes.2,genes[i])
			}
		}
	}
	genes.2 <- unique(genes.2)
	fpkm <- fpkm[genes.2,samples]
	fpkm[fpkm<=fpkm_threshold] <- 0.0
	fpkm <- fpkm[rowSums(fpkm)>0,]

	png(paste(prefix,".png",sep=""),width=1000,height=1000)
	heatmap.2(log(fpkm+0.1),trace="none",keysize=1.0,scale="row",margins=c(8,12),cexRow=0.25,cexCol=1.2,col= redgreen(75),key.title="Gene expression",main=main)
	dev.off()
	write.csv(fpkm,paste(prefix,".csv",sep=""))
}
