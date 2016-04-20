#
#
#

library(gplots)
library("beeswarm")

process_exon_counts <- function(counts,gff,ensembl_gene,ensembl_transcript) {
  gff_gene <- gff[grep(ensembl_gene,gff[,9]),]
  gff_gene <- gff_gene[gff_gene[,3]!="aggregate_gene",]
  gff_transcript <- gff_gene[grep("ENSMUST00000172031",gff_gene[,9]),]

  #collapse by exon


  gene_exon_lengths <- gff_inpp4b[,5] - gff_inpp4b[,4]
}



plot_exon_counts_heatmap <- function(counts_file,out_prefix,ensembl_gene,ensembl_transcript,log_normalize=F) {

  counts <- read_exon_counts(counts_file=counts_file)


  legend_name <- ""
  if (log_normalize) {
    legend_name <- "log2(Normalized counts+1)"
    expression_exons = log2(expression_exons+1)
  } else {
    legend_name <- "Normalized counts"
  }

  png(paste(out_prefix,"-normalized-counts.png",sep=""),width=800,height=800)
  heatmap.2(
    expression_exons,
    trace="none",
    keysize=1.0,
    scale="column",
    margins=c(4,12),
    cexRow=1.5,
    cexCol=0.8,
    col= redgreen(75),
    main=paste("Normalized exon read counts for ",ensembl_gene,sep=""),
    key.title=legend_name,
    key.par=list(cex=0.5),
    Rowv=FALSE,
    dendrogram="column"
  )
  dev.off()

}
