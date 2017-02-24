

library("argparse")

main <- function(args) {

  library("DESeq2")
  library("RColorBrewer")
  library("gplots")
  library("genefilter")

  dir.create(args$outDir, showWarnings = FALSE)

  group1 <- strsplit(args$G1,",")[[1]]
  group2 <- strsplit(args$G2,",")[[1]]
  phenotype <- strsplit(args$phenotype,",")[[1]]

  data <- read.csv(args$counts, header=T, row.names=1, check.names=FALSE)
  data <- data[,c(group1, group2)]
  data[data<args$minCount] <- 0
  data <- data[rowSums(data)>=args$minRowCount,]

  coldata <- data.frame(
    mainFactor=factor(c(
      rep(phenotype[1],length(group1)),
      rep(phenotype[2],length(group2))
      ),
      levels=phenotype
    ),
    row.names=colnames(data)
  )
  countTable <- DESeqDataSetFromMatrix(
    countData=data,
    colData=coldata,
    design=~mainFactor
  )
  result <- DESeq(countTable)
  res <- results(result)
  dd = res[with(res,order(padj)),]

	write.csv(
    dd,
    paste(args$outDir,"/",args$outDir,"_results.csv",sep="")
  )

  pdf(paste(args$outDir,"/",args$outDir,"-pvalue-hist.pdf",sep=""))
	hist(res$pvalue, breaks=20, col="grey")
	dev.off()

  pdf(paste(args$outDir,"/",args$outDir,"-padj-hist.pdf",sep=""))
	hist(res$padj, breaks=20, col="grey")
	dev.off()

  pdf(paste(args$outDir,"/",args$outDir,"-MA-plot.pdf",sep=""))
	plotMA(res)
	dev.off()
}



parser <- ArgumentParser()

parser$add_argument(
  "-counts",
  type="character",
  help="The file where row names are symbol and columns are whatever."
)
parser$add_argument(
  "-G1",
  type="character",
  help="Comma-separate list of sample names for group1."
)
parser$add_argument(
  "-G2",
  type="character",help="Comma-separate list of sample names for group2."
)
parser$add_argument(
  "-phenotype",
  type="character",help="The file where row names are symbol and columns are whatever."
)
parser$add_argument(
  "-minCount",
  type="integer",
  help="Minimum count a gene should have to be counted"
)
parser$add_argument(
  "-minRowCount",
  type="integer",
  help="Minimum sum of the row (gene) counts."
)
parser$add_argument(
  "-outDir",
  type="character",
  help="The new output file name."
)

args <- parser$parse_args()

main(args=args)
