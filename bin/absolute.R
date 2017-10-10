

library("argparse")

main <- function(args) {

  #options(warn=1)

  library(ABSOLUTE)

  dir.create(args$outDir, showWarnings = FALSE)

  tmp <- read.table(args$cns,sep="\t",header=T)
  tmp2 <- tmp[,c("chromosome","start","end","probes","log2")]
  colnames(tmp2) <- c("Chromosome","Start","End","Num_Probes","Segment_Mean")
  tmp <- as.character(tmp2$Chromosome)
  tmp[tmp=="X"] <- "20"
  tmp[tmp=="Y"] <- "21"
  tmp2$Chromosome<-tmp
  write.table(tmp2,file.path(args$outDir,"tmp.tsv"),sep="\t",quote=F,row.names=F)

  RunAbsolute(
    seg.dat.fn=file.path(args$outDir,"tmp.tsv"),
    max.sigma.h=0.02,
    min.ploidy=0.5,
    max.ploidy=3.5,
    platform="SNP_6.0",
    copy_num_type="total",
    sigma.p=0.02,
    results.dir=args$outDir,
    primary.disease="BRCA",
    sample.name=args$sample,
    max.as.seg.count=1500,
    max.non.clonal=0,
    max.neg.genome=0
  )

  obj.name <- "DRAWS_summary"
  results.dir <- file.path(args$outDir, "output", "abs_summary")
  absolute.file <- file.path(args$outDir,paste(args$sample,".ABSOLUTE.RData",sep=""))

  CreateReviewObject(obj.name, absolute.file, results.dir, "total", verbose=TRUE)

  calls.path = file.path(args$outDir, "output", "abs_summary", "DRAWS_summary.PP-calls_tab.txt")
  modes.path = file.path(args$outDir, "output", "abs_summary", "DRAWS_summary.PP-modes.data.RData")
  output.path = file.path(args$outDir, "output", "abs_extract")

  ExtractReviewedResults(calls.path, "test", modes.path, output.path, "absolute", "total")

}

parser <- ArgumentParser()

parser$add_argument(
  "-cns",
  type="character",
  help=".cns file"
)
parser$add_argument(
  "-sample",
  type="character",help="Sample name."
)
parser$add_argument(
  "-outDir",
  type="character",
  help="The new output file name."
)

args <- parser$parse_args()

main(args=args)
