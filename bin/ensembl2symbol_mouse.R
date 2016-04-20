library("argparse")
library("biomaRt")
ensembl=useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl",host="www.ensembl.org")

parser <- ArgumentParser()

parser$add_argument("-fpkm",type="character",help="The FPKM file")
parser$add_argument("-out",type="character",help="The new output FPKM file")

args <- parser$parse_args()


fpkm <- read.csv(args$fpkm,row.names=1)

r <- getBM(attributes=c("mgi_symbol","ensembl_gene_id"),filters=c("ensembl_gene_id"),values=row.names(fpkm),mart=ensembl)

#deal with duplicate gene names

dup <- duplicated(r[,1]) | duplicated(r[,1], fromLast=TRUE)
r[dup,1] <- r[dup,2]
r <- r[!duplicated(r[,1]),]

dup <- duplicated(r[,2]) | duplicated(r[,2], fromLast=TRUE)
r <- r[!dup,]
#r[dup,2] <- r[dup,1]

r.not <- unique(setdiff(row.names(fpkm),r[,2]))

fpkm <- fpkm[c(r[,2],r.not),]

#replace old names and save
row.names(fpkm) <- c(r[,1],r.not)

write.csv(fpkm,args$out)
