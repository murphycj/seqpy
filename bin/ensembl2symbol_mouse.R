library("argparse")
library("biomaRt")
ensembl=useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl",host="useast.ensembl.org")

parser <- ArgumentParser()

parser$add_argument("-data",type="character",help="The csv file where row names are ensembl genes and columns are whatever.")
parser$add_argument("-out",type="character",help="The new output file name.")

args <- parser$parse_args()


data <- read.csv(args$data,row.names=1)

r <- getBM(attributes=c("mgi_symbol","ensembl_gene_id"),filters=c("ensembl_gene_id"),values=row.names(data),mart=ensembl)

#deal with duplicate gene names

dup <- duplicated(r[,1]) | duplicated(r[,1], fromLast=TRUE)
r[dup,1] <- r[dup,2]
r <- r[!duplicated(r[,1]),]

dup <- duplicated(r[,2]) | duplicated(r[,2], fromLast=TRUE)
r <- r[!dup,]
#r[dup,2] <- r[dup,1]

r.not <- unique(setdiff(row.names(data),r[,2]))

data <- data[c(r[,2],r.not),]

#replace old names and save
row.names(data) <- c(r[,1],r.not)

write.csv(data,args$out)
