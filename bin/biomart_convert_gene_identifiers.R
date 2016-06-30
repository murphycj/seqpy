#
# Uses biomart R package to convert between gene identifiers (e.g. entrez to gene symbol)
#
# Mouse
# mmusculus_gene_ensembl
# types: mgi_symbol, ensembl_gene_id
#
# Human
# hsapiens_gene_ensembl
# types: entrezgene, hgnc_symbol
#

library("argparse")

main <- function(args) {

  library("biomaRt")

  ensembl=useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset=args$dataset,host="useast.ensembl.org")

  if (args$d=="tab") {
    data <-read.table(args$data,sep="\t",row.names=1,header=T,check.names=F)
  } else {
    data <-read.csv(args$data,row.names=1,header=T,check.names=F)
  }

  r <- getBM(attributes=c(args$outtype,args$intype),filters=c(args$intype),values=row.names(data),mart=ensembl)

  #deal with non-duplicat genes with no symbol

  r[is.na(r[,1]),1] <- ""

  empty_symbols <- !(duplicated(r[,2]) | duplicated(r[,2], fromLast=TRUE)) & (r[,1]=="")
  r[empty_symbols,1] <- r[empty_symbols,2]

  #deal with duplicate mapping of entrez id to gene symbols

  dups <- duplicated(r[,2]) | duplicated(r[,2], fromLast=TRUE)

  duplicate_genes <- unique(r[dups,2])

  rows_to_remove <- c()

  for (i in duplicate_genes) {
    temp <- r[r[,2]==i,]

    #if one all gene symbols except one is non-empty use that one

    if ((nrow(temp) - sum(temp[,1]==""))==1) {
      new_symbol <- temp[temp[,1]!="",]

      remove <- temp[temp[,1]=="",]

      rows_to_remove <- c(rows_to_remove,row.names(remove))

      r[row.names(new_symbol),1] <- new_symbol[,1]
    } else {

      #else replace gene symbols with entrez gene id

      r[row.names(temp),1] <- temp[,2]

      rows_to_remove <- c(rows_to_remove,row.names(temp)[-1])
    }
  }

  r <- r[!(row.names(r) %in% rows_to_remove),]

  #deal with duplicate gene symbols by replacing symbol with entrez

  dups <- duplicated(r[,1]) | duplicated(r[,1], fromLast=TRUE)

  r[dups,1] <- r[dups,2]

  r.not <- unique(setdiff(row.names(data),r[,2]))

  data <- data[c(r[,2],r.not),]

  #replace old names and save

  row.names(data) <- c(r[,1],r.not)

  write.csv(data,args$out,quote=F)

}


parser <- ArgumentParser()

parser$add_argument(
  "-data",
  type="character",help="The file where row names are symbol and columns are whatever."
)
parser$add_argument(
  "-d",
  type="character",
  help="Delimit character in data. Options: csv, tab"
)
parser$add_argument(
  "-intype",
  type="character",
  help="Biomart input type (e.g. entrezgene, ensembl_gene_id, hgnc_symbol)"
)
parser$add_argument(
  "-outtype",
  type="character",
  help="Biomart output type (e.g. mgi_symbol)"
)
parser$add_argument(
  "-dataset",
  type="character",
  help="Biomart dataset to use (e.g. mmusculus_gene_ensembl, hsapiens_gene_ensembl)"
)
parser$add_argument(
  "-out",
  type="character",
  help="The new output file name."
)

args <- parser$parse_args()

main(args=args)
