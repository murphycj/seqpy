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

  outs = c(args$outtype,args$intype)

  if (!args$no_external_gene_name) {
    outs <- c(outs,"external_gene_name")
  }

  if (args$description) {
    outs <- c(outs,"description")
  }

  r <- getBM(
    attributes=outs,
    filters=c(args$intype),
    values=row.names(data),
    mart=ensembl
  )

  r[is.na(r[,args$outtype]),args$outtype] <- ""

  #for those without a primary output, add the secondary output

  if (!args$no_external_gene_name) {
    r[r[,args$outtype]=="",args$outtype] <- r[r[,args$outtype]=="","external_gene_name"]
  }

  #deal with non-duplicat genes with no symbol

  empty_symbols <- !(duplicated(r[,args$intype]) | duplicated(r[,args$intype], fromLast=TRUE)) & (r[,args$outtype]=="")
  r[empty_symbols,args$outtype] <- r[empty_symbols,args$intype]

  #deal with duplicate mapping of entrez id to gene symbols

  dups <- duplicated(r[,args$intype]) | duplicated(r[,args$intype], fromLast=TRUE)

  duplicate_genes <- unique(r[dups,args$intype])

  rows_to_remove <- c()

  for (i in duplicate_genes) {
    temp <- r[r[,args$intype]==i,]

    #if one all gene symbols except one is non-empty use that one

    if ((nrow(temp) - sum(temp[,args$outtype]==""))==1) {
      new_symbol <- temp[temp[,args$outtype]!="",]

      remove <- temp[temp[,args$outtype]=="",]

      rows_to_remove <- c(rows_to_remove,row.names(remove))

      r[row.names(new_symbol),args$outtype] <- new_symbol[,args$outtype]
    } else {

      #else replace gene symbols with entrez gene id

      r[row.names(temp),args$outtype] <- temp[,args$intype]

      rows_to_remove <- c(rows_to_remove,row.names(temp)[-1])
    }
  }

  r <- r[!(row.names(r) %in% rows_to_remove),]

  #deal with duplicate gene symbols by replacing symbol with entrez

  dups <- duplicated(r[,args$outtype]) | duplicated(r[,args$outtype], fromLast=TRUE)

  r[dups,args$outtype] <- r[dups,args$intype]

  r.not <- unique(setdiff(row.names(data),r[,args$intype]))

  data <- data[c(r[,args$intype],r.not),]

  #replace old names and save

  row.names(data) <- c(r[,args$outtype],r.not)

  #add description if needed

  if (args$description) {

    column_names <- colnames(data)

    r[,"description"] <- gsub(",",";",r[,"description"])

    data["description"] <- c(r[,"description"],rep("",length(r.not)))

    column_names <- c("description",column_names)

    data <- data[,column_names]
  }

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
  help="Biomart output type (e.g. mgi_symbol, hgnc_symbol)"
)
parser$add_argument(
  "-no_external_gene_name",
  action="store_true",
  default=FALSE,
  help="For ensembl ids without gene symbol do not use external_gene_name (default false)."
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
parser$add_argument(
  "-description",
  action="store_true",
  default=FALSE,
  help="Add gene descriptions (default false)"
)

args <- parser$parse_args()

main(args=args)
