library("argparse")
source("/Users/charlesmurphy/Desktop/Research/lib/seqpy/lib/EB_20150806_Functions_TIDE_version201.R")

SGSEQ_DATABASE <- list()
SGSEQ_DATABASE[["9707117"]] = "CTTGTGCTGACTTACCAGAT"
SGSEQ_DATABASE[["9707115"]] = "AAATCTTAGAGTGTCCCATC"
SGSEQ_DATABASE[["9681514"]] = "ACAGATTGTATATCTTGTAA"
SGSEQ_DATABASE[["51049"]] = "CCCCGGACGATATTGAACAA"
SGSEQ_DATABASE[["51048"]] = "CGCTATCTGAGCAGCGCTCA"
SGSEQ_DATABASE[["51047"]] = "CCGGTTCATGCCGCCCATGC"
SGSEQ_DATABASE[["04760"]] = "TGCTAGTCTGGAGTTGATCA"
SGSEQ_DATABASE[["9681583"]] = "ACCGCCAAATTTAATTGCAG"

WT_DATABASE = list()
WT_DATABASE[["9707117"]] = "/Users/charlesmurphy/Desktop/Research/sequential_2015_6/data/genewiz/1.21.16/10-321436138_ab1/9707117-wt-137723972.ab1"
WT_DATABASE[["9707115"]] = "/Users/charlesmurphy/Desktop/Research/sequential_2015_6/data/genewiz/12.23.15/10-319256333_ab1/sgBRCA1-9707115-WT-137723969.ab1"
WT_DATABASE[["9681514"]] = "/Users/charlesmurphy/Desktop/Research/sequential_2015_6/data/genewiz/12.23.15/10-319256333_ab1/sgPTEN-9681514-WT-137723964.ab1"
WT_DATABASE[["51049"]] = "/Users/charlesmurphy/Desktop/Research/sequential_2015_6/data/genewiz/12.23.15/10-319256333_ab1/sgP53-51049-WT-137723986.ab1"
WT_DATABASE[["51048"]] = "/Users/charlesmurphy/Desktop/Research/sequential_2015_6/data/genewiz/12.23.15/10-319256333_ab1/sgP53-51048-WT-137723980.ab1"
WT_DATABASE[["51047"]] = "/Users/charlesmurphy/Desktop/Research/sequential_2015_6/data/genewiz/3.31.16/10-329071847_ab1/sg51047-WT-mcf10a-143975993.ab1"
WT_DATABASE[["9681583"]] = "/Users/charlesmurphy/Desktop/Research/sequential_2015_6/data/genewiz/12.23.15/10-319256333_ab1/sgPTEN-9681583-WT-139932431.ab1"
WT_DATABASE[["04760"]] = "/Users/charlesmurphy/Desktop/Research/sequential_2015_6/data/genewiz/12.23.15/10-319256333_ab1/sgBRCA1-04760-WT-140094010.ab1"

run_TIDE_function <- function(control_file,sg_sequence,rg1,rg2,args) {
  #
  # run the TIDE coad and save results
  #

  result=import(
    control_file=control_file,
    experimental_file=args$test,
    guide=sg_sequence,
    seqstart=args$seqstart,
    seqend=args$seqend,
    maxshift=args$maxshift,
    rg1=rg1,
    rg2=rg2
  )

  png(paste(args$out,"/",args$prefix,"-decomposition.png",sep=""),width=1200,height=600)
  d=decomposition(
    import=result,
    p.threshold=args$pthreshold,
    fout=paste(args$out,"/",args$prefix,"-decomposition.txt",sep="")
  )
  dev.off()

  png(paste(args$out,"/",args$prefix,"-quality.png",sep=""),width=1200,height=600)
  quality(
    import=result,
    fout=paste(args$out,"/",args$prefix,"-quality.txt",sep="")
  )
  dev.off()
}

main <- function(args) {

  dir.create(args$out,showWarnings=F)

  control_trace <- ""
  sg_sequence <- ""
  rg1 = NA
  rg2 = NA

  if (args$rg1==0) {
    rg1 <- NA
  } else {
    rg1 <- args$rg1
  }

  if (args$rg2==0) {
    rg2 <- NA
  } else {
    rg2 <- args$rg2
  }

  if (args$control=="") {
    control_trace <- WT_DATABASE[[args$guide]]
    sg_sequence <- SGSEQ_DATABASE[[args$guide]]
  } else {
    control_trace <- args$control
    sg_sequence <- args$sgseq
  }

  run_TIDE_function(
    control_file=control_trace,
    sg_sequence=sg_sequence,
    rg1=rg1,
    rg2=rg2,
    args=args
  )
}


parser <- ArgumentParser()

parser$add_argument("-test",type="character",help="The test trace file.")
parser$add_argument("-control",type="character",help="The control file (optional).",default="")
parser$add_argument("-prefix",type="character",help="Prefix of output files.",default="")
parser$add_argument("-out",type="character",help="Output dir.",default="")
parser$add_argument("-guide",type="character",help="The sgRNA code.")
parser$add_argument("-sgseq",type="character",help="The sgRNA sequence (optional).",default="")
parser$add_argument("-gene",type="character",help="The gene.")
parser$add_argument("-seqstart",type="integer",help="start of sequence read from where data will be included (because beginning of seq reads tends to be poor quality, default 100)",default=100)
parser$add_argument("-seqend",type="integer",help="last bp to be included in analysis (will be automatically adjusted if reads are shorter, see below, default 700)",default=700)
parser$add_argument("-maxshift",type="integer",help="range of basepair shifts (indels) to be analyzed, both positive and negative (defualt 10)",default=10)
parser$add_argument("-rg1",type="integer",help="[optional] the first (rg1) and last (rg2) base of the sequence region that is used for decomposition; will be automatically set if NA (default 0)",default=0)
parser$add_argument("-rg2",type="integer",help="[optional] the first (rg1) and last (rg2) base of the sequence region that is used for decomposition; will be automatically set if NA (default 0)",default=0)
parser$add_argument("-pthreshold",type="double",help="p.threshold (defualt 0.01)",default=0.01)

args <- parser$parse_args()

main(args=args)
