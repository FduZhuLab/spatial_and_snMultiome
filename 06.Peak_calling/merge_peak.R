# create parser object
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()

parser$add_argument("-i", "--input", required=TRUE, help="list of summit bed")
parser$add_argument("-g", "--genome", default="hg38", help="genome")
parser$add_argument("-blacklist", "--blacklist", required=TRUE, help="blacklist in bed format")
parser$add_argument("-chromSize", "--chromSize", required=TRUE, help="chromsomes and their size in tab-separated format")
parser$add_argument("-o", "--out_dir", required=TRUE, help="output dir")
parser$add_argument("-t", "--threads", default=1, help="number of threads to use")

args <- parser$parse_args()

# load packages
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("BSgenome"))
suppressPackageStartupMessages(library("data.table"))

suppressPackageStartupMessages(library("future"))
suppressPackageStartupMessages(library("future.apply"))
options(future.rng.onMisuse = "ignore")
options(future.globals.maxSize = 1e11)
plan(multicore, workers = as.integer(args$threads))

# prepare input
inF <- args$input
outDir <- args$out_dir

blacklist <- read.table(args$blacklist)
blacklist.gr <- GRanges(blacklist[,1], IRanges(blacklist[,2], blacklist[,3]))

chrom <- read.table(args$chromSize)
chrom.gr <- GRanges(chrom[,1], IRanges(0, chrom[,2]))

genome <- getBSgenome(args$genome)

##########################
# Functions 

# read bed to gr
read2gr <- function(bedF, label){
  df <- as.data.frame(fread(bedF, sep="\t", header=FALSE, nThread=1))
  colnames(df) <- c("chr", "start", "end", "name", "score")
  df$label <- label
  gr <- GRanges(df$chr, IRanges(df$start, df$end))
  mcols(gr)$score <- df$score
  mcols(gr)$name <- df$name
  mcols(gr)$label <- df$label
  return(gr)
}

# extend summit to 500 bp
extendSummit <- function(gr, size=500){
  gr <- resize(gr, width=size, fix="center")
  return(gr)
}

# keep only chosen chromosomes
filter4chrom <- function(gr, chrom.gr){  
  idx <- queryHits(findOverlaps(gr, chrom.gr, type="within"))
  gr <- gr[idx]
  return(gr)
}

# filter blacklist regions
filter4blacklist <- function(gr, blacklist.gr){
  idx <- queryHits(findOverlaps(gr, blacklist.gr))
  if(length(idx) > 0){
    gr <- gr[-idx]
  }
  return(gr)
}

# filter N containing regions
filter4N <- function(gr, genome){
  nucFreq <- BSgenome::alphabetFrequency(getSeq(genome, gr))
  mcols(gr)$GC <- round(rowSums(nucFreq[,c("G","C")]) / rowSums(nucFreq),4)
  mcols(gr)$N <- round(nucFreq[,c("N")] / rowSums(nucFreq),4)
  idx <- which(mcols(gr)$N >= 0.001)
  if(length(idx) > 0){
    gr <- gr[-idx]
  }
  return(gr)
}

# get non-overlapped regions
#' Retreive a non-overlapping set of regions from a Genomic Ranges object
#'
#' This function returns a GRanges object containing a non-overlapping set regions derived from a supplied Genomic Ranges object.
#'
#' @param gr A `GRanges` object.
#' @param by The name of a column in `mcols(gr)` that should be used to determine how overlapping regions should be resolved.
#' The resolution of overlapping regions also depends on `decreasing`. For example, if a column named "score" is used for `by`,
#' `decreasing = TRUE` means that the highest "score" in the overlap will be retained and `decreasing = FALSE` means that the
#' lowest "score" in the overlap will be retained.
#' @param decreasing A boolean value indicating whether the values in the column indicated via `by` should be ordered in decreasing
#' order. If `TRUE`, the higher value in `by` will be retained.
#' @param verbose A boolean value indicating whether the output should include extra reporting.
#' @export
nonOverlappingGR <- function(
    gr = NULL, 
    by = "score", 
    decreasing = TRUE, 
    verbose = FALSE
  ){

  stopifnot(by %in% colnames(mcols(gr)))

  #-----------
  # Cluster GRanges into islands using reduce and then select based on input
  #-----------
  .clusterGRanges <- function(gr = NULL, by = "score", decreasing = TRUE){
    gr <- sort(sortSeqlevels(gr))
    r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
    o <- findOverlaps(gr, r, ignore.strand = TRUE)
    mcols(gr)$cluster <- subjectHits(o)
    gr <- gr[order(mcols(gr)[, by], decreasing = decreasing),]
    gr <- gr[!duplicated(mcols(gr)$cluster),]
    gr <- sort(sortSeqlevels(gr))
    mcols(gr)$cluster <- NULL
    return(gr)
  }

  if(verbose){
    message("Converging", appendLF = FALSE)
  }
  i <-  0
  grConverge <- gr
  while(length(grConverge) > 0){
    if(verbose){
      message(".", appendLF = FALSE)
    }
    i <-  i + 1
    grSelect <- .clusterGRanges(
      gr = grConverge, by = by, decreasing = decreasing
    )

    grConverge <- subsetByOverlaps(
      grConverge, grSelect, invert=TRUE, ignore.strand = TRUE
    ) #blacklist selected gr
    
    if(i == 1){ #if i=1 then set gr_all to clustered
      grAll <- grSelect
    }else{
      grAll <- c(grAll, grSelect)
    } 

  }
  message(sprintf("Converged after %s iterations!", i))

  if(verbose){
    message("\nSelected ", length(grAll), " from ", length(gr))
  }
  grAll <- sort(sortSeqlevels(grAll))

  return(grAll)
}

# normlize to score per million (SPM)
norm2spm <- function(gr, by = "score"){
  mlogp <- mcols(gr)[,by]
  normmlogp <- 10^6 * mlogp / sum(mlogp)
  mcols(gr)$spm <- normmlogp
  return(gr)
}

################################

# load summit
summitF <- read.table(inF, sep="\t", header=FALSE)
label.lst <- as.character(summitF$V1)
file.lst <- as.character(summitF$V2)

message("parse summit peak set")
peak.list <- future_lapply(seq(file.lst), function(i){
  message("working on summit set for ", label.lst[i])
  p.gr <- read2gr(file.lst[i], label=label.lst[i]) %>%
    extendSummit(size=500) %>%
    filter4blacklist(blacklist.gr) %>%
    filter4chrom(chrom.gr) %>%
    filter4N(genome) %>%
    nonOverlappingGR(by = "score", decreasing = TRUE) %>%
    norm2spm(by="score")
  return(p.gr)
})

message("write filtered & fixed peak set")
for(i in seq(file.lst)){
  message("write peak set for ", label.lst[i])
  outPeak <- as.data.frame(peak.list[[i]])
  outfname <- paste0(outDir, "/", label.lst[i], ".filterNfixed.peakset")
  write.table(outPeak, file=outfname, sep="\t", quote = F, col.names = T, row.names = F)
}

message("merge to union peak list")
merged.gr <- do.call(c, peak.list) %>% nonOverlappingGR(by = "spm", decreasing = TRUE)

message("write union peak set")
outUnion <- as.data.frame(merged.gr)
outfname <- paste0(outDir, "/", "filteredNfixed.union.peakSet")
write.table(outUnion, file=outfname, sep="\t", quote = F, col.names = T, row.names = F)
