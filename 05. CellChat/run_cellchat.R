# Load Packages
setwd(Sys.getenv("HOME"))

suppressPackageStartupMessages(library(argparse))

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(qs))

source("R_func/mutate_load.R")
source("R_func/future_load.R")

suppressPackageStartupMessages(library(CellChat))

# load parameters
parser <- ArgumentParser()

parser$add_argument("-i", "--input", required=TRUE)
parser$add_argument("-o", "--out_dir", required=TRUE)
parser$add_argument("-grp", "--group_by", required=TRUE)
parser$add_argument("-keep", "--keep_grp", default=NULL)

parser$add_argument("-db", "--database", default=NULL)

parser$add_argument("-min", "--min_cells", default=10, type="integer")
parser$add_argument("-ppi", "--ppi", default=FALSE, action="store_true")
parser$add_argument("-type", "--type", default="triMean")
parser$add_argument("-t", "--trim", default=NULL, type="double")

args <- parser$parse_args()

dir_create(args$out_dir)

# Database
DB <- args$'database'
if (!is.null(DB)) {
	DB <- readRDS(DB)
} else {
	DB <- CellChatDB.human
}

# Keep selected groups
keep_grp <- args$'keep_grp'
if (!is.null(keep_grp)) {
	keep_grp <- read.csv(keep_grp, header=FALSE)[, 1]
	message(paste("Keeping ", length(keep_grp), "groups"))
}

# Pipeline
my.pipeline <- function(input) {
	# Load Data
	obj <- readRDS(input)
	DefaultAssay(obj) <- "RNA"
	obj %<>% DietSeurat(assays="RNA")
	foo <- gc()

	# Filter
	Idents(obj) <- args$'group_by'
	if (!is.null(keep_grp)) {
		ind <- Idents(obj) %in% keep_grp
		obj$'tmp_keep' <- ind
		obj %<>% subset(tmp_keep)
		obj$'tmp_keep' <- NULL
		message(paste("Remove", sum(!ind), "cells"))
	}
	
	my.level <- Idents(obj) %>% as.character %>% unique %>% sort
	Idents(obj) %<>% factor(level=my.level)

	# Preprocess
	chat.obj <- createCellChat(obj)
	chat.obj@DB <- DB
	chat.obj %<>% subsetData 

	# DEG
	chat.obj %<>% identifyOverExpressedGenes %>% identifyOverExpressedInteractions

	if (args$'ppi') {
		chat.obj %<>% projectData(PPI.human)
	}

	# Run
	chat.obj %<>% computeCommunProb(
		raw.use=!(args$'ppi'), type=args$'type', trim=args$'trim'
	) %>% filterCommunication(min.cells = args$'min_cells')

	return(chat.obj)
}

# Run
chat.obj <- my.pipeline(args$'input')

# postanlys
foo <- gc()
message("Run computeCentrality...")
chat.obj %<>% aggregateNet %>% netAnalysis_computeCentrality(slot.name="net")

message("Saving raw...")
chat.obj %>% qsave(path_join(c(
	args$out_dir, "chat.qs"
)), nthreads=4)

# save interaction
df.net <- subsetCommunication(chat.obj, slot.name = "net")
df.net %>% fwrite(path_join(c(args$out_dir, "chat.net.csv")))

# Save
message("Saving postanlys...")
chat.obj %>% saveRDS(path_join(c(args$out_dir, "chat.rds")))
