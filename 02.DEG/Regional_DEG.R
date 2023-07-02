# packages
setwd(Sys.getenv("HOME"))

suppressPackageStartupMessages(library(argparse))

suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(qs))

suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(fs))

source("R_func/future_load.R")

# parameter
parser <- ArgumentParser()

parser$add_argument("-r", "--region", required=TRUE)
parser$add_argument("-i", "--input", required=TRUE)
parser$add_argument("-o", "--out_dir", required=TRUE)

parser$add_argument("-grp", "--group.label", default="subclass_label")
parser$add_argument("-rgrp", "--region.label", default="Region")
parser$add_argument("-sub", "--sub.ident", default=NULL)

args <- parser$parse_args()

objF <- args$'input'
region <- args$'region'
home_dir <- args$'out_dir'
group.label <- args$'group.label'
region.label <- args$'region.label'
sub.ident <- args$'sub.ident'

# function
source("project/multiomics/GeneAnlys/script/my_runDEG.R")

# load
obj <- qread(objF, nthreads=4)
DefaultAssay(obj) <- "RNA"
obj %<>% DietSeurat(assays="RNA")
obj %<>% NormalizeData(verbose=FALSE)
foo <- gc()

# Run DA
avail_grp <- unique(as.character(obj@meta.data[, group.label]))

if (!is.null(sub.ident)) {
	if (!sub.ident %in% avail_grp) {
		stop(paste0(sub.ident, " is not available in ", group.label))
	} else {
		avail_grp <- sub.ident
	}
}

mdf_list <- list()
for (g in avail_grp) {
    message("\nRunning DEG for ", g)

	dietFun <- function(...) {my_runDEG(
		object=obj, group=g, region=region,
		all.reg = c("dlPFC", "M1C", "ACC", "V1C"),
		group.label=group.label,
		region.label=region.label,
		home_dir=home_dir, ...
	)}

    mdf_list[[g]] <- dietFun(
        min.pct = 0.1, logfc.threshold = 0.25,
		test.use="MAST", latent.vars=c("Sample")
    )
	foo <- gc()
}

# Save
if (length(mdf_list) > 0) {
	saveRDS(mdf_list, paste(home_dir, "agg.cons.markers_list.rds", sep="/"))

	# Merge
	ccols <- c(
		"gene", "cluster", "region", "sign", "pct.1",
		"max.p_val_adj", "min.p_val_adj", "max.avg_log2FC", "min.avg_log2FC",
		"max.pct.2", "min.pct.2", "max.pct.qdiff", "min.pct.qdiff"
	)

	mdf <- do.call(c, mdf_list) %>%
		lapply(function(x){x[, ccols]}) %>% bind_rows()

	data.table::fwrite(
		mdf, paste(home_dir, "agg.cons.mdf.csv", sep="/")
	)
}
