# packages
setwd(Sys.getenv("HOME"))

suppressPackageStartupMessages(library(argparse))

source("R_func/mutate_load.R")
source("R_func/future_load.R")

source("R_func/seurat_load.R")
library(SeuratDisk)

library(MetaNeighbor)

# parameter
parser <- ArgumentParser()

parser$add_argument("-obj_list", "--obj_list", required=TRUE, nargs="+")
parser$add_argument("-out", "--out_dir", required=TRUE)
parser$add_argument("-grp", "--group_by", required=TRUE)
parser$add_argument("-dataset_label", "--dataset_label", default=NULL)

parser$add_argument("-min", "--min_cells", default=5, type="integer")

parser$add_argument("-assay", "--assay", default="RNA")

parser$add_argument("-subset", "--subset", default=NULL)
parser$add_argument("-subidents", "--subidents", default=NULL, nargs="+")

parser$add_argument("-num_hvg", "--num_hvg", default=5000, type="integer")
parser$add_argument("-seurat_hvg", "--seurat_hvg", action="store_true", default=FALSE)
parser$add_argument("-marker_hvg", "--marker_hvg", action="store_true", default=FALSE)
parser$add_argument("-only_share", "--only_share", action="store_true", default=FALSE)

parser$add_argument("-match_size", "--match_size", action="store_true", default=FALSE)
parser$add_argument("-spl_size", "--spl_size", default=10000, type="integer")

parser$add_argument("-one_vs_best", "--one_vs_best", action="store_true", default=FALSE)
parser$add_argument("-fast_version", "--fast_version", action="store_true", default=FALSE)


# arguments
args <- parser$parse_args()

assay <- args$'assay'

subset_label <- args$'subset'
subset_idents <- args$'subidents' %>% unlist

out_dir <- args$'out_dir'
dir_create(out_dir)

group_by <- args$'group_by'

# Load
obj_list <- lapply(args$'obj_list', qread, nthreads=4)

## set assay
message(paste("Setting default assay to", assay))
my.diet <- function(object) {
	DefaultAssay(object) <- assay
	object %>% DietSeurat(assays=assay)
}

obj_list %<>% lapply(my.diet) %>% lapply(NormalizeData)

## subset dataset
my.subset <- function(object, by, idents) {
	object$'tmp.keep' <- object@meta.data[, by] %in% idents
	object %<>% subset(tmp.keep)
	object$'tmp.keep' <- NULL
	return(object)
}

if (!is.null(subset_label)) {
	message(paste("Subsetting all datasets by", subset_label, "to", paste(subset_idents, collapse=" ")))
	obj_list %<>% lapply(
		my.subset, by=subset_label, idents=subset_idents
	)
}

# match size
source("R_func/seurat/Integrate/my.sketchObject.R")

if (args$'match_size') {
	size_list <- obj_list %>% sapply(ncol)
	spl_size <- args$'spl_size'
	message(glue("Downsampling to {spl_size} cells"))
	obj_list %<>% lapply(my.sketchObject, sketch.size=spl_size)
}

# HVG
avail.features <- unique(Reduce(intersect, lapply(obj_list, rownames)))
message(glue("\nThere are {length(avail.features)} avail.features"))

source("R_func/seurat/anlys_pipeline/my.anlys_hvg.R")
raw_method <- FALSE
if (args$'seurat_hvg') {
	message("\nUsing Seurat hvg...")
	obj_list %<>% lapply(
		my.anlys_hvg, hvg_type='hvg', num_hvg=args$'num_hvg'
	)
	source("R_func/seurat/Select.Covar.R")
	# var_genes <- SelectIntegrationFeatures(obj_list, nfeatures=args$'num_hvg')
	cohvg.df <- Select.Covar(obj_list, method='vst')
	var_genes <- rownames(cohvg.df) %>% head(n=args$'num_hvg')
} else if (args$'marker_hvg') {
	message("\nUsing marker hvg...")
	obj_list %<>% lapply(
		my.anlys_hvg, hvg_type=group_by, features=avail.features,
		min.log2fc=0.1, p.thres=0.01, only.pos=FALSE
	)
	var_genes <- obj_list %>% lapply(VariableFeatures) %>% unlist %>% unique
} else {
	raw_method <- TRUE
}

# merge data
if (length(obj_list) == 1) {
	mobj <- obj_list[[1]]
} else {
	if (is.null(args$'dataset_label')) {
		for (i in seq(obj_list)) {
			obj_list[[i]]$'dataset_label' <- paste0("D", i)
		}
	}
	mobj <- merge(x=obj_list[[1]], y=obj_list[2:length(obj_list)])
}

if (!is.null(args$'dataset_label')) {
	mobj$'dataset_label' <- mobj@meta.data[[args$'dataset_label']]
}

if (args$'only_share') {
	mobj <- CreateSeuratObject(
		counts = mobj$RNA@counts[avail.features, ],
		data = mobj$RNA@data[avail.features, ],
		meta.data = mobj@meta.data
	)
}

# remove small groups
clu_cell_num <- mobj@meta.data[[group_by]] %>% table
keep.clusters <- names(clu_cell_num)[clu_cell_num >= args$'min_cells']
ind <- mobj@meta.data[[group_by]] %in% keep.clusters
mobj <- mobj[, ind]

# Format
message("\nHere are the merged object: ")
print(mobj)

message("\nFormat: ")
merge.data <- as.SingleCellExperiment(mobj)
print(merge.data)

rm(obj_list, mobj)
foo <- gc()

# Var genes
if (raw_method) {
	message("\nUsing raw hvg...")
	var_genes <- variableGenes(
		dat = merge.data, exp_labels = merge.data$dataset_label
	)
}

message(glue("\nThere are {length(var_genes)} var genes"))

# run MetaUS
message("\nRunning MetaUS:")
celltype_NV <- MetaNeighborUS(
    var_genes = var_genes,
    dat = merge.data,
    study_id = merge.data$dataset_label,
    cell_type = merge.data@colData[[group_by]],
    one_vs_best = args$'one_vs_best',
    fast_version = args$'fast_version'
)


## save object
celltype_NV %>% write.csv(
	glue("{out_dir}/{group_by}_NV_mat.csv")
)

