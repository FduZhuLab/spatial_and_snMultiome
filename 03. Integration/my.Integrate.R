# packages
setwd(Sys.getenv("HOME"))

suppressPackageStartupMessages(library(argparse))

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(qs))

source("R_func/mutate_load.R")
source("R_func/future_load.R")
suppressPackageStartupMessages(library(data.table))

# parameter
parser <- ArgumentParser()

parser$add_argument("-obj_list", "--obj_list", required=TRUE, nargs="+")

parser$add_argument("-out", "--out_dir", required=TRUE)
parser$add_argument("-obj_name", "--obj_name", default="intg.obj.qs")

parser$add_argument("-assay", "--assay", default="RNA")

parser$add_argument("-subset", "--subset", default=NULL)
parser$add_argument("-subidents", "--subidents", default=NULL, nargs="+")
parser$add_argument("-split_by", "--split_by", default=NULL)

parser$add_argument("-anchor_type", "--anchor_type", default='hvg')
parser$add_argument("-anchor_file", "--anchor_file", default=NULL)
parser$add_argument("-num_hvg", "--num_hvg", default=3000, type="integer")
parser$add_argument("-new_select", "--new_select", action="store_true", default=FALSE)
# parser$add_argument("-only_share", "--only_share", action="store_true", default=FALSE)

parser$add_argument("-anlys_npcs", "--anlys_npcs", default=50, type="integer")
parser$add_argument("-intg_npcs", "--intg_npcs", default=50, type="integer")

parser$add_argument("-do_sct", "--do_sct", action="store_true", default=FALSE)

parser$add_argument("-regress_out", "--regress_out", default=NULL)

parser$add_argument("-ref_id", "--ref_id", default=NULL, type="integer")
parser$add_argument("-use_ref", "--use_ref", action="store_true", default=FALSE,
					help="use reference object in FindIntegrationAnchors")

parser$add_argument("-tf_reduc", "--tf_reduc", default='pcaproject')

parser$add_argument("-match_size", "--match_size", action="store_true", default=FALSE)
parser$add_argument("-spl_size", "--spl_size", default=10000, type="integer")
parser$add_argument("-match_group", "--match_group", default=NULL)
parser$add_argument("-match_group_lst", "--match_group_lst", default=30, type="integer")

# arguments
args <- parser$parse_args()

## subset
split_by <- args$'split_by'
subset_label <- args$'subset'
subset_idents <- args$'subidents' %>% unlist

## anlys
assay <- args$'assay'
do_sct <- args$'do_sct'
anlys.npcs <- args$'anlys_npcs'
intg.npcs <- args$'intg_npcs'

## hvg
num_hvg <- args$'num_hvg'

anchor_type <- args$'anchor_type'
if (anchor_type == "hvg") {
	if (do_sct) {
		hvg.method <- "SCT"
	} else {
		hvg.method <- "vst"
	}
}
anchor_file <- args$'anchor_file'

## which is reference
ref_id <- args$'ref_id'
ref_intg_id <- NULL
if (args$'use_ref') {
	ref_intg_id <- ref_id
}

## save
out_dir <- args$'out_dir'
prc_dir <- paste(out_dir, "process", sep="/")
dir_create(out_dir)
dir_create(prc_dir)

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

## split dataset
if (!is.null(split_by)) {
	message(paste("Spliting all datasets by", split_by))
	obj_list %<>% lapply(
		SplitObject, split.by=split_by
	) %>% unlist
}

foo <- gc()

# match size
source("R_func/seurat/Integrate/my.sketchObject.R")

if (args$'match_size') {
	size_list <- obj_list %>% sapply(ncol)
	spl_size <- args$'spl_size'
	message(glue("Downsampling to {spl_size} cells"))
	obj_list %<>% lapply(
		my.sketchObject, sketch.size=spl_size,
		group.by=args$'match_group',
		lst_num=args$'match_group_lst'
	)
}

# add dataset_id
for (i in 1:length(obj_list)) {
	obj_list[[i]]$'dataset_id' <- i
}

# only shared genes
avail.features <- unique(Reduce(intersect, lapply(obj_list, rownames)))

# if (args$'only_share') {
# 	mobj <- merge(x=obj_list[[1]], y=obj_list[2:length(obj_list)])
# 	mobj <- CreateSeuratObject(
# 		counts = mobj[[assay]]@counts[avail.features, ],
# 		data = mobj[[assay]]@data[avail.features, ],
# 		meta.data = mobj@meta.data
# 	)
# 	obj_list <- SplitObject(mobj, 'dataset_id')

# 	rm(mobj)
# 	foo <- gc()
# }

# anchor features
if (!is.null(anchor_file)) {
	anchor.features <- read.csv(anchor_file)[, 1] %>% unique
	anchor.features <- anchor.features[anchor.features %in% avail.features]
} else {
	source("R_func/seurat/anlys_pipeline/my.anlys_hvg.R")
	obj_list %<>% lapply(
		my.anlys_hvg,
		hvg_type = anchor_type,
		num_hvg = num_hvg,
		do_sct = do_sct,
		min.log2fc=0.1, p.thres=0.05,
		features=avail.features,
		vars.to.regress=args$'regress_out'
	)

	if (anchor_type == "hvg") {
		if (args$'new_select') {
			message("Using new_select method to select hvg from shared gene sets")
			source("R_func/seurat/Select.Covar.R")
			# cohvg.df <- Select.Covar(obj_list, method=hvg.method)
			# cohvg.df %>% write.csv(path_join(c(prc_dir, "cohvg.df.csv")))
			# anchor.features <- rownames(cohvg.df) %>% head(n=num_hvg)

			obj_list %<>% lapply(
				my.find.hvg, method=hvg.method, avail.features=avail.features, nfeatures=num_hvg
			)
			anchor.features <- SelectIntegrationFeatures(obj_list, nfeatures=num_hvg)
		} else {
			anchor.features <- SelectIntegrationFeatures(obj_list, nfeatures=num_hvg)
		}
	} else {
		anchor.features <- obj_list %>% lapply(VariableFeatures) %>% unlist %>% unique
		anchor.features <- anchor.features[anchor.features %in% avail.features]
	}
}

message(paste("\nThere are", length(anchor.features), "anchor.features"))

anchor.features %>% write.csv(
	path_join(c(prc_dir, "anchor.features.csv")), row.names=FALSE
)

# obj_list %<>% lapply(function(object){
# 	VariableFeatures(object) <- anchor.features
# 	return(object)
# })

foo <- gc()

# Integrated
message("\nHere are all the objects for integration: ")
print(obj_list)

set.seed(42)

if (do_sct) {
	message("\nIntegrating with SCT...\n")
	if (anchor_type != "hvg") {
		message("\nRunning SCT with provided anchor.features...\n")
		obj_list %<>% lapply(
			SCTransform, vst.flavor="v2",
			residual.features=anchor.features,
			vars.to.regress=args$'regress_out'
		)
	}

	obj_list %<>% lapply(function(object) {
		DefaultAssay(object) <- "SCT"
		return(object)
	})

	obj_list %<>% PrepSCTIntegration(
		anchor.features = anchor.features
	)

	norm.method <- "SCT"
} else {
	message("\nIntegrating with Default method...\n")
	norm.method <- "LogNormalize"
}

intg.anchors <- FindIntegrationAnchors(
	object.list = obj_list,
	anchor.features = anchor.features,
	normalization.method = norm.method,
	reduction = "cca", dims=1:intg.npcs,
	reference=ref_intg_id
)

message("\nSave anchors...\n")
intg.anchors %>% qsave(
	paste(prc_dir, "intg.anchors.qs", sep="/"), nthreads=4
)

intg.obj <- IntegrateData(
	intg.anchors,
	new.assay.name="intg", dims=1:intg.npcs,
	normalization.method = norm.method
)

## anlys
message("\nAnlys intg...\n")

if (!do_sct) {
	intg.obj %<>% ScaleData(vars.to.regress=args$'regress_out')
}

intg.obj %<>% RunPCA(npcs=anlys.npcs, verbose = FALSE) %>%
    RunUMAP(reduction="pca", dims = 1:anlys.npcs, verbose=FALSE)

## Save embeddings
Embeddings(intg.obj, "pca") %>% write.csv(paste(prc_dir, "pca.csv", sep="/"))
Embeddings(intg.obj, "umap") %>% write.csv(paste(prc_dir, "umap.csv", sep="/"))

## diet seurat
intg.obj[["pca"]]@assay.used <- assay
intg.obj[["umap"]]@assay.used <- assay
DefaultAssay(intg.obj) <- assay
intg.obj %<>% DietSeurat(assays=assay, dimreducs=c('pca', 'umap'))
VariableFeatures(intg.obj) <- anchor.features

## save object
message("\nSave object...\n")
intg.obj %>% qsave(path_join(c(out_dir, args$'obj_name')), nthreads=4)

# Transfer
if (!is.null(ref_id)) {
	message("\nTransfering...\n")

	nl <- length(obj_list)
	ref_id <- min(ref_id, nl)
	qry_ids <- 1:nl
	qry_ids <- qry_ids[qry_ids != ref_id]

	set.seed(42)
	tf.anchors.list <- lapply(obj_list[qry_ids], function(object) {
		FindTransferAnchors(
			reference=obj_list[[ref_id]], query=object,
			features=anchor.features,
			reduction=args$'tf_reduc',
			npcs=intg.npcs, dims=1:intg.npcs,
			mapping.score.k=TRUE,
			normalization.method = norm.method
		)
	})

	tf.anchors.list %>% qsave(
		path_join(c(prc_dir, "tf.anchors.list.qs")), nthreads=4
	)
}
