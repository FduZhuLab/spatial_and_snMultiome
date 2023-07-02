# Load Arguments
args.split <- function(x) {
    x <- stringr::str_replace_all(string = x, pattern = " +", replacement="")
    x <- unique(stringr::str_split(string = x, pattern = ",", simplify = TRUE)[1,])
    return(x)
}

dft.args <- list()
dft.args[["iter"]] <- FALSE
dft.args[["input"]] <- "/home/whe/project/multiomics/samples"
dft.args[["output"]] <- "/home/whe/project/multiomics/preprocess-2"
dft.args[["config"]] <- NULL

arg_list <- R.utils::commandArgs(trailingOnly = T, asValue = T, defaults = dft.args)

stopifnot(length(arg_list$"input") == 1)
stopifnot(length(arg_list$"output") == 1)

if (arg_list$"iter") {
    spl_list <- args.split(arg_list$"sample")
    reg_list <- args.split(arg_list$"region")
    cfg <- data.frame(
        sample = sort(rep(spl_list, length(reg_list))),
        region = rep(reg_list, length(spl_list))
    )
    cfg$input <- paste(arg_list$"input", cfg$sample, cfg$region, sep = "/")
    cfg$output <- paste(arg_list$"output", cfg$sample, cfg$region, sep = "/")
    cfg$out_name <- paste0(cfg$sample, "_", cfg$region, ".rds")
}

if (!is.null(arg_list$"config")) {
    cfg <- read.csv(file = arg_list$"config")
    cols <- colnames(cfg)
    stopifnot("sample" %in% cols, "region" %in% cols)
    if (!"input" %in% cols) {
        cfg$input <- paste(arg_list$"input", cfg$sample, cfg$region, sep = "/")
    }
    if (!"output" %in% cols) {
        cfg$output <- paste(arg_list$"output", cfg$sample, cfg$region, sep = "/")
    }
    if (!"out_name" %in% cols) {
        cfg$out_name <- paste0(cfg$sample, "_", cfg$region, ".rds")
    }
}

cols <- colnames(cfg)
cfg_list <- list()
for (i in 1:nrow(cfg)) {
    cfg_list[[i]] <- as.list(cfg[i, ])
    names(cfg_list[[i]]) <- cols
}

print(cfg_list)

nsamples <- nrow(cfg)

# Load Packages
cat("\n"); print("Loading packages...")

library(dplyr)
library(future)
library(future.apply)
library(Seurat)
library(Signac)

num.cores <- as.integer(Sys.getenv("PBS_NP")) - 1
options(future.rng.onMisuse = "ignore")
options(future.globals.maxSize = 6 * 1024 ^ 3)
plan(multicore, workers = num.cores)

# Define
annotation <- readRDS("/home/whe/project/multiomics/Data/meta_features/EnsDb.Hsapiens.v86.rds")
hg38 <- readRDS("/home/whe/project/multiomics/Data/meta_features/cellranger_chromInfo.rds")

xLoad10X <- function(
    cfg, rna.assay = "RNA", atc.assay = "ATAC"
) {
    Exp <- Read10X(paste(cfg$input, "filtered_feature_bc_matrix", sep = "/"))

    Peak.counts <- Exp$"Peaks"

    obj <- CreateSeuratObject(Exp$"Gene Expression", assay = rna.assay, min.cells = 3)

    obj[["Sample"]] <- cfg$sample
    obj[["Region"]] <- cfg$region
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-", assay = "RNA")

    obj[[atc.assay]] <- CreateChromatinAssay(
        counts = Peak.counts, min.cells = 10, sep = c(":", "-"),
        fragments = paste(cfg$input, "atac_fragments.tsv.gz", sep = "/"),
		genome = hg38, annotation = annotation
    )

    rm(Exp, Peak.counts)
    foo <- gc()

    DefaultAssay(obj) <- atc.assay
    obj <- obj %>% TSSEnrichment() %>% NucleosomeSignal()
    obj$"blacklist_fraction" <- FractionCountsInRegion(object = obj, regions = blacklist_hg38)

#     # Add meta data
#     meta.info <- read.csv(paste(cfg$input, "per_barcode_metrics.csv", sep = "/"))
#     rownames(meta.info) <- meta.info$"barcode"
#     meta.info <- meta.info[colnames(obj), ]

#     obj$"atac_fragments" <- meta.info$"atac_fragments"
#     obj$"atac_peak_region_fragments" <- meta.info$"atac_peak_region_fragments"
#     obj$"atac_pct_reads_in_peaks" <- obj$"atac_peak_region_fragments" / obj$"atac_fragments" * 100

    DefaultAssay(obj) <- rna.assay

    return(list(obj, cfg$output, cfg$out_name))
}

xSave <- function(out) {
    system(paste("mkdir -p", out[[2]]))
    saveRDS(out[[1]], paste(out[[2]], out[[3]], sep = "/"))
    return("Success")
}

# preprocess and save

if (nsamples > 1) {
    plan(list(
        tweak(multicore, workers = nsamples),
        tweak(multicore, workers = as.integer(num.cores/nsamples))
    ))

    setOpt.FUN <- function(in.fun, ...) {
        stopifnot(is.function(in.fun))
        options(future.globals.maxSize = 6 * 1024 ^ 3)
        in.fun(...)
    }

    out_list <- future_lapply(
        cfg_list, setOpt.FUN, in.fun = xLoad10X
    )
    res <- future_lapply(out_list, xSave)
} else {
    out <- xLoad10X(cfg = cfg_list[[1]])
    res <- xSave(out)
}

print(res)
