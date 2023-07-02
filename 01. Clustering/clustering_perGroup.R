library(Seurat)

library(dplyr)
library(stringr)

library(future)
library(future.apply)

num.cores <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE")) - 1
options(future.globals.maxSize = 1e6 ^ 3)
options(future.rng.onMisuse = "ignore")
plan("multicore", workers = 20)

setwd(Sys.getenv("HOME"))
source("R_func/seurat/Select.Covar.R")

home_dir <- "project/multiomics/Consensus/GABAergic"

obj <- readRDS(
    paste(home_dir, "merge_rna_intg.rds", sep="/")
)

DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)

obj.list <- SplitObject(obj, "Sample")

tmp.obj.list <- lapply(
    X = obj.list, FUN = SCTransform,
    ncells=10000, n_genes=3000, min_cells=1,
    method="qpoisson", verbose=F
)
cohvg.df <- Select.Covar(tmp.obj.list)
cohvg <- rownames(cohvg.df)[1:3000]

anchors <- FindIntegrationAnchors(obj.list, anchor.features = cohvg, verbose=F)

obj <- IntegrateData(anchors, new.assay.name= "intg", verbose=F)

obj <- obj %>%
    ScaleData(vars.to.regress=c("log_nCount", "percent.mt")) %>%
    RunPCA(npcs = 50, verbose = FALSE) %>%
    RunUMAP(
        reduction = "pca", dims = 1:50, verbose=FALSE,
        reduction.name="intg.umap", reduction.key="intgUMAP_"
    ) %>%
    FindNeighbors(reduction = "pca", dims = 1:50, verbose=FALSE) %>%
    FindClusters(
        graph.name="intg_snn", resolution = seq(0.6, 2.0, by=0.1), verbose=FALSE
    )

saveRDS(
    obj, paste(home_dir, "merge_rna_intg.rds", sep="/")
)
