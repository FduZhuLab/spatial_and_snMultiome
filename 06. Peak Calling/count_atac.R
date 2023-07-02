# Load packages
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(GenomicRanges))

suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(dplyr))

suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(future.apply))
options(future.rng.onMisuse = "ignore")
options(future.globals.maxSize = 1e11)
num.cores <- as.integer(Sys.getenv("PBS_NP"))
plan(multicore, workers = num.cores)

setwd(Sys.getenv("HOME"))

# Load meta data
full_meta <- read.csv(
    "project/multiomics/Data/meta_data/meta.data.csv", row.names=1
)
full_meta$'barcode' <- rownames(full_meta)
meta_list <- full_meta %>% dplyr::group_by(Sample, Region) %>% dplyr::group_split()

nL <- c()
for (i in seq(meta_list)) {
    df <- meta_list[[i]]
    sample <- df$'Sample'[1]
    region <- df$'Region'[1]
    nL <- c(nL, paste(sample, region))
}
names(meta_list) <- nL

# Load peaks
final.peaks <- read.table(
    "project/multiomics/CallPeak/peaks/merge/final.union.bed",
    sep="\t", header=F,
    col.names = c("seqnames", "start", "end", "name", "score")
) %>% GenomicRanges::makeGRangesFromDataFrame()

# Load fragments
frag.list <- list()
file_base <- "project/multiomics/RenameFrag"
for (donor in c("No_53", "No_61")) {
    for (region in c("M1C", "V1C", "ACC", "dlPFC")) {
        fragF <- paste(file_base, donor, region, "atac_fragments.tsv.gz", sep="/")
        n <- paste(donor, region)
        frag.list[[n]] <- CreateFragmentObject(
            fragF, cells=meta_list[[n]]$'barcode', verbose=F
        )
        message("\nLoaded fragments for ", n)
        message("There are ", length(meta_list[[n]]$'barcode'), " cells")
    }
}

# Count Features
counts <- FeatureMatrix(
    fragments=frag.list, features=final.peaks,
    cells=rownames(full_meta), verbose=F
)
message("\nCounted features for merged peaks")

# Create ATAC
merg.atac <- CreateChromatinAssay(counts, fragments=frag.list)

# Add Annot info
annot.gr <- readRDS("project/multiomics/Data/meta_features/EnsDb.Hsapiens.v98.rds")
gen_info <- readRDS("project/multiomics/Data/meta_features/cellranger_chromInfo.rds")

Annotation(merg.atac) <- annot.gr
genome(merg.atac) <- gen_info

# Save
full.obj <- readRDS(
	"project/multiomics/GeneAnlys/Merge/merge_rna_intg.rds"
)

full.obj[["ATAC"]] <- merg.atac

# Metric
DefaultAssay(full.obj) <- "ATAC"
full.obj %<>% RunTFIDF()

full.obj$'atac_peak_region_fragments' <- colSums(full.obj$ATAC@counts)
full.obj$'log10_nFrag' <- log10(full.obj$'atac_peak_region_fragments')
full.obj$"atac_pct_reads_in_peaks" <- full.obj$"atac_peak_region_fragments" / full.obj$"atac_fragments" * 100

full.obj$'blacklist_fragments' <- CountsInRegion(
    full.obj, assay="ATAC", regions=blacklist_hg38_unified
)
full.obj$"blacklist_ratio" <- full.obj$"blacklist_fragments" / full.obj$"atac_peak_region_fragments" * 100

full.obj %<>% TSSEnrichment()

saveRDS(full.obj, "project/multiomics/PeakAnlys/merge_obj.rds")

message("Saved...")

