setwd(Sys.getenv("HOME"))

source("R_func/mutate_load.R")
library(tibble)
library(arrow)

library(GenomicRanges)
library(ChIPpeakAnno)

# Load TEs
all.human.TE <- read_feather(
    "project/multiomics/PeakAnlys/PeakAnnot/all.hg38.TE.feather"
)
peak.df <- all.human.TE[, c('chr', 'start', 'end')]
peak.df %<>% unique() %>% as.data.frame
all.peak.name <- paste(
    peak.df$'chr', peak.df$'start', peak.df$'end', sep="-"
)
all.human.TE.gr <- makeGRangesFromDataFrame(peak.df)
names(all.human.TE.gr) <- all.peak.name

# overlap
own.gr <- readRDS(
    "project/multiomics/PeakAnlys/Data/peaks.rds"
)

overlap.res <- ChIPpeakAnno::findOverlapsOfPeaks(own.gr, all.human.TE.gr)
overlap.df <- overlap.res$'overlappingPeaks'[[1]]

# To boolean matrix
meta.features <- fread(
    "project/multiomics/PeakAnlys/Data/meta.features.csv"
)
all.peaks <- meta.features[, 1, ] %>% unlist %>% as.vector

bool.list <- list()
for (anot in unique(all.human.TE$'TEfamily2')) {
    ind <- all.human.TE$'TEfamily2' == anot
    tmp.te <- unique(all.human.TE$'peak_name'[ind])

    ind <- overlap.df$'peaks2' %in% tmp.te
    tmp.peak <- unique(overlap.df$'peaks1'[ind])

    ind <- all.peaks %in% tmp.peak
    bool.list[[anot]] <- as.integer(ind)
}

bool.df <- as.data.table(bool.list)
all.columns <- colnames(bool.df)
keep.columns <- all.columns[!(all.columns %>% str_detect("\\?|Unknown"))]
bool.df <- bool.df[, ..keep.columns]
rownames(bool.df) <- all.peaks
bool.df %<>% rownames_to_column("peak_name")

bool.df %>% fwrite(
    "project/multiomics/PeakAnlys/PeakAnnot/bool/TEpeak_bool.csv"
)
