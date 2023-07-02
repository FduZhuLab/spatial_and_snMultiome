setwd(Sys.getenv("HOME"))

source("R_func/mutate_load.R")
library(tibble)

library(GenomicRanges)

# Load peak ages
peak.ages <- fread(
    "project/multiomics/PeakAnlys/PeakAnnot/AgePeak.csv"
)
meta.features <- fread(
    "project/multiomics/PeakAnlys/meta.features.csv"
)
all.peaks <- meta.features[, 1, ] %>% unlist %>% as.vector

# To boolean matrix
bool.list <- list()
for (anot in unique(peak.ages$'Age')) {
    ind <- peak.ages$'Age' == anot
    tmp.peak <- unique(peak.ages$'peak_name'[ind])

    ind <- all.peaks %in% tmp.peak
    bool.list[[anot]] <- as.integer(ind)
}

bool.df <- as.data.table(bool.list)
rownames(bool.df) <- all.peaks
bool.df %<>% rownames_to_column("peak_name")

bool.df %>% fwrite(
    "project/multiomics/PeakAnlys/PeakAnnot/bool/AGEpeak_bool.csv"
)
