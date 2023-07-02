setwd(Sys.getenv("HOME"))

source("R_func/mutate_load.R")
library(tibble)
library(arrow)

# Load background peaks
avail_metric <- read_feather(
    "project/multiomics/SCENIC/Data/Regional_DA/agg.avail.mdf.feather"
)

name_key <- 'cluster'
tmp.list <- avail_metric %>% group_split(!!sym(name_key))
avail.list <- tmp.list %>% lapply(function(x){unique(x$'peak_name')})
names(avail.list) <- tmp.list %>% sapply(function(x){x[[name_key]][1]})

# Load AUPs
AUPs <- read_feather(
    "project/multiomics/SCENIC/Data/Regional_DA/agg.pos.mdf.feather"
)

name_key <- 'DEType'
tmp.list <- AUPs %>% group_split(!!sym(name_key))
deg.list <- tmp.list %>% lapply(function(x){unique(x$'peak_name')})
names(deg.list) <- tmp.list %>% sapply(function(x){x[[name_key]][1]})

# set background list for AUPs list
bkg.list <- list()
for (nm in names(deg.list)) {
    nm2 <- str_split(nm, ":", n=2)[[1]][2]
    bkg.list[[nm]] <- avail.list[[nm2]]
}

# Load Peak-Annotation boolean matrix
base_dir <- "project/multiomics/PeakAnlys/PeakAnnot/bool"
peak_annot <- fread(glue("{base_dir}/TEpeak_bool.csv"))
save_label <- "TE_family"

# counting
my.counting_meet <- function(
    degL, bkgL, peak_df
) {
    nondegL <- bkgL[!(bkgL %in% degL)]

    sum_A <- peak_df %>% filter(peak_name %in% degL) %>%
        select(-peak_name) %>% colSums()

    sum_C <- peak_df %>% filter(peak_name %in% nondegL) %>%
        select(-peak_name) %>% colSums()
    
    meet_sum <- cbind(sum_A, sum_C) %>% as.data.frame %>% rownames_to_column
    colnames(meet_sum) <- c("Annot", "A", "C")
    
    meet_sum$'DEG_num' <- length(degL)
    meet_sum$'BKG_num' <- length(bkgL)
    
    return(meet_sum)
}

meet_stats <- mapply(
    FUN = my.counting_meet, deg.list, bkg.list,
    MoreArgs=list('peak_df'=peak_annot),
    SIMPLIFY=FALSE
) %>% bind_rows(.id="DEType")

meet_stats %<>% mutate(
    B = DEG_num - A,
    D = (BKG_num - DEG_num) - C
)

DEType.mtx <- meet_stats$'DEType' %>% str_split(":", simplify=T)
colnames(DEType.mtx) <- c("type", "cluster")
meet_stats %<>% cbind(DEType.mtx)
meet_stats %<>% select(
    type, cluster, DEType, Annot,
    BKG_num,DEG_num,A,B,C,D
)

# Fisher exact test
meet_test <- meet_stats %>% select(A,B,C,D) %>% apply(1, function(X){
    X %<>% unlist %>% matrix(nrow=2)
    fisher.test(X, alternative="greater")$'p.value'
})
meet_stats$'p.value' <- meet_test

meet_stats %<>% group_split(DEType) %>%
    lapply(function(df){
        df %>% mutate(p.adj.bh = p.adjust(p.value, method="BH"))
    }) %>% bind_rows()

# save
save_base <- "/home/whe/project/multiomics/PlotFig/Peak_annot"
save_dir <- glue("{save_base}/{save_label}")
dir_create(save_dir)
meet_stats %>% fwrite(glue("{save_dir}/meet_stats.csv"))
