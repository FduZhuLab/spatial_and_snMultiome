suppressPackageStartupMessages(require(fs))
suppressPackageStartupMessages(require(magrittr))

my_FindM <- function(
	object, ident.1, ident.2=NULL, features=NULL,
	p.thres=0.05, padj.thres=0.05, ...
) {

	M <- suppressMessages(Seurat::FindMarkers(
		object = object, ident.1 = ident.1, ident.2 = ident.2,
		features = features, verbose = FALSE,
		...
	))

    if (nrow(M) > 0) {
        M %<>% filter((p_val < p.thres) & (p_val_adj < padj.thres)) %>%
            mutate(pct.qdiff = abs(pct.1 - pct.2)/pmax(pct.1, pct.2))

        pos.M <- M %>% filter(avg_log2FC > 0) %>% arrange(desc(abs(avg_log2FC)))
        neg.M <- M %>% filter(avg_log2FC < 0) %>% arrange(desc(abs(avg_log2FC)))

        return(list("pos"=pos.M, "neg"=neg.M))
    } else {
        return(list("pos"=NULL, "neg"=NULL))
    }
}

my.pmax <- function(...) {pmax(..., na.rm=TRUE)}
my.pmin <- function(...) {pmin(..., na.rm=TRUE)}

sum_Mdf <- function(marker_df_list) {
    if (length(marker_df_list) == 0) {return(NULL)}

    glist <- lapply(marker_df_list, row.names)
    gl.df <- as.data.frame(sapply(glist, function(x){paste(x,collapse=" ")}))
    colnames(gl.df) <- "gene"

    f.glist <- glist[!unlist(sapply(glist, is.null))]
    if (length(f.glist) == 0) {return(NULL)}

    C.gl <- Reduce(intersect, f.glist)
    if (length(C.gl) == 0) {return(NULL)}

    C.marker_df_list <- list()
    for (x in names(f.glist)) {
        tmpM <- marker_df_list[[x]][C.gl, ]
        colnames(tmpM) <- paste0(x, ".", colnames(tmpM))
        C.marker_df_list[[x]] <- tmpM
    }

    C.marker_df <- C.marker_df_list %>% dplyr::bind_cols() %>%
        mutate(pct.1 = do.call(my.pmin, c(select(., ends_with(".pct.1"))))) %>%
        select(-ends_with(".pct.1"))

    metricL <- c("p_val_adj", "avg_log2FC", "pct.2", "pct.qdiff")
    for (m in metricL) {
        C.marker_df %<>% mutate(
            !!paste0("max.", m) := do.call(my.pmax, c(abs(select(., ends_with(m))))),
            !!paste0("min.", m) := do.call(my.pmin, c(abs(select(., ends_with(m)))))
        )
    }

    C.marker_df$'gene' <- rownames(C.marker_df)
	
	C.marker_df %<>% arrange(max.p_val_adj)

    return(list(marker_list=gl.df, comb_marker=C.marker_df))
}

if.empty.df <- function(df) {
    if (is.null(df)) {
        return(TRUE)
    } else if (nrow(df) == 0) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

my_runDEG <- function(
    object, group, region, all.reg, group.label, region.label,
	home_dir, silent=FALSE, min.cells.group = 3,
    features=NULL, ...
) {
    object %<>% DietSeurat(assays=DefaultAssay(object)) %>%
    	subset(!!sym(group.label) == group)

	foo <- gc()
    Idents(object) <- region.label

    avail.reg <- table(as.character(object@meta.data[[region.label]]))
    avail.reg <- names(avail.reg)[avail.reg >= min.cells.group]

    if (!region %in% avail.reg) {
		message(region, " is not available in ", group)
		return(NULL)
	}

    all.reg <- all.reg[all.reg %in% avail.reg]
    if (length(all.reg) == 1) {
		message("Not enough regional group in ", group)
		return(NULL)
	}

    other.reg <- all.reg[all.reg != region]

    # DEG
    regM.list <- list()
    regM.list$"pos" <- list()
    regM.list$"neg" <- list()

    for (vs.reg in other.reg) {
        tmpML <- try(my_FindM(
			object, ident.1=region, ident.2=vs.reg, features=features, ...
		), silent=silent)

        if (any(unlist(lapply(tmpML, is.data.frame)))) {
			regM.list$"pos"[[vs.reg]] <- tmpML$"pos"
			regM.list$"neg"[[vs.reg]] <- tmpML$"neg"
		}
    }

	if (!is.null(home_dir)) {
		file_home <- paste(home_dir, group, sep="/")
		dir_create(file_home)
		saveRDS(regM.list, paste(file_home, "merged.markers_list.rds", sep="/"))
	}

    # pos&neg
    comb.markers.list <- list()
    for (p in c("pos", "neg")) {
        M.list <- regM.list[[p]]

        ## saveRDS
		if (!is.null(home_dir)) {
			res_dir <- paste(file_home, p, sep="/")
			dir_create(res_dir)
			saveRDS(M.list, paste(res_dir, "cons.markers_list.rds", sep="/"))
		}

        ## combined markers df
        markers.df <- sum_Mdf(M.list)[[2]]
        if (!if.empty.df(markers.df)) {
            markers.df$'cluster' <- group
            markers.df$'region' <- region
            markers.df$'sign' <- p
			comb.markers.list[[p]] <- markers.df
			if (!is.null(home_dir)) {
				write.csv(
					markers.df, paste(res_dir, "cons.markers_df.csv", sep="/")
				)
			}
        }
    }

    if (length(comb.markers.list) > 0) {
		return(comb.markers.list)
	} else {
		return(NULL)
	}
}

