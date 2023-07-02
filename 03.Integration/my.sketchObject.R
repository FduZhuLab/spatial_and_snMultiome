require(magrittr)
require(dplyr)
require(Seurat)

my.spl_fun <- function(df, group.label, map_key, seed=42L) {
    set.seed(seed)
    df %>% filter(!!sym(group.label) %in% names(map_key)) %>%
		group_by(!!sym(group.label)) %>%
		sample_n(size=map_key[!!sym(group.label)])
}


my.sketchObject <- function(
	object, sketch.size, hvg.num = 2000, dimPC = 30,
	idents = NULL, group.by=NULL, lst_num=30
) {
	# pip install geosketch

	sketch.size <- as.integer(sketch.size)
	if (1.2*sketch.size >= ncol(object)) {
	return(object)
	}

	object.2 <- object

	geosketch <- reticulate::import('geosketch')

	if (!is.null(idents)) {
		object.2 %<>% subset(idents = idents)
	}

	object.2 %<>% FindVariableFeatures(
		selection.method = "vst", nfeatures = hvg.num, verbose=FALSE
	) %>% ScaleData(verbose=FALSE) %>%
		RunPCA(npcs = dimPC, verbose = FALSE)

	X.pcs <- object.2@reductions$pca@cell.embeddings

	cells.all <- Cells(object.2)

	rm(object.2)
	foo <- gc()

	sketch.index <- geosketch$gs(X.pcs, sketch.size, seed=42L)
	sketch.index <- unlist(sketch.index) + 1
	sketch.cells <- cells.all[sketch.index]
	
	if (!is.null(group.by)) {
		# add more cells
		lst_num_list <- object@meta.data[[group.by]] %>% table
		lst_num_list %<>% sapply(function(x){min(lst_num, x)})
		
		spl_num_list <- object@meta.data[sketch.cells, group.by] %>% unlist %>% table
		
		add_num_list <- lst_num_list
		add_num_list[names(spl_num_list)] <- add_num_list[names(spl_num_list)] - spl_num_list
		add_num_list <- add_num_list[add_num_list > 0]

		if (length(add_num_list) > 0) {
			all_meta_data <- object@meta.data
			all_meta_data$'cell_names' <- rownames(all_meta_data)
			rest_meta_data <- all_meta_data %>%
				filter(!(cell_names %in% sketch.cells))
			
			add_meta_data <- rest_meta_data %>%
				my.spl_fun(group.by, add_num_list)
			
			sketch.cells <- unique(c(sketch.cells, add_meta_data$'cell_names'))
		}
	}

	object %<>% subset(cells = sketch.cells)

	return(object)
}
