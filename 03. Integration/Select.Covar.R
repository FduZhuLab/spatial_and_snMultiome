require(magrittr)

my.find.hvg <- function(object, method, avail.features, nfeatures=3000) {
    if (method == "SCT") {
        hvg.df <- object[[DefaultAssay(object)]]@SCTModel.list[[1]]@feature.attributes
        avail.features <- intersect(avail.features, rownames(hvg.df))
        hvg.df <- hvg.df[avail.features, , drop=FALSE] %>% arrange(desc(residual_variance))
    } else if (method == "vst") {
        hvg.df <- object[[DefaultAssay(object)]]@meta.features
        avail.features <- intersect(avail.features, rownames(hvg.df))
        hvg.df <- hvg.df[avail.features, , drop=FALSE] %>% arrange(desc(vst.variance.standardized))
    } else {
        stop("Please specific right method.")
    }

    VariableFeatures(object) <- hvg.df %>% head(n=nfeatures) %>% row.names
    return(object)
}


Select.Covar <- function(obj.list, method, hvg.mask=NULL) {
    stopifnot(length(obj.list) > 1)

    if (method == "SCT") {
        obj.var.list <- lapply(
            obj.list, function(x) {x[[DefaultAssay(x)]]@SCTModel.list[[1]]@feature.attributes}
        )
        co.features <- Reduce(
            intersect, lapply(obj.var.list, rownames)
        )
        if (!is.null(hvg.mask)) {
            co.features <- co.features[!co.features %in% hvg.mask]
        }
        var.rank.list <- lapply(
            obj.var.list, function(x) {rank(-x[co.features, 'residual_variance'])}
        )
    } else if (method == "vst") {
        obj.var.list <- lapply(
            obj.list, function(x) {x[[DefaultAssay(x)]]@meta.features}
        )
        co.features <- Reduce(
            intersect, lapply(obj.var.list, rownames)
        )
        if (!is.null(hvg.mask)) {
            co.features <- co.features[!co.features %in% hvg.mask]
        }
        var.rank.list <- lapply(
            obj.var.list, function(x) {rank(-x[co.features, 'vst.variance.standardized'])}
        )
    } else {
        stop("Please specific right method.")
    }

    names(var.rank.list) <- paste("rank", names(var.rank.list), sep=".")
    var.rank.df <- data.frame(var.rank.list)
    rm(obj.var.list, var.rank.list)
    rownames(var.rank.df) <- co.features
    var.rank.df$'mean.rank' <- apply(var.rank.df, 1, mean)
    var.rank.df <- var.rank.df[
        order(var.rank.df$'mean.rank', decreasing=F), 
    ]
    return(var.rank.df)
}
