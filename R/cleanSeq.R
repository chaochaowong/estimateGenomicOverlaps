cleanSeq <- function(query, subject = NULL, rename = NULL, keep = NULL, ...)
{
    if (is.null(query) && is.null(rename))
        stop("rename values are present but query is not")
    ## rename 
    if (!is.null(rename)) {
        old <- rename[["old"]]
        new <- rename[["new"]]
        if (any(!old %in% seqlevels(query)))
            stop("all values in 'old' should be present in seqlevels(query)")
        seqlevels(query)[seqlevels(query) %in% old] <- new 
    }
    ## subset
    if (!is.null(keep)) {
       query <- query[seqnames(query) %in% keep]
       seqlevels(query) <-
           seqlevels(query)[seqlevels(query) %in% keep]
       if (!is.null(subject)) {
           subject <- subject[seqnames(subject) %in% keep, ]
           ## slow for lists but retains partitioning
           seqlevels(subject) <-
               seqlevels(subject)[seqlevels(subject) %in% keep]
          return(list(query=query, subject=subject))
       }
    }
   query
}

#matchSeqLevels <- function(target, template)
#{
#    ## match target to template, return target
#    target <- target[seqnames(target) %in% runValue(seqnames(template))]
#    keeplevels <- seqlevels(target) %in% seqlevels(template)
#    cat("removing the following seqlevels from target: ", 
#        seqlevels(target)[!keeplevels])
#    seqlevels(target) <- seqlevels(target)[keeplevels]
#    target
#}


setGeneric("keepSeqLevels", signature = c("target", "keep"),
           function(target, keep, ...)
           standardGeneric("keepSeqLevels")
)

setMethod("keepSeqLevels",  c("GenomicRanges", "character"),
            function(target, keep, ...)
{
    callGeneric(target, keep, ...)
})

setMethod("keepSeqLevels",  c("GappedAlignments", "GenomicRanges"),
            function(target, keep, ...)
{
    keep <- seqlevels(keep)
    callGeneric(target, keep, ...)
})

keepSeqLevels <- function(target, keep)
{
    ## keep is a character vector of seqlevels
    if (!all(keep %in% seqlevels(target)))
        warning("not all values in 'keep' are present in seqlevels(target)")
    target <- target[!seqnames(target) %in% keep] 
    seqlevels(target) <- seqlevels(target)[!seqlevels(target) %in% keep]
    target
}

renameSeqLevels <- function(target, rename)
{
    ## rename is a named list, the name is old, the element is new
    old <- names(rename) 
    new <- unlist(rename, use.names=FALSE)
    if (any(!old %in% seqlevels(target)))
        warning("not all values in names(rename) are present in seqlevels(target)")
    seqlevels(target)[seqlevels(target) %in% old] <- new
    target 
}
