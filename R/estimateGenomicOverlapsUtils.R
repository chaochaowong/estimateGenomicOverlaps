### preprocessExons performs several tasks including:
### return pseudoExons with counts

## NOTE: make 'which arguement more clear in the front?
.internal_preprocessExons <- function(query, subject, rename = NULL, 
    which = NULL, ignore.strand = FALSE, ...)
{
    cleanExons <- exons(subject, vals=list(exon_chrom=seqlevels(which)),
        columns = c("exon_id", "tx_id", "gene_id"))

    ## composite genes 
    exonByGene <- exonsBy(subject, "gene")
    originalRange <- unlist(range(exonByGene))
    pseudoRange <- reduce(originalRange)
    go <- findOverlaps(cleanExons, pseudoRange)
    values(cleanExons) <- append(values(cleanExons),
        DataFrame(compositeGenes=subjectHits(go)))

    ## pseduo exons
    pseudoExons <- disjoin(cleanExons)
    eo <- findOverlaps(cleanExons, pseudoExons)
    originalExons <-  IntegerList(split(queryHits(eo), subjectHits(eo)))
    expandtx <- rep.int(subjectHits(eo),
        elementLengths(values(cleanExons)[["tx_id"]][queryHits(eo)]))
    tx <- IntegerList(split(unlist(
        values(cleanExons)[["tx_id"]][queryHits(eo)]), expandtx))
    compositeGenes <- split(values(cleanExons)[["compositeGenes"]][
        queryHits(eo)], subjectHits(eo))
    compositeGenes <- IntegerList(lapply(compositeGenes, unique))
    values(pseudoExons) <- DataFrame(originalExons, tx, compositeGenes)

    ## clean and count reads
    co <- countOverlaps(pseudoExons, cleanExons)
    mapOriginalToPseudo <-
        findOverlaps(cleanExons, pseudoExons[co == 1])
    fls <- query[["files"]]
    lst <- lapply(fls, .dispatch, pseudoExons, cleanExons, rename,
        which, ignore.strand, mapOriginalToPseudo)
    counts <- DataFrame(do.call(cbind, lst))
    values(pseudoExons) <-
        append(values(pseudoExons), counts)

    ## FIXME : can be avoided for glm if not used
    totalReads <- lapply(fls, function(x, which) {
        sum(countBam(x, param=ScanBamParam(which=which))$records)
        }, which)
 
    metadata <- list(groups=query[["group"]], counts=names(counts),
        totalReads=unlist(totalReads), genes="compositeGenes", tx="tx")

    return(list(Exons=pseudoExons, metadata=metadata))

}


.checkValidQuery <- function(DF, checkColName="file", checkFile.exist=FALSE)
{
    for (i in seq_len(length(checkColName)))
        if (!any(names(DF) %in% checkColName[i]))
            stop(paste("the query much have column named '",
                       checkColName[i], "'.", sep=""))
    ## check if the file exists
    if (any(names(DF) %in% "file")) {
        exist_flag <- file.exists(as.character(DF[["file"]]))
        if (!all(exist_flag))
            stop(paste("file does not exist:",
                       paste(DF[["file"]][!exist_flag], collapse=" and ")))
    }
}

#.getGeneModels <- function(eBygene, tBygene)
#{
#    ## exon-tx relationship
#    tx <- values(eBygene@unlistData)[["tx_id"]]
#    txFactor <- new("CompressedIntegerList", 
#                    unlistData = factor(unlist(tx)),
#                    partitioning = tx@partitioning)
#    exonMaptx <- rep(seq_along(values(eBygene@unlistData)[["exon_id"]]), 
#                     elementLengths(txFactor))
#
#    ## matrix of gene models
#    df <- data.frame(exon=exonMaptx, tx=unlist(txFactor))
#    geneModels <- 
#        matrix(0L, nrow=length(values(eBygene@unlistData)[["exon_id"]]), 
#               ncol=length(unique(unlist(txFactor))))
#    colnames(geneModels) <- sort(unique(unlist(txFactor)))
#    geneModels[cbind(as.integer(df$exon), as.integer(df$tx))] <- 1L
#
#    ## indices to subset model matrix by gene
#    exonMapgene <- 
#        split(seq_along(values(eBygene@unlistData)[["exon_id"]]),
#              rep.int(seq_len(length(eBygene)), elementLengths(eBygene)))
#              #rep(names(eBygene), elementLengths(eBygene)))
#    names(exonMapgene) <- names(eBygene)
#
#    txMapgene <- 
#        split(values(tBygene@unlistData)[["tx_id"]],
#              rep.int(seq_len(length(tBygene)), elementLengths(tBygene)))
#    names(txMapgene) <- names(tBygene)
#              #rep(names(tBygene), elementLengths(tBygene)))
#              #factor(as.vector(rep(names(tBygene@partitioning),
#              #       width(tBygene@partitioning)))))
#
#    return(list(geneModels=geneModels, exonMapgene=exonMapgene,
#                txMapgene=txMapgene))
#}

.dispatch <- 
    function(fls, pseudoSubject, originalSubject, rename, 
             which, ignore.strand, mapOriginalToPseudo, ...)
{
 ga <- readGappedAlignments(fls, which=which)
    ## rename
    if (!is.null(rename))
            ga <- renameSeqlevels(ga, rename)
    ## subset 
    if (!is.null(which))
            ga <- keepSeqlevels(ga, seqlevels(which))
    grlclean <- as(ga, "GRangesList")

    ## count 
    .countPseudoExons(grlclean, pseudoSubject, originalSubject,
                      ignore.strand, mapOriginalToPseudo)
}

.makeGeneModels <- function(Exons, metadata, ...)
{
    if (class(Exons) != "GRanges")
        stop("'Exons' argument must be a GRanges object")
    exExpand <- rep.int(seq_len(length(Exons)),
        elementLengths(values(Exons)[ , metadata[["tx"]]]))
    txExpand <- unlist(values(Exons)[ , metadata[["tx"]]])
    txUnique <- unique(txExpand)
    txIdx <- match(txExpand, txUnique)
    gm <- matrix(0L, nrow=length(Exons), ncol=length(txUnique))
    index <- cbind(exExpand, txIdx)
    gm[index] <- 1L
    gm
}

.countPseudoExons <- 
    function(query, pseudoSubject, originalSubject,
             ignore.strand, mapOriginalToPseudo, ...)
{
    ## pseudo regions
    co <- countOverlaps(query, pseudoSubject, type="within")
    pseudo <- countGenomicOverlaps(query[co > 0], pseudoSubject,
        type="within", resolution="none", ignore.strand)
    pseudoHits <- values(pseudo)[["hits"]]

    ## original regions 
    original <- countGenomicOverlaps(query[co == 0], originalSubject,
        type="any", resolution="uniqueDisjoint", ignore.strand)
    originalHits <- rep.int(0L, length(pseudoSubject))
    h <- values(original)[["hits"]][queryHits(mapOriginalToPseudo)]
    originalHits[subjectHits(mapOriginalToPseudo)] <- h

    pseudoHits + originalHits
}

setGeneric("keepSeqlevels", signature = c("target", "keep"),
           function(target, keep, ...)
           standardGeneric("keepSeqlevels")
)

setMethod("keepSeqlevels",  c("ANY", "GenomicRanges"),
            function(target, keep, ...)
{
    keep <- seqlevels(keep)
    callGeneric(target, keep, ...)
})

setMethod("keepSeqlevels",  c("ANY", "character"),
            function(target, keep, ...)
{
    if (!any(keep %in% seqlevels(target)))
        warning("no values in 'keep' are present in seqlevels(target)")
    target <- target[seqnames(target) %in% keep]
    seqlevels(target) <- seqlevels(target)[seqlevels(target) %in% keep]
    target
})

renameSeqlevels <- function(target, rename)
{
    old <- names(rename)
    new <- unlist(rename, use.names=FALSE)
    if (!any(old %in% seqlevels(target)))
        warning("no values in names(rename) are present in ",
                "seqlevels(target)")
    seqlevels(target)[seqlevels(target) %in% old] <- new
    target
}

#cleanSeq <- function(query, subject = NULL, rename = NULL, keep = NULL, ...)
#{
#    if (is.null(query) && is.null(rename))
#        stop("rename values are present but query is not")
#    ## rename 
#    if (!is.null(rename)) {
#        old <- rename[["old"]]
#        new <- rename[["new"]]
#        if (any(!old %in% seqlevels(query)))
#            stop("all values in 'old' should be present in seqlevels(query)")
#        seqlevels(query)[seqlevels(query) %in% old] <- new
#    }
#    ## subset
#    if (!is.null(keep)) {
#       query <- query[seqnames(query) %in% keep]
#       seqlevels(query) <-
#           seqlevels(query)[seqlevels(query) %in% keep]
#       if (!is.null(subject)) {
#           subject <- subject[seqnames(subject) %in% keep, ]
#           ## slow for lists but retains partitioning
#           seqlevels(subject) <-
#               seqlevels(subject)[seqlevels(subject) %in% keep]
#          return(list(query=query, subject=subject))
#       }
#    }
#   query
#}
#

