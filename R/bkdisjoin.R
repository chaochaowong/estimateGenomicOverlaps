setMethod("disjoin", "GRangesList",
    function(x, ...)
    {
        if (length(list(...)) != 0L)
            stop("\"disjoin\" method for GRangeList objects only ",
                 "takes a single object")
        gr <- GenomicRanges:::deconstructGRLintoGR(x)
        gr <- disjoin(gr)
        GenomicRanges:::reconstructGRLfromGR(gr, x)
    }
)          


disjoinKeepValues <- function(x)
{
    disj <- isDisjoint(x) 
    ## call .bkdisjoin
    bkdisj <- .bkdisjoin(x[!disj])
    nd <- .processMetadata(x[disj]) # remove gene_id, change exon_id type
    newx <- c(nd, bkdisj$grl)
    newx <- newx[names(x)]
    list(grl=newx, disj=disj, disModel=bkdisj$disModel)
}




.bkdisjoin <- function(x, keep=c("tx_id", "exon_id")) {
    ## x : query, a return value from countGenomicOverlaps,
    ## x and ga can be a list of GRangesList
    ## ga: gappedalilgment, subject
    gr <- GenomicRanges:::deconstructGRLintoGR(x)
    
    dgr <- disjoin(gr)
    olaps <-  findOverlaps(query=dgr, subject=gr)
    disModel <- xtabs(~subjectHits(olaps) + queryHits(olaps))
    psuedox <- as.logical(colSums(disModel)-1)
    ## only need to fix
    if (any(psuedox)) {
        res <- .getTxId(olaps, gr, psuedox, disModel)
        values(dgr) <-
                DataFrame(exon_id = res$exid,
                          tx_id = res$txid)
         #             gene_id = CharacterList(as.list(gid)))

        gr <- dgr
    }
    
    list(grl=GenomicRanges:::reconstructGRLfromGR(gr, x), disModel=disModel)
}

## .restore_seqnames() is copied from GenomicRanges:::reconstructGRLfromGR
.restore_seqnames <- function(gr, x) {
    ## map the seqnames back to normal
    snames <- strsplit(as.character(seqnames(gr)), "|", fixed=TRUE)
    m12  <- matrix(as.integer(unlist(snames)), ncol=2, byrow=TRUE)
    f2 <- m12[, 2L]
    x_seqlevels <- seqlevels(x)
    gr@seqnames <- Rle(factor(x_seqlevels[f2], x_seqlevels))
    gr@seqinfo <- seqinfo(x)
    gr
}

.getTxId <- function(olaps, gr, psuedox, disModel)
{
    ## original code: takes too long
    #txid <- apply(disModel, 2, function(x) 
    #                  unique(unlist(values(gr)$tx_id[as.logical(x)])))
    #exonid <- apply(disModel, 2, function(x)
    #                paste(values(gr)$exon_id[as.logical(x)], collapse="|")
    mo <- as.matrix(olaps)
    np <- which(psuedox==FALSE)
    m  <- mo[mo[,1] %in% np, ]
    f1 <- values(gr)$tx_id[m[,2]]
    names(f1) <- m[,1]
    f2 <- apply(as.matrix(disModel[, psuedox]), 2, function(y) 
                      unique(unlist(values(gr)$tx_id[as.logical(y)])))
    f2 <- IntegerList(f2)
    if (length(f2)==1L) names(f2)=as.character(which(psuedox==TRUE))
    f <- c(f1, f2)
    f <- f[as.character(sort(as.numeric(names(f))))]

    ## for exonid
    e1 <- values(gr)$exon_id[m[,2]]
    names(e1) <- m[,1]
    e2 <- apply(as.matrix(disModel[, psuedox]), 2, function(x)
                    paste(values(gr)$exon_id[as.logical(x)], collapse="|"))
    ee <- c(e1,e2)
    ee <- ee[as.character(sort(as.numeric(names(ee))))]
    return(list(exid=ee, txid=f))
}


    
.processMetadata <- function(x)
{
    gr <- GenomicRanges:::deconstructGRLintoGR(x)
    df <- values(gr)
    df <- df[, !(colnames(df) %in% "gene_id")]
    df[, "exon_id"] <- as.character(df[, "exon_id"])
    values(gr) <- df
    GenomicRanges:::reconstructGRLfromGR(gr, x)
}

.getGeneModelsOneTx <- function(grl)
{   ## build gene model only for genes with one transcript
    ## x: GRangesList
    df <- values(grl@unlistData)
    df[, "tx_id"] <- factor(unlist(df[, "tx_id"]))
    df <- as.data.frame(df)
    gM<- model.matrix(~ tx_id + 0, df)
    colnames(gM) <- levels(df$tx_id)
    gM   
}
