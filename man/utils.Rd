\name{utils}
\docType{methods}
\alias{keepSeqlevels}
\alias{renameSeqlevels}

\title{
  Rename or subset a \code{GRanges} or \code{gappedAlignments} object
}

\description{
  Rename or subset a \code{GRanges} or \code{gappedAlignments} object
}

\usage{
  keepSeqlevels(target, keep, ...)
  renameSeqlevels(target, rename, ...)
}

\arguments{
  \item{target}{The \code{GRanges}, \code{GRangesList} or \code{GappedAlignments}
    object to be subset or renamed by seqlevels. 
  }
  \item{keep}{A character vector of seqlevels to keep. 
  }
  \item{rename}{A named list where the names are the old seqlevels and the
    values are the new seqlevels. 
  }
  \item{...}{Arguments to be passed to other functions.
  }
}

\value{
  The \code{target} with seqlevels renamed or removed.
}

\details{
  Many operations on \code{GRanges} objects require the seqlevels 
  to match before a comparison can be made (e.g., \code{findOverlaps}).
  \code{keepSeqlevels} and \code{renameSeqlevels} are designed to
  subset objects that may have unwanted seqlevels or rename the
  existing seqlevels. 

  The \code{keepSeqlevels} function subsets the \code{target} on 
  the seqlevels provided in a character list. \code{renameSeqlevels}
  renames the seqlevels in the \code{target}.
}

\examples{
gr1 <- GRanges(seqnames=Rle(c('Chr1', 'Chr2', 'Chr3'), c(3, 3, 4)),
    IRanges(105:114, width=5), strand='-')
gr2 <- GRanges(seqnames=Rle(c('Chrom1', 'Chrom2', 'Chrom3'), c(3, 4, 3)),
    IRanges(101:110, width=10), strand='-')

rename <- list(Chr1="Chrom1", Chr2="Chrom2", Chr3="Chrom3") 

gr1_rename <- renameSeqlevels(gr1, rename)
keepSeqlevels(gr1_rename, "Chrom2")

gr2_subset <- keepSeqlevels(gr2, c("Chrom2", "Chrom3"))
keepSeqlevels(gr1_rename, gr2_subset)
}


\author{Valerie Obenchain \email{vobencha@fhcrc.org}}

\keyword{methods}
\keyword{utilities}

