\name{cleanSeq}
\docType{methods}
\alias{cleanSeq}

\title{Rename or subset the seqlevels of a GRanges object}

\description{
  Rename and or subject the seqlevels of a GRanges object or a pair of objects.
}

\usage{
  cleanSeq(query, subject = NULL, rename = NULL, keep = NULL, ...)
}

\arguments{
  \item{query}{A \code{GRanges} or \code{GRangesList} object to subset 
     and/or rename by seqlevels.
  }

  \item{subject}{A \code{GRanges} or \code{GRangesList} object to subset 
     by seqlevels.
  }
  \item{rename}{A \code{DataFrame} object with columns 
    \code{old} and \code{new} that are character vectors
    indicating which seqlevels in the \code{query} should be
    renamed. The seqlevels in \code{old} must match existing levels
    in the \code{query}.
  }
  \item{keep}{A character vector of seqlevels to keep.
  }
}

\value{
  When only a \code{query} is provided, the \code{query} object is 
  returned with seqlevels renamed as specified in \code{rename} and 
  subset on the seqlevels in the \code{keep}. If both a \code{query} 
  and \code{subject} are provided, both objects are returned with the
  \code{query} subset and renamed but the \code{subject} will only
  be subset. 
}

\details{
  \code{cleanSeq} prepares \code{GRanges} objects
  for comparison by subsetting and renaming the seqlevels.
  Many operations on \code{GRanges} objects require the seqlevels 
  to match before a comparison can be made (e.g., \code{findOverlaps}).

  The \code{cleanSeq} function accepts both a \code{query} and
  a \code{subject} or just a \code{query}. The \code{query} can be
  both subset on seqlevels and existing seqlevels can be renamed.
  The \code{subject} can only be subset on seqlevels. 

  If a single object needs seqlevels renamed or subset, only the
  \code{query} argument along with the appropriate \code{rename} and
  or \code{keep} arguments need be used (i.e., no \code{subject}).
}

\examples{

gr1 <- GRanges(seqnames=Rle(c('Chr1', 'Chr2', 'Chr3'), c(3, 3, 4)),
                      IRanges(105:114, width=5), strand='-')
gr2 <- GRanges(seqnames=Rle(c('Chrom1', 'Chrom2', 'Chrom3'), c(3, 4, 3)),
               IRanges(101:110, width=10), strand='-')

newNames <- DataFrame(old = c("Chr1", "Chr2", "Chr3"), 
    new = c("Chrom1", "Chrom2", "Chrom3"))
subsetChrom <- c("Chrom2", "Chrom3")

res <- cleanSeq(gr1, gr2, newNames, subsetChrom)

res <- cleanSeq(gr1, rename=newNames)
}


\author{Valerie Obenchain \email{vobencha@fhcrc.org}}

\keyword{methods}
\keyword{utilities}

