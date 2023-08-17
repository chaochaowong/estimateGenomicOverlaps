\name{estimateGenomicOverlaps_glm-methods}
\docType{methods}
\alias{estimateGenomicOverlaps_glm,DataFrame,TranscriptDb-method}
\alias{estimateGenomicOverlaps_glm,DataFrame,character-method}


\title{Estimate Transcript Expression in RNA-Seq with Generalized Linear Models}

\description{
    Estimate transcript expression using the 'glm' function.
}

\usage{
  \S4method{estimateGenomicOverlaps_glm}{DataFrame, TranscriptDb}(query,
  subject, rename=NULL, which=NULL, ignore.strand=FALSE, 
  control=glm.control(), ...)

  \S4method{estimateGenomicOverlaps_glm}{DataFrame, character}(query,
  subject, rename=NULL, which=NULL, ignore.strand=FALSE, 
  control=controlglm.control(), ...)
}

\arguments{
  \item{query}{
    A DataFrame with \code{file} column specifying the full path to the
    bam files and a \code{groups} column specifying what groups the files
    correspond to. 
  }
  \item{subject}{
    A \code{\link{TranscriptDb}} object or character string specifying the full 
    path to an sqlite database that can be made into a \code{\link{TranscriptDb}}
  }
  \item{rename}{This argument is passed to the \code{cleanSeq} function.
    A \code{DataFrame} object with columns \code{old} and \code{new} that 
    are character vectors indicating which seqlevels in the \code{query} 
    should be renamed. The seqlevels in \code{old} must match existing levels
    in the \code{query}.
  }
  \item{which}{This argument is the 'which' argument in \code{ScanBamParam}.
    See \code{?ScanBamParam}. Both the bam file specified in the \code{query}
    and the \code{TranscriptDb} object created from the \code{subject}
    will be subset on the sequences specified in \code{which}.
  }
  \item{ignore.strand}{This argument is passed to the 
    \code{countGenomicOverlaps} function. A logical value indicating if strand 
    should be considered when matching.
  }
  \item{control}{The the 'control' argument from \code{glm}. See
    \code{?glm.control}. 
  }
  \item{...}{Additional arguments, perhaps used by methods defined on
    this generic.
  }
}

\value{
  When a single bam file is supplied in the \code{query}
  the result is a \code{list} of 'expression', 'standard_error',
  and 'p-value' for each transcript.

  If multiple files are supplied, and the files represent
  different groups, the results are provided for
  each transcript as well as for the group. If multiple
  files of the same group type are provided, the results are
  an average over all files. 
}

\details{
  \code{estimateGenomicOverlaps_glm} counts reads stored in 
  bam files against a reference \code{\link{TranscriptDb}} using 
  \code{countGenomicOverlaps}. Transcript expression is estimated
  with \code{\link{glm}} using with the negative binomial or
  poisson family distribution. 

  A two-pass strategy is used for the counting. The first pass involves
  dividing the exons into disjoint regions and counting reads that fall
  entirely within these regions. The second pass counts reads that hit
  the original exon regions and resolves double-hits with \code{resolution}
  'uniqueDisjoint'. Reads are tallied according to the pseudo exon regions.

  Overlapping genes are merged into composite gene groups and transcript
  estimation is performed for each of these composite gene groups using
  the counts from the pseudo exon regions.

  The 'expression' value in the output is provided as a log value.
  The exponential of the 'expression' is representative of the counts
  for each transcript.

  FIXME : Explain p-value for group 
}


\author{Valerie Obenchain \email{vobencha@fhcrc.org}}

\examples{
\dontrun{
## data on orca4
txdb <- makeTranscriptDbFromUCSC(genome = "hg18", tablename = "refGene")
path <- "/mnt/cpl/data/Solexa/SOC/101101/Samples/"
files <- c(
    paste(path, "tophat_Normal_3018_246639/accepted_hits.bam", sep=""),
    paste(path, "tophat_Normal_5309_266588/accepted_hits.bam", sep=""),
    paste(path, "tophat_SOC_5622_266562/accepted_hits.bam", sep=""),
    paste(path, "tophat_SOC_5936_294158/accepted_hits.bam", sep=""))

metadata <- DataFrame(files=files, group=c(0,0,1,1))
which <- GRanges(c("chr1", "chr2", "chr3", "chr4", "chr5"), IRanges(1, 1e+07))

res <- estimateGenomicOverlaps(query=DF, subject=txdb, method="glm", rename=NULL, keep=keep, which=which,
    control=glm.control(maxit=50))
}
}

\keyword{methods}
\keyword{utilities}

