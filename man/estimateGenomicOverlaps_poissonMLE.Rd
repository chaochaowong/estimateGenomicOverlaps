\name{poissonMLE-methods}
\docType{methods}
\alias{poissionMLE}
\alias{estimateGenomicOverlaps_poisson}
\alias{estimateGenomicOverlaps_poisson-methods}
\alias{estimateGenomicOverlaps,DataFrame,TranscriptDb-method}
\alias{estimateGenomicOverlaps,DataFrame,character-method}
  
\title{Estimation Transcript Expression in RNA-Seq with Poisson Models}
\description{
  Statistical inference for transcript expression in RNA-Seq using
  Poisson maximum likelihood functions. The expression level is
  quantified as read per kilobase of the transcript per million mapped
  read (RPKM).
}

\usage{
  \S4method{estimateGenomicOverlaps_poisson}{DataFrame, TranscriptDb}(query,
  subject, split.overlappingExons=TRUE, ignore.strand=FALSE)

  \S4method{estimateGenomicOverlaps_poissonx}{DataFrame, character}(query,
  subject, split.overlappingExons=TRUE, ignore.strand=FALSE)
}

\arguments{
  \item{query}{A DataFrame with \code{file} column specifing the full
    path and file names of the bam files.}
  
  \item{subject}{A \code{TranscriptDb} instance or a character string
    specifying the path to an sqlite database.}
  
  \item{split.overlappingExons}{A logical value indicating whether to
    split overlapping exons with genes. See details.}
  
  \item{ignore.strand}{A logical value indicating whether to ignore stand.}
}

\value{
  A \code{GRangesList} representing the transcript expression grouped by
  genes. The expression level is quantified as read per kilobase of the
  transcript per million mapped reads to the transcript (RPKM). 
  If the expression level is \code{NA}, it means either an error
  from the optimization method used to estimate expression or the
  iteration limit had been reached.
}

\details{
  We use Poisson models to estimate transcript expression, especially
  for transcripts that share the same exons with the others in the same
  gene. Paragraphs below are the high-level description of the
  mathematical models used to model the sequencing process. See H. Jiang (2009)
  for more details.

  
  Suppose a gene \eqn{g} has \eqn{n} transcripts with expression
  \eqn{\Theta=[\theta_1, \theta_2, \cdots, \theta_i, \cdots, \theta_n]}
  and \eqn{m} exons with length \eqn{L = [l_1, l_2, \cdots, l_j, \cdots,
    l_m]}. Denote the observed reads falling into the exons to be \eqn{X
    = [x_1, x_2, \cdots, x_j, \cdots, x_m]}. Assuming each read is sampled
  independently, \eqn{x_j} are i.i.d. random variables followed by
  Poisson distributing where its parameter \eqn{\lambda_j} is formulated
  by 
  \deqn{\lambda_j = l_j w \Sigma_{i=1}^{n} C_{ij} \theta_i,}
  where \eqn{w} is the total number of reads in a sample, and \eqn{C_{ij}}
  is \eqn{1} if transcript \eqn{i} contains exon \eqn{j} and \eqn{0}
  otherwise. Note thatthe probability mass function of the Poisson
  distribution for \eqn{x_i} is
  \deqn{P(x_j=x|\Theta) = \frac{e^{\lambda_j} \lambda_j^x}{x!}.}
  If the corresponding regions of \eqn{x_i} are disjoint, we can write
  the joint likelihood function as

  \deqn{L(\Theta|x) = P(X=x|\Theta) = \Pi_{j=1}^{m}P(x_j=x|\Theta).}

  The estimation of \eqn{\Theta} is obtained by

  \deqn{\hat{\Theta} = \arg \max \log L(\Theta|x).}

  To solve this maximum likelihood estimation problem, we employ an
  optimization function \code{\link[stats]{optim}} using \code{L-BFGS-B}
  method.

}



\seealso{\code{\link[GenomicRanges]{countGenomicOverlaps}},
  \code{\link[stats]{optim}}
  % superdisjoin
}
\examples{
\dontrun{
txdbFile <-
  "/home/cwon2/Project/cwon2/TransExpression/extdata/sacCer2_sgdGene.sqlite"

bamFile <-
  "/home/cwon2/Project/cwon2/TransExpression/extdata/SRR002051.chrI-V.bam"
query <- DataFrame(file=bamFile)
subject <- loadFeatures(txdbFile)
trs <- estimateGenomicOverlaps_poisson(query, subject,
           split.overlappingExons=TRUE)
trs
  } 
}

\references{
  Hui Jiang and Wing Hung Wong, \emph{Statistical inferences for isoform
  expression in RNA-Seq}, Bioinformatics, vol. 25, no. 8 2009, pages
  1026-1032.  
}

\author{Chao-Jen Wong \email{cwon2@fhcrc.org}}

\keyword{methods}
\keyword{utilities}

