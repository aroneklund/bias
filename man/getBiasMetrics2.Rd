\name{getBiasMetrics2}
\alias{getBiasMetrics2}
\title{ Calculate an extended set of bias metrics }
\description{ Calculate an extended set of bias metrics from an AffyBatch object.  }
\usage{
getBiasMetrics2(x.batch, x.rma = rma(x.batch))
}
\arguments{
  \item{x.batch}{ An \code{\link[affy:AffyBatch-class]{AffyBatch}}.}
  \item{x.rma}{ An \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} calculated from \code{x.batch}.}
}
\details{
  This function calculates several bias metrics for each microarray in a data set.  The first four bias metrics are equivalent to those calculated by \code{\link{getBiasMetrics}}, and are those described in the publication.  The other \dQuote{extended} bias metrics are other numbers that might or might not be interesting.
  
  Most of these bias metrics only make sense for Affymetrix expression arrays. Furthermore, four of the bias metrics are based on specific probes that may be missing from some Affymetrix platforms (in which case they will be \code{NA}).
}
\value{
  A data frame containing bias metrics for each array.
  \item{pm.median}{Median intensity of all perfect-match probes on the array.}
  \item{pm.IQR}{Interquartile range of all perfect-match probes on the array.}
  \item{rma.IQR}{Interquartile range of all expression values on the array.}
  \item{degradation}{The slope returned by \code{\link[affy]{AffyRNAdeg}}; indicates the relative preference for 3'-biased probes in a probe set.}
  \item{present.calls}{The fraction of probe sets called present by \code{\link[affy]{mas5calls}}.}
  \item{rma.spikes.mRNA}{Mean expression value of mRNA-level control bacterial spikes.}
  \item{rma.spikes.cRNA}{Mean expression value of cRNA-level control bacterial spikes.}
  \item{rma.rRNA}{Mean expression value of five rRNA probesets.}
  \item{rma.alu}{ Expression value of the Alu probeset.}
  \item{border.plus}{Median intensity of the \dQuote{plus} elements of the border checkerboard (see \code{\link[bias]{borders}}).}
  \item{border.minus}{Median intensity of the \dQuote{minus} elements of the border checkerboard (see \code{\link[bias]{borders}}).}
  
  \emph{The exact values returned here may change in future versions of this package!}
  
}
\references{
  Eklund, A.C. and Szallasi, Z. (2008) Correction of technical bias
  in clinical microarray data improves concordance with known
  biological information.  Genome Biology, 9:R26.
}
\seealso{ \code{\link{getBiasMetrics}} }
\examples{
  library(affydata)
  data(Dilution)
  bm2 <- getBiasMetrics2(Dilution)
  pairs(bm2)
}
\keyword{ misc }
