\name{getBiasMetrics}
\alias{getBiasMetrics}
\title{ Calculate bias metrics }
\description{
  Calculate bias metrics from an AffyBatch object. These are the original four bias metrics as described in the Genome Biology article.
}
\usage{
  getBiasMetrics(x.batch, x.rma = rma(x.batch))
}
\arguments{
  \item{x.batch}{ An \code{\link[affy:AffyBatch-class]{AffyBatch}}.}
  \item{x.rma}{ An \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} calculated from \code{x.batch} (optional).}
}
\details{
  This function calculates four bias metrics for each microarray in a data set.  
  
  Three of the bias metrics are derived from the raw intensities in \code{x.batch} and one is derived from the probeset-summarized expression values in \code{x.rma}.  If \code{x.rma} is missing, the expression values will be computed using \code{\link[affy]{rma}}. 

}
\value{
  A data frame containing bias metrics for each array.
  \item{pm.median}{Median intensity of all perfect-match probes on the array.}
  \item{pm.IQR}{Interquartile range of all perfect-match probes on the array.}
  \item{rma.IQR}{Interquartile range of all expression values on the array.}
  \item{degradation}{The slope returned by \code{\link[affy]{AffyRNAdeg}}; indicates the relative preference for 3'-biased probes in a probe set.}
     
}
\references{
  Eklund, A.C. and Szallasi, Z. (2008) Correction of technical bias
  in clinical microarray data improves concordance with known
  biological information.  Genome Biology, 9:R26.
}
\seealso{ \code{\link{biasCorrection}}; also \code{\link[simpleaffy]{qc}} in the \pkg{simpleaffy} package. 
}
\examples{
  library(affydata)
  data(Dilution)
  bm <- getBiasMetrics(Dilution)
  pairs(bm)
}
\keyword{ misc }
