\name{biasCorrection}
\alias{biasCorrection}
\title{ Linear bias correction for microarray data }
\description{
  Remove bias from an expression matrix along specified components.
}
\usage{
  biasCorrection(x, metrics)
}
\arguments{
  \item{x}{ An \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} object, or a matrix. }
  \item{metrics}{  A data frame with bias metrics for each array. }
}
\details{
  Each individual gene is fit to a multivariate linear model with \code{metrics} as the independent variables.
  The residuals from the fit are retained as the \dQuote{bias corrected} values.
}
\value{
  Either an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} or a matrix (corresponding to the class of \code{x}).
  
}
\references{
  Eklund, A.C. and Szallasi, Z. (2008) Correction of technical bias
  in clinical microarray data improves concordance with known
  biological information.  Genome Biology, 9:R26.
}
\seealso{ \code{\link{getBiasMetrics}} }
\examples{

  ## For testing/example purposes only:
  ## This data set has too few samples to do meaningful bias correction
  library(affydata)
  data(Dilution)
  eset <- rma(Dilution)
  bm <- getBiasMetrics(Dilution, eset)
  esetCorrected <- biasCorrection(eset, bm)

}
\keyword{ misc }
