\name{batchCorrection}
\alias{batchCorrection}
\title{ Batch correction for microarray data }
\description{
  Removes batch differences from an expression matrix.
}
\usage{
  batchCorrection(x, b)
}
\arguments{
  \item{x}{ An \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} object, or a matrix. }
  \item{b}{ A vector indicating batch membership for each column/sample in \code{x}. }
}
\details{
  The idea is to correct for \dQuote{batch-effect} systematic bias caused by inconsistent sample processing.  In this context a \dQuote{batch} may include groups of samples that were processed together, or scanned on the same date, etc.
  
  This function adjusts each row (gene) of the expression matrix independently:  For each row, a constant is added to all measurements within a batch, such that the mean within each batch is equal to the mean of the entire row.
  
}
\value{
  Either an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} or a matrix (corresponding to \code{x}).    
}
\seealso{ \code{\link{biasCorrection}} }
\examples{

  ## The Dilution data set was measured on two different scanners.
  library(affydata)
  data(Dilution)
  eset <- rma(Dilution)
  esetCorrected <- batchCorrection(eset, b = Dilution$scanner)
}
\keyword{ misc }
