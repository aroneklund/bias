\name{quantileNormalize}
\alias{quantileNormalize}
\title{ Quantile normalize a matrix }
\description{
  Quantile normalize a matrix.
}
\usage{
quantileNormalize(x, margin = 2)
}
\arguments{
  \item{x}{ A numeric matrix. }
  \item{margin}{ Margin over which to normalize (1 = rows, 2 = columns). }
}
\details{
  Quantile normalization gives all columns (or rows, if \code{margin=1}) the same distribution.
}
\value{
  A numeric matrix of the same dimension as the input.
}
\note{ 
  Tied values within each column are placed in random order.

\code{margin} should be 2 for \dQuote{standard} microarray matrices in which rows = genes.

}
\seealso{ \code{\link[affy]{normalize.quantiles}} in the \pkg{affy} package, 
          \code{normalizeQuantiles} in the \pkg{limma} package.
}
\examples{
  \dontrun{
  library(affy)
  library(affydata)
  data(Dilution)
  d <- exprs(rma(Dilution))

  par(mfrow = c(1,2))
  plotDensity(d, main = '4 arrays: before...')
  plotDensity(quantileNormalize(d), 
    main = 'and after quantile normalization')

  }
}
\keyword{ misc }
