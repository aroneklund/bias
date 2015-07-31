\name{calcCM}
\alias{calcCM}
\title{ Calculate a binned correlation matrix }
\description{
  Calculate the bin-wise averages of the correlations between rows of a matrix.  This may be useful for visualizing gross trends in the correlation matrix.
}
\usage{
calcCM(x, n = 50, 
       use = 'all', method = 'pearson', 
       FUN = median, ...) 
}
\arguments{
  \item{x}{ A matrix, or an exprSet. }
  \item{n}{ Number of bins. }
  \item{use, method}{ Passed directly to \code{\link{cor}}. }
  \item{FUN}{ Function to summarize values in each bin. }
  \item{\dots}{ Additional argments are passed to \code{FUN}. }
}
\details{
The rows of the input matrix \code{x} are split into \code{n} bins.  For each possible pair of bins, all possible correlation coefficients between rows in the two bins are calculated.  The median (or other \code{FUN}) correlation for each pair of bins is returned as a matrix.  

For the diagonals of the returned matrix, which correspond to the correlation between a bin and itself, only the correlations between non-identical probe sets are included.
}
\value{
  A list with the following components:
  \item{x, y }{Numeric vectors of length \code{n} + 1; the breakpoints defining the bins.  \code{x} and \code{y} are equivalent; both are included so the result can be passed directly to plotting functions.}
  \item{z }{\code{n} x \code{n} numeric matrix; the median (or other \code{FUN}) correlation for each bin.}
  \item{count }{\code{n} x \code{n} numeric matrix; the number of observations in each bin.}
  \item{call }{ Character; the function call.}
}
\references{
  Eklund, A.C. and Szallasi, Z. (2008) Correction of technical bias
  in clinical microarray data improves concordance with known
  biological information.  Genome Biology, 9:R26.
}
\note{ 
  The results will be more informative if the input matrix \code{x} is sorted by some criterion first.  For example, try sorting by (row-wise) mean, variance, probe sequence composition...
}
\seealso{ \code{colorgram} in the in the \pkg{squash} package for visualization}
\examples{

  library(affydata)
  data(Dilution)
  d <- exprs(rma(Dilution))

  ## Sort the rows (probe sets) by mean expression, and calculate correlation matrix
  dsort <- d[order(rowMeans(d)),]
  cm <- calcCM(dsort)
  
  ## Visualize the correlation matrix
  image(cm, col = rainbow(20),
    xlab = 'Expression rank', ylab = 'Expression rank')
  
  ## Maybe use the "squash" package for more control of the plot
  \dontrun{
  par(mar = c(7, 4, 4, 1) + 0.1)
  map <- makecmap(cm$z, colFn = blueorange, 
            symm = TRUE, n = 20)
  colorgram(cm, map = map, 
        xlab = 'Expression rank', 
        ylab = 'Expression rank',
        zlab = 'Median Pearson correlation')
  }

}
\keyword{ misc }
