\name{borders}
\alias{borders}
\title{ Get border intensity values from an AffyBatch object }
\description{ Get border intensity values from an AffyBatch object }
\usage{
borders(x.batch)
}
\arguments{
  \item{x.batch}{ An \code{\link[affy:AffyBatch-class]{AffyBatch}} object.}
}
\details{ The borders of Affy arrays form a plus/minus checkerboard; this function returns their intensities. }
\value{ A list of lists:
  \item{plus}{A list of four vectors, containing the bottom, left, top, and right \dQuote{plus} values.}
  \item{minus}{A list of four vectors, containing the bottom, left, top, and right \dQuote{minus} values.}
}
\seealso{ The package \pkg{affyQCReport}. }
\examples{  

  library(affydata)
  data(Dilution)
  
  bord <- borders(Dilution)

  bord.all <- sapply(bord, unlist)

  plotDensity(log2(bord.all),
    xlab = 'log2 intensity',
    main = 'Array border intensity distribution')
  legend('topright', legend = sampleNames(Dilution),
    lty = 1:4, col = 1:4)


}
\keyword{ misc }
