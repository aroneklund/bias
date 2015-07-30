\name{getScanDate}
\alias{getScanDate}
\title{ Extract scan dates from CEL files }
\description{
  Extract the scan dates from CEL file headers.
}
\usage{
getScanDate(filenames = list.celfiles())
}
\arguments{
  \item{filenames}{ Character vector of file names. }
}
\value{
  A vector of class \code{\link{Date}}.
}
\note{ 
  This function seems to work with both binary and text CEL files.
}
\seealso{ \code{\link[affyio]{read.celfile.header}} in the \pkg{affyio} package }
\examples{
\dontrun{
     library(affydata)
     cel.path <- system.file("celfiles", package="affydata")
     filenames <- list.celfiles(path = cel.path, full.names = TRUE) 
     getScanDate(filenames)
}
}
\keyword{ misc }
