\name{summary.gee}
\alias{summary.gee}
\title{Summary of GEE Models}
\description{
Function to summarise fit of GEE model.
}
\usage{
\method{summary}{gee}(object, working.correlation = FALSE, \dots)
}
\arguments{
  \item{object}{fitted model}
  \item{working.correlation}{logical; if \code{TRUE}, the working correlation matrix is returned and printed}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
Default printing does not show working correlation matrix, which can be large for ordinal data.
}
\keyword{models}

 