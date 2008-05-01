\name{HHSpain}
\alias{HHSpain}
\docType{data}
\title{Harris Hip Pain Scores}
\description{
Harris hip pain scores, measured on 4 a point ordinal score scale, post-surgery for 58 patients at 3 occasions. 
}
\usage{data(HHSpain)}
\format{
  A data frame with 174 observations on the following 4 variables.
  \describe{
    \item{\code{Patient}}{a patient identifier variable.}
    \item{\code{Sex}}{a factor with levels \code{F} \code{M}}
    \item{\code{Time}}{recorded at baseline (1), two years (2) and five years (5) post-surgery.}
    \item{\code{HHSpain}}{a pain score coded as none (1), slight (2), mild (3) and moderate or marked (4)}
  }
}
\source{
Insert reference when paper is published.
}
\examples{
data(HHSpain)
str(HHSpain)
}
\keyword{datasets}
