\name{QIC}
\alias{QIC}
\title{Quasilikelihood Information Criterion}
\description{
The quasilikelihood information criterion (QIC) developed by Pan (2001) is a modification of the Akaike information criterion (AIC) for models fit by GEE. 
}
\usage{
QIC(x)
}
\arguments{
  \item{x}{is a fitted model using \code{grepolr}}
}
\value{
The quasilikelihood for the fitted model.
}
\references{
Pan, W. (2001). Akaikes information criterion in generalized estimating equations. \emph{Biometrics}, 57, 120-125.\cr
}
\keyword{models}

\examples{
## linear polynomial
data(QoL)
space <- 1:10
mod1 <- grepolr(QoL~factor(Time)+factor(Treat),subjects="Patient",
         data=QoL, categories=10, times=c(1,2,3),poly=1,space=space)
## cubic polynomial
mod3 <- grepolr(QoL~factor(Time)+factor(Treat),subjects="Patient",
         data=QoL, categories=10, times=c(1,2,3),poly=3,space=space)
QIC(mod1)
QIC(mod3)

}

