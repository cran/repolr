\name{polycuts}
\alias{polycuts}
\title{
Provide Estimates of Cut-point Parameters for grepolr
}
\description{
After fitting model using \code{grepolr}, function \code{polycuts} gives estimates of \emph{K}-1 cut-point parameters, based on polynomial model from fit of \code{grepolr}, that provide useful way of comparing alternate polynomial models.
}
\usage{
polycuts(x)
}
\arguments{
  \item{x}{is a fitted model using \code{grepolr}}
}
\value{
A vector of \emph{K}-1 cut-point parameters.
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
polycuts(mod1)
polycuts(mod3)

}