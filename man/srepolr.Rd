\name{srepolr}
\alias{srepolr}
\title{Explore Sensitivity of repolr}
\description{
Determines the logarithm of the determinant of the covariance matrix of the regression parameters (the logarithm of the generalized variance) from a previous fit of \code{repolr} for an arbitrary correlation parameter.
}
\usage{
srepolr(data, mod.gee, alpha = 0.5)
}
\arguments{
  \item{data}{a dataframe in which to interpret the model fitted by \code{repolr}}
  \item{mod.gee}{a fitted model from \code{repolr}}
  \item{alpha}{a value for the correlation parameter}
}
\value{
Returns logarithm of the determinant of the covariance matrix.
}
\examples{

## fit ar1 correlation model to hip resurfacing pain data
data(HHSpain)
fitted.mod <- repolr(HHSpain~factor(Sex)*factor(Time), subjects="Patient",
data=HHSpain, categories=4, times=c(1,2,5), corstr="ar1", alpha=0.5)

## plot sensitivity
sensit <- vector(length=81); alphas <- seq(0.1,0.9,0.01)
for (j in 1:81){
sensit[j] <- srepolr(mod.gee=fitted.mod,alpha=alphas[j],data=HHSpain)
}
Nsensit <-sensit/
srepolr(mod.gee=fitted.mod, alpha=fitted.mod[["corr"]]$alpha, data=HHSpain)
plot(x=alphas,y=Nsensit,type="l",lwd=2)
}
\keyword{models}
