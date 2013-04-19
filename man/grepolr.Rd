\name{grepolr}
\alias{grepolr}
\title{
Polynomial Cut-point Models for repolr
}
\description{
Regression models are fitted to repeated ordinal scores, for the proportional odds model, using a modified version of the generalized estimating equation (GEE) method. The algorithm used estimates the correlation parameter by minimizing the generalized variance of the regression parameters at each step of the fitting algorithm. Parameter estimation is available for the uniform and first-order autoregressive correlation models, for data potentially recorded at irregularly spaced time intervals.
}
\usage{
grepolr(formula, subjects, data, times, categories, scalevalue = 1, corstr = "independence", maxiter = 10, tol = 0.001, alpha = 0.5, po.test = FALSE, fixed = FALSE, poly = 3, space = NULL)
}
\arguments{
  \item{formula}{a formula, as for other regression models.}
  \item{subjects}{a character string specifying the name of the subject variable.}
  \item{data}{a dataframe in which to interpret the variables occurring in the formula.}
  \item{times}{a vector of times which occur within subject clusters; e.g. for four evenly spaced times c(1,2,3,4).}
  \item{categories}{a numeric variable indicating the number of ordinal score categories.}
  \item{scalevalue}{a numeric variable giving the value to which the scale parameter should be fixed.}
  \item{corstr}{a character string specifying the correlation structure. The following are permitted: \dQuote{\code{ar1}}, \dQuote{\code{uniform}} and \dQuote{\code{independence}}.}
  \item{maxiter}{the maximum number of iterations.}
  \item{tol}{the tolerance used in the fitting algorithm.}
  \item{alpha}{an initial value for the correlation parameter.}
  \item{po.test}{a logical variable; if true a score test for proportional odds is reported.}
  \item{fixed}{a logical variable; if true the correlation is fixed at the initial value (\code{alpha}) during model fitting.}
  \item{poly}{a numeric variable indicating the order of the polynomial contrasts used for the cut-point model.}
  \item{space}{a vector indicating the category spacing when fitting the polynomial model; can generally be set to 1:{\code{categories}}}
}
\details{
The \code{grepolr} function fits generalized estimating equations using the same algorithm as \code{repolr}, and the majority of arguments are the same for both functions. However, rather than estimate \emph{K}-1 cut-point parameters, as function \code{repolr} does, function \code{grepolr} fits orthogonal polynomial contrasts to model the cut-point parameters. This is particularly useful for long ordinal scores, where \emph{K} is large, or where a particular form for the cut-point parameters is preferred; e.g. increasing  uniformaly at the extremes of the score scale. The function requires two additional arguments, \code{poly} and \code{space}, that indicate the order of the polynomial and the spacing between categories. The latter is normally set to 1:{\code{categories}} when all categories are observed, but can be modified to provide more realsitic models if for instance some categories are not observed.}
\value{
The fitted model is an object of class \dQuote{\code{gee}} and the correlation model has the following values in addition to the model set-up variables:
  \item{crit}{convergence criterion.}
  \item{score.test}{score test statistic.}
  \item{iter}{number of iterations.}
  \item{alpha}{estimate of the correlation parameter.}
  \item{grad1}{first derivative of generalized variance at convergence.}
  \item{grad2}{second derivative of generalized variance at convergence.}
}
\references{
Parsons, N.R. (2013). Proportional-odds models for repeated composite and long ordinal outcome scales. \emph{Statistics in Medicine} In Press.\cr
}

\keyword{models}
\seealso{
\code{\link{repolr}}, \code{\link{srepolr}}
}

\examples{
## conventional equal spacing model
data(QoL)
base.mod <- repolr(QoL~factor(Time)+factor(Treat),subjects="Patient",
         data=QoL, categories=10, times=c(1,2,3))

## linear polynomial
space <- 1:10
mod1 <- grepolr(QoL~factor(Time)+factor(Treat),subjects="Patient",
         data=QoL, categories=10, times=c(1,2,3),poly=1,space=space)
## cubic polynomial
mod3 <- grepolr(QoL~factor(Time)+factor(Treat),subjects="Patient",
         data=QoL, categories=10, times=c(1,2,3),poly=3,space=space)
summary(mod3[["gee"]])
print(paste(round(QIC(base.mod),2),round(QIC(mod3),2),round(QIC(mod1),2)))
## plot cut-point models
plot(x=space[2:10],y=base.mod$gee$coeff[1:9],type="p",xlim=c(1,10),ylab="Cut-points",
          xlab="Category Divisions")
lines(x=space[2:length(space)],y=polycuts(mod1),type="l",lty=2)
lines(x=space[2:length(space)],y=polycuts(mod3),type="l",lty=1)

## alternative spacing model
space2 <- c(1,1.5,3,5,5.5,7,7.25,8,9,10)
base.modx <- grepolr(QoL~factor(Time)+factor(Treat),subjects="Patient",
   data=QoL, categories=10, times=c(1,2,3),corstr="uniform",poly=8,space=space2)
## cubic polynomial
modx <- grepolr(QoL~factor(Time)+factor(Treat),subjects="Patient",
   data=QoL, categories=10, times=c(1,2,3),corstr="uniform",poly=3,space=space2)
summary(modx[["gee"]])
lines(x=space2[2:length(space2)],y=polycuts(modx),type="l",lty=1,col=2)
points(x=space2[2:length(space2)],y=polycuts(base.modx),type="p",col=2)

}