\name{repolr-package}
\alias{repolr-package}
\alias{repolr}
\docType{package}
\title{
Repeated Measures Proportional Odds Logistic Regression using GEE}
\description{
The package allows regression models to be fitted to repeated ordinal scores, for the proportional odds model, using a modified version of the generalized estimating equation (GEE) method. The algorithm used estimates the correlation parameter by minimizing the generalized variance of the regression parameters at each step of the fitting algorithm. Parameter estimation is available for the uniform and first-order autoregressive correlation models, for data potentially recorded at irregularly spaced time intervals and with arbitrary patterns of missingness. A sensitivity analysis is used to assess the importance of the correlation model and a score test is available that tests the assumption of proportional odds.}
\usage{
repolr(formula, subjects, data, times, categories, scalevalue = 1, corstr = "independence", maxiter = 10, tol = 0.001, alpha = 0.5, po.test = FALSE, fixed = FALSE)
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
}
\details{
The \code{repolr} function fits models for repeated ordinal scores using GEE and \code{srepolr} can be used to assess the significance of the estimated correlation parameter.

The user is required to specify: (i) a data set name (\code{data}), (ii) a model formula (\code{formula}), (iii) a cluster identification variable (\code{subjects}), (iv) a time variable (\code{time}) and (v) the number of categories used for the response variable (\code{categories}). The data set may contain records with missing data for either the response variable or the explanatory variables. The response variable must have at least three ordered categories (\emph{K} greater than or equal to 3) and, as \emph{K}-1 cut-point parameters are estimated, an additional intercept term can not be explicitly included in the model formula. A subject variable, which takes integer values from 1 to \emph{N} with no missing values allowed, indicates the data clusters (patients or experimental units) and a time variable indicates the within cluster ordering; times must be ordered integers starting from one and spaced to indicate the relative distance between successive times. For instance, four observations at equally spaced times would be entered as 1, 2, 3 and 4, whereas if the first two observations were separated by half the time interval of the other observations then coding would be 1, 2, 4 and 6. The data must be pre-sorted by time clusters within each subject, and complete, i.e. where data is missing for a particular time it must be indicated as such. The available options for the correlation model (\code{corstr}) are AR1, uniform, fixed and independence, with default setting independence. 

Additionally there are a number of other algorithm related options. A fixed scale parameter (\code{scalevalue}) is assumed to be 1, but otherwise it can be set to any positive value. The algorithm is generally robust to the initial value for alpha (default setting = 0.5), however a starting value for alpha can be set. The maximum number of iterations allowed (\code{maxiter}; default setting = 10) can be specified and convergence (\code{tol}) is monitored using the average relative change in the regression parameters at each iteration of the algorithm; the default setting (0.001) is such that the algorithm continues until this drops to below 0.1 percent. The algorithm will fail if the working correlation matrix is not positive definite and, in order to speed convergence for the AR1 and uniform correlation models, is constrained to be between 0.05 and 0.95. The proportional odds assumption can be tested, using a score test, by setting the option \code{po.test} to TRUE.

The functions require the user to have pre-installed the latest version of the R GEE solver \code{gee} and the fitted model \code{fitted.mod[[gee]]} inherits the class attributes of this package; correlation model ouput is available from \code{fitted.mod[[corr]]}. }

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
Parsons, N.R., Costa, M.L., Achten, J. and Stallard, N. (2009). Repeated measures proportional odds logistic regression analysis of ordinal score data in the statistical software package R. \emph{Computational Statistics and Data Analysis}, 53, 632-641.\cr

Parsons, N.R., Edmondson, R.N. and Gilmour, S.G. (2006). A generalized estimating equation method for fitting autocorrelated ordinal score data with an application in horticultural research. \emph{Applied Statistics}, 55, 507-524.\cr

Stiger, T.R., Barnhart, H.X., Williamson, J.M. (1999). Testing proportionality in the proportional odds model fitted with GEE. \emph{Statistics in Medicine} 18, 1419-1433.\cr
}
\keyword{models}
\seealso{
\code{\link{grepolr}}, \code{\link{srepolr}}
}


\examples{
## fit uniform correlation model to achilles tendon rupture data
data(achilles)
fitted.mod <- repolr(formula=Activity~factor(Treat)*Time, subjects="Patient",
data=achilles, categories=3, times=c(1,2,3), corstr="uniform", tol=0.001,
scalevalue=1, alpha=0.5, po.test=TRUE, fixed=FALSE)
summary(fitted.mod[["gee"]])
fitted.mod[["corr"]]

## fit ar1 correlation model to hip resurfacing pain data
data(HHSpain)
fitted.mod <- repolr(HHSpain~factor(Sex)*factor(Time), subjects="Patient",
data=HHSpain, categories=4, times=c(1,2,5), corstr="ar1", alpha=0.5)
summary(fitted.mod[["gee"]])

## plot sensitivity
sensit <- vector(length=81); alphas <- seq(0.1,0.9,0.01)
for (j in 1:81){
sensit[j] <- srepolr(mod.gee=fitted.mod,alpha=alphas[j],data=HHSpain)
}
Nsensit <-sensit/
srepolr(mod.gee=fitted.mod, alpha=fitted.mod[["corr"]]$alpha, data=HHSpain)
plot(x=alphas,y=Nsensit,type="l",lwd=2)

}
