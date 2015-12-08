\name{begonia}
\alias{begonia}
\docType{data}
\title{Begonia pot plant quality scores}
\description{
Begonia pot plant quality scores (for 2 varieties and 3 transport chains), during 5 weeks in simulated shelf-life conditions (temperature, light and irrigation). Quality scores were originally made on an ordinal scale from 1 to 10 (highest quality). However, only categories 2 to 9 were used, so these are re-coded to scale from 1 to 8. In addition to overall quality scores, a range of plant physiological characteristics were all observed. 
}
\usage{data(begonia)}
\format{
A data frame with 720 observations on the following variables.
  \describe{
    \item{\code{Pot}}{Pot plant identifier; 1 to 144.}
    \item{\code{Plot}}{Location in growing compartment; 1 to 48.}
    \item{\code{Week}}{Week number in simulated shelf-life; 1 to 5.}
    \item{\code{Temp}}{Temperature in compartment; \code{16C} or \code{21C}.}
    \item{\code{Light}}{Light level in compartment; \code{High} or \code{Low}.}
    \item{\code{Chain}}{Transport chain; \code{Cold}, \code{Comm} (commercial) or \code{Optm} (optimum).}
    \item{\code{Irrig}}{Irrigation in compartment; \code{Cont} (control) or \code{Fluct} (fluctuating).}    
    \item{\code{Variety}}{Variety; \code{Balli} or \code{Batik}.}
    \item{\code{Qual}}{Quality score; 1 to 8.}
    \item{\code{FDrop}}{Count of dropped flowers.}
    \item{\code{CBDrop}}{Count of dropped coloured buds.}
    \item{\code{GBDrop}}{Count of dropped green buds.}
    \item{\code{FDam}}{Count of damaged flowers.}    
    \item{\code{FSing}}{Count of single flowers.}
    \item{\code{FDoub}}{Count of double flowers.}
    \item{\code{LYell}}{Count of yellow leaves.}
  }
}
\source{
Parsons NR, Edmondson RN, Gilmour SG. A generalized estimating equation method for fitting autocorrelated ordinal score data with an application in horticultural research. \emph{Applied Statistics} 2006; 55:507-524.\cr
}
\examples{
data(begonia)
head(begonia)
}
\keyword{datasets}
