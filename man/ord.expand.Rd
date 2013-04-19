\name{ord.expand}
\alias{ord.expand}
\title{Expand Ordinal Data for GEE Modeling}
\description{
Expands ordinal score data into appropriate form, \emph{K}-1 new binary variables for ordinal score with \emph{K} categories, for model fitting using GEE. 
}
\usage{
ord.expand(scores, data, subjects, categories)
}
\arguments{
  \item{scores}{a character string identifying the response variable.}
  \item{data}{a dataframe in which to interpret the variables occurring in the formula.}
  \item{subjects}{a character string specifying the name of the subject variable.}
  \item{categories}{a numeric variable indicating the number of ordinal score categories.}
}
\details{
Mainly for internal use with \code{repolr} and \code{grepolr}.
}
\value{
An expanded dataframe.
}
\keyword{models}
