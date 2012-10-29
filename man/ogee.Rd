\name{ogee}
\alias{ogee}
\title{Internal Function to fit Proportional Odds Model}
\description{
An internal function called by \code{repolr}, that would not normally be called by user. This is a modified version of \code{gee} for fitting the proportional odds model. 
}
\usage{
ogee(formula = formula(data), id = id, data = parent.frame(), subset, na.action, R = NULL, b = NULL, tol = 0.001, maxiter = 25, family = gaussian, corstr = "independence", Mv = 1, silent = TRUE, contrasts = NULL, scale.fix = FALSE, scale.value = 1, v4.4compat = FALSE)
}

\arguments{
  \item{formula}{as \code{gee}}
  \item{id}{as \code{gee}}
  \item{data}{as \code{gee}}
  \item{subset}{as \code{gee}}
  \item{na.action}{as \code{gee}}
  \item{R}{as \code{gee}}
  \item{b}{as \code{gee}}
  \item{tol}{as \code{gee}}
  \item{maxiter}{as \code{gee}}
  \item{family}{as \code{gee}}
  \item{corstr}{as \code{gee}}
  \item{Mv}{as \code{gee}}
  \item{silent}{as \code{gee}}
  \item{contrasts}{as \code{gee}}
  \item{scale.fix}{as \code{gee}}
  \item{scale.value}{as \code{gee}}
  \item{v4.4compat}{as \code{gee}}
}
\details{
Not generally to be called by user.
}
\keyword{models}
