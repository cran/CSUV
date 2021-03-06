% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/csuv.r
\name{get.csuv.unique.fit}
\alias{get.csuv.unique.fit}
\title{Helper function, please do not use it}
\usage{
get.csuv.unique.fit(
  X,
  Y,
  intercept,
  method.names,
  B,
  fit.percent,
  current.fit = NULL,
  num.core = 1
)
}
\arguments{
\item{X}{covariates (n times p matrix, n: number of entries, p: number of covariates)}

\item{Y}{response (vector with n entries)}

\item{intercept}{TRUE to fit the data with an intercept, FALSE to fit the data without an intercept}

\item{method.names}{vector of method names to be used in CSUV. Choose among "lasso", "elastic", "relaxo", "mcp" and "scad". Default is to use all methods listed above}

\item{B}{number of subsampling. Default is 100}

\item{fit.percent}{percentage of observations used in fitting in CSUV}

\item{current.fit}{(optional) all fitted models}

\item{num.core}{number of cores to use. Default is 1 (i.e. no parallel running)}
}
\value{
a list of current fit
}
\description{
Helper function, please do not use it
}
