% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_gamma_distribution.R
\name{fit_gamma_distribution}
\alias{fit_gamma_distribution}
\title{Find gamma distribution parameters to match a point estimate and confidence interval}
\usage{
fit_gamma_distribution(
  mean,
  lower,
  upper,
  confidence_level = 0.95,
  gamma_shape_range = c(0, 1000)
)
}
\arguments{
\item{mean}{The point estimate}

\item{lower}{The lower confidence limit}

\item{upper}{The upper confidence limit}

\item{confidence_level}{The coverage probability of the confidence interval}

\item{gamma_shape_range}{The range of gamma (shape) parameter values to search over (if insufficient, an error will be produced).}
}
\value{
a vector containing the estimated shape and scale parameters
}
\description{
Find gamma distribution parameters to match a point estimate and confidence interval
}
