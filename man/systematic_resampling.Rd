% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/systematic_resampling.R
\name{systematic_resampling}
\alias{systematic_resampling}
\title{Compute the number of cases for a given population and incidence rate}
\usage{
systematic_resampling(weights)
}
\arguments{
\item{weights}{A vector of weights}
}
\description{
The \code{systematic_resampling()} does SE1E2IR model calculation for a given time period. This is used
population for both sexes and incidence rate
}
\examples{
pop <- setup_population(disease, country); calculate_cases(disease, country, pop)
}
