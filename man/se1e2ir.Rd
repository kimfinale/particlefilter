% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/se1e2ir.R
\name{se1e2ir}
\alias{se1e2ir}
\title{Implements the rate of change of the SE1E2IR model for integration using deSolve}
\usage{
se1e2ir(t, y, param)
}
\arguments{
\item{t}{A vector of times for integration}

\item{y}{A vector of state variables}

\item{param}{A vector of parameter values (i.e., sigma, gamma, R)}
}
\description{
The \code{se1e2ir} returns the rate of change of the SE1E2IR model that can be integrated using deSolve package
}
\examples{
pop <- setup_population(disease, country); calculate_cases(disease, country, pop)
}
