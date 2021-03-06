% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_model.R
\name{process_model}
\alias{process_model}
\title{Models the process (infection transmission, recovery, etc.).
State variables and beta are given as a vector in which the number of elements are the number of particles}
\usage{
process_model(
  theta = NULL,
  y = NULL,
  tbegin = 0,
  tend = 1,
  dt = 0.2,
  beta = NULL
)
}
\arguments{
\item{theta}{Parameters}

\item{y}{A vector of a vector of state variables}

\item{tend}{end time with default value of 1 (i.e., one day process of SEI1I2R)}

\item{dt}{A step size}

\item{beta}{A vector of beta}

\item{tstart}{start time for simulation with default value of 0}
}
\description{
The \code{process_model()} does SE1E2IR model integration using Euler method
}
\examples{
pop <- process_model(theta, y, sdisease, country); calculate_cases(disease, country, pop)
}
