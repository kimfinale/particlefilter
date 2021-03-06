% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/particle_filter.R
\name{particle_filter}
\alias{particle_filter}
\title{Estimate posterior distribution of Compute the number of cases for a given population and incidence rate}
\usage{
particle_filter(
  theta = NULL,
  y = NULL,
  data = NULL,
  data_type = c("infection", "symptom onset", "confirmation"),
  nparticle = 100,
  tend,
  dt = 0.2
)
}
\arguments{
\item{theta}{Parameters}

\item{y}{A vector of state variables (e.g., c(S,E1,E2,I,X))}

\item{data}{data to calculate the likelihoods against}

\item{data_type}{data to calculate the likelihoods against}

\item{nparticle}{number of particles}
}
\description{
The \code{particle_filter()} estimates posterior distribution using the particle filtering method
population for both sexes and incidence rate
}
\examples{
sample <- particle_filter(theta, y = y0, nparticle = 100, dt = 1)
}
