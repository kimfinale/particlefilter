% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smc_model.R
\name{smc_model}
\alias{smc_model}
\title{Compute the number of cases for a given population and incidence rate}
\usage{
smc_model(theta, nn, dt = 1)
}
\arguments{
\item{theta}{Parameters}

\item{y}{A vector of state variables (e.g., c(S,E1,E2,I,X))}

\item{simzetaA}{A vector of state variables (e.g., c(S,E1,E2,I,X))}
}
\description{
The \code{smc_model()} sequential Monte Carlo model
population for both sexes and incidence rate
}
\examples{
pop <- setup_population(disease, country); calculate_cases(disease, country, pop)
}
