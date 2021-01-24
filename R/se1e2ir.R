#' Implements the rate of change of the SE1E2IR model for integration using deSolve
#'
#' The \code{se1e2ir} returns the rate of change of the SE1E2IR model that can be integrated using deSolve package
#' @param t A vector of times for integration
#' @param y A vector of state variables
#' @param param A vector of parameter values (i.e., sigma, gamma, R)
#' @export
#' @examples
#' pop <- setup_population(disease, country); calculate_cases(disease, country, pop)
# Likelihood calc for SMC --------------------------------------------
se1e2ir <- function (t, y, param) {
  sigma <- param$sigma
  gamma <- param$gamma
  beta <- param$R[floor(t)+1]*gamma
  presymp_infect <- param$presymp_infect
  if (presymp_infect) {
    beta <- param$R[floor(t)+1]*(1/(1/(2*sigma) + 1/gamma))
  }

  S <- y["S"]
  E1 <- y["E1"]
  E2 <- y["E2"]
  I <- y["I"]
  R <- y["R"]
  ## cumulative values
  CE1 <- y["CE1"]
  CE2 <- y["CE2"]
  CI <- y["CI"]

  N <- S + E1 + E2 + I + R

  foi <- beta * I / N
  if (presymp_infect) {
    foi <- beta*(E2 + I)/N #force of infection
  }
  dSdt = - S*foi
  dE1dt = + S*foi - 2*sigma*E1
  dE2dt = + 2*sigma*E1 - 2*sigma*E2
  dIdt = + 2*sigma*E2 - gamma*I
  dRdt = + gamma*I

  dCE1dt = + S*foi
  dCE2dt = + 2*sigma*E1
  dCIdt = + 2*sigma*E2


  dydt <- c(dSdt, dE1dt, dE2dt, dIdt, dRdt, dCE1dt, dCE2dt, dCIdt)

  list(dydt)
}
