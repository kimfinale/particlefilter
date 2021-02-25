#' Models the process (infection transmission, recovery, etc.).
#' State variables and beta are given as a vector in which the number of elements are the number of particles
#'
#' The \code{process_model()} does SE1E2IR model integration using Euler method
#' @param tstart start time for simulation with default value of 0
#' @param tend end time with default value of 1 (i.e., one day process of SEI1I2R)
#' @param dt A step size
#' @param params Parameters
#' @param y A vector of a vector of state variables
#' @param beta A vector of beta
#' @export
#' @import data.table
#' @examples
#' pop <- process_model(params, y, tbegin = 0, tend = 1,)
# Process model for simulation --------------------------------------------
process_model <- function (params = NULL, y = NULL, tbegin = 0, tend = 1, dt = 0.2, beta = NULL) {
  S <- y[, "S"]
  E1 <- y[, "E1"]
  E2 <- y[, "E2"]
  I <- y[, "I"]
  R <- y[, "R"]
  CE1 <- y[, "CE1"]
  CE2 <- y[, "CE2"]
  CI <- y[, "CI"]

  pop <- S + E1 + E2 + I + R
  infect_rate <- beta / pop * dt
  incub_rate <- params[["sigma"]] * 2 * dt # as two compartments
  removal_rate <- params[["gamma"]] * dt

  for (ii in seq((tbegin + dt), tend, dt)) {
    foi <- I * infect_rate
    S_to_E1 <- S * foi
    E1_to_E2 <- E1 * incub_rate
    E2_to_I <- E2 * incub_rate
    I_to_R <- I * removal_rate

    # Process model for SEIR
    S <- S - S_to_E1
    E1 <- E1 + S_to_E1 - E1_to_E2
    E2 <- E2 + E1_to_E2 - E2_to_I
    I <- I + E2_to_I - I_to_R
    R <- R + I_to_R
    CE1 <- CE1 + S_to_E1
    CE2 <- CE2 + E1_to_E2
    CI <- CI + E2_to_I
  }

  nm <- c("S", "E1", "E2", "I", "R", "CE1", "CE2", "CI")
  for (s in nm) y[, s] <- eval(parse(text = s))

  return(y)
}
