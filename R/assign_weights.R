#' Assign weights (i.e., likelihoods) to each of the value
#'
#' The \code{assign_weights()} calculate the likelihoods assuming case incidences are Poisson-distributed random numbers
#' @param var A vector of state variables from simulation
#' @param t current time
#' @param data data on the times of the infection transmission process (e.g., infection, symptom onset, or confirmation)
#' @param data_type c("infection", "symptom onset", "confirmation")
#' @export
#' @examples
#'  wt <- assign_weights(var = latent_var, t = t, data = data, data_type = type)
assign_weights <- function (var,
                            t,
                            data,
                            data_type = c("infection", "symptom onset", "confirmation")) {

  type <- match.arg(data_type)
  if (type == "infection") {
    case_expected <- var[, t, "CE1"] - var[, t-1, "CE1"]
    case_data <- data[t, "daily_infect"]
  }
  else if (type == "symptom onset") {
    case_expected <- var[, t, "CI"] - var[, t-1, "CI"]
    case_data <- data[t, "daily_onset"]
  }
  else if (type == "confirmation") {
    case_expected <- var[, t, "R"] - var[, t-1, "R"]
    case_data <- data[t, "daily_confirm"]
  }

  if (!is.na(case_data)) {
    expected_val <- pmax(0, case_expected)
    log_lik <- dpois(as.integer(case_data), lambda = expected_val, log = T)
  }
  else {
    log_lik <- -Inf
  }
  return (exp(log_lik)) # convert to normal probability
}

