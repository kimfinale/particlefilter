#' Assign weights (i.e., likelihoods) to each of the valuee
#'
#' The \code{assign_weights()} calculate the likelihoods assuming case incidences are Poisson-distributed random numbers
#' @param var A vector of state variables from simulation
#' @param t current time
#' @param data data on the times of the infection transmisison process (e.g., infection, symptom onset, or confirmation)
#' @param data_type c("infection", "symptom onset", "confirmation")
#' @export
#' @examples
#' pop <- setup_population(disease, country); calculate_cases(disease, country, pop)
assign_weights <- function (var,
                            t,
                            data,
                            data_type = c("infection", "symptom onset", "confirmation")) {
  #  # Calculation likelihood
  # loglik <- dpois(y_lam,lambda=x_lam,log=T)
  # loglikSum <- rowSums(matrix(loglik, nrow = nn, byrow = T))
  # exp(loglikSum) # convert to normal probability
  #
  # Gather data
  # case_data_t <- data_list$case_confirmed_daily[ t ]
  # # Gather variables
  type <- match.arg(data_type)
  # case_data_t <- data_list$case_local[ t ] + data_list$case_imported[ t ]
  if (type == "infection") {
    case_diff <- var[, t, "CE1"] - var[, t-1, "CE1"]
    case_data <- data[t, "daily_infect"]
  }
  else if (type == "symptom onset") {
    case_diff <- var[, t, "CI"] - var[, t-1, "CI"]
    case_data <- data[t, "daily_onset"]
  }
  else if (type == "confirmation") {
    case_diff <- var[, t, "R"] - var[, t-1, "R"]
    case_data <- data[t, "daily_confirm"]
  }

  if (!is.na(case_data)) {
    expected_val <- pmax(0, case_diff)
    log_lik <- dpois(as.integer(case_data), lambda = expected_val, log = T)
  }
  else {
    log_lik <- -Inf
  }
  return (exp(log_lik)) # convert to normal probability
}

