#' Assign weights (i.e., likelihoods) to each of the value
#'
#' The \code{negloglik()} calculate the negative log likelihoods assuming case incidences are Poisson-distributed random numbers
#' @param var A vector of state variables from simulation
#' @param t current time
#' @param data data on the times of the infection transmission process (e.g., infection, symptom onset, or confirmation)
#' @param data_type c("infection", "symptom onset", "confirmation")
#' @export
#' @examples
#'  nll <- negloglik(var = latent_var, t = t, data = data, data_type = type)
negloglik <- function (pars = NULL,
                       theta = NULL,
                       y = NULL,
                       tbegin = NULL,
                       tend = NULL,
                       dt = NULL,
                       data = NULL,
                       data_type = c("infection", "symptom onset", "confirmation")) {


  beta <- pars # numeric
  y_next <- process_model_smpl(theta = theta,
                y = y,
                tbegin = tbegin,
                tend = tend,
                dt = dt,
                beta = beta)

  type <- match.arg(data_type)

  if (type == "infection") {
    case_expected <- y_next[1, "CE1"] - y[1, "CE1"]
    case_data <- data[tend, "daily_infect"]
  }
  else if (type == "symptom onset") {
    case_expected <- y_next[1, "CI"] - y[t-1, "CI"]
    case_data <- data[tend, "daily_onset"]
  }
  else if (type == "confirmation") {
    case_expected <- y_next[1, "R"] - y[t-1, "R"]
    case_data <- data[tend, "daily_confirm"]
  }

  if (!is.na(case_data)) {
    expected_val <- pmax(0, case_expected)
    log_lik <- dpois(as.integer(case_data), lambda = expected_val, log = T)
  }
  else {
    log_lik <- -Inf
  }
  return (-(log_lik)) # negative log likelihood
}

