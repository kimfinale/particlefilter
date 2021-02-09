#' Estimate parameter values using the maximum likelihood approach
#'
#' The \code{maxlik()} estimates the parameter values via the maximum likelihood approach
#' @param theta Parameters
#' @param y A vector of state variables (e.g., c(S,E1,E2,I,X))
#' @param data data to calculate the likelihoods against
#' @param data_type data to calculate the likelihoods against
#' @export
#' @import data.table
#' @examples
#' sample <- maxlik(theta = theta, y = y0, data = Rt_data, data_type = "infection", tend = 200, dt = 0.2)
maxlik <- function (theta = theta0,
                    y = y0,
                    data = Rt_data,
                    data_type = c("infection", "symptom onset", "confirmation"),
                    tend = 200,
                    dt = 0.2) {

  type <- match.arg(data_type)
  latent_var <- data.frame(matrix(NA, nrow = tend, ncol = nrow(y)))
  names(latent_var) <- y[["name"]]
  # Add initial values
  for (i in seq_along(y$name)) {
    latent_var[1, y$name[i]] <- y[name == y$name[i], val]
  }
  beta <- theta[name == "R0", val] * theta[name == "gamma", val]
  lik_traj <- rep(NA, tend) # likelihood
  lik_traj[1] <- 0 # An arbitrary value
  beta_traj <- data.frame(matrix(NA, nrow = tend, ncol = 3)) # 3 columns = MLE, lwr, and upr
  names(beta_traj) <- c("mle", "lwr", "upr")
  beta_traj[1, "mle"] <- beta
  beta_traj[1, "lwr"] <- beta - 0.1 # An arbitrary choice not to use zero or NA
  beta_traj[1, "upr"] <- beta + 0.1 # An arbitrary choice not to use zero or NA
  # begin parameter estimation using maximum likelihood
  for (t in 2:tend) {
    # DEBUG  t=2
    beta_init <- beta_traj[t-1, "mle"]
    # run optim to run observation model (negloglik), which runs the process model
    fit <- optim(par = beta_init, fn = negloglik, hessian = TRUE, theta = theta,
                 y = latent_var[t-1,], tbegin = t-1, tend = t, dt = dt, data = data)

    lik_traj[t] <- fit$value
    beta_traj[t, "mle"] <- fit$par
    stderr <- sqrt(abs(diag(solve(fit$hessian))))
    beta_traj[t, "lwr"] <- fit$par - qnorm(0.975) * stderr
    beta_traj[t, "upr"] <- fit$par + qnorm(0.975) * stderr

    latent_var[t,] <- process_model_smpl(theta = theta, y = latent_var[t-1,], tbegin = 0, tend = 1, dt = 0.2, beta = fit$par)

  }

  return (list(beta_trace = beta_traj, lik = lik_traj, latent_var = latent_var))
}

