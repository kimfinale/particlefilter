#' Estimate posterior distribution of Compute the number of cases for a given population and incidence rate
#'
#' The \code{particle_filter()} estimates posterior distribution using the particle filtering method
#' population for both sexes and incidence rate
#' @param theta Parameters
#' @param y A vector of state variables (e.g., c(S,E1,E2,I,X))
#' @param data data to calculate the likelihoods against
#' @param data_type data to calculate the likelihoods against
#' @param nparticle number of particles
#' @export
#' @import data.table
#' @examples
#' sample <- particle_filter(theta, y = y0, nparticle = 100, dt = 1)
particle_filter <- function(theta = NULL,
                            y = NULL,
                            data = NULL,
                            data_type = c("infection", "symptom onset", "confirmation"),
                            nparticle = 100,
                            tend,
                            dt = 0.2) {

  # data_type <- match.arg(data_type)
  # Assumptions - using daily growth rate
  type <- match.arg(data_type)
  latent_var <- array(0,
                  dim = c(nparticle, tend, nrow(y)),
                  dimnames = list(NULL, NULL, y$name))
  # Add initial values
  latent_var[, 1, "S"] <- y[name == "S", val]
  latent_var[, 1, "I"] <- y[name == "I", val]
  latent_var[, 1, "E1"] <- y[name == "E1", val]
  latent_var[, 1, "E2"] <- y[name == "E2", val]
  latent_var[, 1, "R"] <- y[name == "R", val]
  latent_var[, 1, "CE1"] <- y[name == "CE1", val]
  latent_var[, 1, "CE2"] <- y[name == "CE2", val]
  latent_var[, 1, "CI"] <- y[name == "CI", val]

  beta <- theta[name =="R0", val] * theta[name =="gamma", val]
  #beta_vol <- matrix(rlnorm(nparticle * tend, mean = -theta[name == "betavol", val]^2/2, sd = theta[name == "betavol", val]), ncol = tend)
  beta_vol <- matrix(rnorm(nparticle * tend, mean = 0, sd = theta[name == "betavol", val]), nrow = tend)
  beta_vol[1,] <- exp(beta_vol[1,]) * beta # initial beta

  # Latent variables sampled from the filtered distribution
  S_traj = matrix(NA, ncol = 1, nrow = tend)
  E1_traj = matrix(NA, ncol = 1, nrow = tend)
  E2_traj = matrix(NA, ncol = 1, nrow = tend)
  I_traj = matrix(NA, ncol = 1, nrow = tend)
  R_traj = matrix(NA, ncol = 1, nrow = tend)
  CE1_traj = matrix(NA, ncol = 1, nrow = tend)
  CE2_traj = matrix(NA, ncol = 1, nrow = tend)
  CI_traj = matrix(NA, ncol = 1, nrow = tend)
  beta_traj = matrix(NA, ncol = 1, nrow = tend)

  wt <- matrix(NA, nrow = nparticle, ncol = tend) # weight (likelihood)
  wt[, 1] <- 1  # initial weights
  W <- matrix(NA, nrow = nparticle, ncol = tend) # normalized weights
  A <- matrix(NA, nrow = nparticle, ncol = tend) # particle parent matrix
  loc_sample <- rep(NA, tend)
  lik_values <- rep(NA, tend)

  # begin particle loop
  for (t in 2:tend) {
    # DEBUG  t=2
    # Add random walk on transmission
    beta_vol[t, ] <- beta_vol[t-1, ] * exp(beta_vol[t,])
    # run process model
    latent_var[, t, ] <- process_model(theta = theta,
                                       y = latent_var[, t-1,],
                                       tbegin = t-1,
                                       tend = t,
                                       dt = dt,
                                       beta = beta_vol[t,])
    # calculate weights (likelihood)
    wt[, t] <- assign_weights(var = latent_var, t = t, data = data, data_type = type)
    # normalize particle weights
    W[1:nparticle, t] <- wt[1:nparticle, t] / sum(wt[1:nparticle, t])
    # resample particles by sampling parent particles according to weights:
    A[, t] <- sample(1:nparticle, prob = W[1:nparticle, t], replace = T)
    # Resample particles for corresponding variables
    latent_var[, t,] <- latent_var[ A[, t], t,]
    beta_vol[t,] <- beta_vol[t, A[, t]] #- needed for random walk on beta
  } # end particle loop

  # Estimate likelihood:
  for (t in 1:tend) {
    lik_values[t] <- log(sum(wt[1:nparticle, t])) # log-likelihoods
  }
  loglik <- - tend * log(nparticle) + sum(lik_values) # log-likelihoods

  # Sample latent variables:
  # A[, t] <- sample(1:nparticle, prob = W[1:nparticle, t], replace = T)
  loc <-  sample(1:nparticle, prob = W[1:nparticle, tend], replace = T)
  loc_sample[tend] <- loc[1]

  beta_traj[tend,] <- beta_vol[tend, loc_sample[tend]]
  S_traj[tend,] <- latent_var[loc_sample[tend], tend, "S"]
  E1_traj[tend,] <- latent_var[loc_sample[tend], tend, "E1"]
  E2_traj[tend,] <- latent_var[loc_sample[tend], tend, "E2"]
  I_traj[tend,] <- latent_var[loc_sample[tend], tend, "I"]
  R_traj[tend,] <- latent_var[loc_sample[tend], tend, "R"]
  CE1_traj[tend,] <- latent_var[loc_sample[tend], tend, "CE1"]
  CE2_traj[tend,] <- latent_var[loc_sample[tend], tend, "CE2"]
  CI_traj[tend,] <- latent_var[loc_sample[tend], tend, "CI"]

  for (i in seq(tend, 2, -1)) {
    loc_sample[i-1] <- A[loc_sample[i], i] # updated indexing
    # loc_sample[i-1] <- A[1, i-1] # this works just fine
    beta_traj[i-1,] <- beta_vol[i-1, loc_sample[i-1]]
    S_traj[i-1,] <- latent_var[loc_sample[i-1], i-1, "S"]
    E1_traj[i-1,] <- latent_var[loc_sample[i-1], i-1, "E1"]
    E2_traj[i-1,] <- latent_var[loc_sample[i-1], i-1, "E2"]
    I_traj[i-1,] <- latent_var[loc_sample[i-1], i-1, "I"]
    R_traj[i-1,] <- latent_var[loc_sample[i-1], i-1, "R"]
    CE1_traj[i-1,] <- latent_var[loc_sample[i-1], i-1, "CE1"]
    CE2_traj[i-1,] <- latent_var[loc_sample[i-1], i-1, "CE2"]
    CI_traj[i-1,] <- latent_var[loc_sample[i-1], i-1, "CI"]
  }
  # DEBUG  plot(Rep_traj[,1]-C_traj[,1])
  return(list(S_trace = S_traj, E1_trace = E1_traj, E2_trace = E2_traj,
              I_trace = I_traj, R_trace = R_traj, CE1_trace = CE1_traj, CE2_trace = CE2_traj,
              CI_trace = CI_traj, beta_trace = beta_traj, lik = loglik))
}

