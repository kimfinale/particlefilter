#' Estimate posterior distribution using particle filter
#'
#' The \code{particle_filter()} estimates posterior distribution using the particle filtering method
#' population for both sexes and incidence rate
#' @param theta Parameters
#' @param y A vector of state variables (e.g., c(S,E1,E2,I,X))
#' @param data data to calculate the likelihoods against
#' @param data_type data to calculate the likelihoods against
#' @param npart number of particles
#' @export
#' @import data.table
#' @examples
#' sample <- particle_filter(theta, y = y0, npart = 100, dt = 1)
particle_filter <- function (theta = NULL,
                            y = NULL,
                            data = NULL,
                            data_type = c("infection", "symptom onset", "confirmation"),
                            npart = 100,
                            tend,
                            dt = 0.2) {

  # Assumptions - using daily growth rate
  type <- match.arg(data_type)
  latent_var <- array(0,
                  dim = c(npart, tend, nrow(y)),
                  dimnames = list(NULL, NULL, y$name))
  # Add initial values
  # lapply(y$name, function(z) latent_var[, 1, z] <- y[name == z, val])

  for (i in seq_along(y$name)) {
    latent_var[, 1, y$name[i]] <- y[name == y$name[i], val]
  }
  beta <- theta[name == "R0", val] * theta[name == "gamma", val]

  # latent_var[, 1, "S"] <- y[name == "S", val]
  # latent_var[, 1, "E1"] <- y[name == "E1", val]
  # latent_var[, 1, "E2"] <- y[name == "E2", val]
  # latent_var[, 1, "I"] <- y[name == "I", val]
  # latent_var[, 1, "R"] <- y[name == "R", val]
  # latent_var[, 1, "CE1"] <- y[name == "CE1", val]
  # latent_var[, 1, "CE2"] <- y[name == "CE2", val]
  # latent_var[, 1, "CI"] <- y[name == "CI", val]

  #beta_vol <- matrix(rlnorm(npart * tend, mean = -theta[name == "betavol", val]^2/2, sd = theta[name == "betavol", val]), ncol = tend)
  beta_vol <- matrix(rnorm(npart * tend, mean = 0, sd = theta[name == "betavol", val]), nrow = tend)
  beta_vol[1,] <- exp(beta_vol[1,]) * beta # initial beta

  # Latent variables sampled from the filtered distribution
  traj <- replicate(9, matrix(NA, ncol = 1, nrow = tend), simplify = F)
  names(traj) <- c(y0[["name"]],"beta")

  # S_traj = matrix(NA, ncol = 1, nrow = tend) # npart <- 1
  # E1_traj = matrix(NA, ncol = 1, nrow = tend)
  # E2_traj = matrix(NA, ncol = 1, nrow = tend)
  # I_traj = matrix(NA, ncol = 1, nrow = tend)
  # R_traj = matrix(NA, ncol = 1, nrow = tend)
  # CE1_traj = matrix(NA, ncol = 1, nrow = tend)
  # CE2_traj = matrix(NA, ncol = 1, nrow = tend)
  # CI_traj = matrix(NA, ncol = 1, nrow = tend)
  # beta_traj = matrix(NA, ncol = 1, nrow = tend)

  wt <- matrix(NA, nrow = npart, ncol = tend) # weight (likelihood)
  wt[, 1] <- 1  # initial weights
  W <- matrix(NA, nrow = npart, ncol = tend) # normalized weights
  A <- matrix(NA, nrow = npart, ncol = tend) # Resample according to the normalized weight
  loc_sample <- rep(NA, tend)
  lik_values <- rep(NA, tend)

  # beta_rs <- wt[, 1]
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
    # W[1:npart, t] <- wt[1:npart, t] / sum(wt[1:npart, t])
    W[, t] <- wt[, t] / sum(wt[, t])
    # resample particles by sampling parent particles according to weights
    A[, t] <- sample(1:npart, prob = W[1:npart, t], replace = T)
    # Resample particles for corresponding variables
    latent_var[, t,] <- latent_var[A[, t], t,]
    beta_vol[t,] <- beta_vol[t, A[, t]] #- needed for random walk on beta
    # systematic resampling
    # beta_rs <- systematic_resampling(W[1:npart, t])
  } # end particle loop

  # Estimate likelihood:
  for (t in 1:tend) {
    lik_values[t] <- log(sum(wt[1:npart, t])) # log-likelihoods
  }
  loglik <- - tend * log(npart) + sum(lik_values) # averaged log likelihoods log(L/(npart^tend))

  # Sample latent variables:
  # A[, t] <- sample(1:npart, prob = W[1:npart, t], replace = T)
  loc <-  sample(1:npart, prob = W[1:npart, tend], replace = T)
  loc_sample[tend] <- loc[1]


  traj$beta[tend,] <- beta_vol[tend, loc_sample[tend]]
  lapply(y$name, function(z) traj[[z]][tend,] <- latent_var[loc_sample[tend], tend, z])

  # traj$S[tend,] <- latent_var[loc_sample[tend], tend, "S"]
  # traj$E1[tend,] <- latent_var[loc_sample[tend], tend, "E1"]
  # traj$E2[tend,] <- latent_var[loc_sample[tend], tend, "E2"]
  # traj$I[tend,] <- latent_var[loc_sample[tend], tend, "I"]
  # traj$R[tend,] <- latent_var[loc_sample[tend], tend, "R"]
  # traj$CE1[tend,] <- latent_var[loc_sample[tend], tend, "CE1"]
  # traj$CE2[tend,] <- latent_var[loc_sample[tend], tend, "CE2"]
  # traj$CI[tend,] <- latent_var[loc_sample[tend], tend, "CI"]

  #
  for (i in seq(tend, 2, -1)) {
    loc_sample[i-1] <- A[loc_sample[i], i] # updated indexing
    # loc_sample[i-1] <- A[1, i-1] # this works just fine
    traj$beta[i-1,] <- beta_vol[i-1, loc_sample[i-1]]
    lapply(y$name, function(z) traj[[z]][i-1,] <- latent_var[loc_sample[i-1], i-1, z])
    # traj$S[i-1,] <- latent_var[loc_sample[i-1], i-1, "S"]
    # traj$E1[i-1,] <- latent_var[loc_sample[i-1], i-1, "E1"]
    # traj$E2[i-1,] <- latent_var[loc_sample[i-1], i-1, "E2"]
    # traj$I[i-1,] <- latent_var[loc_sample[i-1], i-1, "I"]
    # traj$R[i-1,] <- latent_var[loc_sample[i-1], i-1, "R"]
    # traj$CE1[i-1,] <- latent_var[loc_sample[i-1], i-1, "CE1"]
    # traj$CE2[i-1,] <- latent_var[loc_sample[i-1], i-1, "CE2"]
    # traj$CI[i-1,] <- latent_var[loc_sample[i-1], i-1, "CI"]
  }
  # DEBUG  plot(Rep_traj[,1]-C_traj[,1])
  return (list(trace = traj, lik = loglik,
              latent_var_filtered = latent_var,
              beta_filtered = beta_vol))
}

