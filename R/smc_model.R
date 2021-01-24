#' Compute the number of cases for a given population and incidence rate
#'
#' The \code{smc_model()} sequential Monte Carlo model
#' population for both sexes and incidence rate
#' @param theta Parameters
#' @param y A vector of state variables (e.g., c(S,E1,E2,I,X))
#' @param simzetaA A vector of state variables (e.g., c(S,E1,E2,I,X))
#' @export
#' @import data.table
#' @examples
#' pop <- setup_population(disease, country); calculate_cases(disease, country, pop)
# SMC function --------------------------------------------
smc_model <- function(theta, nn, dt=1) {
  # nn = 100;   dt <- 0.25
  # Assumptions - using daily growth rate
  t_period <- ttotal <- as.numeric(sum(diff(DAT$date))) + 1
  t_length <- ttotal

  storeL <- array(0, dim = c(nn, t_length, length(theta_init_names)),
                  dimnames = list(NULL, NULL, theta_init_names))
  # Add initial condition
  storeL[, 1, "I"] <- theta[["init_cases"]]
  storeL[, 1, "S"] <- theta[["pop"]] - theta[["init_cases"]]

  #simzeta <- matrix(rlnorm(nn*t_length, mean = -theta[["betavol"]]^2/2, sd = theta[["betavol"]]),ncol=ttotal)
  simzeta <- matrix(rnorm(nn*t_length, mean = 0, sd = theta[["betavol"]]), nrow = ttotal)
  simzeta[1,] <- exp(simzeta[1,])*theta[["beta"]] # define IC

  # Latent variables
  S_traj = matrix(NA, ncol=1, nrow=ttotal)
  E1_traj = matrix(NA, ncol=1, nrow=ttotal)
  E2_traj = matrix(NA, ncol=1, nrow=ttotal)
  I_traj = matrix(NA, ncol=1, nrow=ttotal)
  # I_imported_traj = matrix(NA,ncol=1,nrow=ttotal)

  CI_traj = matrix(NA, ncol=1, nrow = ttotal)
  X_traj = matrix(NA, ncol=1, nrow = ttotal)
  beta_traj = matrix(NA, ncol=1, nrow = ttotal);
  w <- matrix(NA, nrow = nn, ncol = ttotal);
  w[,1] <- 1  # weights
  W <- matrix(NA, nrow = nn, ncol = ttotal) # normialized weights
  A <- matrix(NA, nrow = nn, ncol = ttotal) # particle parent matrix
  l_sample <- rep(NA, ttotal)
  lik_values <- rep(NA, ttotal)

  # Iterate through steps
  for (tt in 2:ttotal) {
    # DEBUG  tt=2
    # Add random walk on transmission ?
    simzeta[tt, ] <- simzeta[tt-1, ]*exp(simzeta[tt,])
    # run process model
    storeL[, tt, ] <- process_model(tt-1, tt, dt, theta, storeL[,tt-1,], simzeta[tt,])
    # cat( "tt =", tt,  ", median I at t =", median(storeL[ , tt, "I" ]) , "\n")
    # calculate weights
    w[ , tt] <- assign_weights(data_list = DAT, storeL, nn, theta, tt)

    # normalise particle weights
    sum_weights <- sum(w[1:nn, tt])
    W[1:nn, tt] <- w[1:nn, tt]/sum_weights

    # resample particles by sampling parent particles according to weights:
    A[, tt] <- sample(1:nn, prob = W[1:nn,tt], replace = T)
    # Resample particles for corresponding variables
    storeL[,tt,] <- storeL[ A[, tt] , tt,]
    simzeta[tt,] <- simzeta[tt, A[, tt]] #- needed for random walk on beta

  } # END PARTICLE LOOP

  # Estimate likelihood:
  for (tt in 1:ttotal) {
    lik_values[tt] = log(sum(w[1:nn,tt])) # log-likelihoods
  }

  loglik = -ttotal*log(nn) + sum(lik_values) # log-likelihoods

  # Sample latent variables:
  locs <-  sample(1:nn, prob = W[1:nn,tt], replace = T)
  l_sample[ttotal] <- locs[1]
  X_traj[ttotal,] <- storeL[l_sample[ttotal], ttotal, "X"]
  S_traj[ttotal,] <- storeL[l_sample[ttotal], ttotal, "S"]
  I_traj[ttotal,] <- storeL[l_sample[ttotal], ttotal, "I"]
  # I_imported_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"I_imported"]
  CI_traj[ttotal,] <- storeL[l_sample[ttotal], ttotal, "CI"]
  beta_traj[ttotal,] <- simzeta[ttotal, l_sample[ttotal]]

  for (ii in seq(ttotal, 2, -1)) {
    l_sample[ii-1] <- A[l_sample[ii], ii] # have updated indexing
    X_traj[ii-1,] <- storeL[l_sample[ii-1], ii-1, "X"]
    S_traj[ii-1,] <- storeL[l_sample[ii-1], ii-1, "S"]
    I_traj[ii-1,] <- storeL[l_sample[ii-1], ii-1, "I"]
    # I_imported_traj[ii-1,] <- storeL[l_sample[ii-1],ii-1,"I_imported"]
    CI_traj[ii-1,] <- storeL[l_sample[ii-1], ii-1, "CI"]
    beta_traj[ii-1,] <- simzeta[ii-1, l_sample[ii-1]]
  }
  # DEBUG  plot(Rep_traj[,1]-C_traj[,1])
  # return( list( S_trace=S_traj, I_trace=I_traj, I_imported_trace=I_imported_traj, CI_trace=CI_traj, X_trace=X_traj, beta_trace=beta_traj, lik=likelihood0 ) )
  return(list(S_trace = S_traj, I_trace = I_traj, CI_trace = CI_traj, X_trace = X_traj, beta_trace = beta_traj, lik = loglik))
}

