#' This function implements forward filtering-backward smoothing algorithm
#'
#' The \code{particle_filter()} estimates posterior distribution using the particle filtering method
#' population for both sexes and incidence rate
#' @param params Parameters
#' @param y A vector of state variables (e.g., c(S,E1,E2,I,X))
#' @param data data to calculate the likelihoods against
#' @param data_type data to calculate the likelihoods against
#' @param npart number of particles
#' @export
#' @import data.table
#' @examples
#' sample <- particle_filter()
particle_filter <- function (params = theta,
                            y = y0,
                            data = Rt_data,
                            data_type = c("infection", "symptom onset", "confirmation"),
                            npart = 100,
                            tend = 200,
                            dt = 0.2,
                            systematic_resampling = FALSE,
                            filter_traj = TRUE) {

  if (missing(npart) || is.null(npart)) {
    stop("Number of particles (npart) must be specified.")
  }
   # Assumptions - using daily growth rate
  type <- match.arg(data_type)
  latent_var <- array(0,
                  dim = c(npart, tend, length(y)),
                  dimnames = list(NULL, NULL, names(y)))
  # Add initial values
  for (nm in names(y)) {
    latent_var[, 1, nm] <- y[[nm]]
  }
  beta <- params[["R0"]] * params[["gamma"]]
  #beta_vol <- matrix(rlnorm(npart * tend, mean = -params[name == "betavol", val]^2/2, sd = params[name == "betavol", val]), ncol = tend)
  beta_vol <- matrix(rnorm(npart * tend, mean = 0, sd = params[["betavol"]]), nrow = tend)
  beta_vol[1,] <- beta * exp(beta_vol[1,]) # initial beta

  wt <- matrix(NA, nrow = npart, ncol = tend) # weight (likelihood)
  wt[, 1] <- 1 / npart  # initial weights
  W <- matrix(NA, nrow = npart, ncol = tend) # normalized weights
  A <- matrix(NA, nrow = npart, ncol = tend) # Resample according to the normalized weight

  # begin particle loop
  for (t in 2:tend) {
    # DEBUG  t=2
     beta_vol[t, ] <- beta_vol[t-1, ] * exp(beta_vol[t, ])
     # beta_vol[t, ] <- beta * exp(rnorm(npart, mean = 0, sd = params[["betavol"]]))
     # run process model
     latent_var[, t, ] <- process_model(params = params,
                                         y = latent_var[, t-1, ],
                                         tbegin = t-1,
                                         tend = t,
                                         dt = dt,
                                         beta = beta_vol[t,])
    # calculate weights (likelihood)
    wt[, t] <- assign_weights(var = latent_var, t = t, data = data, data_type = type)
    # normalize particle weights
    W[, t] <- wt[, t] / sum(wt[, t])
    # resample particles by sampling parent particles according to weights
    A[, t] <- sample(1:npart, prob = W[1:npart, t], replace = T)
    if (systematic_resampling) {
      A[, t] <- systematic_resampling(W[1:npart, t])
    }
    # Resample particles for corresponding variables
    latent_var[, t,] <- latent_var[A[, t], t,]
    beta_vol[t,] <- beta_vol[t, A[, t]] #- needed for random walk on beta
  } # end particle loop

  # Marginal likelihoods:
  lik_values <- rep(NA, tend)
  for (t in 1:tend) {
    lik_values[t] <- log(sum(wt[1:npart, t])) # log-likelihoods
  }
  loglik <- - tend * log(npart) + sum(lik_values) # averaged log likelihoods log(L/(npart^tend))

  traj <- replicate(length(y)+1, matrix(NA, ncol = 1, nrow = tend), simplify = F)
  names(traj) <- c(names(y), "beta")
  if (filter_traj) {
    # Latent variables sampled from the filtered distribution (backward smoothing)
    loc <- rep(NA, tend)
    # loc_tend <-  sample(1:npart, prob = W[1:npart, tend], replace = T)
    loc[tend] <-  sample.int(npart, size = 1, prob = W[, tend], replace = T)
    # particle for the last time step
    traj$beta[tend,] <- beta_vol[tend, loc[tend]]
    traj[["S"]][tend,] <- latent_var[loc[tend], tend, "S"]
    traj[["E1"]][tend,] <- latent_var[loc[tend], tend, "E1"]
    traj[["E2"]][tend,] <- latent_var[loc[tend], tend, "E2"]
    traj[["I"]][tend,] <- latent_var[loc[tend], tend, "I"]
    traj[["R"]][tend,] <- latent_var[loc[tend], tend, "R"]
    traj[["CE1"]][tend,] <- latent_var[loc[tend], tend, "CE1"]
    traj[["CE2"]][tend,] <- latent_var[loc[tend], tend, "CE2"]
    traj[["CI"]][tend,] <- latent_var[loc[tend], tend, "CI"]

    lapply(names(y), function(x) traj[[x]][tend,] <- latent_var[loc[tend], tend, x])
    # update backward
    for (i in seq(tend, 2, -1)) {
      loc[i-1] <- A[loc[i], i]
      traj$beta[i-1,] <- beta_vol[i-1, loc[i-1]]
      # invisible(lapply(names(y), function(x) traj[[x]][i-1,] <- latent_var[loc[i-1], i-1, x]))
      traj[["S"]][i-1,] = latent_var[loc[i-1], i-1, "S"]
      traj[["E1"]][i-1,] = latent_var[loc[i-1], i-1, "E1"]
      traj[["E2"]][i-1,] = latent_var[loc[i-1], i-1, "E2"]
      traj[["I"]][i-1,] = latent_var[loc[i-1], i-1, "I"]
      traj[["R"]][i-1,] = latent_var[loc[i-1], i-1, "R"]
      traj[["CE1"]][i-1,] = latent_var[loc[i-1], i-1, "CE1"]
      traj[["CE2"]][i-1,] = latent_var[loc[i-1], i-1, "CE2"]
      traj[["CI"]][i-1,] = latent_var[loc[i-1], i-1, "CI"]
    }
  }
  ## code snippet from pomp package
  # if (filter.traj) { ## select a single trajectory
  #   b <- sample.int(n=length(weights),size=1L,replace=TRUE)
  #   filt.t[,1L,ntimes+1] <- xparticles[[ntimes]][,b]
  #   for (nt in seq.int(from=ntimes-1,to=1L,by=-1L)) {
  #     b <- pedigree[[nt+1]][b]
  #     filt.t[,1L,nt+1] <- xparticles[[nt]][,b]
  #   }
  #   if (times[2L] > times[1L]) {
  #     b <- pedigree[[1L]][b]
  #     filt.t[,1L,1L] <- init.x[,b]
  #   } else {
  #     filt.t <- filt.t[,,-1L,drop=FALSE]
  #   }
  # }
  # DEBUG  plot(Rep_traj[,1]-C_traj[,1])
  return (list(trace = traj, lik_marginal = lik_values,
               lik_overall_average = loglik, latent_var_filtered = latent_var,
              beta_filtered = beta_vol))
}

