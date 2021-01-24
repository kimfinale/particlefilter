# SMC function --------------------------------------------
smc_model <- function(theta, nn, dt=1) {
  # nn = 100;   dt <- 0.25
  # Assumptions - using daily growth rate
  t_period <- ttotal <- as.numeric(sum(diff(DAT$date))) + 1
  t_length <- ttotal

  storeL <- array(0, dim = c(nn, t_length, length(theta_init_names)),
                  dimnames = list(NULL, NULL, theta_init_names))
  # Add initial condition
  storeL[,1,"I"] <- theta[["init_cases"]]
  storeL[,1,"S"] <- theta[["pop"]] - theta[["init_cases"]]

  #simzeta <- matrix(rlnorm(nn*t_length, mean = -theta[["betavol"]]^2/2, sd = theta[["betavol"]]),ncol=ttotal)
  simzeta <- matrix(rnorm(nn*t_length, mean = 0, sd = theta[["betavol"]]), nrow = ttotal)
  simzeta[1,] <- exp(simzeta[1,])*theta[["beta"]] # define IC

  # Latent variables
  S_traj = matrix(NA, ncol=1, nrow=ttotal)
  E1_traj = matrix(NA, ncol=1, nrow=ttotal)
  E2_traj = matrix(NA, ncol=1, nrow=ttotal)
  I_traj = matrix(NA, ncol=1, nrow=ttotal)
  # I_imported_traj = matrix(NA,ncol=1,nrow=ttotal)

  CI_traj = matrix(NA, ncol=1, nrow=ttotal)
  X_traj = matrix(NA, ncol=1, nrow=ttotal)
  beta_traj = matrix(NA, ncol=1, nrow=ttotal);
  w <- matrix(NA, nrow = nn, ncol = ttotal);
  w[,1] <- 1  # weights
  W <- matrix(NA, nrow = nn, ncol = ttotal) # normialized weights
  A <- matrix(NA, nrow = nn, ncol = ttotal) # particle parent matrix
  l_sample <- rep(NA, ttotal)
  lik_values <- rep(NA, ttotal)

  # Iterate through steps
  for(tt in 2:ttotal){
    # DEBUG  tt=2
    # Add random walk on transmission ?
    simzeta[ tt, ] <- simzeta[ tt-1, ]*exp( simzeta[tt,] )
    # run process model
    storeL[ , tt, ] <- process_model(tt-1, tt, dt, theta, storeL[,tt-1,], simzeta[tt,])
    # cat( "tt =", tt,  ", median I at t =", median(storeL[ , tt, "I" ]) , "\n")
    # calculate weights
    w[ ,tt ] <- AssignWeights(data_list = DAT, storeL, nn, theta, tt)

    # normalise particle weights
    sum_weights <- sum(w[1:nn, tt])
    W[1:nn, tt] <- w[1:nn, tt]/sum_weights

    # resample particles by sampling parent particles according to weights:
    A[, tt] <- sample(1:nn, prob = W[1:nn,tt], replace = T)
    # Resample particles for corresponding variables
    storeL[,tt,] <- storeL[ A[, tt] ,tt,]
    simzeta[tt,] <- simzeta[tt, A[, tt]] #- needed for random walk on beta

  } # END PARTICLE LOOP

  # Estimate likelihood:
  for (tt in 1:ttotal) {
    lik_values[tt] = log(sum(w[1:nn,tt])) # log-likelihoods
  }

  likelihood0 = -ttotal*log(nn)+ sum(lik_values) # log-likelihoods

  # Sample latent variables:
  locs <-  sample(1:nn,prob = W[1:nn,tt],replace = T)
  l_sample[ttotal] <- locs[1]
  X_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"X"]
  S_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"S"]
  I_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"I"]
  # I_imported_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"I_imported"]
  CI_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"CI"]
  beta_traj[ttotal,] <- simzeta[ttotal,l_sample[ttotal]]

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
  return( list( S_trace=S_traj, I_trace=I_traj, CI_trace=CI_traj, X_trace=X_traj, beta_trace=beta_traj, lik=likelihood0 ) )
}
# Likelihood calc for SMC --------------------------------------------
AssignWeights <- function( data_list, storeL, nn, theta, tt ){
  #  # Calculation likelihood
  # loglik <- dpois(y_lam,lambda=x_lam,log=T)
  # loglikSum <- rowSums(matrix(loglik,nrow=nn,byrow=T))
  # exp(loglikSum) # convert to normal probability
  #
  # Gather data
  # case_data_tt <- data_list$case_confirmed_daily[ tt ]
  # # Gather variables


  # case_data_tt <- data_list$case_local[ tt ] + data_list$case_imported[ tt ]

  case_diff <- storeL[,tt,"CI"] - storeL[,tt-1,"CI"]
  case_data_tt <- data_list$case_local[ tt ]

  case_val <- pmax( 0, case_diff )

  # Local confirmed cases (by onset)

  if( !is.na(case_data_tt) ){
    expected_val <- case_val # scale by reporting proportion and known onsets
    log_lik_sum <- dpois( case_data_tt, lambda = expected_val, log=T )
  } else {
    log_lik_sum <- 0
  }

  exp(log_lik_sum) # convert to normal probability
}



# SMC function --------------------------------------------
smc_model2 <- function( theta, nn, dt=1 ){
  library(truncnorm)
  # nn = 100;   dt <- 0.25
  # Assumptions - using daily growth rate
  t_period <- ttotal <- as.numeric(sum(diff(DAT$date)))+1
  t_length <- ttotal

  storeL <- array( 0, dim=c( nn, t_length, length(theta_init_names) ), dimnames = list(NULL,NULL,theta_init_names) )

  # Add initial condition
  storeL[,1,"I"] <- theta[["init_cases"]]
  storeL[,1,"S"] <- theta[["pop"]] - theta[["init_cases"]]

  simzeta <- matrix( NA, nrow=ttotal, ncol=nn )
  # simzeta[1,] <- runif( nn, min=0.5, max=30 ) # define IC
  simzeta[1,] <- rtruncnorm( nn, a=1e-2, b=Inf, mean=1.1, sd = 5)

  # Latent variables
  S_traj = matrix(NA,ncol=1,nrow=ttotal)
  E1_traj = matrix(NA,ncol=1,nrow=ttotal)
  E2_traj = matrix(NA,ncol=1,nrow=ttotal)
  I_traj = matrix(NA,ncol=1,nrow=ttotal)
  X_traj = matrix(NA,ncol=1,nrow=ttotal)
  beta_traj = matrix(NA,ncol=1,nrow=ttotal);
  w <- matrix(NA,nrow=nn,ncol=ttotal);
  w[,1] <- 1  # weights
  W <- matrix(NA,nrow=nn,ncol=ttotal) # normialized weights
  A <- matrix(NA,nrow=nn,ncol=ttotal) # particle parent matrix
  l_sample <- rep(NA,ttotal)
  lik_values <- rep(NA,ttotal)

  # Iterate through steps
  for(tt in 2:ttotal){
    # cat("t " = tt, "\n")
    # DEBUG  tt=2
    # Add random walk on transmission ?
    simzeta[ tt, ] <- rtruncnorm(nn, a=1e-2, b=Inf, mean = simzeta[ tt-1, ], sd=1)
    # simzeta[ tt, ] <- rtruncnorm( nn, a=1e-2, b=Inf, mean = 2, sd = 15)
    # run process model
    storeL[ , tt, ] <- process_model( tt-1, tt, dt, theta, storeL[,tt-1,], simzeta[tt,] )
    # calculate weights
    w[ ,tt ] <- AssignWeights( data_list=DAT, storeL, nn, theta, tt )
    # normalise particle weights
    sum_weights <- sum(w[1:nn,tt])
    W[1:nn,tt] <- w[1:nn,tt]/sum_weights
    # resample particles by sampling parent particles according to weights:
    A[, tt] <- resample_systematic( W[1:nn,tt] )

    # Resample particles for corresponding variables
    storeL[,tt,] <- storeL[ A[, tt] ,tt,]
    simzeta[tt,] <- simzeta[tt, A[, tt]] #- needed for random walk on beta

  } # END PARTICLE LOOP

  # Estimate likelihood:
  for(tt in 1:ttotal){
    lik_values[tt] = log(sum(w[1:nn,tt])) # log-likelihoods
  }

  likelihood0 = -ttotal*log(nn)+ sum(lik_values) # log-likelihoods

  # Sample latent variables:
  locs <-  sample(1:nn,prob = W[1:nn,tt],replace = T)
  l_sample[ttotal] <- locs[1]
  X_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"X"]
  S_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"S"]
  I_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"I"]
  beta_traj[ttotal,] <- simzeta[ttotal,l_sample[ttotal]]

  for(ii in seq(ttotal,2,-1)){
    l_sample[ii-1] <- A[l_sample[ii],ii] # have updated indexing
    X_traj[ii-1,] <- storeL[l_sample[ii-1],ii-1,"X"]
    S_traj[ii-1,] <- storeL[l_sample[ii-1],ii-1,"S"]
    I_traj[ii-1,] <- storeL[l_sample[ii-1],ii-1,"I"]
    beta_traj[ii-1,] <- simzeta[ii-1,l_sample[ii-1]]
  }

  # DEBUG  plot(Rep_traj[,1]-C_traj[,1])
  return( list( S_trace=S_traj, I_trace=I_traj, X_trace=X_traj, beta_trace=beta_traj, lik=likelihood0 ) )

}

resample_systematic <- function(weights) {
  # Systematic resampling
  # input:
  ### weights: a vector of length N with (unnormalized) importance weights
  # output:
  ### a vector of length N with indices of the replicated particles
  N <- length(weights)
  weights <- weights/sum(weights)# normalize weights
  csum <- cumsum(weights)
  u1 <- runif(1,min=0,max=1/N) # draw a single uniform number
  u <- (1:N - 1)/N + u1
  idx <- vector("integer",length=length(weights)) # allocate a vector for the results
  j <- 1
  for(i in 1:N) {
    while (u[i] > csum[j]) {
      j <- j + 1
    }
    idx[i] <- j
  }
  return(idx)
}
