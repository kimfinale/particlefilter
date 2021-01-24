#' Run particle filtering models and sample trajectories of variables
#'
#' The \code{run_model()} runs particle filtering models
#' @param theta Parameter values
#' @param y lstate variables
#' @param data data to calculate likelihood against
#' @param data_type infection, symptom onset, or confirmation
#' @export
#' @examples
#' pop <- setup_population(disease, country); calculate_cases(disease, country, pop)
# Likelihood calc for SMC --------------------------------------------
run_model <- function (theta = NULL,
                       y = NULL,
                       data = NULL,
                       data_type = c("infection", "symptom onset", "confirmation"),
                       rep = 1,
                       nparticle = 100,
                       tend = 2,
                       dt = 0.2) {
  # rep <- 10; nparticle <- 100; cut_off=0
  ## create a output container using a list
  type <- match.arg(data_type)

  output <- replicate((nrow(y)+1), matrix(NA, nrow = tend, ncol = rep), simplify = FALSE)
  names(output) <-  c(y0[["name"]], "Rt")

  for (i in 1:rep) {
    cat( "i = ", i, "\n" )
    res <- particle_filter(theta = theta, y = y, data = data, data_type = type, nparticle = nparticle, tend = tend, dt = dt)
    output$S[, i] <- res$S_trace
    output$E1[, i] <- res$E1_trace
    output$E2[, i] <- res$E2_trace
    output$I[, i] <- res$I_trace
    output$R[, i] <- res$R_trace
    output$CE1[, i] <- res$CE1_trace
    output$CE2[, i] <- res$CE2_trace
    output$CI[, i] <- res$CI_trace
    output$Rt[, i] <- res$beta_trace/theta[name == "gamma", val]
  }
  return (output)
  # saveRDS(output, file = paste0("outputs/fit_", filename, ".rds"))
}
