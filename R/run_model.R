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
                       npart = 100,
                       tend = 2,
                       dt = 0.2) {
  # rep <- 10; nparticle <- 100; cut_off=0
  ## create a output container using a list
  type <- match.arg(data_type)

  output <- replicate((nrow(y)+1), matrix(NA, nrow = tend, ncol = rep), simplify = FALSE)
  names(output) <-  c(y0[["name"]], "Rt")

  for (i in 1:rep) {
    cat( "i = ", i, "\n" )
    res <- particle_filter(theta = theta, y = y, data = data, data_type = type, npart = npart, tend = tend, dt = dt)
    lapply(y0[["name"]], function(z) output[[z]][,i] <- res$trace[[z]])
    # output$S[, i] <- res$trace$S
    # output$E1[, i] <- res$trace$E1
    # output$E2[, i] <- res$trace$E2
    # output$I[, i] <- res$trace$I
    # output$R[, i] <- res$trace$R
    # output$CE1[, i] <- res$trace$CE1
    # output$CE2[, i] <- res$trace$CE2
    # output$CI[, i] <- res$trace$CI
    #
    output$Rt[, i] <- res$trace$beta/theta[name == "gamma", val]
  }
  return (output)
  # saveRDS(output, file = paste0("outputs/fit_", filename, ".rds"))
}
