#' Run particle filtering models and sample trajectories of variables
#'
#' The \code{run_model()} runs particle filtering models
#' @param theta Parameter values
#' @param y lstate variables
#' @param data data to calculate likelihood against
#' @param data_type infection, symptom onset, or confirmation
#' @export
#' @examples
#' res <- run_model()
# Likelihood calc for SMC --------------------------------------------
run_model <- function (params = theta,
                       y = y0,
                       data = Rt_data,
                       data_type = c("infection", "symptom onset", "confirmation"),
                       rep = 1,
                       npart = 100,
                       tend = 200,
                       dt = 0.2,
                       systematic_resampling = FALSE) {

  type <- match.arg(data_type)
  output <- replicate((length(y)+1), matrix(NA, nrow = tend, ncol = rep), simplify = FALSE)
  names(output) <-  c(names(y), "Rt")

  for (i in seq_len(rep)) {
    # cat( "i =", i, "\n" )
    res <- particle_filter(params = params, y = y, data = data, data_type = type, npart = npart, tend = tend, dt = dt,
                           systematic_resampling = systematic_resampling)
    # invisible(lapply(names(y), function(x) output[[x]][,i] <- res$trace[[x]]))
    output[["S"]][, i] <- res$trace[["S"]]
    output[["E1"]][, i] <- res$trace[["E1"]]
    output[["E2"]][, i] <- res$trace[["E2"]]
    output[["I"]][, i] <- res$trace[["I"]]
    output[["R"]][, i] <- res$trace[["R"]]
    output[["CE1"]][, i] <- res$trace[["CE1"]]
    output[["CE2"]][, i] <- res$trace[["CE2"]]
    output[["CI"]][, i] <- res$trace[["CI"]]
    output[["Rt"]][, i] <- res$trace[["beta"]]/params[["gamma"]]
  }
  return (output)
  # saveRDS(output, file = paste0("outputs/fit_", filename, ".rds"))
}
