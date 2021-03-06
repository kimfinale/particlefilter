#' Re-sample systematically (by comparing to the cumulative distribution)
#'
#' The \code{systematic_resampling()} a
#' population for both sexes and incidence rate
#' @param weights A vector of weights
#' @export
#' @examples
#' id <- systematic_resampling

systematic_resampling <- function(weights) {
  N <- length(weights)
  weights <- weights/sum(weights)
  csum <- cumsum(weights)
  u1 <- runif(1, min = 0, max = 1/N )
  u <- (1:N - 1)/N + u1
  idx <- vector("integer", length = length(weights)) # allocate a vector for the results
  j <- 1
  for (i in 1:N) {
    while (u[i] > csum[j]) {
      j <- j + 1
    }
    idx[i] <- j
  }
  return (idx)
}
