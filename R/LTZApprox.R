#' Liu-Tang-Zhang method
#'
#' Quantiles of a weighted sum of chi-squared distributions
#' using Liu-Tang-Zhang approximation method.
#'
#' The `q_mixchisq_liu()` function computes the quantiles of
#' weighted sum of chi-squared distributions using Liu-Tang-Zhang approximation method.
#' See details in R fuction liu() in the package CompQuadForm.

#' @param p A numeric value whose sample quantile is wanted
#' @param lambda The weights of the weighted sum of chi-squared distributions
#' to be evaluated
#'
#'@return A value is returned.
#'
#'@export
q_mixchisq_liu <- function(p, lambda){
  # Define function to find quantile by numerical inversion
  f <- function(q) {
    res <- liu(q, lambda=lambda)
    return(res - (1-p))  # difference between CDF and target p
  }
  # Root finding (bisection)
  q <- uniroot(f, c(-max(lambda)*10000, max(lambda)*10000))$root
  return(q)
}



