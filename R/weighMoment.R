#' Weighted sample moments
#'
#' Weighted moments of two matrix.
#'
#' The `Swab()` function computes the weighted moments of two matrix
#' with the same number of rows: \eqn{n^{-1}\sum_{i=1}^{n}w_ia_ib_i^{'}}

#' @param a A matrix of \eqn{n\times L_1}
#' @param b A matrix of \eqn{n\times L_2}
#' @param w A numeric vector of length \eqn{n}

#'@return A value is returned.
#'
#'@export
Swab <- function(a, b, w){
  n = length(w)
  L1 = ncol(a)
  L2 = ncol(b)
  resmat = matrix(0, nrow=L1, ncol=L2)
  for(i in 1:n){
    resmat = resmat + w[i]*tcrossprod(x=a[i,], y=b[i,])
  }
  return(resmat/n)
}


