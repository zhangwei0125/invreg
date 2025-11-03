#' Inverse Regression for Stratified Randomized Experiments with Multiple Outcomes
#'
#' Implements an inverse regression approach to estimating the average treatment effect (ATE)
#' for multiple outcomes in stratified randomized experiments based on the stratification strategy.
#'
#' The `compATE_SRES()` function computes the ATE estimator, its standard error, and confidence interval
#' for a composite outcome based on the inverse regression framework. The function supports two settings:
#' estimation without covariates and estimation with covariates included in the inverse regression.
#' In both settings, the function reports results when the vector of true marginal treatment effects
#' is either zero or nonzero.
#'
#' @param y A numeric matrix of observed multiple outcomes, where rows represent individuals and
#' columns represent distinct outcomes.
#' @param z A numeric vector of assigned treatments, where each element corresponds to an individual.
#' @param c A numeric vector of observed stratum indicators.
#' @param cf A factor of possible strata indicators.
#' @param x A numeric vector or matrix of covariates to be included in the inverse regression.
#' Rows correspond to individuals and columns to covariates.
#' @param ctau A logical vector indicating whether the stratum-specific vectors of marginal
#'  treatment effects are zero.
#' @param alpha A numeric value specifying the confidence level for the confidence interval.
#'
#' @return A list of class `'res'` containing the following components:
#' \describe{
#'   \item{ATEc}{Estimated average treatment effect (ATE) for the composite outcome.}
#'   \item{SE}{Standard error of the ATE estimator for the composite outcome.}
#'   \item{CI}{Confidence interval for the ATE of the composite outcome.}
#'   \item{lambda}{Estimated weights for the asymptotic distribution of the ATE estimator for the composite outcome
#'   when all stratum-specific vectors of marginal treatment effects are zero.}
#' }
#'
#' @examples
#' n <- 1000
#' L <- 9
#' beta <- rep(0.1, L)
#' K <- 5
#' eta <- matrix(0.01, nrow = K, ncol = L)
#' S <- 3
#' cf <- as.factor(c(1:S))
#' ctau <- c("Nonzero", "Nonzero", "Nonzero")
#' p <- c(0.3, 0.5, 0.7)
#' n.stra <- rmultinom(n = 1, size = n, prob = c(0.3, 0.3, 0.4))
#' z <- NULL
#' x <- NULL
#' y <- NULL
#' c <- NULL
#' for(s in 1:S){
#' ns <- n.stra[s,1]
#' tempz <- rbinom(ns, size=1, prob=p[s])
#' z <- c(z, tempz)
#' tempx <- matrix(rnorm(ns*K), ncol=K, nrow=ns)
#' x <- rbind(x, tempx)
#' epsilons <- matrix(rnorm(ns*L), ncol=L, nrow=ns)
#' tempy <- matrix(tempz, nrow=ns, ncol=L)*matrix(beta, nrow=ns, ncol=L, byrow=T) + tempx%*%eta + epsilons
#' y <- rbind(y, tempy)
#' tempc <- rep(s, ns)
#' c <- c(c, tempc)
#' }
#' compATE_SRES(y=y, z=z, c=c, x=x, cf=cf, ctau=ctau, alpha=0.05)
#'
#' @export
compATE_SRES <- function(y, z, c, cf, x, ctau, alpha){
  n = length(z)
  L = ncol(y)
  S = length(cf)

  ### evalaute the asymptotic distribution by the value of tau
  if(any(ctau == "Nonzero")){
    V.stra = rep(NA, S)
    n.stra = rep(NA, S)
    tauc.stra = rep(NA, S)
    for(s in 1:S){
      tempLoc = which(c==cf[s])
      n.stra[s] = length(tempLoc)
      ys = y[tempLoc,]
      zs = z[tempLoc]
      if(is.null(x)){
        xs = NULL
      }else{
        if(is.vector(x)) x = matrix(x, ncol=1)
        xs = x[tempLoc,]
      }
      runCRE_s = compATE_CRE(y=ys, z=zs, x=xs, tau="Nonzero", alpha=alpha)
      tauc.stra[s] = runCRE_s$ATEc
      V.stra[s] = n.stra[s]*(runCRE_s$SE)^2
    }
    tauSRc.hat = sum((n.stra/n)*tauc.stra)
    asyvar = sum((n.stra/n)*V.stra*as.numeric(ctau=="Nonzero"))
    tauSRc.se = sqrt(asyvar/n)
    ### confidence interval for the composite treatment effect
    lowerCI = tauSRc.hat-qnorm(1-alpha/2)*tauSRc.se
    upperCI = tauSRc.hat+qnorm(1-alpha/2)*tauSRc.se
    ## result to be returned
    res = list(ATEc=tauSRc.hat, CI=c(lowerCI, upperCI), SE=tauSRc.se)
  }

  if(all(ctau=="Zero")){
    lambda.stra = matrix(NA, nrow=S, ncol=L)
    n.stra = rep(NA, S)
    tauc.stra = rep(NA, S)
    for(s in 1:S){
      tempLoc = which(c==cf[s])
      n.stra[s] = length(tempLoc)
      ys = y[tempLoc,]
      zs = z[tempLoc]
      if(is.null(x)){
        xs = NULL
      }else{
        if(is.vector(x)) x = matrix(x, ncol=1)
        xs = x[tempLoc,]
      }
      runCRE_s = compATE_CRE(y=ys, z=zs, x=xs, tau="Zero", alpha=alpha)
      tauc.stra[s] = runCRE_s$ATEc
      lambda.stra[s,] = runCRE_s$lambda
    }
    tauSRc.hat = sum((n.stra/n)*tauc.stra)
    quan.lower = q_mixchisq_liu(p=alpha/2, lambda=c(lambda.stra))
    quan.upper = q_mixchisq_liu(p=1-alpha/2, lambda=c(lambda.stra))
    lowerCI = tauSRc.hat-quan.upper/n
    upperCI = tauSRc.hat-quan.lower/n
    tauSRc.se = sqrt(2*sum(c(lambda.stra)^2)/(n^2))
    res = list(ATEc=tauSRc.hat, CI=c(lowerCI, upperCI), lambda=c(lambda.stra),
               SE=tauSRc.se)
  }

  ####
  return(res)
}







