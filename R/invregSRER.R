#' Inverse Regression for Stratified Randomized Experiments with Multiple Outcomes
#'
#' Implements an inverse regression approach to estimating the average treatment effect (ATE)
#' for multiple outcomes in stratified randomized experiments based on the regression strategy.
#'
#' The `compATE_SRE()` function computes the ATE estimator, its standard error, and confidence interval
#' for a composite outcome based on the inverse regression framework. The function supports two settings:
#' estimation without covariates and estimation with covariates included in the inverse regression.
#' In both settings, the function reports results when the vector of variance-weighted stratum-specific
#'  treatment effects is either zero or nonzero.
#'
#' @param y A numeric matrix of observed multiple outcomes, where rows represent individuals and
#' columns represent distinct outcomes.
#' @param z A numeric vector of assigned treatments, where each element corresponds to an individual.
#' @param c A numeric vector of observed stratum indicators #'
#' @param cf A factor of possible strata indicators
#' @param x A numeric vector or matrix of covariates to be included in the inverse regression.
#' Rows correspond to individuals and columns to covariates.
#' @param r An adjusting coefficient in the covariate-adjusted ATE estimators.
#' @param tauSR A logical value indicating whether the vector of variance-weighted stratum-specific
#'  treatment effects is equal zero.
#' @param alpha A numeric value specifying the confidence level for the confidence interval.
#'
#' @return A list of class `'res'` containing the following components:
#' \describe{
#'   \item{ATEc}{Estimated average treatment effect (ATE) for the composite outcome.}
#'   \item{SE}{Standard error of the ATE estimator for the composite outcome.}
#'   \item{CI}{Confidence interval for the ATE of the composite outcome.}
#'   \item{lambda}{Estimated weights for the asymptotic distribution of the ATE estimator for the composite outcome
#'   when \code{tauSR = 0}.}
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
#' compATE_SRER(y=y, z=z, c=c, cf=cf, x=NULL, r=-1, tauSR="Nonzero", alpha=0.05)
#'
#' @export
compATE_SRER <- function(y, z, c, cf, x, r, tauSR, alpha){
  n = length(z)
  L = ncol(y)
  S = length(cf)

  ### convert the stratification variable into dummy variables
  g = matrix(0, nrow=n, ncol=S)
  for(s in 1:S){
    g[c==cf[s], s] = 1
  }

  ### SRE without covariates based on the inverse regression strategy
  if(is.null(x)){
    ### run the inverse regression to obtain the composite outcome and its ATE
    invReg = lm(z~0+y+g)
    invregCoef = coef(invReg)
    betaSR.hat = invregCoef[1:L]
    invResid = resid(invReg)
    fwdReg = lm(y~0+z+g)
    fwdregCoef = coef(fwdReg)
    tauSR.hat = fwdregCoef[1,]
    tauSRc.hat = sum(betaSR.hat*tauSR.hat)

    ### calculate sample-based conditional mean and variance
    Szz = 0
    Syy = matrix(0, nrow=L, ncol=L)
    for(s in 1:S){
      tempLoc = which(c==cf[s])
      ns = length(tempLoc)
      zs = z[tempLoc]
      Szzs = var(zs)
      Szz = Szz + (ns/n)*Szzs
      ys = y[tempLoc,]
      Syys = cov(ys)
      Syy = Syy + (ns/n)*Syys
    }

    ### Inference about tauc.hat based on the true value of tauSR
    # when tauSR =! 0
    if(tauSR == "Nonzero"){
      Rmat = matrix(0, nrow=L, ncol=L)
      for(s in 1:S){
        tempLoc = which(c==cf[s])
        ns = length(tempLoc)
        zs = z[tempLoc]
        ys = y[tempLoc,]
        invResids = invResid[tempLoc]
        Szzs = var(zs)
        Syys = cov(ys)
        meanzs = mean(zs)
        meanys = colMeans(ys)
        for(i in 1:ns){
          yi.center = matrix(ys[i,]-meanys, ncol=1)
          zi.center = zs[i]-meanzs
          Term1 = (yi.center%*%t(yi.center)-Syys)%*%matrix(betaSR.hat, ncol=1)
          Term2 = Szz^(-1)*(zi.center^2-Szzs)*Syy%*%matrix(betaSR.hat, ncol=1)
          Term3 = 2*invResids[i]*yi.center
          rsi = Term1-Term2+Term3
          Rmat = Rmat + n^(-1)*matrix(rsi, ncol=1)%*%matrix(rsi, nrow=1)
        }
      }
      asyvar = Szz^(-2)*matrix(betaSR.hat, nrow=1)%*%Rmat%*%matrix(betaSR.hat, ncol=1)
      tauSRc.se = sqrt(asyvar/n)
      ### confidence interval for the composite treatment effect
      lowerCI = tauSRc.hat-qnorm(1-alpha/2)*tauSRc.se
      upperCI = tauSRc.hat+qnorm(1-alpha/2)*tauSRc.se
      ## result to be returned
      res = list(ATEc=tauSRc.hat, CI=c(lowerCI, upperCI), SE=tauSRc.se)
    }

    # when tauSR == 0
    if(tauSR == "Zero"){
      InsideMat = matrix(0, nrow=L, ncol=L)
      for(s in 1:S){
        tempLoc = which(c==cf[s])
        ns = length(tempLoc)
        zs = z[tempLoc]
        ys = y[tempLoc,]
        invResids = invResid[tempLoc]
        meanzs = mean(zs)
        meanys = colMeans(ys)
        for(i in 1:ns){
          yi.center = matrix(ys[i,]-meanys, ncol=1)
          InsideMat = InsideMat + n^(-1)*(invResids[i])^2*(yi.center%*%t(yi.center))
        }
      }
      OuterMat = Syy
      gammaSR.hat = Szz^(-1)*(solve(OuterMat)%*%InsideMat)
      eigval = eigen(gammaSR.hat)$values
      quan.lower = q_mixchisq_liu(p=alpha/2, lambda=eigval)
      quan.upper = q_mixchisq_liu(p=1-alpha/2, lambda=eigval)
      lowerCI = tauSRc.hat-quan.upper/n
      upperCI = tauSRc.hat-quan.lower/n
      tauSRc.se = sqrt(2*sum(eigval^2)/(n^2))
      ## result to be returned
      res = list(ATEc=tauSRc.hat, CI=c(lowerCI, upperCI), lambda=eigval,
                 SE=tauSRc.se)
    }
  }

  ### SRE with covariates based on the inverse regression strategy
  if(!is.null(x)){
    if(is.vector(x)){
      K = 1
      x = matrix(x, ncol=1)
    }else{
      K = ncol(x)
    }
    u = cbind(x, y)

    ### run the inverse regression to obtain the composite outcome and its ATE
    invReg = lm(z~0+y+x+g)
    invregCoef = coef(invReg)
    betaYSR.hat = invregCoef[1:L]
    betaXSR.hat = invregCoef[(L+1):(L+K)]
    invResid = resid(invReg)
    fwdYReg = lm(y~0+z+g)
    tauYSR.hat = coef(fwdYReg)[1,]
    #tauYSRc.hat = sum(betaYSR.hat*tauYSR.hat)
    fwdXReg = lm(x~0+z+g)
    tauXSR.hat = coef(fwdXReg)[1,]
    #tauXSRc.hat = sum(betaXSR.hat*tauXSR.hat)
    #tauSRc.hat = tauYSRc.hat-r*tauXSRc.hat
    betaUSR.hat = c(betaXSR.hat, betaYSR.hat)
    tauUSR.hat = c(tauXSR.hat, tauYSR.hat)
    D = matrix(0 , nrow=K+L, ncol=K+L)
    D[1:K, 1:K] = diag(rep(-r, K))
    D[(K+1):(K+L), (K+1):(K+L)] = diag(rep(1, L))
    tauSRc.hat = matrix(betaUSR.hat, nrow=1)%*%D%*%matrix(tauUSR.hat,ncol=1)

    ### calculate sample-based conditional mean and variance
    Szz = 0
    Suu = matrix(0, nrow=K+L, ncol=K+L)
    for(s in 1:S){
      tempLoc = which(c==cf[s])
      ns = length(tempLoc)
      zs = z[tempLoc]
      Szzs = var(zs)
      Szz = Szz + (ns/n)*Szzs
      us = u[tempLoc,]
      Suus = cov(us)
      Suu = Suu + (ns/n)*Suus
    }

    ###
    if(tauSR == "Nonzero"){
      Rmat = matrix(0, nrow=K+L, ncol=K+L)
      for(s in 1:S){
        tempLoc = which(c==cf[s])
        ns = length(tempLoc)
        zs = z[tempLoc]
        us = u[tempLoc,]
        invResids = invResid[tempLoc]
        Szzs = var(zs)
        Suus = cov(us)
        meanzs = mean(zs)
        meanus = colMeans(us)
        for(i in 1:ns){
          ui.center = matrix(us[i,]-meanus, ncol=1)
          zi.center = zs[i]-meanzs
          Term1 = D%*%(ui.center%*%t(ui.center)-Suus)%*%matrix(betaUSR.hat, ncol=1)
          Term2 = Szz^(-1)*(zi.center^2-Szzs)*D%*%Suu%*%matrix(betaUSR.hat, ncol=1)
          Term3 = invResids[i]*(Suu%*%D%*%solve(Suu)%*%ui.center)
          Term4 = invResids[i]*(D%*%ui.center)
          rsi = Term1-Term2+Term3+Term4
          Rmat = Rmat + n^(-1)*matrix(rsi, ncol=1)%*%matrix(rsi, nrow=1)
        }
      }
      asyvar = Szz^(-2)*matrix(betaUSR.hat, nrow=1)%*%Rmat%*%matrix(betaUSR.hat, ncol=1)
      tauSRc.se = sqrt(asyvar/n)
      ### confidence interval for the composite treatment effect
      lowerCI = tauSRc.hat-qnorm(1-alpha/2)*tauSRc.se
      upperCI = tauSRc.hat+qnorm(1-alpha/2)*tauSRc.se
      ## result to be returned
      res = list(ATEc=tauSRc.hat, CI=c(lowerCI, upperCI), SE=tauSRc.se)
    }
    ###
    if(tauSR == "Zero"){
      V.uresid = matrix(0, nrow=K+L, ncol=K+L)
      for(s in 1:S){
        tempLoc = which(c==cf[s])
        ns = length(tempLoc)
        zs = z[tempLoc]
        us = u[tempLoc,]
        invResids = invResid[tempLoc]
        meanzs = mean(zs)
        meanus = colMeans(us)
        for(i in 1:ns){
          ui.center = matrix(us[i,]-meanus, ncol=1)
          V.uresid = V.uresid + n^(-1)*(invResids[i])^2*(ui.center%*%t(ui.center))
        }
      }
      gammaSR.hat = (1/2)*Szz^(-1)*V.uresid%*%(solve(Suu)%*%D+D%*%solve(Suu))
      eigval = eigen(gammaSR.hat)$values
      quan.lower = q_mixchisq_liu(p=alpha/2, lambda=eigval)
      quan.upper = q_mixchisq_liu(p=1-alpha/2, lambda=eigval)
      lowerCI = tauSRc.hat-quan.upper/n
      upperCI = tauSRc.hat-quan.lower/n
      tauSRc.se = sqrt(2*sum(eigval^2)/(n^2))
      ## result to be returned
      res = list(ATEc=tauSRc.hat, CI=c(lowerCI, upperCI), lambda=eigval,
                 SE=tauSRc.se)
    }
  }
  ####
  return(res)
}




