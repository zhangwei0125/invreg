
#' Inverse Regression for Completely Randomized Experiments with Multiple Outcomes
#'
#' Implements an inverse regression approach to estimating the average treatment effect (ATE)
#' for multiple outcomes in completely randomized experiments.
#'
#' The `compATE_CRE()` function computes the ATE estimator, its standard error, and confidence interval
#' for a composite outcome based on the inverse regression framework. The function supports two settings:
#' estimation without covariates and estimation with covariates included in the inverse regression.
#' In both settings, the function reports results when the vector of true marginal treatment effects
#' is either zero or nonzero.
#'
#' @param y A numeric matrix of observed multiple outcomes, where rows represent individuals and
#' columns represent distinct outcomes.
#' @param z A numeric vector of assigned treatments, where each element corresponds to an individual.
#' @param x A numeric vector or matrix of covariates to be included in the inverse regression.
#' Rows correspond to individuals and columns to covariates.
#' @param tau A logical value indicating whether the vector of true marginal treatment effects is zero.
#' @param alpha A numeric value specifying the confidence level for the confidence interval.
#'
#' @return A list of class `'res'` containing the following components:
#' \describe{
#'   \item{ATEc}{Estimated average treatment effect (ATE) for the composite outcome.}
#'   \item{SE}{Standard error of the ATE estimator for the composite outcome.}
#'   \item{CI}{Confidence interval for the ATE of the composite outcome.}
#'   \item{lambda}{Estimated weights of the asymptotic distribution of the ATE estimator
#'    for the composite outcome when \code{tau = 0}.}
#' }
#'
#' @examples
#' n <- 1000
#' L <- 9
#' beta <- matrix(rep(0.1, L), nrow = n, ncol = L, byrow = TRUE)
#' K <- 5
#' eta <- matrix(0.01, nrow = K, ncol = L)
#' p <- 0.5
#' z <- rbinom(n, size = 1, prob = p)
#' x <- matrix(rnorm(n * K), ncol = K, nrow = n)
#' epsilon <- matrix(rnorm(n * L), ncol = L, nrow = n)
#' y <- matrix(z, nrow = n, ncol = L) * beta + x %*% eta + epsilon
#'
#' compATE_CRE(y = y, z = z, x = x, tau = "Nonzero", alpha = 0.05)
#'
#' @export
compATE_CRE <- function(y, z, x, tau, alpha){

  n = length(z)
  L = ncol(y)

  ### completely randomized experiments without covariates
  if(is.null(x)){
    ### run the inverse regression to obtain the composite outcome and its ATE
    invReg = lm(z~1+y)
    invregCoef = coef(invReg)
    beta.hat = matrix(invregCoef[-1],ncol=1)
    invResid = resid(invReg)
    tau.hat = coef(lm(y~1+z))[2,]
    tauc.hat = sum(beta.hat*tau.hat)

    ### sample mean and variance
    meanz = mean(z)
    Szz = meanz*(1-meanz)
    meany = colMeans(y)
    Syy = (n-1)/n*cov(y)

    ### Inference about tauc.hat based on the true value of tau
    # when tau = 0
    if(tau == "Nonzero"){
      ### calculate the asymptotic variance of tauc.hat
      Rmat = matrix(0, nrow=L, ncol=L)
      for(i in 1:n){
        yi.center = matrix(y[i,]-meany, ncol=1)
        Term1 = (yi.center%*%t(yi.center)-Syy)%*%beta.hat
        Term2 = -Szz^(-1)*((z[i]-meanz)^2-Szz)*Syy%*%beta.hat
        Term3 = 2*invResid[i]*yi.center
        ri = Term1+Term2+Term3
        Rmat = Rmat+n^(-1)*ri%*%t(ri)
      }
      asyvar = Szz^(-2)*t(beta.hat)%*%Rmat%*%beta.hat
      tauc.se = sqrt(asyvar/n)
      ### confidence interval for the composite treatment effect
      lowerCI = tauc.hat-qnorm(1-alpha/2)*tauc.se
      upperCI = tauc.hat+qnorm(1-alpha/2)*tauc.se

      ## result to be returned
      res = list(ATEc=tauc.hat, CI=c(lowerCI, upperCI), SE=tauc.se)
    }
    # when tau != 0
    if(tau == "Zero"){
      ### calculate the asymptotic variance of tauc.hat
      InsideMat = matrix(0, nrow=L, ncol=L)
      for(i in 1:n){
        yi.center = matrix(y[i,]-meany, ncol=1)
        InsideMat = InsideMat + n^(-1)*(invResid[i]^2)*(yi.center%*%t(yi.center))
      }
      OuterMat = Syy
      gamma.hat = Szz^(-1)*(solve(OuterMat)%*%InsideMat)
      eigval = eigen(gamma.hat)$values
      quan.lower = q_mixchisq_liu(p=alpha/2, lambda=eigval)
      quan.upper = q_mixchisq_liu(p=1-alpha/2, lambda=eigval)
      lowerCI = tauc.hat-quan.upper/n
      upperCI = tauc.hat-quan.lower/n
      tauc.se = sqrt(2*sum(eigval^2)/(n^2))
      ## result to be returned
      res = list(ATEc=tauc.hat, CI=c(lowerCI, upperCI), lambda=eigval,
                 SE=tauc.se)
    }
  }

  ### completely randomized experiments with covariates
  if(!is.null(x)){
    if(is.vector(x)){
      K = 1
      x = matrix(x, ncol=1)
    }else{
      K = ncol(x)
   }

    ### run the inverse regression to obtain the weights
    invReg = lm(z~1+y+x)
    invregCoef = coef(invReg)
    beta.hat = matrix(invregCoef[2:(1+L)],ncol=1)
    invResid = resid(invReg)
    tau.hat = coef(lm(y~1+z))[2,]
    tauc.hat = sum(beta.hat*tau.hat)

    ### sample mean and variance
    meanz = mean(z)
    Szz = meanz*(1-meanz)
    meany = colMeans(y)
    Syy = (n-1)/n*cov(y)
    meanx = colMeans(x)
    Sxx = (n-1)/n*cov(x)
    invSxx = solve(Sxx)
    Syx = (n-1)/n*cov(y, x)
    Sxy = t(Syx)

    ### Inference about tauc.hat based on the true value of tau
    # when tau != 0
    if(tau == "Nonzero"){
      ### calculate the asymptotic variance of tauc.hat
      Rmat = matrix(0, nrow=L, ncol=L)
      for(i in 1:n){
        yi.center = matrix(y[i,]-meany, ncol=1)
        xi.center = matrix(x[i,]-meanx, ncol=1)
        Term1 = (yi.center%*%t(yi.center)-Syy)%*%beta.hat
        Term2 = -Syx%*%invSxx%*%(xi.center%*%t(yi.center)-Sxy)%*%beta.hat
        Term3 = -Szz^(-1)*((z[i]-meanz)^2-Szz)*(Syy-Syx%*%invSxx%*%Sxy)%*%beta.hat
        Term4 = 2*invResid[i]*yi.center
        Term5 = -2*invResid[i]*Syx%*%invSxx%*%xi.center
        ri = Term1+Term2+Term3+Term4+Term5
        Rmat = Rmat+n^(-1)*ri%*%t(ri)
      }
      asyvar = Szz^(-2)*t(beta.hat)%*%Rmat%*%beta.hat
      tauc.se = sqrt(asyvar/n)
      ### confidence interval for the composite treatment effect
      lowerCI = tauc.hat-qnorm(1-alpha/2)*tauc.se
      upperCI = tauc.hat+qnorm(1-alpha/2)*tauc.se

      ## result to be returned
      res = list(ATEc=tauc.hat, CI=c(lowerCI, upperCI), SE=tauc.se)
    }
    # when tau = 0
    if(tau == "Zero"){
      ### calculate the asymptotic variance of tauc.hat
      Mat1 = matrix(0, nrow=L, ncol=L)
      Mat2 = matrix(0, nrow=L, ncol=K)
      Mat3 = matrix(0, nrow=K, ncol=K)
      for(i in 1:n){
        yi.center = matrix(y[i,]-meany, ncol=1)
        xi.center = matrix(x[i,]-meanx, ncol=1)
        Mat1 = Mat1 + n^(-1)*(invResid[i]^2)*(yi.center%*%t(yi.center))
        Mat2 = Mat2 + n^(-1)*(invResid[i]^2)*(yi.center%*%t(xi.center))
        Mat3 = Mat3 + n^(-1)*(invResid[i]^2)*(xi.center%*%t(xi.center))
      }
      InsideMat = Mat1-Mat2%*%invSxx%*%Sxy-Syx%*%invSxx%*%t(Mat2)+Syx%*%invSxx%*%Mat3%*%invSxx%*%Sxy
      OuterMat = Syy-Syx%*%invSxx%*%Sxy
      gamma.hat = Szz^(-1)*(solve(OuterMat)%*%InsideMat)
      eigval = eigen(gamma.hat)$values
      quan.lower = q_mixchisq_liu(p=alpha/2, lambda=eigval)
      quan.upper = q_mixchisq_liu(p=1-alpha/2, lambda=eigval)
      lowerCI = tauc.hat-quan.upper/n
      upperCI = tauc.hat-quan.lower/n
      tauc.se = sqrt(2*sum(eigval^2)/(n^2))
      ## result to be returned
      res = list(ATEc=tauc.hat, CI=c(lowerCI, upperCI), lambda=eigval,
                 SE=tauc.se)
    }
  }

  ####
  return(res)
}



