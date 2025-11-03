#' Inverse Regression for Observational Studies with Multiple Outcomes
#'
#' Implements a weighted inverse regression approach to estimating the average treatment effect (ATE)
#' for multiple outcomes in observational studies.
#'
#' The `compATE_OS()` function computes the ATE estimator, its standard error, and confidence interval
#' for a composite outcome based on the weighted inverse regression framework. The function supports two settings:
#' estimation without covariates and estimation with covariates included in the weighted inverse regression.
#' In both settings, the function reports results when the vector of true marginal treatment effects
#' is either zero or nonzero.
#'
#' @param y A numeric matrix of observed multiple outcomes, where rows represent individuals and
#' columns represent distinct outcomes.
#' @param z A numeric vector of assigned treatments, where each element corresponds to an individual.
#' @param x A numeric vector or matrix of pretreatment covariates. Rows correspond to individuals and columns to covariates.
#' @param xinInvreg A logical value indicating whether the covariates are included in the weighted inverse regression
#' @param tau A logical value indicating whether the vector of true marginal treatment effects is zero.
#' @param alpha A numeric value specifying the confidence level for the confidence interval.
#'
#' @return A list of class `'res'` containing the following components:
#' \describe{
#'   \item{ATEc}{Estimated average treatment effect (ATE) for the composite outcome.}
#'   \item{SE}{Standard error of the ATE estimator for the composite outcome.}
#'   \item{CI}{Confidence interval for the ATE for the composite outcome.}
#'   \item{lambda}{Estimated weights for the asymptotic distribution of the ATE estimator for the composite outcome
#'   when \code{tau = 0}.}
#' }
#'
#' @examples
#' n <- 1000
#' L <- 9
#' beta <- matrix(rep(0.1, L), nrow = n, ncol = L, byrow = TRUE)
#' K <- 5
#' eta <- matrix(0.01, nrow = K, ncol = L)
#' alpha <- rep(0.1, K+1)
#' x <- matrix(rnorm(n * K), ncol = K, nrow = n)
#' xring.alpha <- cbind(rep(1,n), x) %*% matrix(alpha, ncol=1)
#' px <- exp(xring.alpha)/(1+exp(xring.alpha))
#' z <- apply(matrix(px, ncol=1), MARGIN=1, function(u) rbinom(1, size=1, prob=u))
#' epsilon <- matrix(rnorm(n*L), ncol=L, nrow=n)
#' y <- matrix(z, nrow=n, ncol=L)*beta + x %*% eta + epsilon
#'
#' compATE_OS(y=y, z=z, x=x, xinInvreg=FALSE, tau="Nonzero", alpha=0.05)
#'
#' @export
compATE_OS <- function(y, z, x, xinInvreg, tau, alpha){

  n = length(z)
  L = ncol(y)
  if(is.vector(x)){
    K = 1
    x = matrix(x, ncol=1)
  }else{
    K = ncol(x)
  }

  ### estimate the propensity score
  fitpsmodel = glm(z~x, family=binomial(link="logit"))
  pscore = fitpsmodel$fitted.values
  w = z/pscore+(1-z)/(1-pscore)
  ## observed score function
  alpha_hat = coef(fitpsmodel)
  p = length(alpha_hat)
  x_mat = model.matrix(fitpsmodel)
  eta_hat = x_mat%*%alpha_hat
  score = x_mat*matrix(z-pscore, nrow=n, ncol=p, byrow=F)
  ## observed Fisher Information matrix
  Imat_inv = n*(summary(fitpsmodel)$cov.unscaled)  # inverse of Fisher info
  Imat = solve(Imat_inv)


  ### observational studies without covariates in weighted inverse regression
  if(!xinInvreg){
    ### run the inverse weighted least square regression
    invReg = lm(z~1+y, weights=w)
    invregCoef = coef(invReg)
    betaOS.hat = invregCoef[2:(L+1)]
    invResid = resid(invReg)
    fwdReg = lm(y~1+z, weights=w)
    fwdregCoef = coef(fwdReg)
    tauOS.hat = fwdregCoef[2,]
    tauOSc.hat = sum(betaOS.hat*tauOS.hat)

    ##### calculate the asymptotic variance of the composite treatment effect under tau_SR!=0
    # calculate S_ab^w
    Swyy = Swab(a=y, b=y, w=w)
    Swy1 = Swab(a=y, b=matrix(rep(1,n),ncol=1), w=w)
    Sw11 = mean(w)
    Sw1y = t(Swy1)
    Swzz = mean(w*z*z)
    Swz1 = mean(w*z)
    Sw1z = Swz1
    phi_wyy = Swyy-Sw11^(-1)*Swy1%*%Sw1y
    ## calculate the first derivative of the weight function
    ##  w=z/e(x;alpha)+(1-z)/{1-e(x;alpha)}
    der1w = matrix(NA, nrow=n, ncol=p)
    for(i in 1:n){
      Term1 = -z[i]/(pscore[i])^2+(1-z[i])/(1-pscore[i])^2
      Term2 = (exp(eta_hat[i])/(1+exp(eta_hat[i]))^2)*x_mat[i,]
      der1w[i,] = Term1*Term2
    }
    tempMat1 = matrix(0, nrow=L, ncol=p)
    tempMat2 = matrix(0, nrow=1, ncol=p)
    tempMat3 = matrix(0, nrow=L, ncol=p)
    for(i in 1:n){
      yi = y[i,]
      zi = z[i]
      yi.center = yi-Swy1/Sw11
      tempMat1 = tempMat1+n^(-1)*tcrossprod(x=yi.center, y=yi.center)%*%betaOS.hat%*%matrix(der1w[i,], nrow=1)
      tempMat2 = tempMat2+n^(-1)*matrix(der1w[i,], nrow=1)
      tempMat3 = tempMat3+n^(-1)*invResid[i]*tcrossprod(x=yi.center, y=der1w[i,])
    }

    ### Inference about tauc.hat based on the true value of tau
    # when tau = 0
    if(tau == "Nonzero"){
      Rmat = matrix(0, nrow=L, ncol=L)
      for(i in 1:n){
        wi = w[i]
        yi = y[i,]
        zi = z[i]
        yi.center = yi-Swy1/Sw11
        Term11 = (wi*tcrossprod(x=yi.center, y=yi.center)-Swyy)%*%betaOS.hat
        Term12 = tempMat1%*%Imat_inv%*%score[i,]
        Term1 = Term11+Term12
        Term21 = wi-Sw11
        Term22 = tempMat2%*%Imat_inv%*%score[i,]
        Term2 = as.numeric(Term21+Term22)*phi_wyy%*%betaOS.hat
        Term31 = invResid[i]*wi*yi.center
        Term32 = tempMat3%*%Imat_inv%*%score[i,]
        Term3 = Term31+Term32
        ri = 2*Term1-Term2+4*Term3
        Rmat = Rmat + tcrossprod(x=ri, y=ri)
      }
      asyvar = t(betaOS.hat)%*%(Rmat/n)%*%betaOS.hat
      tauOSc.se = sqrt(asyvar/n)

      ### confidence interval for the composite treatment effect
      lowerCI = tauOSc.hat-qnorm(1-alpha/2)*tauOSc.se
      upperCI = tauOSc.hat+qnorm(1-alpha/2)*tauOSc.se
      ## result to be returned
      res = list(ATEc=tauOSc.hat, CI=c(lowerCI, upperCI), SE=tauOSc.se)
    }
    # when tau != 0
    if(tau == "Zero"){
      ### calculate the asymptotic variance of tauc.hat
      InsideMat = matrix(0, nrow=L, ncol=L)
      for(i in 1:n){
        wi = w[i]
        yi = y[i,]
        zi = z[i]
        yi.center = yi-Swy1/Sw11
        Term1 = invResid[i]*wi*yi.center
        Term2 = tempMat3%*%Imat_inv%*%score[i,]
        InsideMat = InsideMat + n^(-1)*tcrossprod(x=Term1+Term2, y=Term1+Term2)
      }
      OuTermat = phi_wyy
      gamma.hat = 2*(solve(OuTermat)%*%InsideMat)
      eigval = eigen(gamma.hat)$values
      quan.lower = q_mixchisq_liu(p=alpha/2, lambda=eigval)
      quan.upper = q_mixchisq_liu(p=1-alpha/2, lambda=eigval)
      lowerCI = tauOSc.hat-quan.upper/n
      upperCI = tauOSc.hat-quan.lower/n
      tauOSc.se = sqrt(2*sum(eigval^2)/(n^2))
      ## result to be returned
      res = list(ATEc=tauOSc.hat, CI=c(lowerCI, upperCI), lambda=eigval,
                 SE=tauOSc.se)
    }
  }

  ### observational studies without covariates in weighted inverse regression
  if(xinInvreg){
    ### run the inverse weighted least square regression
    invReg = lm(z~1+y+x, weights=w)
    invregCoef = coef(invReg)
    betaOS.hat = invregCoef[2:(L+1)]
    invResid = resid(invReg)
    fwdReg = lm(y~1+z+x, weights=w)
    fwdregCoef = coef(fwdReg)
    tauOS.hat = fwdregCoef[2,]
    tauOSc.hat = sum(betaOS.hat*tauOS.hat)
    xring = cbind(rep(1,n), x)

    # calculate weighted sample moments S_ab^w
    Swyy = Swab(a=y, b=y, w=w)
    Swy1 = Swab(a=y, b=matrix(rep(1,n),ncol=1), w=w)
    Sw11 = mean(w)
    Sw1y = t(Swy1)
    Swzz = mean(w*z*z)
    Swz1 = mean(w*z)
    Sw1z = Swz1
    Swyx = Swab(a=y, b=xring, w=w)
    Swxy = t(Swyx)
    Swxx = Swab(a=xring, b=xring, w=w)
    Swzx = Swab(a=matrix(z, ncol=1), b=xring, w=w)
    Swxz = t(Swzx)
    phi_wyyx = Swyy-Swyx%*%solve(Swxx)%*%Swxy

    ## calculate the first derivative of the weight function
    der1w = matrix(NA, nrow=n, ncol=p)
    for(i in 1:n){
      Term1 = -z[i]/(pscore[i])^2+(1-z[i])/(1-pscore[i])^2
      Term2 = (exp(eta_hat[i])/(1+exp(eta_hat[i]))^2)*x_mat[i,]
      der1w[i,] = Term1*Term2
    }
    tempMat1 = matrix(0, nrow=L, ncol=p)
    tempMat2 = matrix(0, nrow=1, ncol=p)
    tempMat3 = matrix(0, nrow=L, ncol=p)
    for(i in 1:n){
      wi = w[i]
      yi = y[i,]
      zi = z[i]
      xringi = xring[i,]
      yi.center = yi-Swyx%*%solve(Swxx)%*%matrix(xring[i,], ncol=1)
      tempMat1 = tempMat1+n^(-1)*tcrossprod(x=yi.center, y=yi.center)%*%betaOS.hat%*%matrix(der1w[i,], nrow=1)
      tempMat2 = tempMat2+n^(-1)*matrix(der1w[i,], nrow=1)
      tempMat3 = tempMat3+n^(-1)*invResid[i]*tcrossprod(x=yi.center, y=der1w[i,])
    }

    ### Inference about tauc.hat based on the true value of tau
    # when tau != 0
    if(tau == "Nonzero"){
      Rmat = matrix(0, nrow=L, ncol=L)
      for(i in 1:n){
        wi = w[i]
        yi = y[i,]
        zi = z[i]
        xringi = xring[i,]
        yi.center = yi-Swyx%*%solve(Swxx)%*%matrix(xring[i,], ncol=1)
        Term11 = (wi*tcrossprod(x=yi.center, y=yi.center)-phi_wyyx)%*%betaOS.hat
        Term12 = tempMat1%*%Imat_inv%*%score[i,]
        Term1 = Term11+Term12
        Term21 = wi-Sw11
        Term22 = tempMat2%*%Imat_inv%*%score[i,]
        Term2 = as.numeric(Term21+Term22)*phi_wyyx%*%betaOS.hat
        Term31 = invResid[i]*wi*yi.center
        Term32 = tempMat3%*%Imat_inv%*%score[i,]
        Term3 = Term31+Term32
        ri = 2*Term1-Term2+4*Term3
        Rmat = Rmat + tcrossprod(x=ri, y=ri)
      }
      asyvar = t(betaOS.hat)%*%(Rmat/n)%*%betaOS.hat
      tauOSc.se = sqrt(asyvar/n)

      ### confidence interval for the composite treatment effect
      lowerCI = tauOSc.hat-qnorm(1-alpha/2)*tauOSc.se
      upperCI = tauOSc.hat+qnorm(1-alpha/2)*tauOSc.se
      ## result to be returned
      res = list(ATEc=tauOSc.hat, CI=c(lowerCI, upperCI), SE=tauOSc.se)
    }
    # when tau = 0
    if(tau == "Zero"){
      ### calculate the asymptotic variance of tauc.hat
      InsideMat = matrix(0, nrow=L, ncol=L)
      for(i in 1:n){
        wi = w[i]
        yi = y[i,]
        zi = z[i]
        xringi = xring[i,]
        yi.center = yi-Swyx%*%solve(Swxx)%*%matrix(xring[i,], ncol=1)
        Term1 = invResid[i]*wi*yi.center
        Term2 = tempMat3%*%Imat_inv%*%score[i,]
        InsideMat = InsideMat + n^(-1)*tcrossprod(x=Term1+Term2, y=Term1+Term2)
      }
      OuTermat = phi_wyyx
      gamma.hat = 2*(solve(OuTermat)%*%InsideMat)
      eigval = eigen(gamma.hat)$values
      quan.lower = q_mixchisq_liu(p=alpha/2, lambda=eigval)
      quan.upper = q_mixchisq_liu(p=1-alpha/2, lambda=eigval)
      lowerCI = tauOSc.hat-quan.upper/n
      upperCI = tauOSc.hat-quan.lower/n
      tauOSc.se = sqrt(2*sum(eigval^2)/(n^2))
      ## result to be returned
      res = list(ATEc=tauOSc.hat, CI=c(lowerCI, upperCI), lambda=eigval,
                 SE=tauOSc.se)
    }
  }

  ####
  return(res)
}





