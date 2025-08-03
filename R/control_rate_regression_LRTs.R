#  Copyright 2016 Annamaria Guolo (University of Padova)
#  Permission to use, copy, modify and distribute this software and
#  its documentation, for any purpose and without fee, is hereby granted,
#  provided that:
#  1) this copyright notice appears in all copies and
#  2) the source is acknowledged through a citation to the paper
#     Guolo A. (2016). Improving likelihood-based inference in control rate regression. Submitted.
#  The Authors make no representation about the suitability of this software
#  for any purpose.  It is provided "as is", without express or implied warranty


lik <- function(theta, beta1.null){
  lik.single <- function(xx, theta, beta1.null){
    p <- length(theta)
    beta0 <- theta[1]
    if(!is.null(beta1.null)) ## under H0, searching for constrained MLE
      beta1 <- beta1.null
    else
      beta1 <- theta[2]    ## searching for MLE
    mu <- theta[p-2]
    sigma2 <- theta[p-1]
    tau2 <- theta[p]
    if(any(theta[(p-1):p]<=0))
      return(NA)
    else{
      yi <- xx[1:2]
      fi <- c(beta0+beta1*mu, mu)
      Vi <- matrix(xx[3:6], ncol=2)+ matrix(c(tau2+(beta1^2)*sigma2, beta1*sigma2, beta1*sigma2, sigma2), ncol=2)
      return( mvtnorm::dmvnorm(yi, mean=fi, sigma=Vi, log=TRUE) )
    }
  }
  values <- apply(data, 1, lik.single, theta=theta, beta1.null=beta1.null)
  return( sum(values) )
}
## Mean vector for a single study
f.single <- function(xx, theta){
  p <- length(theta)
  beta0 <- theta[1]
  beta1 <- theta[2]
  mu <- theta[3]
  sigma2 <- theta[p-1]
  tau2 <- theta[p]
  fi <- matrix(c(beta0+beta1*mu, mu), ncol=1)
  return( fi )
}
## Variance/covariance matrix for a single study
V.single <- function(xx, theta){
  p <- length(theta)
  beta0 <- theta[1]
  beta1 <- theta[2]
  mu <- theta[3]
  sigma2 <- theta[p-1]
  tau2 <- theta[p]
  Vi <- matrix(xx[3:6], ncol=2)+ matrix(c(tau2+(beta1^2)*sigma2, beta1*sigma2, beta1*sigma2, sigma2), ncol=2)
  return( Vi )
}
## Gradient of the mean vector for a single study
f.grad.single <- function(xx, theta, idx){
  p <- length(theta)
  beta0 <- theta[1]
  beta1 <- theta[2]
  mu <- theta[3]
  sigma2 <- theta[p-1]
  tau2 <- theta[p]
  if(idx==1) ##beta0
    return(matrix(c(1,0), ncol=1))
  if(idx==2) ##beta1
    return(matrix(c(mu,0), ncol=1))
  if(idx==3) ##mu
    return(matrix(c(beta1,1), ncol=1))
  if(idx==4 | idx==5) ## variance components
    return( matrix(c(0,0), ncol=1) )
}

## Gradient of the variance/covariance matrix for a single study
V.grad.single <- function(xx, theta, idx){
  p <- length(theta)
  beta0 <- theta[1]
  beta1 <- theta[2]
  mu <- theta[3]
  sigma2 <- theta[p-1]
  tau2 <- theta[p]
  if(idx==1 | idx==3) ##beta0 o mu
    return(matrix(0.0, ncol=2, nrow=2))
  if(idx==2) ##beta1
    return(matrix(c(2*beta1*sigma2, sigma2, sigma2 ,0), ncol=2))
  if(idx==4) ##sigma2
    return(matrix(c(beta1^2, beta1, beta1, 1), ncol=2))
  if(idx==5) ##tau2
    return( matrix(c(1, 0, 0, 0), ncol=2) )
}

## Hessian of the mean vector for a single study
f.hess.single <- function(xx, theta, idx1, idx2){
  p <- length(theta)
  beta0 <- theta[1]
  beta1 <- theta[2]
  mu <- theta[3]
  sigma2 <- theta[p-1]
  tau2 <- theta[p]
  m <- matrix(0.0, ncol=1, nrow=2)
  if( (idx1==2 & idx2==3) | (idx1==3 & idx2==2) ) ## (beta1, mu)
    m <- matrix(c(1,0), ncol=1, nrow=2)
  return( m )
}

## Hessian of the variance/covariance matrix for a single study
V.hess.single <- function(xx, theta, idx1, idx2){
  p <- length(theta)
  beta0 <- theta[1]
  beta1 <- theta[2]
  mu <- theta[3]
  sigma2 <- theta[p-1]
  tau2 <- theta[p]
  m <- diag(0, 2)
  if( (idx1==2 & idx2==2) | (idx2==2 & idx1==2) ) ## (beta1, beta1)
    m <- matrix(c(2*sigma2, 0, 0, 0), ncol=2, nrow=2)
  if( (idx1==2 & idx2==4) | (idx2==2 & idx1==4) ) ## (beta1, sigma2)
    m <- matrix(c(2*beta1, 1, 1, 0), ncol=2, nrow=2)
  return( m )
}

## Inverse of the derivative of the variance/covariance matrix with respect to idx
V.ginv.single <- function(xx, theta, idx){
  return( -solve(V.single(xx, theta))%*%V.grad.single(xx, theta, idx=idx)%*%solve(V.single(xx, theta)) )
}

## Inverse of the Hessian of the variance/covariance matrix with respect to (idx1, idx2)
V.hess.inv.single <- function(xx, theta, idx1, idx2){
  V <- V.single(xx, theta)
  V.idx1 <- V.grad.single(xx, theta, idx1)
  V.idx2 <- V.grad.single(xx, theta, idx2)
  V.idx1.idx2 <- V.hess.single(xx, theta, idx1, idx2)
  m <- solve(V) %*% (V.idx1%*%solve(V)%*%V.idx2 - V.idx1.idx2 + V.idx2%*%solve(V)%*%V.idx1) %*% solve(V)
  return( m )
}

S.matrix <- function(theta.hat, theta.tilde){
  p <- length(theta.hat)
  S <- matrix(0.0, ncol=p, nrow=p)
  for(j in 1:2)
    for(k in 1:2)
      S[j,k] <- sum( apply(data, 1, function(x) 0.5*sum(diag(V.ginv.single(x, theta.hat, j)%*%V.single(x, theta.hat)%*%V.ginv.single(x, theta.tilde, k)%*%V.single(x, theta.hat))) + t(f.grad.single(x, theta.hat, j))%*%V.ginv.single(x, theta.tilde, k)%*%(f.single(x, theta.tilde)-f.single(x, theta.hat)) + t(f.grad.single(x, theta.hat, j))%*%solve(V.single(x, theta.tilde))%*%f.grad.single(x, theta.tilde, k)) )
  for(j in 1:2)
    S[j,3] <- sum( apply(data, 1, function(x) t(f.grad.single(x, theta.hat, j))%*%solve(V.single(x, theta.tilde))%*%f.grad.single(x, theta.tilde, 3)) )
  for(j in 1:2)
    for(k in 4:5)
      S[j,k] <- sum( apply(data, 1, function(x) 0.5*sum(diag(V.ginv.single(x, theta.hat, j)%*%V.single(x, theta.hat)%*%V.ginv.single(x, theta.tilde, k)%*%V.single(x, theta.hat))) + t(f.grad.single(x, theta.hat, j))%*%V.ginv.single(x, theta.tilde, k)%*%(f.single(x, theta.tilde)-f.single(x, theta.hat))) )
  S[3,3] <- sum( apply(data, 1, function(x) t(f.grad.single(x, theta.hat, 3))%*%solve(V.single(x, theta.tilde))%*%f.grad.single(x, theta.tilde, 3)) )
  for(k in 4:5)
    S[3,k] <- sum( apply(data, 1, function(x) t(f.grad.single(x, theta.hat, 3))%*%V.ginv.single(x, theta.tilde, k)%*%f.single(x, theta.tilde) - t(f.single(x, theta.hat))%*%solve(V.single(x, theta.tilde))%*%f.grad.single(x, theta.hat, 3)) )
  for(j in 4:5)
    for(k in 4:5)
      S[j,k] <- sum( apply(data, 1, function(x) 0.5*sum(diag(V.ginv.single(x, theta.hat, j)%*%V.single(x, theta.hat)%*%V.ginv.single(x, theta.tilde, k)%*%V.single(x, theta.hat)))) )
  for(k in 1:2)
    S[3,k] <- sum( apply(data, 1, function(x) t(f.grad.single(x, theta.hat, 3))%*%V.ginv.single(x, theta.tilde, k)%*%(f.single(x, theta.tilde)-f.single(x, theta.hat)) + t(f.grad.single(x, theta.hat, 3))%*%solve(V.single(x, theta.tilde))%*%f.grad.single(x, theta.tilde, k)) )
  for(j in 4:5)
    for(k in 1:2)
      S[j,k] <- sum( apply(data, 1, function(x) 0.5*sum(diag(V.ginv.single(x, theta.hat, j)%*%V.single(x, theta.hat)%*%V.ginv.single(x, theta.tilde, k)%*%V.single(x, theta.hat)))) )
  return(S)
}

q.vector <-  function(theta.hat, theta.tilde){
  p <- length(theta)
  q <- matrix(0.0, ncol=1, nrow=p)
  for(j in 1:p)
    q[j] <- sum(apply(data, 1, function(x) 0.5*sum(diag(V.ginv.single(x, theta.hat, j)%*%V.single(x, theta.hat))) - 0.5*sum(diag(V.ginv.single(x, theta.hat, j)%*%V.single(x, theta.hat) %*%solve(V.single(x, theta.tilde))%*%V.single(x, theta.hat))) + t(f.grad.single(x, theta.hat, j))%*%solve(V.single(x, theta.tilde))%*%(f.single(x, theta.hat)-f.single(x, theta.tilde)) ))
  return(q)
}

## expected information matrix
i.matrix <- function(theta){
  p <- length(theta)
  i.mat <- matrix(0.0, ncol=p, nrow=p)
  for(j in 1:p)
    for(k in 1:p)
      i.mat[j,k] <- sum( apply(data, 1, function(x) 0.5*sum(diag(V.ginv.single(x, theta, k)%*%V.grad.single(x, theta, j) + solve(V.single(x, theta))%*%V.hess.single(x, theta, j, k) + V.hess.inv.single(x, theta, j, k)%*%V.single(x, theta) )) + t(f.grad.single(x, theta, j))%*%solve(V.single(x, theta))%*%f.grad.single(x, theta, k)) )
  return(i.mat)
}

## observed information matrix
j.matrix <- function(theta){
  p <- length(theta)
  j.mat <- matrix(0.0, ncol=p, nrow=p)
  for(j in 1:p)
    for(k in 1:p)
      j.mat[j,k] <- sum( apply(data, 1, function(x) 0.5*sum(diag( V.ginv.single(x, theta, k)%*%V.grad.single(x, theta, j) + solve(V.single(x, theta))%*%V.hess.single(x, theta, j, k) )) + 0.5*t(x[1:2])%*%V.hess.inv.single(x, theta, j, k)%*%x[1:2] - t(f.hess.single(x, theta, j, k))%*%solve(V.single(x, theta))%*%x[1:2] - t(f.grad.single(x, theta, j))%*%V.ginv.single(x, theta, k)%*%x[1:2] - t(f.grad.single(x, theta, k))%*%V.ginv.single(x, theta, j)%*%x[1:2] - t(f.single(x, theta))%*%V.hess.inv.single(x, theta, j, k)%*%x[1:2] + t(f.hess.single(x, theta, j, k))%*%solve(V.single(x, theta))%*%f.single(x, theta) + t(f.grad.single(x, theta, j))%*%V.ginv.single(x, theta, k)%*%f.single(x, theta) + t(f.grad.single(x, theta, j))%*%solve(V.single(x, theta))%*%f.grad.single(x, theta, k) + t(f.grad.single(x, theta, k))%*%V.ginv.single(x, theta, j)%*%f.single(x, theta) + 0.5*t(f.single(x, theta))%*%V.hess.inv.single(x, theta, j, k)%*%f.single(x, theta)) )
  return(j.mat)
}

## correction term u
u.stat <- function(theta.hat, theta.tilde){
  S <- S.matrix(theta.hat, theta.tilde)
  q <- q.vector(theta.hat, theta.tilde)
  j.hat <- j.matrix(theta.hat)
  j.tilde <- j.matrix(theta.tilde)
  i.hat <- i.matrix(theta.hat)
  if(det(j.hat)<0){
    j.hat <- i.hat
    print('j.hat substituted by i.hat: check the MLEs')
  }
  if(det(j.tilde[-2,-2])<0){
    j.tilde <- i.matrix(theta.tilde)
    print('j.tilde substituted by i.tilde: check the MLEs')
  }
  return( (solve(S)%*%q)[2]*sqrt( (det(j.hat)))*solve(det(i.hat))*det(S)*( (det(j.tilde[-2,-2])))^(-1/2) )
}


# per il modello naive dovrei considerare tutti i baseline vs tutti i trattamenti presi uno per volta?




## parameter vector theta:c(beta0, beta1, mu, sigma2, tau2)
## beta1.null= values of beta1 under H0
## vector of information xx = c(eta, xi, var.eta, cov.etaxi, cov.etaxi, var.xi)


crr.test <- function(data, beta1.null, alternative = c("two.sided",
    "less", "greater"), maxit=1000){
    ans <- list()
    ans$value <- beta1.null
    alternative <- match.arg(alternative)
    ans$alternative <- alternative

    w <- 1/data$var.eta
    ## naive model, WLS
    model.naive <- lm(eta.obs~xi.obs, data=data, weights=w) ##naive model
    ## starting value for the evaluation of the MLE
    theta <- c(coef(model.naive),                 ## beta0, beta1
               mean(data[,2]),               ## mux
               var(data[,2]),        ## sigmax^2
               (mean(resid(model.naive)^2))   ## tau^2
               )
    ans$theta.wls <- theta
    ans$se.theta.wls <- sqrt(diag(vcov(model.naive)))
    ## Wald statistic
    wald <- (coef(model.naive)[2]-beta1.null)/sqrt(vcov(model.naive)[2,2])
    ans$wald <- wald
    ## MLE
    model.mle <- try(optim(theta, lik, control=list(fnscale=-1, maxit=maxit), beta1.null=NULL), silent=TRUE)
    if(class(model.mle)=='try-error' | model.mle$convergence!=0)
        print('Possible convergence problem when searching for the MLE')
    ans$mle <- model.mle$par
    theta.hat <- model.mle$par
    se <- try(sqrt(diag(solve(i.matrix(theta.hat)))), silent=TRUE)
    ans$se.mle <- se
    model.mle.constrained <- optim(theta[-2], lik, beta1.null=beta1.null, control=list(fnscale=-1, maxit=maxit))
    theta.constrained <- c( model.mle.constrained$par[1], beta1.null,  model.mle.constrained$par[-1])
    ## first-order statistic
    r <- sign(theta.hat[2]-beta1.null)*sqrt(2*(lik(theta.hat, beta1.null=NULL)-lik(theta.constrained, beta1.null=NULL)))
    u <- try(u.stat(theta.hat, theta.constrained), silent=TRUE)
    ## Skovgaard's statistic
    r.skovgaard <- r + log( (u/r) )/r
    ans$r <- r
    ans$r.skovgaard <- r.skovgaard
    if (alternative == "less") {
      ans$pvalue.wald <- pnorm(wald)
      ans$pvalue.r <- pnorm(r)
      ans$pvalue.r.skovgaard <- pnorm(r.skovgaard)
    }
    else if (alternative == "greater") {
      ans$pvalue.wald <- pnorm(wald, lower.tail = FALSE)
      ans$pvalue.r <- pnorm(r, lower.tail = FALSE)
      ans$pvalue.r.skovgaard <- pnorm(r.skovgaard, lower.tail = FALSE)
    }
    else {
      ans$pvalue.wald <- 2 * pnorm(-abs(wald))
      ans$pvalue.r <- 2 * pnorm(-abs(r))
      ans$pvalue.r.skovgaard <- 2 * pnorm(-abs(r.skovgaard))
    }
    class(ans) <- "crr.test"
    return(ans)
  }

#' @export
print.crr.test <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

    cat("\nEstimate of beta1:\n")
    tab <- matrix(NA, nrow=2, ncol=2)
    tab[,1] <- c(x$theta.wls[2], x$mle[2])
    tab[,2] <- c(x$se.theta.wls[2], x$se.mle[2])
    rownames(tab) <- c('WLS', 'MLE')
    colnames(tab) <- c('Estimate', 'Std.Err.')
    print.default(format(tab, digits = digits), print.gap = 2L, quote = FALSE)

    cat("\nHypothesis test for beta1:\n" )
    tab <- matrix(NA, nrow=3, ncol=2)
    tab[,1] <- c(x$wald, x$r, x$r.skovgaard)
    tab[,2] <- c(x$pvalue.wald, x$pvalue.r, x$pvalue.r.skovgaard)
    rownames(tab) <- c('Wald statistic', 'Signed profile log-likelihood ratio statistic', 'Skovgaard statistic')
    colnames(tab) <- c('Value','P-value')
    print.default(format(tab, digits = digits), print.gap = 2L, quote = FALSE)
    if (x$alternative == "two.sided")
      cat("\nalternative hypothesis: parameter is different from ",
          round(x$value, digits), sep = "", "\n")
    else cat("\nalternative hypothesis: parameter is ", x$alternative,
             " than ", round(x$value, digits), sep = "", "\n")
  }
