library(gglasso)
library(expm)

CenterData <- function(X, W) {
  # input X: n by p matrix of gene expression levels
  # input W: n by q matrix of confounding variables (No intercept!)
  n <- nrow(W)
  W.new <- cbind(rep(1, n), W)
  X.center <- X - W.new %*% solve(t(W.new) %*% W.new) %*% t(W.new) %*% X
  return(X.center)
}

InteractionMatrix <- function(X, W) {
  # input X: n by p matrix of gene expression levels
  # input W: n by q matrix of confounding variables (No intercept!)
  p <- ncol(X)
  q <- ncol(W)
  n <- nrow(X)
  W.new <- cbind(rep(1, n), W)
  
  XW <- sweep(W.new, 1, X[,1], `*`)
  for(j in 2:p) {
    XW <- cbind(XW, sweep(W.new, 1, X[,j], `*`))
  }
  
  return(XW)
}

EstAlpha <- function(X, W, j, nfolds = 10) {
  # input X: n by p matrix of gene expression levels
  # input W: n by q matrix of confounding variables (No intercept!)
  # input j: index of which gene (column of X) is response vector
  q <- ncol(W)
  p <- ncol(X)
  n <- nrow(X)
  
  X.cent <- CenterData(X, W)
  XW <- InteractionMatrix(X.cent, W)
  J <- ((q + 1) * (j - 1) + 1):((q + 1) * j) # columns of XW pertaining to gene j
  
  XW.scale <- apply(XW, 2, sd) # to scale columns of design matrix

  ### cross validation
  cv.fit <- cv.gglasso(x = XW[,-J] %*% diag(1/XW.scale[-J]), y = X.cent[,j],
                       nfolds = nfolds, group = rep(1:(p-1), each = q + 1),
                       intercept = FALSE)
  opt.lam <- cv.fit$lambda.min
  opt.est <- cv.fit$gglasso.fit$beta[,which(cv.fit$lambda == opt.lam)] /
             XW.scale[-J] # return estimates to original scale
  
  sigma2.hat <- sum((X.cent[,j] - (XW[,-J] %*% opt.est))^2)/(n - sum(opt.est != 0))
  out <- list(alpha.hat = opt.est,
              sigma2.hat = sigma2.hat)
  return(out)
}

EstGamma <- function(X, W, j, k, nfolds = 5) {
  # input X: n by p matrix of gene expression levels
  # input W: n by q matrix of confounding variables (No intercept!)
  # input j: index of which gene (column of X) is response vector
  # input k: index of gene's coefficients we're de-biasing
  q <- ncol(W)
  p <- ncol(X)
  n <- nrow(X)
  
  X.cent <- CenterData(X, W)
  XW <- InteractionMatrix(X.cent, W)
  
  XW.scale <- apply(XW, 2, sd) # to scale columns of design matrix
  
  J <- ((q + 1) * (j - 1) + 1):((q + 1) * j) # columns of XW pertaining to gene j
  K <- ((q + 1) * (k - 1) + 1):((q + 1) * k) # columns of XW pertaining to gene k
  
  Gamma.hat <- matrix(0, nrow = (p - 2) * (q + 1), ncol = q + 1)
  for(t in 1:(q + 1)) {
    ### cross validation version 
    cv.fit <- cv.gglasso(x = XW[,-c(J,K)] %*% diag(1/XW.scale[-c(J,K)]), y = XW[,K[t]],
                         nfolds = nfolds, group = rep(1:(p-2), each = q + 1),
                         intercept = FALSE)
    opt.lam <- cv.fit$lambda.min
    opt.est <- cv.fit$gglasso.fit$beta[, which(cv.fit$lambda == opt.lam)] /
               XW.scale[-c(J, K)]
    
    Gamma.hat[,t] <- opt.est
  }
  
  return(Gamma.hat)
}

DebiasAlpha <- function(X, W, alpha.hat, Gamma.hat, sigma2.hat, j, k) {
  # input X: n by p matrix of gene expression levels
  # input W: n by q matrix of confounding variables (No intercept!)
  # input alpha.hat: (p - 1) * (q + 1) dimensional vector with group lasso estimate of alpha
  # input Gamma.hat: p - 2 by q + 1 dimensional matrix output from EstGamma
  # input j: index of which gene (column of X) is response vector
  # input k: index of gene's coefficients we're de-biasing
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(W)
  
  X.cent <- CenterData(X, W)
  XW <- InteractionMatrix(X.cent, W)
  
  J <- ((q + 1) * (j - 1) + 1):((q + 1) * j) # columns of XW pertaining to gene j
  K <- ((q + 1) * (k - 1) + 1):((q + 1) * k) # columns of XW pertaining to gene k
  
  
  ### Calcuate de-biased estimate of alpha
  bias <- solve(t(XW[,K] - (XW[,-c(J,K)] %*% Gamma.hat)) %*% XW[,K]) %*% 
          t(XW[,K] - XW[,-c(J,K)] %*% Gamma.hat) %*%
          (X.cent[,j] - XW[,-J] %*% alpha.hat)
  
  alpha.hat.jk <- numeric(q + 1)
  if(k < j){
    alpha.hat.jk <- alpha.hat[K]
  } else{
    alpha.hat.jk <- alpha.hat[(K - q - 1)]
  }
  
  alpha.check <- alpha.hat.jk + bias
  
  ### Compute asymptotic variance of debiased estimate
  Omega <- solve(sqrtm(t(XW[,K] - (XW[,-c(J,K)] %*% Gamma.hat)) %*%
                        (XW[,K] - (XW[,-c(J,K)] %*% Gamma.hat)))) %*%
           t(XW[,K] - (XW[,-c(J,K)] %*% Gamma.hat)) %*% XW[,K]
  Omega.inv <- solve(Omega)
  
  var.est <- sigma2.hat * Omega.inv %*% t(Omega.inv)
  
  ### return output
  out <- list(alpha.check = alpha.check,
              var.est = var.est)
  
  return(out)
}

NeighborhoodInference <- function(X, W, j, nfolds = 5) {
  # input X: n by p matrix of gene expression levels
  # input W: n by q matrix of confounding variables (No intercept!)
  # input j: index of which gene (column of X) is response vector
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(W)
  
  
  ### Step 1 - Group Lasso Estimate
  grpLassoEst <- EstAlpha(X, W, j, nfolds = nfolds)
  alpha.hat <- grpLassoEst$alpha.hat
  sigma2.hat <- grpLassoEst$sigma2.hat
  
  ### Step 2 - Debias Group Lasso Estimate
  alpha.check.j <- matrix(0, nrow = p - 1, ncol = q + 1)
  var.est.j <- array(0, dim = c(p - 1, q + 1, q + 1))
  
  ks <- (1:p)[-j]
  for(t in 1:(p - 1)) {
    k <- ks[t]
    Gamma.hat <- EstGamma(X, W, j, k)
    debias.grpLassoEst <- DebiasAlpha(X, W, alpha.hat, Gamma.hat, sigma2.hat, j, k)
    alpha.check.j[t,] <- debias.grpLassoEst$alpha.check
    var.est.j[t,,] <- debias.grpLassoEst$var.est
  }
  
  ### Return output
  alpha.init <- matrix(alpha.hat, nrow = p - 1, ncol = q + 1, byrow = T)
  alpha.debias <- alpha.check.j
  var.est <- var.est.j
  
  # rownames(alpha.init) <- (1:p)[-j]
  # rownames(alpha.debias) <- (1:p)[-j]
  # dimnames(var.est) <- list((1:p)[-j], NULL, NULL)

  out = list(alpha.init = alpha.init,
             alpha.debias = alpha.debias,
             var.est = var.est,
             sigma2.hat = sigma2.hat)
  
  return(out)
}

FullNetworkInference <- function(X, W, nfolds = 5) {
  # input X: n by p matrix of gene expression levels
  # input W: n by q matrix of confounding variables (No intercept!)
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(W)
  
  alpha.init <- array(0, dim = c(p, p, q + 1))
  alpha.debias <- array(0, dim = c(p, p, q + 1))
  var.est <- array(0, dim = c(p, p, q + 1, q + 1))
  p.vals <- matrix(NA, p, p)
  test.stats <- matrix(NA, p, p)
  
  for(j in 1:p) {
    #obtain de-biased estimates of alpha
    print(j)
    neighborhood.j <- NeighborhoodInference(X, W, j, nfolds)
    alpha.init[j, (1:p)[-j],] <- neighborhood.j$alpha.init
    alpha.debias[j,(1:p)[-j],] <- neighborhood.j$alpha.debias
    var.est[j,(1:p)[-j],,] <- neighborhood.j$var.est
    # test the null that edge j and k are disconnected for all w
    for(k in (1:p)[-j]) {
      test.stat <- t(alpha.debias[j,k,]) %*% solve(var.est[j,k,,]) %*% alpha.debias[j,k,]
      test.stats[j,k] <- test.stat
      p.vals[j,k] <- pchisq(test.stat, df = q + 1, lower.tail = FALSE)
    }
  }
  
  out = list(alpha.init = alpha.init,
             alpha.debias = alpha.debias,
             var.est = var.est,
             test.stats = test.stats,
             p.vals = p.vals)
  
  return(out)
}

ConfAdjDiffNetwork <- function(NetInf1, NetInf2) {
  # input: NetInf1 - a network inference object from FullNetworkInference for population 1
  # input: NetInf2 - a network inference object from FullNetworkInference for population 2
  p <- ncol(NetInf1$p.vals)
  q <- dim(NetInf1$alpha.init)[3] - 1
  
  p.vals <- matrix(NA, p, p)
  for(j in 1:p) {
    for(k in (1:p)[-j]) {
      alpha.dif <- NetInf1$alpha.debias[j,k,] - NetInf2$alpha.debias[j,k,]
      var.sum <- NetInf1$var.est[j,k,,] + NetInf2$var.est[j,k,,]
      test.stat <- t(alpha.dif) %*% solve(var.sum) %*% alpha.dif
      p.vals[j,k] <- pchisq(test.stat, df = q + 1, lower.tail = FALSE)
    }
  }
  
  return(p.vals)
}

############################################################
### TOY EXAMPLE
############################################################

library(mvtnorm)

set.seed(206)
n <- 40
p <- 5
q <- 1

Sigma.1 <- matrix(0, p, p)
diag(Sigma.1) <- 1
Omega.2 <- matrix(0, p, p)
diag(Omega.2) <- 1
Sigma.2 <- solve(Omega.2)

X.1 <- rmvnorm(n, rep(0,p), Sigma.1)
W.1 <- rmvnorm(n, rep(0,q), diag(q))

X.2 <- rmvnorm(n, rep(0, p), Sigma.2)
W.2 <- rmvnorm(n, rep(0, q), diag(q))

out.1 <- FullNetworkInference(X.1, W.1)
out.2 <- FullNetworkInference(X.2, W.2)
out <- ConfAdjDiffNetwork(out.1, out.2)

