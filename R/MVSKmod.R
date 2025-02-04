pacman::p_load(
  MatrixMixtures,
  tidyverse,
  Rcpp,
  sn,
  clusterGeneration,
  rmutil,
  MixMatrix,
  matlib,
  rootSolve,
  DistributionUtils,
  Bessel,
  pracma,
  maxLik,
  GeneralizedHyperbolic,
  truncnorm,
  SSLASSO,
  mvtnorm,
  parallel
)

# turns unlisted vector into formatted list according to how i do distributions
reconstitute <- function(unlisted, nparams, nresp, extra_names = NULL){
  lst <- unlisted
  thetalen <- 1:(nparams*nresp)
  alen <- 1:nresp
  rholen <- 1:1
  psilen <- 1:(nresp^2)
  Theta <- matrix(lst[thetalen], nparams, nresp)
  lst <- lst[-c(thetalen)]
  A <- matrix(lst[alen], nresp, 1)
  lst <- lst[-c(alen)]
  rho <- lst[rholen]
  lst <- lst[-c(rholen)]
  Psi <- matrix(lst[psilen], nresp, nresp)
  lst <- lst[-c(psilen)]
  ret <- list("Theta" = Theta, "A" = A, "rho" =  rho, "Psi" = Psi)
  # browser()
  ret2 <- as.list(lst)
  if(!is.null(extra_names)){
    names(ret2) <- extra_names
  }
  c(ret, ret2)
}

## bessel related functions

# take abs of nu since nuAsym is weird
besselK.adj <- function(x, nu, k.max = 5, exp = FALSE) {

  if(!is.numeric(abs(nu))){
    print(paste("ASD", abs(nu)))
  }
  if (abs(nu) >= 30 & x >= 300) {
    return(besselK.nuAsym(x, abs(nu), k.max, exp))
  }
  else if (x >= 100) {
    return(besselKasym(x, abs(nu), k.max, exp))
  }
  else{
    return(besselK(x, nu, exp))
  }
}

# closed form derivative of besselK wrt x
dbesselK_omega <- function(x, nu) {
  -besselK.adj(x, nu - 1, 5) - (nu / x) * besselK.adj(x, nu, 5)
}

# find log modified bessel function #2
# employ scaling to be able to find large x
lbesselK <- function(x, nu) {
  log(besselK.adj(x, nu, 5, TRUE)) - x
}


# order derivative of log-besselK
ldbesselK_lambda <- function(x, nu) {
  # helper just to get nu to be first argument

  rev <- function(nu,
                  x) {
    lbesselK(x, nu)
  }

  # order derivative of besselK
  numer <-
    tryCatch(
      numderiv(rev, x0 = nu, x = x),
      error = function(e) {
        print(c(x, nu))
      }
    )

  return(numer$df)

  # delete below if this works

  # besselK
  denom <- besselK.adj(x, nu, 5, TRUE)


  # in case derivative was not computable
  if (numer$df) {
    numer$df / denom
  }
  else{
    numderiv(rev,
             x0 = nu,
             x = x,
             exp = TRUE,)$df / denom
  }

}

# derivative of log besselK wrt x
ldbesselK_omega <- function(x, nu) {
  brat <- besselK.adj(x, nu - 1, 5, TRUE) / besselK.adj(x, nu, 5, TRUE)
  -besselRatio(x, nu, -1) - nu/x
  # -brat - nu / x

  #dbesselK_omega(x, nu)/besselK.adj(x, nu)
}

## functions related to mathing out updating/e-step

# construct n_i x n_i CS matrix from rho
csform <- function(ni, rho, inv = FALSE) {
  i <- diag(ni)
  j <- matrix(rep(1, ni ^ 2), ni, ni)
  ret <- c()
  if (inv) {
    ret <- 1 / (1 - rho) * (i - rho / (1 - rho + ni * rho) * j)
  }
  else{
    ret <- rho * j + (1 - rho) * i
  }
  ret
}

# closed form determinant of compound symmetric matrix
sigma_det <- function(ni, rho) {
  ret <- (1 - rho) ^ ni * (1 + ni * rho / (1 - rho))
  ret
}

# construct skewness matrix from given rowvector
# outer product between vector of ones and given skewness vector
a.mat <- function(ni, a) {
  tryCatch(
    rep(1, ni) %*% t(a),
    error = function(e) {
      print(ni)
    }
  )
}

# Gallaugher's Rho given in his paper
# Invert covariance matrices before passing
# A should be a matrix
rhofun <- function(A, Sigma_inv, Psi_inv) {
  tr(Sigma_inv %*% A %*% Psi_inv %*% t(A))
}

# Gallaugher's Delta given in his paper
# Invert cov matrices before passing
deltafun <- function(Y, X, Theta, Sigma_inv, Psi_inv) {
  loc <- (Y - X %*% Theta)
  tr(Sigma_inv %*% loc %*% Psi_inv %*% t(loc))
  #t(vec(loc)) %*% (Psi_inv %x% Sigma_inv)  %*% vec(loc)
}

# That other quadratic form thing Gallaugher had in his paper
# Invert cov matrices before passing
# A should be matrix
bilinform <- function(Y, X, A, Theta, Sigma_inv, Psi_inv) {
  loc <- (Y - X %*% Theta)
  tr(Sigma_inv %*% loc %*% Psi_inv %*% t(A))
  #t(vec(A)) %*% (Psi_inv %x% Sigma_inv) %*% vec(loc)
}


# generate random theta for guessing
# p -> number of response columns
# q -> number of covariates in coefficient matrix
genRandomTheta <- function(p, q, dist = "mvgh") {
  Theta <- matrix(rnorm(q * p), q, p) * runif(1, 0, 3)
  A <- matrix(rnorm(p) * runif(1, 0, 5), ncol =1)
  rho <- rtruncnorm(1, 0, 0.98, .5, .5)
  Psi <- genPositiveDefMat(p)$Sigma
  Psi <- Psi/(det(Psi)^(1/p))

  baseparams <- list(
    Theta = Theta,
    A = A,
    rho = rho,
    Psi = Psi

  )
  mix <- list()
  if(dist == "mvgh"){
    mix <- list(lambda = rgamma(1,5), omega = rgamma(1,3,3))
  }
  else if(dist == "mvvg"){
    mix <- list(gamma = rgamma(1,3,3))
  }
  else if(dist == "mvnig"){
    mix <- list(tgamma = rgamma(1,3,3))
  }
  else if(dist == "mvcn"){

  }

  c(baseparams, mix)
}


## latent update functions

# these three are obsolete

# update latent expectation for a single subject
# u -> expected value of mixing variable
mvgh.ui_update <- function(Yi, Xi, theta) {
  # indexing
  ni <- nrow(Yi)
  p <- ncol(Yi)

  # steps to produce rho and delta quadratic forms for computation
  amat <- a.mat(ni, theta$A)
  Sigma_inv <- csform(ni, theta$rho, TRUE)
  Psi_inv <- solve(theta$Psi)
  rho <- rhofun(amat, Sigma_inv, Psi_inv)
  delta <- deltafun(Yi, Xi, theta$Theta, Sigma_inv, Psi_inv)

  # return

  order <- theta$lambda - ni * p / 2

  ret <-
    1 / sqrt(delta + theta$omega) * besselRatio(rho + theta$omega, order, 1)

  ret

}

# update latent expectation for 1 subject
# zeta -> expected inverse of mixing variable
mvgh.zetai_update <- function(Yi, Xi, theta) {
  # same setup as ui update
  ni <- nrow(Yi)
  p <- ncol(Yi)
  amat <- a.mat(ni, theta$A)
  Sigma_inv <- csform(ni, theta$rho, TRUE)
  Psi_inv <- solve(theta$Psi)
  rho <- rhofun(amat, Sigma_inv, Psi_inv)
  delta <- deltafun(Yi, Xi, theta$Theta, Sigma_inv, Psi_inv)
  order <- theta$lambda - ni * p / 2

  1 / mvgh.ui_update(Yi, Xi, theta) - 2 * (theta$lambda - ni * p / 2) / ((rho + theta$omega) /
                                                                           (delta + theta$omega))
}

# update latent expectation for 1 subject
# xi -> expected log of mixing variable
mvgh.xii_update <- function(Yi, Xi, theta) {
  # same setup as ui update
  ni <- nrow(Yi)
  p <- ncol(Yi)
  amat <- a.mat(ni, theta$A)
  Sigma_inv <- csform(ni, theta$rho, TRUE)
  Psi_inv <- solve(theta$Psi)
  rho <- rhofun(amat, Sigma_inv, Psi_inv)
  delta <- deltafun(Yi, Xi, theta$Theta, Sigma_inv, Psi_inv)
  order <- theta$lambda - ni * p / 2

  # numer <- besselKAsym(rho + theta$omega, order, 10, TRUE)
  # denom <- numderiv(function(nu, x){besselKasym(x, nu, 10, TRUE)}, x0 = order, x= rho + theta$omega)-0.5 *
  log(delta + theta$omega) + ldbesselK_lambda(rho + theta$omega, order)
  # -0.5*log(delta + theta$omega) + numer/denom

}

## log-likelihood updating functions (Q-function updates)

# raw log likelihood mvgh update for single subject
q1i.mvgh.update <- function(Yi, Xi, theta){
  ni <- nrow(Yi)
  p <- ncol(Yi)
  Sigma_inv <- csform(ni, theta$rho, TRUE)
  Psi_inv <- solve(theta$Psi)
  amat <- a.mat(ni, theta$A)
  rho <- rhofun(amat, Sigma_inv, Psi_inv)
  delta <- deltafun(Yi, Xi, theta$Theta, Sigma_inv, Psi_inv)
  rhopw <- rho + theta$omega
  deltapw <- delta + theta$omega
  bln <- bilinform(Yi, Xi, amat, theta$Theta, Sigma_inv, Psi_inv)
  order <- theta$lambda - ni*p/2
  dt <- sigma_det(ni, theta$rho)

  bln -
    ni * p/2 * log(2*pi) -
    (p/2)*log(dt) -
    (ni/2)*log(det(theta$Psi)) -
    log(besselK.adj(theta$omega, theta$lambda)) +
    (order/2)*(log(deltapw) - log(rhopw)) +
    lbesselK(sqrt(rhopw*deltapw), order)


}

# raw log likelihood mvvg update for single subject
q1i.mvvg.update <- function(Yi, Xi, theta){
  ni <- nrow(Yi)
  p <- ncol(Yi)
  Sigma_inv <- csform(ni, theta$rho, TRUE)
  Psi_inv <- solve(theta$Psi)
  amat <- a.mat(ni, theta$A)
  rho <- rhofun(amat, Sigma_inv, Psi_inv)
  delta <- deltafun(Yi, Xi, theta$Theta, Sigma_inv, Psi_inv)
  rhopw <- rho + 2*theta$gamma
  deltapw <- delta
  bln <- bilinform(Yi, Xi, amat, theta$Theta, Sigma_inv, Psi_inv)
  order <- (theta$gamma - ni*p/2)
  dt <- sigma_det(ni, theta$rho)

  log(2) + theta$gamma*log(theta$gamma) + bln -
    ni * p/2 * log(2*pi) -
    (p/2)*log(dt) -
    (ni/2)*log(det(theta$Psi)) -
    lgamma(theta$gamma) +
    (order/2)*(log(deltapw) - log(rhopw)) +
    lbesselK(sqrt(rhopw*deltapw), order)


}

# raw log likelihood mvnig update for single subject
q1i.mvnig.update <- function(Yi, Xi, theta){
  ni <- nrow(Yi)
  p <- ncol(Yi)
  Sigma_inv <- csform(ni, theta$rho, TRUE)
  Psi_inv <- solve(theta$Psi)
  amat <- a.mat(ni, theta$A)
  rho <- rhofun(amat, Sigma_inv, Psi_inv)
  delta <- deltafun(Yi, Xi, theta$Theta, Sigma_inv, Psi_inv)
  rhopw <- rho + theta$tgamma^2
  deltapw <- delta + 1
  bln <- bilinform(Yi, Xi, amat, theta$Theta, Sigma_inv, Psi_inv) + theta$tgamma
  order <- -(1 + ni*p)/2
  dt <- sigma_det(ni, theta$rho)

  log(2/(2*pi)) + bln -
    ni * p/2 * log(2*pi) -
    (p/2)*log(dt) -
    (ni/2)*log(det(theta$Psi)) +
    (order/2)*(log(deltapw) - log(rhopw)) +
    lbesselK(sqrt(rhopw*deltapw), order)


}

# update Q1 over all subjects
q1.mvgh.update <- function(Y, X, theta){
  ret <- 0
  for (i in 1:length(Y)) {
    ret <-
      ret + q1i.mvgh.update(Y[[i]],
                            X[[i]],
                            theta)
  }

  ret
}

q1.mvvg.update <- function(Y,X,theta){
  ret <- 0
  for (i in 1:length(Y)) {
    ret <-
      ret + q1i.mvvg.update(Y[[i]],
                            X[[i]],
                            theta)
  }

  ret
}

q1.mvnig.update <- function(Y,X,theta){
  ret <- 0
  for (i in 1:length(Y)) {
    ret <-
      ret + q1i.mvnig.update(Y[[i]],
                             X[[i]],
                             theta)
  }

  ret
}

# complete data Q-update for 1 subject
q2i.mvgh.update <- function(Yi, Xi, ui, zetai, xii, theta) {
  ni <- nrow(Yi)
  p <- ncol(Yi)
  Sigma_inv <- csform(ni, theta$rho, TRUE)
  Psi_inv <- solve(theta$Psi)
  amat <- a.mat(ni, theta$A)
  rho <- rhofun(amat, Sigma_inv, Psi_inv)
  delta <- deltafun(Yi, Xi, theta$Theta, Sigma_inv, Psi_inv)

  t1 <- theta$lambda * xii
  t2 <- ni * p / 2 * log(1 - theta$rho)
  t3 <- p / 2 * log(1 - ni + ni / (1 - theta$rho))
  t4 <- ni / 2 * log(det(theta$Psi))
  t5 <- lbesselK(theta$omega, theta$lambda)
  t6 <- bilinform(Yi, Xi, amat, theta$Theta, Sigma_inv, Psi_inv)
  t7a <- zetai * (delta + theta$omega)
  t7b <- ui * (rho + theta$omega)
  t7 <- 0.5 * (t7a + t7b)
  # if (sum(is.nan(c(t1, t2, t3, t4, t5, t6, t7)))) {
  #   browser()
  # }
  t1 - t2 - t3 - t4 - t5 + t6 - t7
}




# update Q2 over all subjects
q2.mvgh.update <- function(Y, X, latentvars, theta) {
  ret <- 0
  for (i in 1:length(Y)) {
    ret <-
      ret + q2i.mvgh.update(Y[[i]],
                            X[[i]],
                            latentvars$u[i],
                            latentvars$zeta[i],
                            latentvars$xi[i],
                            theta)
  }

  ret
}

# proposed stopping condition:
# Compute L-infinity norm of theta
linf <- function(theta, prev){
  #max(abs(unlist(theta)))
  max(abs(unlist(theta) - unlist(prev)))
}

# update all latent variables for all subjects (mvgh)

# check over this,
lat.update <- function(Y, X, theta, dist) {
  N <- length(Y)
  ret <- list(u = c(),
              zeta = c(),
              xi = c())

  # get expected latent variables for each subject and add them to ret
  for (i in 1:N) {
    # init helpful variables
    Yi <- Y[[i]]
    Xi <- X[[i]]
    ni <- nrow(Yi)
    p <- ncol(Yi)
    amat <- a.mat(ni, theta$A)
    Sigma_inv <- csform(ni, theta$rho, TRUE)
    Psi_inv <- solve(theta$Psi)
    rho <- rhofun(amat, Sigma_inv, Psi_inv)
    delta <- deltafun(Yi, Xi, theta$Theta, Sigma_inv, Psi_inv)


    a <- switch(dist,
                "mvgh" = rho + theta$omega,
                "mvvg" = rho + 2*theta$gamma,
                "mvnig" = rho + (theta$tgamma)^2)
    b <- switch(dist,
                "mvgh" = delta + theta$omega,
                "mvvg" = delta,
                "mvnig" = delta + 1)
    order <- switch(dist,
                    "mvgh" = theta$lambda - ni * p / 2,
                    "mvvg" = theta$gamma - ni*p/2,
                    "mvnig" = -(1 + ni * p)/2)
    x <- sqrt(a * b)

    # latent vars for subject i
    ui <- sqrt(b / a) * besselRatio(x, order, 1)
    zetai <- sqrt(a / b) * besselRatio(x, order, 1) - 2 * order / b
    xii <- log(sqrt(b / a)) + ldbesselK_lambda(x, order)

    ret$u <- c(ret$u, ui)
    ret$zeta <- c(ret$zeta, zetai)
    ret$xi <- c(ret$xi, xii)
    # ret$u <- c(ret$u, mvgh.ui_update(Y[[i]], X[[i]], theta))
    # ret$zeta <- c(ret$zeta, mvgh.zetai_update(Y[[i]], X[[i]], theta))
    # ret$xi <- c(ret$xi, mvgh.xii_update(Y[[i]], X[[i]], theta))
  }
  #browser()
  ret
}

## CM-step updating functions

# Update coefficient matrix Theta
# closed form, complete-data likelihood update
theta.update <- function(Y, X, latentvars, theta) {
  N <- length(Y)
  # get the inverse sum thing

  # term to be inverted
  sum1 <- 0
  # second sum
  sum2 <- 0

  # get the summation terms
  for (i in 1:N) {
    ni <- nrow(Y[[i]])
    zetai <- latentvars$zeta[i]
    amat <- a.mat(ni, theta$A)
    Sigma_inv <- csform(ni, theta$rho, TRUE)

    sum1 <- sum1 + zetai * t(X[[i]]) %*% Sigma_inv %*% X[[i]]
    sum2 <-
      sum2 + t(X[[i]]) %*% Sigma_inv %*% (zetai * Y[[i]] - amat)
  }


  #browser()
  # update
  ret <- theta
  ret$Theta <- solve(sum1) %*% sum2
  ret
}


# update skewness vector a
# closed form, complete-data likelihood update
a.update <- function(Y, X, latentvars, theta) {
  N <- length(Y)

  # term to be inverted
  sum1 <- 0
  # second sum
  sum2 <- 0

  # get the summation terms
  for (i in 1:N) {
    ni <- nrow(Y[[i]])
    ui <- latentvars$u[i]
    Theta <- theta$Theta
    Sigma_inv <- csform(ni, theta$rho, TRUE)
    ones <- rep(1, ni)

    sum1 <- sum1 + ui * t(ones) %*% Sigma_inv %*% ones

    sum2 <-
      sum2 + t(Y[[i]] - X[[i]] %*% Theta) %*% Sigma_inv %*% ones
  }

  ret <- theta
  ret$A <- 1 / as.numeric(sum1) * sum2
  ret
}

# compute the raw log-likelihood for rho, the compound symmetry parameter
# used in iterative search update
rho.mvgh.logl <- function(rho, Y, X, params) {
  theta <- params
  N <- length(Y)
  p <- ncol(Y[[1]])

  s1 <- 0
  s2 <- 0
  s3 <- 0
  s4 <- 0

  for (i in 1:N) {
    Yi <- Y[[i]]
    Xi <- X[[i]]
    ni <- nrow(Yi)

    Sigma_inv <- csform(ni, rho, TRUE)
    Psi_inv <- solve(theta$Psi)
    amat <- a.mat(ni, theta$A)
    rho2 <- rhofun(amat, Sigma_inv, Psi_inv)
    delta <- deltafun(Yi, Xi, theta$Theta, Sigma_inv, Psi_inv)
    bln <- bilinform(Yi, Xi, amat, theta$Theta, Sigma_inv, Psi_inv)
    order <- theta$lambda - (ni * p / 2)

    s1 <- s1 + log(sigma_det(ni, rho))
    s2 <- s2 + bln
    s3 <- s3 + order / 2 * (log(delta + theta$omega) - log(rho2 + theta$omega))

    # if(is.nan(s1)){
    #   browser()
    # }

    ## trying to exp and then adjust on the log scale
    x <- sqrt((rho2 + theta$omega) * (delta + theta$omega))
    s4 <- s4 + lbesselK(x, order)
    if(sum(is.infinite(c(s1,s2,s3,s4)))){
      browser()
    }
  }

  s1 <- p / 2 * s1
  - s1 + s2 + s3 + s4
}

# update rho, the compound symmetry parameter
# uses raw likelihood
rho.mvgh.update <- function(Y, X, latentvars, theta) {

  # NM sends warnings
  # rho is constrained to prevent singularities
  theta$rho <-
    suppressWarnings(
      maxLik::maxLik(
        rho.mvgh.logl,
        start = theta$rho,
        method = "NM",
        Y = Y,
        X = X,
        params = theta,
        constraints = list(ineqA = matrix(c(-1, 1), 2, 1), ineqB = matrix(c(.988, .988), 2, 1))
      )$estimate
    )
  theta
}


# update Psi, response covariance matrix
# closed form, complete-data likelihood
psi.update <- function(Y, X, latentvars, theta) {
  N <- length(Y)

  # term to be inverted
  sum1 <- 0
  # second sum
  sum2 <- 0

  # get the summation terms
  for (i in 1:N) {
    ni <- nrow(Y[[i]])
    ui <- latentvars$u[i]
    zetai <- latentvars$zeta[i]
    Theta <- theta$Theta
    amat <- a.mat(ni, theta$A)
    Sigma_inv <- csform(ni, theta$rho, TRUE)
    loc <- Y[[i]] - X[[i]] %*% Theta
    # check this
    sum1 <- sum1 + ni

    sum2a <- zetai * t(loc) %*% Sigma_inv %*% loc
    sum2b <- 2 * t(loc) %*% Sigma_inv %*% amat
    sum2c <- ui * t(amat) %*% Sigma_inv %*% amat

    sum2 <- sum2 + sum2a - sum2b + sum2c
  }

  ret <- theta
  ret$Psi <- (1 / sum1) * sum2
  ret
}

# compute raw log-likelihood for lambda
lambda.logl <- function(mix, parms) {
  # inits
  lambda <- mix[1]
  Y <- parms[[1]]
  X <- parms[[2]]
  theta <- parms[[3]]
  N <- length(Y)
  p <- ncol(Y[[1]])

  # sum components
  s1 <- 0
  s2 <- N * lbesselK(theta$omega, lambda)
  s3 <- 0

  # get sum
  for (i in 1:N) {
    ni <- nrow(Y[[i]])
    amat <- a.mat(ni, theta$A)
    Sigma_inv <- csform(ni, theta$rho, TRUE)
    Psi_inv <- solve(theta$Psi)
    rho <- rhofun(amat, Sigma_inv, Psi_inv)
    delta <-
      deltafun(Y[[i]], X[[i]], theta$Theta, Sigma_inv, Psi_inv)
    rhopw <- rho + theta$omega
    deltapw <- delta + theta$omega
    order <- lambda - ni * p / 2
    x <- sqrt(rhopw * deltapw)

    s1 <- s1 + lbesselK(x, order)
    s3 <- s3 + (lambda / 2) * (log(deltapw) - log(rhopw))


    if (sum(is.nan(c(s1, s2, s3)))) {
      browser()
    }
  }

  s1 - s2 + s3
}

# same as lambda.logl but for omega
omega.logl <- function(mix, parms) {
  # inits
  omega <- mix[1]
  Y <- parms[[1]]
  X <- parms[[2]]
  theta <- parms[[3]]
  N <- length(Y)
  p <- ncol(Y[[1]])

  # components of sum
  s1 <- 0
  s2 <- N * lbesselK(omega, theta$lambda)
  s3 <- 0

  # get sum
  for (i in 1:N) {
    ni <- nrow(Y[[i]])
    amat <- a.mat(ni, theta$A)
    Sigma_inv <- csform(ni, theta$rho, TRUE)
    Psi_inv <- solve(theta$Psi)
    rho <- rhofun(amat, Sigma_inv, Psi_inv)
    delta <-
      deltafun(Y[[i]], X[[i]], theta$Theta, Sigma_inv, Psi_inv)
    rhopw <- rho + omega
    deltapw <- delta + omega
    order <- theta$lambda - ni * p / 2
    x <- sqrt(rhopw * deltapw)

    s1 <- s1 + lbesselK(x, order)
    s3 <- s3 + (order / 2) * (log(deltapw) - log(rhopw))


    # if (sum(is.nan(c(s1, s2, s3)))) {
    #   browser()
    # }
  }

  s1 - s2 + s3
}


lambda.score <- function(mix, parms){
  # inits
  lambda <- mix[1]
  Y <- parms[[1]]
  X <- parms[[2]]
  theta <- parms[[3]]
  N <- length(Y)
  p <- ncol(Y[[1]])

  # sum components
  s1 <- 0
  s2 <- N * ldbesselK_lambda(theta$omega, lambda)
  s3 <- 0

  # get sum
  for (i in 1:N) {
    ni <- nrow(Y[[i]])
    amat <- a.mat(ni, theta$A)
    Sigma_inv <- csform(ni, theta$rho, TRUE)
    Psi_inv <- solve(theta$Psi)
    rho <- rhofun(amat, Sigma_inv, Psi_inv)
    delta <-
      deltafun(Y[[i]], X[[i]], theta$Theta, Sigma_inv, Psi_inv)
    rhopw <- rho + theta$omega
    deltapw <- delta + theta$omega
    order <- lambda - ni * p / 2
    x <- sqrt(rhopw * deltapw)

    s1 <- s1 + ldbesselK_lambda(x, order)
    s3 <- s3 + (1 / 2) * (log(deltapw) - log(rhopw))


    if (sum(is.nan(c(s1, s2, s3)))) {
      browser()
    }
  }

  s1 - s2 + s3
}

omega.score <- function(mix, parms) {

  #return(numderiv(omega.logl, mix[1], parms = parms)$df)
  # inits
  omega <- mix[1]
  Y <- parms[[1]]
  X <- parms[[2]]
  theta <- parms[[3]]
  N <- length(Y)
  p <- ncol(Y[[1]])

  # components of sum
  s1 <- 0
  s2 <- N * ldbesselK_omega(omega, theta$lambda)
  s3 <- 0

  # get sum
  for (i in 1:N) {
    ni <- nrow(Y[[i]])
    amat <- a.mat(ni, theta$A)
    Sigma_inv <- csform(ni, theta$rho, TRUE)
    Psi_inv <- solve(theta$Psi)
    rho <- rhofun(amat, Sigma_inv, Psi_inv)
    delta <-
      deltafun(Y[[i]], X[[i]], theta$Theta, Sigma_inv, Psi_inv)
    rhopw <- rho + omega
    deltapw <- delta + omega
    order <- theta$lambda - ni * p / 2
    x <- sqrt(rhopw * deltapw)

    s1 <- s1 + ldbesselK_omega(x, order) * (rho + delta + 2*omega) / (2*x)
    s3 <- s3 + (order / 2) * (1/deltapw - 1/rhopw)


    # if (sum(is.nan(c(s1, s2, s3)))) {
    #   browser()
    # }
  }

  s1 - s2 + s3
}




# update lambda using raw likelihood
lambda.update <- function(Y, X, latentvars, theta) {
  #return(theta)

  parms <- list(Y,X,theta)
  est <-
    suppressWarnings(maxLik::maxLik(
      lambda.logl,
      grad = lambda.score,
      start = theta$lambda,
      method = "BFGS",
      parms = parms,
      constraints=list(ineqA=matrix(c(-1,1),2,1), ineqB=c(10,10))
    )$estimate)
  #browser()

  #print(theta$lambda)

  theta$lambda <- est
  theta
}


# update omega using raw likelihood
omega.update <- function(Y, X, latentvars, theta) {
  #return(theta)
  parms <- list(Y,X,theta)
  theta$omega <-

    suppressWarnings(maxLik::maxLik(
      omega.logl,
      grad = omega.score,
      start = theta$omega,
      method = "BFGS",
      parms = parms,
      constraints=list(ineqA=matrix(c(1,-1), 2,1), ineqB=c(0,10))
    )$estimate)

  theta
}

## mvvg specific functions

# test both update functions then run aecm then sims

rho.mvvg.logl <- function(rho, Y, X, params){
  theta <- params
  N <- length(Y)
  p <- ncol(Y[[1]])


  s1 <- 0
  s2 <- 0
  s3 <- 0
  s4 <- 0

  for(i in 1:N) {
    Yi <- Y[[i]]
    Xi <- X[[i]]
    ni <- nrow(Yi)

    Sigma_inv <- csform(ni, rho, TRUE)
    Psi_inv <- solve(theta$Psi)
    amat <- a.mat(ni, theta$A)
    rho2 <- rhofun(amat, Sigma_inv, Psi_inv)
    delta <- deltafun(Yi, Xi, theta$Theta, Sigma_inv, Psi_inv)
    bln <- bilinform(Yi, Xi, amat, theta$Theta, Sigma_inv, Psi_inv)
    order <-  theta$gamma - ni * p/2
    rhopg <- rho2 + 2*(theta$gamma)

    s1 <- s1 + bln
    s2 <- s2 + p/2 * log(sigma_det(ni, rho))
    s3 <- s3 + order/2 * (log(delta) - log(rhopg))
    s4 <- s4 + lbesselK(sqrt(delta * rhopg),order)

  }
  # s4 plus since order is negative
  s1 - s2 + s3 + s4
}

rho.mvvg.update <- function(Y, X, latentvars, theta) {
  # rho is constrained to prevent singularities
  theta$rho <-
    suppressWarnings(
      maxLik::maxLik(
        rho.mvvg.logl,
        start = theta$rho,
        method = "NM",
        Y = Y,
        X = X,
        params = theta,
        constraints = list(ineqA = matrix(c(-1, 1), 2, 1), ineqB = matrix(c(.988, .988), 2, 1))
      )$estimate
    )
  theta
}

gamma_logl <- function(mix, parms){
  # inits
  gamma <- mix[1]
  Y <- parms[[1]]
  X <- parms[[2]]
  theta <- parms[[3]]
  N <- length(Y)
  p <- ncol(Y[[1]])

  # components of sum
  s1 <- N* gamma * log(gamma)
  s2 <- N*lgamma(gamma)
  s3 <- 0
  s4 <- 0

  # get sum
  for (i in 1:N) {
    ni <- nrow(Y[[i]])
    amat <- a.mat(ni, theta$A)
    Sigma_inv <- csform(ni, theta$rho, TRUE)
    Psi_inv <- solve(theta$Psi)
    rho <- rhofun(amat, Sigma_inv, Psi_inv)
    delta <-
      deltafun(Y[[i]], X[[i]], theta$Theta, Sigma_inv, Psi_inv)
    rhopg <- rho + 2*gamma
    order <- gamma - ni * p / 2

    s3 <- s3 + order/2 * (log(delta) - log(rhopg))
    s4 <- s4 + lbesselK(sqrt(delta * rhopg),order)


    # if (sum(is.nan(c(s1, s2, s3)))) {
    #   browser()
    # }
  }

  s1 - s2 + s3 + s4
}

gamma_update <- function(Y, X, latentvars, theta){
  #return(theta)
  theta$gamma <-

    suppressWarnings(
      maxLik::maxLik(
        gamma_logl,
        start = theta$gamma,
        method = "NM",
        parms = list(Y, X, theta),
        constraints = list(ineqA = matrix(c(1, -1), 2, 1), ineqB = c(0, 10))
      )$estimate
    )

  theta
}

## mvnig specific functions

# rho looking good, do latent variables next

rho.mvnig.logl <- function(rho, Y, X, params) {
  #rho <- exp(logrho)
  #print(rho)
  theta <- params
  N <- length(Y)
  p <- ncol(Y[[1]])


  s1 <- 0
  s2 <- theta$tgamma * N
  s3 <- 0
  s4 <- 0
  s5 <- 0

  for(i in 1:N) {
    Yi <- Y[[i]]
    Xi <- X[[i]]
    ni <- nrow(Yi)

    Sigma_inv <- csform(ni, rho, TRUE)
    Psi_inv <- solve(theta$Psi)
    amat <- a.mat(ni, theta$A)
    rho2 <- rhofun(amat, Sigma_inv, Psi_inv)
    delta <- deltafun(Yi, Xi, theta$Theta, Sigma_inv, Psi_inv)
    bln <- bilinform(Yi, Xi, amat, theta$Theta, Sigma_inv, Psi_inv)
    order <- - (1 + ni * p)/2
    deltap1 <- delta + 1
    rhopg <- rho2 + (theta$tgamma)^2

    s1 <- s1 + bln
    s3 <- s3 + p/2 * log(sigma_det(ni, rho))
    s4 <- s4 + order/2 * (log(deltap1) - log(rhopg))
    s5 <- s5 + lbesselK(sqrt(deltap1 * rhopg),order)

  }
  # s4 plus since order is negative
  s1 + s2 - s3 + s4 + s5
}

rho.mvnig.update <- function(Y, X, latentvars, theta){
  # NM sends warnings
  # rho is constrained to prevent singularities

  #return(theta)
  theta$rho <-
    suppressWarnings(
      maxLik::maxLik(
        rho.mvnig.logl,
        start = theta$rho,
        method = "NM",
        Y = Y,
        X = X,
        params = theta,
        constraints = list(ineqA = matrix(c(-1, 1), 2, 1), ineqB = matrix(c(.988, .988), 2, 1))
      )$estimate
    )
  theta
}

tgamma_update <- function(Y, X, latentvars, theta){
  #return(theta)
  N <- length(Y)

  theta$tgamma <- N/sum(latentvars$u)
  theta
}

# Handles MVGH CM step. Returns updated parameter space.
mvgh.cmupdate <- function(cycle, Y, X, latentvars, theta) {
  ret <- switch(
    cycle,
    theta.update(Y, X, latentvars, theta),
    a.update(Y, X, latentvars, theta),
    rho.mvgh.update(Y, X, latentvars, theta),
    psi.update(Y, X, latentvars, theta),
    #mvgh.mix.update(Y, X, latentvars, theta),
    lambda.update(Y, X, latentvars, theta),
    omega.update(Y, X, latentvars, theta)
  )
  ret
}

# Handles mvnig CM step. Returns updated parameter space.
mvvg.cmupdate <- function(cycle, Y, X, latentvars, theta) {
  ret <- switch(
    cycle,
    theta.update(Y, X, latentvars, theta),
    a.update(Y, X, latentvars, theta),
    rho.mvvg.update(Y, X, latentvars, theta),
    psi.update(Y, X, latentvars, theta),
    gamma_update(Y, X, latentvars, theta)
  )
  ret
}


# Handles mvnig CM step. Returns updated parameter space.
mvnig.cmupdate <- function(cycle, Y, X, latentvars, theta) {
  ret <- switch(
    cycle,
    theta.update(Y, X, latentvars, theta),
    a.update(Y, X, latentvars, theta),
    rho.mvnig.update(Y, X, latentvars, theta),
    psi.update(Y, X, latentvars, theta),
    tgamma_update(Y, X, latentvars, theta)
  )
  ret
}


# handles general cm step
cmupdate <- function(cycle, Y, X, latentvars, theta, dist) {
  ret <- 0
  if (dist == "mvgh") {
    ret <- mvgh.cmupdate(cycle, Y, X, latentvars, theta)
  }
  else if(dist == "mvvg"){
    ret <- mvvg.cmupdate(cycle, Y, X, latentvars, theta)
  }
  else if(dist == "mvnig"){
    ret <- mvnig.cmupdate(cycle, Y, X, latentvars, theta)
  }


  ret
}




# return function, just for aggregating output
aecmreport <- function(i,theta_g, theta, qfinal) {
  list("Iteration" = i, "Starting Value" = theta_g, "Final Value" = theta, "Stopping Criteria" = qfinal)
}


aecm_mvgh <-
  function(Y,
           X,
           theta_g = NULL, dist = "mvgh", stopping = 1e-3, thresh = Inf, iter = 50) {
    # Model Parameter Dimensions:
    # Y_i is ni x p
    # X_i is n_i x q
    # Theta is q x p
    # A_i is n_i x p (assuming structure of 1_n %*% a^T)
    # Sigma_i is n_i x n_i (assuming Sigma has CS structure)
    # Psi is p x p
    # lambda, omega are scalar mixing parameters

    # vector of sites per patient
    ni <- unlist(lapply(Y, nrow))

    # num covariates
    q <- ncol(X[[1]])

    # num responses
    p <- ncol(Y[[1]])



    # initialize first guesses if not given by user
    if (is.null(theta_g)) {
      theta_g <- genRandomTheta(p, q, dist)
    }

    theta <- theta_g

    # Start EM iterations

    # list of vectors of latent vars for each subject
    latentvars <- list(u = c(),
                       zeta = c(),
                       xi = c())

    # store previous theta estimates in case browser opens
    thetahist <- list()

    # stores l-infinity of each iteration starting at iteration 0
    linfhist <- c(Inf)


    # main loop
    for (i in seq(iter)) {

      # save the previous guess in history in case
      thetahist <- append(thetahist, list(theta))

      # 5-6 cycles per iteration, corresponding to the following partitions:
      # 1. Theta
      # 2. A
      # 3. rho
      # 4. Psi
      # 5+: mixing variable
      cat("\r",(paste0("Completed ", i - 1, "/", iter, " iterations")))

      for (j in 1:length(theta_g)) {
        # thetahist <- c(thetahist, theta)

        ## E-step ##

        # Update latent vars
        #print(j)
        latentvars <- lat.update(Y, X, theta, dist)
        #print(latentvars)
        #print(latentvars)
        #print(j)



        #q_current <- q1.mvgh.update(Y, X, theta)
        #qhistrow <- c(qhistrow, q_current)

        # if (is.nan(q_current)) {
        #   browser()
        # }
        #

        ## CM-Step ##

        theta <-   tryCatch(
          cmupdate(j, Y, X, latentvars, theta, dist),
          error = function(e) {
            print(e)
            # browser()
            return(aecmreport(-1,theta_g, theta, linfhist))

          }
        )
        # note that this does not break program because Iteration not assigned until return statement
        if(!is.null(theta$Iteration)){
          return(theta)
        }

        #print(theta[[j]])



      }

      #qdiff <- qhistrow - tail(qhist, 1)
      #print(matrix(qdiff, 1, 6))

      # compute l-infinity norm and check stopping condition

      tryCatch(
        linf_i <- max(abs(unlist(theta) - unlist(tail(
          thetahist, 1
        )))),
        error = function(e) {
          print(e)
          browser()
        }
      )

      #linf_i <-  linf(theta)
      linfhist <- c(linfhist, linf_i)
      #diff <- abs(c(-1, 1) %*% tail(linfhist, 2))
      # print(linf_i)
      if (linf_i < stopping || ( i >10 && linf_i > thresh)) {
        #browser()
        return(aecmreport(i, theta_g, theta, linfhist))
      }


      # qhist <- rbind(qhist, qhistrow)
      # maxdiff <- max(qdiff)
      # if (maxdiff < stopping & maxdiff > 0) {
      #   return(aecmreport(i, theta_g,theta, qhist))
      # }

    }
    cat("\r",(paste0("Completed ", iter, "/", iter, " iterations")))
    cat("\n")
    # worst case: no convergence
    aecmreport(Inf, theta_g, theta, linfhist)
  }


#' AECM Estimation for Matrix-Variate Skew Models
#'
#' This function fits linear models for matrix data with varying row sizes,
#'  according to the three matrix-variate skew distributions described in
#'   (Gallaugher & McNicholas, 2019). Exchangeable observation row correlation
#'   and skewness structures are imposed to accommodate the varying row counts
#'   across matrices.
#'
#'
#' @param Y List of \eqn{n_i \times p} response matrices. Matrices must have same number of columns.
#' @param X List of \eqn{n_i \times q} design matrices. Matrices must have same number of columns.
#' @param theta_g List of parameters to pass as initial values in the AECM algorithm.
#'  If NULL, will be randomly generated. See Details for an in-depth explanation.
#' @param dist Modeling distribution. Select "mvgh" for matrix-variate generalized
#'  hyperbolic, "mvvg" for variance-gamma, and "mvnig" for normal-inverse Gaussian.
#' @param stopping Stopping threshold for the L-infinity norm of differences in consecutive parameter space, evaluated at iteration \eqn{t+1} as
#' \eqn{|\hat\theta^{t+1} - \hat\theta^{t}|_\infty}. Default is 0.001
#' @param iter Maximum number of iterations, default is 50.
#' @details
#' Fits the matrix-variate skew regression model
#'
#'  \deqn{Y_i = X_i \Theta + E_i,}
#'
#' where each response \eqn{Y_i} is a \eqn{n_i \times p} matrix that indexes \eqn{n_i} observations and \eqn{p} response variables. \eqn{X_i} corresponds to a \eqn{n_i \times q} design matrix, and \eqn{\Theta} corresponds to a \eqn{q \times p} coefficient matrix. \eqn{E_i} corresponds to a \eqn{n_i \times p} error matrix, following a matrix-variate skew distribution.
#'
#' Supported distributions include the matrix-variate generalized hyperbolic (MVGH), matrix-variate variance-gamma (MVVG), and matrix-variate normal-inverse Gaussian (MVNIG). Certain reparametrizations are implemented to avoid identifiability issues (Gallaugher & McNicholas, 2019). Additionally, the row correlation matrix \eqn{\Sigma_i} is reduced to a correlation matrix with exchangeable structure based on correlation parameter \eqn{\rho}, and skewness parameter \eqn{A_i} is represented as a function of skewness vector \eqn{a}: \eqn{A_i = \underline{1}_{n_i} a^T}.
#'
#' The structure of `theta_g` and parameter estimates returned by the function must be in the form of a list with the following named elements:
#'
#' \itemize{
#'  \item{`Theta`: } {\eqn{q \times p} coefficient matrix}
#'  \item{`a`: } {\eqn{p \times 1} skewness vector}
#'  \item{`rho`: } {Compound symmetry parameter for row correlation matrix}
#'  \item{`Psi`: } {\eqn{p \times p} column covariance matrix}
#'  \item{`gamma`: } {Univariate mixing parameter, only included when running MVVG model}
#'  \item{`tgamma`: } {Univariate mixing parameter, only included when running MVVG model}
#' }
#'
#' @return MVSKmod returns a list with the following elements:
#'
#' \itemize{
#'  \item{`Iteration`: } {Number of iterations taken to convergence. `Inf` if convergence not reached.}
#'  \item{`Starting Value`: } {List of initial parameter values.}
#'  \item{`Final Value`: } {List of final parameter estimates.}
#'  \item{`Stopping Criteria`: } {Vector of \eqn{|\hat\theta^{t+1} - \hat\theta^{t}|_\infty} at each iteration.}
#' }

#' @export

#' @author Samuel Soon
#' @author Dipankar Bandyopadhyay
#' @author Qingyang Liu
#' @examples
#' set.seed(1234)
#' # num response variables
#' p <- ncol(gaad_res[[1]])
#' # num covariates
#' q <- ncol(gaad_cov[[1]])
#' # generate initial value to input, then run AECM with MVVG distribution
#' initial_mvvg_theta <- list(Theta = matrix(rnorm(p*q), nrow = q, ncol = p),
#'                       A = rep(1,p),
#'                      rho = 0.3,
#'                      Psi = diag(p),
#'                      gamma = 4)
#' MVSKmod(gaad_res[1:50], gaad_cov[1:50], initial_mvvg_theta, "mvvg")
#'
#' p <- ncol(gaad_res[[1]])
#' q <- ncol(gaad_cov[[1]])
#' initial_mvnig_theta <- list(Theta = matrix(rnorm(p*q), nrow = q, ncol = p),
#'                       A = rep(1,p),
#'                      rho = 0.3,
#'                      Psi = diag(p),
#'                      tgamma = 3)
#' MVSKmod(gaad_res[1:50], gaad_cov[1:50], initial_mvnig_theta, "mvnig")
MVSKmod <- function(Y,
                    X,
                    theta_g = NULL,
                    dist,
                    stopping = 1e-3,
                    max_iter = 50) {
  # add mvsk funciton here with error checks, likelihood/aic summary, etc


  # You stopped here
  stop("Add error checks")
  aecm_mvgh(Y,X,theta_g, dist, stopping = stopping, thresh = Inf, iter = max_iter)

}
