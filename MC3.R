library(abind)
library(rTensor)
library(GIGrvg)
library(gtools)
library(marima)
library(doSNOW)
library(parallel)
library(foreach)
library(Rcpp)
sourceCpp('MC3.cpp')

## jocabian matrix
jocobian_omega <- function(omega) {
  H <- dim(omega)[1]
  P <- dim(omega)[1]
  jocobian_matrix <- array(0, c(H, P-1, P-1))
  for (hh in 1:H) {
    tmp <- 1
    for (pp in 1:(P-1)) {
      jocobian_matrix[hh,pp,pp] <- 1/tmp
      if (pp>1) {
        jocobian_matrix[hh,pp,1:(pp-1)] <- omega[hh,pp]/tmp^2
      }
      tmp <- tmp - omega[hh,pp]
    }
  }
  return(sapply(1:H, function(hh) determinant(jocobian_matrix[hh,,], logarithm = TRUE)$modulus))
}

jocobian_Ising <- function(theta, kappa) {
  H <- length(theta)
  jocobian_matrix <- array(0, c(H, 2, 2))
  jocobian_matrix[,1,1] <- -(exp(kappa)-1)*exp(theta)*(exp(2*theta+kappa)-1)/
    ((exp(theta)+1)^2*(exp(theta+kappa)+1)^2)
  jocobian_matrix[,1,2] <- exp(theta+kappa)/(exp(theta+kappa)+1)^2
  jocobian_matrix[,2,1] <- exp(theta)*(exp(2*theta+kappa)+2*exp(theta+kappa)+1)/(exp(2*theta+kappa)+2*exp(theta)+1)^2
  jocobian_matrix[,2,2] <- (exp(theta)+1)*exp(2*theta+kappa)/(exp(2*theta+kappa)+2*exp(theta)+1)^2
  return(sapply(1:H, function(hh) determinant(jocobian_matrix[hh,,], logarithm = TRUE)$modulus))
}

## Log prior
log_prior <- function(a_lambda, b_lambda, a, a_tau, b_tau, beta1, beta2, a_w, b_w, a_sigma, b_sigma,
                      a_p1, b_p1, a_p2, b_p2, lambda1, lambda2, phi, tau, W1, W2, W3, z, nu, omega, sigma, theta,
                      kappa, p1, p2, alpha1, alpha2, alpha3, gamma) {
  H <- dim(W3)[1]
  P <- dim(W3)[2]
  prob <- t(apply(omega, 1, cumsum))
  indicator <- z<=rep(1:P, each = H)
  return(sum(dgamma(lambda1, a_lambda, b_lambda, log = TRUE)) + 
           sum(dgamma(lambda2, a_lambda, b_lambda, log = TRUE)) +
           sum(LaplacesDemon::ddirichlet(phi, rep(a, H), log = TRUE)) +
           sum(dgamma(tau, a_tau, b_tau, log = TRUE)) +
           sum(dexp(W1, 0.5*lambda1^2, log = TRUE)) +
           sum(dexp(W2, 0.5*lambda2^2, log = TRUE)) +
           sum(dbeta(nu[,-P], beta1, beta2, log = TRUE)) + sum(jocobian_omega(omega)) +
           sum(indicator[,-P]*log(prob[,-P]) + (1-indicator[,-P])*log(1-prob[,-P])) + 
           sum((1-indicator)*invgamma::dinvgamma(W3, a_w, b_w, log = TRUE)) +
           sum(dgamma(sigma, a_sigma, b_sigma, log = TRUE)) +
           sum(dbeta(p1, a_p1, b_p1, log = TRUE)) +
           sum(dbeta(p2, a_p2, b_p2, log = TRUE)) + sum(jocobian_Ising(theta, kappa)) +
           sum(dnorm(alpha1, 0, sqrt(tau*phi*W1), log = TRUE)) +
           sum(dnorm(alpha2, 0, sqrt(tau*phi*W2), log = TRUE)) +
           sum(dnorm(alpha3, 0, sqrt(tau*phi*W3), log = TRUE)) +
           sum(theta*gamma) + sum(kappa*gamma[,-1]*gamma[,-dim(gamma)[2]]))
}

## log likelihood
log_likelihood <- function(obs, alpha1, alpha2, alpha3, gamma, sigma, rho) {
  Time <- dim(obs)[1]
  N <- dim(obs)[2]
  H <- dim(alpha1)[1]
  P <- dim(alpha3)[2]
  A_matrix <- array(0, c(N, N*P, Time-P))
  est <- matrix(NA, Time-P, N)
  err <- matrix(NA, Time-P, N)
  log_likelihood <- 0
  for (tt in 1:(Time-P)) {
    for (hh in 1:H) {
      A_matrix[,,tt] <- A_matrix[,,tt] +  gamma[hh,tt] * kronecker(t(alpha3[hh,]), alpha1[hh,] %*% t(alpha2[hh,]))
    }
    est[tt,] <- A_matrix[,,tt] %*% c(t(obs[(tt+P-1):tt,]))
    err[tt,] <- obs[tt+P,] - est[tt,]
    log_likelihood <- log_likelihood -1/2 * t(err[tt,]) %*% diag(rho*1/sigma) %*% err[tt,]
  }
  return(log_likelihood)
}

## log posterior
log_posterior <- function(obs, alpha1, alpha2, alpha3, gamma, sigma, rho, a_lambda, b_lambda, a, a_tau, b_tau,
                          beta1, beta2, a_w, b_w, a_sigma, b_sigma, a_p1, b_p1, a_p2, b_p2, lambda1, lambda2, phi,
                          tau, W1, W2, W3, z, nu, omega, theta, kappa, p1, p2) {
  return(log_prior(a_lambda, b_lambda, a, a_tau, b_tau, beta1, beta2, a_w, b_w, a_sigma, b_sigma,
                   a_p1, b_p1, a_p2, b_p2, lambda1, lambda2, phi, tau, W1, W2, W3, z, nu, omega, sigma, theta,
                   kappa, p1, p2, alpha1, alpha2, alpha3, gamma) +
           log_likelihood(obs, alpha1, alpha2, alpha3, gamma, sigma, rho))
}

### FUNCTION TO UPDATE Z IN STEP 2.1
draw_z <- function(x, y, omega, alpha3, tau, phi, W_infinite, a_w, b_w) {
  P <- dim(alpha3)[2]
  prob <- rep(NA, P)
  if (x == P) {
    prob <- omega[y,]/sum(omega[y,])
  }
  else {
    prob[1:x] <- omega[y,1:x]*dnorm(alpha3[y,x], 0, sqrt(tau*phi[y]*W_infinite))
    prob[(x+1):P] <- omega[y,(x+1):P]*dt(alpha3[y,x]/sqrt(tau*phi[y]*b_w/a_w), 2*a_w)/sqrt(tau*phi[y]*b_w/a_w)
    prob <- prob/sum(prob)
  }
  return(sample(1:P,1,prob = prob))
}

### FUNCTION RELATED TO NDARMA
fun_theta <- function(p1, p2) {
  return(log(p2)+log(1-p1)-log(1-p2*(1-p1)))
}

fun_kappa <- function(p1, p2) {
  return(log(p1+p2*(1-p2)*(1-p1)^2) - log(p2) - log(1-p2) - 2*log(1-p1))
}

theta_star <- function(theta, kappa) {
  return(log(exp(theta)*(exp(theta)+1)/(exp(theta+kappa)+1)))
}

BTVTVAR <- function(seed = NULL, iter, obs, P, H, a_lambda, b_lambda, a, a_tau, b_tau, beta1, beta2, a_w, b_w, W_infinite,
                    a_sigma, b_sigma, theta_upper, theta_lower, kappa_upper, a_p1, b_p1, a_p2, b_p2, rho, A_marima,
                    lambda1_init = NULL, lambda2_init = NULL,
                    phi_init = NULL, tau_init = NULL, W1_init = NULL, W2_init = NULL, W3_init = NULL, z_init = NULL,
                    nu_init = NULL, omega_init = NULL, sigma_init = NULL, theta_init = NULL, kappa_init = NULL, auxiliary_init = NULL, alpha1_init = NULL,
                    alpha2_init = NULL, alpha3_init = NULL, gamma_init = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  } 
  Time <- dim(obs)[1]
  N <- dim(obs)[2]
  lambda1 <- matrix(NA, iter, H)
  lambda2 <- matrix(NA, iter, H)
  phi <- array(NA, dim = c(iter, H))
  tau <- rep(NA, iter)
  W1 <- array(NA, dim = c(iter, H, N))
  W2 <- array(NA, dim = c(iter, H, N))
  W3 <- array(NA, dim = c(iter, H, P))
  z <- array(NA, dim = c(iter, H, P))
  nu <- array(NA, dim = c(iter, H, P))
  omega <- array(NA, dim = c(iter, H, P))
  sigma <- array(NA, dim = c(iter, N))
  # p1 <- matrix(NA, iter, H)
  # p2 <- matrix(NA, iter, H)
  theta <- matrix(NA, iter, H)
  kappa <- matrix(NA, iter, H)
  auxiliary <- array(NA, c(iter, H, Time-P))
  
  ### INITIALIZATION
  if (is.null(lambda1_init)) {
    lambda1[1,] <- rgamma(H, a_lambda, b_lambda)
  } else {
    lambda1[1,] <- lambda1_init
  }
  if (is.null(lambda2_init)) {
    lambda2[1,] <- rgamma(H, a_lambda, b_lambda)
  } else {
    lambda2[1,] <- lambda2_init
  }
  if (is.null(phi_init)) {
    phi[1,] <- gtools::rdirichlet(1, rep(a, H))
  } else {
    phi[1,] <- phi_init
  }
  if (is.null(tau_init)) {
    tau[1] <- rgamma(1, a_tau, b_tau)
  } else {
    tau[1] <- tau_init
  }
  if (is.null(W1_init)) {
    W1[1,,] <- matrix(rexp(H*N, 0.5*lambda1[1,]^2), H, N)
  } else {
    W1[1,,] <- W1_init
  }
  if (is.null(W2_init)) {
    W2[1,,] <- matrix(rexp(H*N, 0.5*lambda2[1,]^2), H, N)
  } else {
    W2[1,,] <- W2_init
  }
  if (is.null(nu_init)) {
    nu[1,,] <- cbind(matrix(rbeta(H*(P-1), beta1, beta2), H, P-1), rep(1,H))
  } else {
    nu[1,,] <- nu_init
  }
  if (is.null(omega_init)) {
    omega[1,,] <- abind::adrop(nu[1,,,drop=FALSE], drop = 1)*
      t(apply(cbind(rep(1,H), (1-abind::adrop(nu[1,,,drop=FALSE], drop = 1)))[,-(P+1),drop=FALSE], 1, cumprod))
  } else {
    omega[1,,] <- omega_init
  }
  if (is.null(z_init)) {
    z[1,,] <- t(apply(abind::adrop(omega[1,,,drop=FALSE], drop = 1), 1, function(x) sample(1:P, P, replace = TRUE, prob = x)))
  } else {
    z[1,,] <- z_init
  }
  if (is.null(W3_init)) {
    W3[1,,] <- t(apply(abind::adrop(z[1,,,drop=FALSE], drop = 1), 1, function(x) return(x <= 1:P)))
    W3[1,,] <- apply(abind::adrop(W3[1,,,drop=FALSE], drop = 1), 1:2, function(x) (1-x)/rgamma(1,a_w,b_w) + x*W_infinite)
  } else {
    W3[1,,] <- W3_init
  }
  if (is.null(sigma_init)) {
    sigma[1,] <- 1/rgamma(N, a_sigma, b_sigma)
  } else {
    sigma[1,] <- sigma_init
  }
  # if (is.null(p1_init)) {
  #   p1[1,] <- rbeta(H, a_p1, b_p1)
  # } else {
  #   p1[1,] <- p1_init
  # }
  # if (is.null(p2_init)) {
  #   p2[1,] <- rbeta(H, a_p2, b_p2)
  # } else {
  #   p2[1,] <- p2_init
  # }
  if (is.null(theta_init)) {
    theta[1,] <- runif(H, theta_lower, theta_upper)
  } else {
    theta[1,] <- theta_init
  }
  if (is.null(kappa_init)) {
    kappa[1,] <- runif(H, 0, kappa_upper)
  } else {
    kappa[1,] <- kappa_init
  }
  if (is.null(auxiliary_init)) {
    for (h in 1:H) {
      auxiliary[1,h,] <-
        NDARMA1(Time, P,
                exp(theta[1,h])*(exp(kappa[1,h])-1)/((exp(theta[1,h]+kappa[1,h])+1)*(exp(theta[1,h])+1)),
                exp(theta[1,h])*(exp(theta[1,h]+kappa[1,h])+1)/(exp(2*theta[1,h]+kappa[1,h])+2*exp(theta[1,h])+1))
    }
  } else {
    auxiliary[1,,] <- auxiliary_init
  }
  
  ### A FUNCTION OF DATA THAT WILL BE USED IN MCMC STEP
  vectX <- matrix(NA, Time-P, N*P)
  for (t in 1:P) {
    vectX[,(N*t-N+1):(N*t)] <- obs[(P-t+1):(Time-t),]
  }
  
  ### VARIABLES THAT WILL BE SAVED
  alpha1 <- array(NA, dim = c(iter, H, N))
  alpha2 <- array(NA, dim = c(iter, H, N))
  alpha3 <- array(NA, dim = c(iter, H, P))
  gamma <- array(NA, dim = c(iter, H, Time-P))
  
  ### INITIALIZATION
  if (is.null(alpha1_init)) {
    A_marima_tensor <- rTensor::as.tensor(array(A_marima, dim = c(N,N,P)))
    cp_A_marima <- rTensor::cp(A_marima_tensor, num_components = H, max_iter = 500, tol = 1e-05)
    alpha1[1,,] <- t(rTensor::khatri_rao(cp_A_marima$U[[1]], t(cp_A_marima$lambdas^(1/3))))
    alpha2[1,,] <- t(rTensor::khatri_rao(cp_A_marima$U[[2]], t(cp_A_marima$lambdas^(1/3))))
    alpha3[1,,] <- t(rTensor::khatri_rao(cp_A_marima$U[[3]], t(cp_A_marima$lambdas^(1/3)))) #initial values from CP decomposition after lasso
  } else {
    alpha1[1,,] <- alpha1_init
    alpha2[1,,] <- alpha2_init
    alpha3[1,,] <- alpha3_init
  }
  if (is.null(gamma_init)) {
    for (h in 1:H) {
      gamma[1,h,] <-
        NDARMA1(Time, P,
                exp(theta[1,h])*(exp(kappa[1,h])-1)/((exp(theta[1,h]+kappa[1,h])+1)*(exp(theta[1,h])+1)),
                exp(theta[1,h])*(exp(theta[1,h]+kappa[1,h])+1)/(exp(2*theta[1,h]+kappa[1,h])+2*exp(theta[1,h])+1))
    }
  } else {
    gamma[1,,] <- gamma_init
  }
  
  ### START MCMC
  for (it in 2:iter) {
    ## 1. update  phi, tau given A and W
    C <- sapply(1:H, function(x) {sum(crossprod(alpha1[it-1,x,],diag(1/W1[it-1,x,])) %*% alpha1[it-1,x,],
                                      crossprod(alpha2[it-1,x,],diag(1/W2[it-1,x,])) %*% alpha2[it-1,x,],
                                      crossprod(alpha3[it-1,x,],diag(1/W3[it-1,x,])) %*% alpha3[it-1,x,])})
    phi[it,] <- sapply(C, function(x) GIGrvg::rgig(1, lambda=a-N-P/2, chi=x, psi=2*b_tau))
    phi[it,] <- phi[it,]/sum(phi[it,])
    tau[it] <- GIGrvg::rgig(1, lambda=H*a-H*(N+P/2), chi=sum(C/phi[it,]), psi=2*b_tau)
    
    ## 2. update A, W given phi, tau, sigma, y
    ## 2.1. update W given A, phi, tau, sigma and y
    lambda1[it,] <- sapply(1:H, function(x) return(
      rgamma(1, a_lambda+N, b_lambda+sum(abs(alpha1[it-1,x,]))/sqrt(tau[it]*phi[it,x]))))
    lambda2[it,] <- sapply(1:H, function(x) return(
      rgamma(1, a_lambda+N, b_lambda+sum(abs(alpha2[it-1,x,]))/sqrt(tau[it]*phi[it,x]))))
    W1[it,,] <- sapply(1:N, function(x) return(sapply(1:H, function(y) return(
      GIGrvg::rgig(1, lambda=1/2, chi=alpha1[it-1,y,x]^2/(tau[it]*phi[it,y]), psi=lambda1[it,y]^2)))))
    W2[it,,] <- sapply(1:N, function(x) return(sapply(1:H, function(y) return(
      GIGrvg::rgig(1, lambda=1/2, chi=alpha2[it-1,y,x]^2/(tau[it]*phi[it,y]), psi=lambda2[it,y]^2)))))
    nu[it,,] <- cbind(matrix(sapply(1:(P-1), function(x) return(sapply(1:H, function(y) return(
      rbeta(1, beta1+sum(z[it-1,y,]==x), beta2+sum(z[it-1,y,]>x)))))), H, P-1),rep(1,H))
    omega[it,,] <- abind::adrop(nu[it,,,drop=FALSE], drop = 1)*
      t(apply(cbind(rep(1,H), (1-abind::adrop(nu[it,,,drop=FALSE], drop = 1)))[,-(P+1),drop=FALSE], 1, cumprod))
    z[it,,] <- matrix(sapply(1:P, function(x) return(sapply(1:H, function(y) return(
      draw_z(x, y, abind::adrop(omega[it,,,drop=F], 1), abind::adrop(alpha3[it-1,,,drop=F], 1), tau[it], phi[it,], W_infinite, a_w, b_w))))), H, P)
    W3[it,,] <- t(apply(abind::adrop(z[it,,,drop=F], 1), 1, function(x) return(x <= 1:P)))
    W <- sapply(1:P, function(x) return(sapply(1:H, function(y) return(
      1/rgamma(1, a_w+1/2, b_w+alpha3[it-1,y,x]^2/(2*tau[it]*phi[it,y]))))))
    W3[it,,] <- (1-W3[it,,])*W + W3[it,,]*W_infinite
    
    ## 2.2 update A given W, phi, tau, sigma and y
    if (it==2) {Y <- vectX%*%rTensor::khatri_rao(t(abind::adrop(alpha3[1,,,drop=F],1)), t(abind::adrop(alpha2[1,,,drop=F],1)))}
    alpha1[it,,] <- sample_alpha1(H, N, Time, P, obs, Y, abind::adrop(gamma[it-1,,,drop=F],1), abind::adrop(alpha1[it-1,,,drop=F], 1),
                                  abind::adrop(W1[it,,,drop=F],1), tau[it], phi[it,], 1/rho*sigma[it-1,])
    alpha2[it,,] <- sample_alpha2(H, N, Time, P, obs, abind::adrop(gamma[it-1,,,drop=F],1), abind::adrop(alpha1[it,,,drop=F], 1),
                                  abind::adrop(alpha2[it-1,,,drop=F],1), abind::adrop(alpha3[it-1,,,drop=F],1),
                                  abind::adrop(W2[it,,,drop=F],1), tau[it], phi[it,], 1/rho*sigma[it-1,])
    alpha3[it,,] <- sample_alpha3(H, N, Time, P, obs, abind::adrop(gamma[it-1,,,drop=F],1), abind::adrop(alpha1[it,,,drop=F], 1),
                                  abind::adrop(alpha2[it,,,drop=F],1), abind::adrop(alpha3[it-1,,,drop=F],1),
                                  abind::adrop(W3[it,,,drop=F],1), tau[it], phi[it,], 1/rho*sigma[it-1,])
    
    ## 3. update gamma given A, y and sigma
    Y <- vectX %*% rTensor::khatri_rao(t(abind::adrop(alpha3[it,,,drop=F],1)), t(abind::adrop(alpha2[it,,,drop=F],1)))
    y_hat <- rTensor::khatri_rao(t(abind::adrop(alpha1[it,,,drop=F],1)), Y)
    y_hat <- array(y_hat, c(Time-P, N, H))
    gamma[it,,] <- sample_gamma(H, N, Time, P, obs, abind::adrop(gamma[it-1,,,drop=F],1), y_hat,
                                matrix(c(theta[it-1,], rep(theta_star(theta[it-1,], kappa[it-1,]),
                                                           times=Time-P-2), theta[it-1,]), H, Time-P),
                                matrix(kappa[it-1,], H, Time-P-1), 1/rho*sigma[it-1,])
    y_hat <- aperm(y_hat, c(2,3,1))
    
    mple <- t(sapply(1:H, function(h) optim(c(0,0), function(theta)
      psudo_likelihood(Time, P, c(theta[1], rep(theta_star(theta[1], theta[2]), times=Time-P-2), theta[1]),
                       rep(theta[2], Time-P-1), gamma[it,h,]))$par))
    theta_kappa_auxiliary <- sample_theta(H, N, Time, P, abind::adrop(gamma[it,,,drop=F],1), theta[it-1,], kappa[it-1,],
                                          mple, theta_upper, theta_lower, kappa_upper, a_p1, b_p1, a_p2, b_p2,
                                          abind::adrop(auxiliary[it-1,,,drop=F],1))
    theta[it,] <- theta_kappa_auxiliary[[1]]
    kappa[it,] <- theta_kappa_auxiliary[[2]]
    auxiliary[it,,] <- theta_kappa_auxiliary[[3]]
    
    ## 4. update sigma given A, gamma and y
    obs_tilde <- obs[(P+1):Time,]-matrix(apply(sapply(1:H, function(x) return(gamma[it,x,]*t(y_hat[,x,]))), 1,
                                               sum), Time-P, N)
    obs_tilde <- obs_tilde^2
    sigma[it,] <- 1/rgamma(N, a_sigma+rho*(Time-P)/2, b_sigma+rho*0.5*apply(obs_tilde, 2, sum))
  }
  return(list('rho' = rho, 'lambda1' = lambda1, 'lambda2' = lambda2, 'phi' = phi, 'tau' = tau, 'W1' = W1,
              'W2' = W2, 'W3' = W3, 'z' = z, 'nu' = nu, 'omega' = omega, 'sigma' = sigma, 'theta' = theta,
              'kappa' = kappa, 'auxiliary' = auxiliary, 'alpha1' = alpha1,
              'alpha2' = alpha2, 'alpha3' = alpha3, 'gamma' = gamma, 'lambda1_last' = lambda1[iter,],
              'lambda2_last' = lambda2[iter,], 'phi_last' = phi[iter,], 'tau_last' = tau[iter], 'W1_last' = W1[iter,,],
              'W2_last' = W2[iter,,], 'W3_last' = W3[iter,,], 'z_last' = z[iter,,], 'nu_last' = nu[iter,,],
              'omega_last' = omega[iter,,], 'sigma_last' = sigma[iter,], 'theta_last' = theta[iter,],
              'kappa_last' = kappa[iter,], 'auxiliary_last' = auxiliary[iter,,], 'alpha1_last' = alpha1[iter,,],
              'alpha2_last' = alpha2[iter,,], 'alpha3_last' = alpha3[iter,,], 'gamma_last' = gamma[iter,,]))
}

## tempering
reps <- 10 
steps <- 100
iter <- 300
Time <- dim(obs)[1]
H <- 8
P <- 4
Mod <- define.model(kvar=N, ar=1:P)
A_marima <- matrix(NA, N, N*P)
obs.fit <- marima(obs, Mod$ar.pattern)
A_marima <- -matrix(obs.fit$ar.estimates[,,-1], N, N*P)

### HYPERPARAMETERS
a_lambda <- 3
b_lambda <- 1            # hyperparameters for lambda
a <- 1/H
a_tau <- H*1/H
b_tau <- H^2                          # hyperparameters for tau
beta1 <- 1
beta2 <- 5                            # hyperparameters for nu
a_w <- 2
b_w <- 2                            
W_infinite <- 0.01                    # hyperparameters for W                
a_sigma <- 10
b_sigma <- 2                          # hyperparameters for sigma
# a_p1 <- 100
# b_p1 <- 10
# a_p2 <- 50
# b_p2 <- 500
theta_upper <- 0
theta_lower <- -8
kappa_upper <- 10
rho <- seq(1,0.1,-0.1)
chains <- length(rho)
seed <- 2022
init <- list()
set.seed(seed)

for (rr in 1:reps) {
  lambda1 <- array(NA, c(chains, steps*iter, H))
  lambda2 <- array(NA, c(chains, steps*iter, H))
  phi <- array(NA, dim = c(chains, steps*iter, H))
  tau <- matrix(NA, chains, steps*iter)
  W1 <- array(NA, dim = c(chains, steps*iter, H, N))
  W2 <- array(NA, dim = c(chains, steps*iter, H, N))
  W3 <- array(NA, dim = c(chains, steps*iter, H, P))
  z <- array(NA, dim = c(chains, steps*iter, H, P))
  nu <- array(NA, dim = c(chains, steps*iter, H, P))
  omega <- array(NA, dim = c(chains, steps*iter, H, P))
  sigma <- array(NA, dim = c(chains, steps*iter, N))
  # p1 <- array(NA, c(chains, steps*iter, H))
  # p2 <- array(NA, c(chains, steps*iter, H))
  theta <- array(NA, c(chains, steps*iter, H))
  kappa <- array(NA, c(chains, steps*iter, H))
  auxiliary <- array(NA, c(chains, steps*iter, H, Time-P))
  alpha1 <- array(NA, dim = c(chains, steps*iter, H, N))
  alpha2 <- array(NA, dim = c(chains, steps*iter, H, N))
  alpha3 <- array(NA, dim = c(chains, steps*iter, H, P))
  gamma <- array(NA, dim = c(chains, steps*iter, H, Time-P))
  swap_acceptance <- rep(0, steps)
  
  for (jj in 1:steps) {
    cl <- makePSOCKcluster(chains, outfile="", master=nsl(Sys.info()['nodename']))
    registerDoSNOW(cl)
    if (rr==1 & jj==1) {
      res <- foreach(ii=rho, .multicombine=TRUE, .packages = 'Rcpp',
                     .noexport = c('sample_alpha1', 'sample_alpha2', 'sample_alpha3', 'sample_gamma',
                                   'sample_theta', 'NDARMA1', 'psudo_likelihood')) %dopar% {
                                     source('tempering_fun.R')
                                     BTVTVAR(NULL, iter, obs, P, H, a_lambda, b_lambda, a, a_tau, b_tau, beta1, beta2, a_w, b_w,
                                             W_infinite, a_sigma, b_sigma, theta_upper, theta_lower, kappa_upper, NULL, NULL, NULL, NULL,
                                             ii, A_marima)
                                   }
    } else {
      res <- foreach(ii=1:chains, .multicombine=TRUE, .packages = 'Rcpp',
                     .noexport = c('sample_alpha1', 'sample_alpha2', 'sample_alpha3', 'sample_gamma',
                                   'sample_theta', 'NDARMA1', 'psudo_likelihood')) %dopar% {
                                     source('tempering_fun.R')
                                     BTVTVAR(NULL, iter, obs, P, H, a_lambda, b_lambda, a, a_tau, b_tau, beta1, beta2, a_w, b_w,
                                             W_infinite, a_sigma, b_sigma, theta_upper, theta_lower, kappa_upper, NULL, NULL, NULL, NULL,
                                             rho[ii], A_marima, init[[ii]]$lambda1_init, init[[ii]]$lambda2_init,
                                             init[[ii]]$phi_init, init[[ii]]$tau_init, init[[ii]]$W1_init,
                                             init[[ii]]$W2_init, init[[ii]]$W3_init, init[[ii]]$z_init,
                                             init[[ii]]$nu_init, init[[ii]]$omega_init, init[[ii]]$sigma_init,
                                             init[[ii]]$theta_init, init[[ii]]$kappa_init, init[[ii]]$auxiliary_init,
                                             init[[ii]]$alpha1_init, init[[ii]]$alpha2_init, init[[ii]]$alpha3_init,
                                             init[[ii]]$gamma_init)
                                   }
    }
    stopCluster(cl)
    for (ii in 1:chains) {
      lambda1[ii,((jj-1)*iter+1):(jj*iter),] <- res[[ii]]$lambda1
      lambda2[ii,((jj-1)*iter+1):(jj*iter),] <- res[[ii]]$lambda2
      phi[ii,((jj-1)*iter+1):(jj*iter),] <- res[[ii]]$phi
      tau[ii,((jj-1)*iter+1):(jj*iter)] <- res[[ii]]$tau
      W1[ii,((jj-1)*iter+1):(jj*iter),,] <- res[[ii]]$W1
      W2[ii,((jj-1)*iter+1):(jj*iter),,] <- res[[ii]]$W2
      W3[ii,((jj-1)*iter+1):(jj*iter),,] <- res[[ii]]$W3
      z[ii,((jj-1)*iter+1):(jj*iter),,] <- res[[ii]]$z
      nu[ii,((jj-1)*iter+1):(jj*iter),,] <- res[[ii]]$nu
      omega[ii,((jj-1)*iter+1):(jj*iter),,] <- res[[ii]]$omega
      sigma[ii,((jj-1)*iter+1):(jj*iter),] <- res[[ii]]$sigma
      theta[ii,((jj-1)*iter+1):(jj*iter),] <- res[[ii]]$theta
      kappa[ii,((jj-1)*iter+1):(jj*iter),] <- res[[ii]]$kappa
      auxiliary[ii,((jj-1)*iter+1):(jj*iter),,] <- res[[ii]]$auxiliary
      alpha1[ii,((jj-1)*iter+1):(jj*iter),,] <- res[[ii]]$alpha1
      alpha2[ii,((jj-1)*iter+1):(jj*iter),,] <- res[[ii]]$alpha2
      alpha3[ii,((jj-1)*iter+1):(jj*iter),,] <- res[[ii]]$alpha3
      gamma[ii,((jj-1)*iter+1):(jj*iter),,] <- res[[ii]]$gamma
    }
    swap <- sort(sample(1:chains, 2, replace = FALSE))
    for (ii in (1:chains)) {
      init[[ii]] <- list('lambda1_init' = res[[ii]]$lambda1_last, 'lambda2_init' = res[[ii]]$lambda2_last,
                         'phi_init' = res[[ii]]$phi_last, 'tau_init' = res[[ii]]$tau_last, 'W1_init' = res[[ii]]$W1_last,
                         'W2_init' = res[[ii]]$W2_last, 'W3_init' = res[[ii]]$W3_last, 'z_init' = res[[ii]]$z_last,
                         'nu_init' = res[[ii]]$nu_last, 'omega_init' = res[[ii]]$omega_last,
                         'sigma_init' = res[[ii]]$sigma_last, 'theta_init' = res[[ii]]$theta_last,
                         'kappa_init' = res[[ii]]$kappa_last, 'auxiliary_init' = res[[ii]]$auxiliary_last,
                         'alpha1_init' = res[[ii]]$alpha1_last, 'alpha2_init' = res[[ii]]$alpha2_last,
                         'alpha3_init' = res[[ii]]$alpha3_last, 'gamma_init' = res[[ii]]$gamma_last)
    }
    MH_ratio <- log_likelihood(obs, res[[swap[1]]]$alpha1_last, res[[swap[1]]]$alpha2_last, res[[swap[1]]]$alpha3_last,
                               res[[swap[1]]]$gamma_last, res[[swap[1]]]$sigma_last, res[[swap[2]]]$rho) +
      log_likelihood(obs, res[[swap[2]]]$alpha1_last, res[[swap[2]]]$alpha2_last, res[[swap[2]]]$alpha3_last,
                     res[[swap[2]]]$gamma_last, res[[swap[2]]]$sigma_last, res[[swap[1]]]$rho) -
      log_likelihood(obs, res[[swap[1]]]$alpha1_last, res[[swap[1]]]$alpha2_last, res[[swap[1]]]$alpha3_last,
                     res[[swap[1]]]$gamma_last, res[[swap[1]]]$sigma_last, res[[swap[1]]]$rho) -
      log_likelihood(obs, res[[swap[2]]]$alpha1_last, res[[swap[2]]]$alpha2_last, res[[swap[2]]]$alpha3_last,
                     res[[swap[2]]]$gamma_last, res[[swap[2]]]$sigma_last, res[[swap[2]]]$rho)
    if (exp(MH_ratio)>runif(1)) {
      init[[swap[1]]] <- list('lambda1_init' = res[[swap[2]]]$lambda1_last, 'lambda2_init' = res[[swap[2]]]$lambda2_last,
                              'phi_init' = res[[swap[2]]]$phi_last, 'tau_init' = res[[swap[2]]]$tau_last, 'W1_init' = res[[swap[2]]]$W1_last,
                              'W2_init' = res[[swap[2]]]$W2_last, 'W3_init' = res[[swap[2]]]$W3_last, 'z_init' = res[[swap[2]]]$z_last,
                              'nu_init' = res[[swap[2]]]$nu_last, 'omega_init' = res[[swap[2]]]$omega_last,
                              'sigma_init' = res[[swap[2]]]$sigma_last, 'theta_init' = res[[swap[2]]]$theta_last,
                              'kappa_init' = res[[swap[2]]]$kappa_last, 'auxiliary_init' = res[[swap[2]]]$auxiliary_last,
                              'alpha1_init' = res[[swap[2]]]$alpha1_last, 'alpha2_init' = res[[swap[2]]]$alpha2_last,
                              'alpha3_init' = res[[swap[2]]]$alpha3_last, 'gamma_init' = res[[swap[2]]]$gamma_last)
      init[[swap[2]]] <- list('lambda1_init' = res[[swap[1]]]$lambda1_last, 'lambda2_init' = res[[swap[1]]]$lambda2_last,
                              'phi_init' = res[[swap[1]]]$phi_last, 'tau_init' = res[[swap[1]]]$tau_last, 'W1_init' = res[[swap[1]]]$W1_last,
                              'W2_init' = res[[swap[1]]]$W2_last, 'W3_init' = res[[swap[1]]]$W3_last, 'z_init' = res[[swap[1]]]$z_last,
                              'nu_init' = res[[swap[1]]]$nu_last, 'omega_init' = res[[swap[1]]]$omega_last,
                              'sigma_init' = res[[swap[1]]]$sigma_last, 'theta_init' = res[[swap[1]]]$theta_last,
                              'kappa_init' = res[[swap[1]]]$kappa_last, 'auxiliary_init' = res[[swap[1]]]$auxiliary_last,
                              'alpha1_init' = res[[swap[1]]]$alpha1_last, 'alpha2_init' = res[[swap[1]]]$alpha2_last,
                              'alpha3_init' = res[[swap[1]]]$alpha3_last, 'gamma_init' = res[[swap[1]]]$gamma_last)
      swap_acceptance[jj] <- 1
    } 
  }
  
  mcmc_index <- seq(round(steps*iter/2), steps*iter, 10)
  A_BTVTVAR_cmp <- array(NA, c(H, N, N*P))
  for (hh in 1:H) {
    A_BTVTVAR_cmp[hh,,] <- dynamic_coefficients(Time, N, P, 1, array(1, c(length(mcmc_index),1,1)),
                                                adrop(alpha1[1,mcmc_index,hh,,drop=F], 1),
                                                adrop(alpha2[1,mcmc_index,hh,,drop=F], 1),
                                                adrop(alpha3[1,mcmc_index,hh,,drop=F], 1))
  }
  A_BTVTVAR <- dynamic_coefficients(Time, N, P, H, adrop(gamma[1,mcmc_index,,,drop=F], 1),
                                    adrop(alpha1[1,mcmc_index,,,drop=F], 1),
                                    adrop(alpha2[1,mcmc_index,,,drop=F], 1),
                                    adrop(alpha3[1,mcmc_index,,,drop=F], 1))
  # save.image(file = sprintf('BTVTVAR_tempering_seed%i_rep%i.RData', seed, rr))
  alpha1_samples <- alpha1[1,,,]
  alpha2_samples <- alpha2[1,,,]
  alpha3_samples <- alpha3[1,,,]
  gamma_samples <- gamma[1,,,]
  save(obs, s, jocobian_omega, jocobian_Ising, log_likelihood, draw_z, fun_theta, fun_kappa, theta_star, BTVTVAR,
       reps, steps, iter, Time, H, P, A_marima, a_lambda, b_lambda, a, a_tau, b_tau, beta1, beta2, a_w, b_w,
       W_infinite, a_sigma, b_sigma, theta_upper, theta_lower, kappa_upper, rho, chains, seed, alpha1_samples,
       alpha2_samples, alpha3_samples, gamma_samples, mcmc_index, A_BTVTVAR_cmp, A_BTVTVAR,
       file = sprintf('BTVTVAR_tempering_sub%i_seed%i_rep%i.RData', s, seed, rr))
}
