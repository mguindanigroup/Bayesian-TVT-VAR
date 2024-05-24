library(abind)
library(rTensor)
library(GIGrvg)
library(gtools)

install.packages("BTVTVAR_1.0.tar.gz", repos = NULL, type = "source")
library(BTVTVAR)

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

BTVT_VAR <- function(iter, H, obs, A_init, burn_in, thin) {
  Time <- dim(obs)[1]
  N <- dim(obs)[2]
  P <- dim(A_init)[2]/dim(A_init)[1]
  
  a_lambda <- 3
  b_lambda <- 1                         # hyperparameters for lambda
  a <- 1/H
  a_tau <- H*1/H
  b_tau <- H^4                          # hyperparameters for tau
  beta1 <- 1
  beta2 <- 5                            # hyperparameters for nu
  a_w <- 2
  b_w <- 2                            
  W_infinite <- 0.01                    # hyperparameters for W                
  a_sigma <- 10
  b_sigma <- 2                          # hyperparameters for sigma
  theta_upper <- 8
  theta_lower <- -8
  kappa_upper <- 10
  
  lambda1 <- matrix(NA, iter, H)
  lambda2 <- matrix(NA, iter, H)
  phi <- matrix(NA, iter, H)
  tau <- rep(NA, iter)
  W1 <- array(NA, dim = c(iter, H, N))
  W2 <- array(NA, dim = c(iter, H, N))
  W3 <- array(NA, dim = c(iter, H, P))
  z <- array(NA, dim = c(iter, H, P))
  nu <- array(NA, dim = c(iter, H, P))
  omega <- array(NA, dim = c(iter, H, P))
  sigma <- matrix(NA, iter, N)
  theta <- matrix(NA, iter, H)
  kappa <- matrix(NA, iter, H)
  auxiliary <- array(NA, c(iter, H, Time-P))
  
  ### INITIALIZATION
  lambda1[1,] <- rgamma(H, a_lambda, b_lambda)
  lambda2[1,] <- rgamma(H, a_lambda, b_lambda)
  phi[1,] <- rdirichlet(1, rep(a, H))
  tau[1] <- rgamma(1, a_tau, b_tau)
  W1[1,,] <- matrix(rexp(H*N, 0.5*lambda1[1,]^2), H, N)
  W2[1,,] <- matrix(rexp(H*N, 0.5*lambda2[1,]^2), H, N)
  nu[1,,] <- cbind(matrix(rbeta(H*(P-1), beta1, beta2), H, P-1), rep(1,H))
  omega[1,,] <- adrop(nu[1,,,drop=FALSE], drop = 1)*
    t(apply(cbind(rep(1,H), (1-adrop(nu[1,,,drop=FALSE], drop = 1)))[,-(P+1),drop=FALSE], 1, cumprod))
  z[1,,] <- t(apply(adrop(omega[1,,,drop=FALSE], drop = 1), 1,
                    function(x) sample(1:P, P, replace = TRUE, prob = x)))
  W3[1,,] <- t(apply(adrop(z[1,,,drop=FALSE], drop = 1), 1, function(x) return(x <= 1:P)))
  W3[1,,] <- apply(adrop(W3[1,,,drop=FALSE], drop = 1), 1:2, function(x) (1-x)/rgamma(1,a_w,b_w) + x*W_infinite)
  sigma[1,] <- 1/rgamma(N, a_sigma, b_sigma)
  theta[1,] <- runif(H, theta_lower, theta_upper)
  int_theta[1,] <- runif(H, 0, int_theta_upper)
  for (h in 1:H) {
    auxiliary[1,h,] <-
      NDARMA1(Time, P,
              exp(theta[1,h])*(exp(kappa[1,h])-1)/((exp(theta[1,h]+kappa[1,h])+1)*(exp(theta[1,h])+1)),
              exp(theta[1,h])*(exp(theta[1,h]+kappa[1,h])+1)/(exp(2*theta[1,h]+kappa[1,h])+2*exp(theta[1,h])+1))
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
  A_init_tensor <- as.tensor(array(A_init, dim = c(N,N,P)))
  cp_A_init <- cp(A_init_tensor, num_components = H, max_iter = 500, tol = 1e-05)
  alpha1[1,,] <- t(khatri_rao(cp_A_init$U[[1]], t(cp_A_init$lambdas^(1/3))))
  alpha2[1,,] <- t(khatri_rao(cp_A_init$U[[2]], t(cp_A_init$lambdas^(1/3))))
  alpha3[1,,] <- t(khatri_rao(cp_A_init$U[[3]], t(cp_A_init$lambdas^(1/3)))) #initial values from CP decomposition after lasso
  for (h in 1:H) {
    gamma[1,h,] <-
      NDARMA1(Time, P,
              exp(theta[1,h])*(exp(kappa[1,h])-1)/((exp(theta[1,h]+kappa[1,h])+1)*(exp(theta[1,h])+1)),
              exp(theta[1,h])*(exp(theta[1,h]+kappa[1,h])+1)/(exp(2*theta[1,h]+kappa[1,h])+2*exp(theta[1,h])+1))
  }
  
  ### START MCMC
  for (it in 2:iter) {
    ## 1. update  phi, tau given A and W
    C <- sapply(1:H, function(x) {sum(crossprod(alpha1[it-1,x,],diag(1/W1[it-1,x,])) %*% alpha1[it-1,x,],
                                      crossprod(alpha2[it-1,x,],diag(1/W2[it-1,x,])) %*% alpha2[it-1,x,],
                                      crossprod(alpha3[it-1,x,],diag(1/W3[it-1,x,])) %*% alpha3[it-1,x,])})
    phi[it,] <- sapply(C, function(x) rgig(1, lambda=a-N-P/2, chi=x, psi=2*b_tau))
    phi[it,] <- phi[it,]/sum(phi[it,])
    tau[it] <- rgig(1, lambda=H*a-H*(N+P/2), chi=sum(C/phi[it,]), psi=2*b_tau)
    
    ## 2. update A, W given phi, tau, sigma, y
    ## 2.1. update W given A, phi, tau, sigma and y
    lambda1[it,] <- sapply(1:H, function(x) return(
      rgamma(1, a_lambda+N, b_lambda+sum(abs(alpha1[it-1,x,]))/sqrt(tau[it]*phi[it,x]))))
    lambda2[it,] <- sapply(1:H, function(x) return(
      rgamma(1, a_lambda+N, b_lambda+sum(abs(alpha2[it-1,x,]))/sqrt(tau[it]*phi[it,x]))))
    W1[it,,] <- sapply(1:N, function(x) return(sapply(1:H, function(y) return(
      rgig(1, lambda=1/2, chi=alpha1[it-1,y,x]^2/(tau[it]*phi[it,y]), psi=lambda1[it,y]^2)))))
    W2[it,,] <- sapply(1:N, function(x) return(sapply(1:H, function(y) return(
      rgig(1, lambda=1/2, chi=alpha2[it-1,y,x]^2/(tau[it]*phi[it,y]), psi=lambda2[it,y]^2)))))
    nu[it,,] <- cbind(matrix(sapply(1:(P-1), function(x) return(sapply(1:H, function(y) return(
      rbeta(1, beta1+sum(z[it-1,y,]==x), beta2+sum(z[it-1,y,]>x)))))), H, P-1),rep(1,H))
    omega[it,,] <- adrop(nu[it,,,drop=FALSE], drop = 1)*
      t(apply(cbind(rep(1,H), (1-adrop(nu[it,,,drop=FALSE], drop = 1)))[,-(P+1),drop=FALSE], 1, cumprod))
    z[it,,] <- matrix(sapply(1:P, function(x) return(sapply(1:H, function(y) return(
      draw_z(x, y, adrop(omega[it,,,drop=F], 1), adrop(alpha3[it-1,,,drop=F], 1), tau[it], phi[it,], W_infinite, a_w, b_w))))), H, P)
    W3[it,,] <- t(apply(adrop(z[it,,,drop=F], 1), 1, function(x) return(x <= 1:P)))
    W <- sapply(1:P, function(x) return(sapply(1:H, function(y) return(
      1/rgamma(1, a_w+1/2, b_w+alpha3[it-1,y,x]^2/(2*tau[it]*phi[it,y]))))))
    W3[it,,] <- (1-W3[it,,])*W + W3[it,,]*W_infinite
    
    ## 2.2 update A given W, phi, tau, sigma and y
    if (it==2) {Y <- vectX%*%khatri_rao(t(adrop(alpha3[1,,,drop=F],1)), t(adrop(alpha2[1,,,drop=F],1)))}
    alpha1[it,,] <- sample_alpha1(H, N, Time, P, obs, Y, adrop(gamma[it-1,,,drop=F],1), adrop(alpha1[it-1,,,drop=F], 1),
                                  adrop(W1[it,,,drop=F],1), tau[it], phi[it,], sigma[it-1,])
    alpha2[it,,] <- sample_alpha2(H, N, Time, P, obs, adrop(gamma[it-1,,,drop=F],1), adrop(alpha1[it,,,drop=F], 1),
                                  adrop(alpha2[it-1,,,drop=F],1), adrop(alpha3[it-1,,,drop=F],1),
                                  adrop(W2[it,,,drop=F],1), tau[it], phi[it,], sigma[it-1,])
    alpha3[it,,] <- sample_alpha3(H, N, Time, P, obs, adrop(gamma[it-1,,,drop=F],1), adrop(alpha1[it,,,drop=F], 1),
                                  adrop(alpha2[it,,,drop=F],1), adrop(alpha3[it-1,,,drop=F],1),
                                  adrop(W3[it,,,drop=F],1), tau[it], phi[it,], sigma[it-1,])
    
    ## 3. update gamma given A, y and sigma
    Y <- vectX %*% khatri_rao(t(adrop(alpha3[it,,,drop=F],1)), t(adrop(alpha2[it,,,drop=F],1)))
    y_hat <- khatri_rao(t(adrop(alpha1[it,,,drop=F],1)), Y)
    y_hat <- array(y_hat, c(Time-P, N, H))
    gamma[it,,] <- sample_gamma(H, N, Time, P, obs, adrop(gamma[it-1,,,drop=F],1), y_hat,
                                matrix(c(theta[it-1,], rep(theta_star(theta[it-1,], kappa[it-1,]),
                                                           times=Time-P-2), theta[it-1,]), H, Time-P),
                                matrix(kappa[it-1,], H, Time-P-1), sigma[it-1,])
    y_hat <- aperm(y_hat, c(2,3,1))
    
    mple <- t(sapply(1:H, function(h) optim(c(0,0), function(theta)
      psudo_likelihood(Time, P, c(theta[1], rep(theta_star(theta[1], theta[2]), times=Time-P-2), theta[1]),
                       rep(theta[2], Time-P-1), gamma[it,h,]))$par))
    theta_kappa_auxiliary <- sample_theta(H, N, Time, P, adrop(gamma[it,,,drop=F],1), theta[it-1,], kappa[it-1,],
                                          mple, theta_upper, theta_lower, kappa_upper, NULL, NULL, NULL, NULL,
                                          adrop(auxiliary[it-1,,,drop=F],1))
    theta[it,] <- theta_kappa_auxiliary[[1]]
    kappa[it,] <- theta_kappa_auxiliary[[2]]
    auxiliary[it,,] <- theta_kappa_auxiliary[[3]]
    
    ## 4. update sigma given A, gamma and y
    obs_tilde <- obs[(P+1):Time,]-matrix(apply(sapply(1:H, function(x) return(gamma[it,x,]*t(y_hat[,x,]))), 1,
                                               sum), Time-P, N)
    obs_tilde <- obs_tilde^2
    sigma[it,] <- 1/rgamma(N, a_sigma+(Time-P)/2, b_sigma+0.5*apply(obs_tilde, 2, sum))
  }
  mcmc_index <- seq(burn_in, iter, thin)
  est_gamma <- apply(gamma[mcmc_index,,], 2:3, median)
  A_BTVTVAR_cmp <- array(NA, c(H, N, N*P))
  for (hh in 1:H) {
    A_BTVTVAR_cmp[hh,,] <- dynamic_coefficients(Time, N, P, 1, array(1, c(length(mcmc_index),1,1)),
                                                alpha1[mcmc_index,hh,,drop=F], alpha2[mcmc_index,hh,,drop=F],
                                                alpha3[mcmc_index,hh,,drop=F])
  }
  A_BTVTVAR <- dynamic_coefficients(Time, N, P, H, gamma[mcmc_index,,], alpha1[mcmc_index,,],
                                    alpha2[mcmc_index,,], alpha3[mcmc_index,,])
  return(list(alpha1=alpha1, alpha2=alpha2, alpha3=alpha3, gamma=gamma, est_gamma=est_gamma, A_BTVTVAR=A_BTVTVAR,
              A_BTVTVAR_cmp=A_BTVTVAR_cmp))
}
