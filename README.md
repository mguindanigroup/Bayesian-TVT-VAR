# Bayesian Time-Varying Tensor Var Model 

This repository contains the raw code of the manuscript

"BAYESIAN  TIME-VARYING  TENSOR  VECTOR  AUTOREGRESSIVEmMODELS FOR DYNAMIC  EFFECTIVE  CONNECTIVITY"

by Wei Zhang, Ivor Cribben, Sonia Petrone and Michele Guindani

# SIMULATION1

SIMULATION1_DATA_GENERATION.R and SIMULATION1.R replicate the first simulation in the manuscript. They demonstrate an example of using the BTVT-VAR package.

```
library(gtools)
library(marima)
library(tvReg)
library(doSNOW)
library(parallel)
library(foreach)
library(mgm)

load('SIMULATION1_SNR2.RData')
source('BTVT-VAR_fun.R')

Rep <- 100
seeds_list <- 1:Rep
P <- 3
H <- 4
iter <- 10000
burn_in <- 3000
thin <- 3
mcmc_index <- seq(burn_in, iter, thin)

### SET UP PARALLEL BACKEND
cl <- makePSOCKcluster(Rep, outfile="", master=nsl(Sys.info()['nodename']))
clusterExport(cl, c('define.model','marima', 'tvVAR', 'bwSelect', 'tvmvar', 'rdirichlet', 'adrop', 'BTVT_VAR',
                    'sample_alpha1', 'sample_alpha2', 'sample_alpha3', 'sample_gamma',
                    'psudo_likelihood', 'sample_theta', 'dynamic_coefficients', 'as.tensor', 'cp',
                    'khatri_rao', 'rgig'))

registerDoSNOW(cl)
### SIMULATE COEFFICIENTS WITH STATIONARITY CHECK
output <- foreach (ii=1:Rep, .multicombine=TRUE) %dopar% {
  set.seed(seeds_list[ii])
  
  ### MARIMA
  Mod <- define.model(kvar=N, ar=1:P)
  A_marima <- matrix(NA, N, N*P)
  obs.fit <- marima(obs[ii,,], Mod$ar.pattern)
  A_marima <- -matrix(obs.fit$ar.estimates[,,-1], N, N*P)
  
  ### TVREG PACKAGE
  tvp_var <- tvVAR(obs[ii,,], P)
  A_tvp <- array(NA, c(Time-P, N, N*P))
  for (i in 1:length(tvp_var$coefficients)) {
    A_tvp[,i,] <- tvp_var$coefficients[[i]][,-(N*P+1)]
  }
  
  ## MGM PACKAGE
  bwSeq <- seq(0.1, 1, by = 0.1)
  bw_tvmvar <- bwSelect(data = obs[ii,,], type = rep('g', N), level = rep(1, N), lags = 1:P,
                        bwSeq = bwSeq, bwFolds = 5, bwFoldsize = 5, modeltype = 'mvar')
  tvmvar_fit <- tvmvar(obs[ii,,], rep('g', N), level = rep(1, N), lags = 1:P,
                       estpoints = seq(1/(Time-1), 1, length.out = Time-P), 
                       bandwidth = bwSeq[which.min(bw_tvmvar$meanError)])
  A_mgm <- array(NA, c(Time-P, N, N*P))
  for (tt in 1:(Time-P)) {
    A_mgm[tt,,] <- tvmvar_fit$signs[,,,tt]*tvmvar_fit$wadj[,,,tt]
  }
  A_mgm[is.na(A_mgm)] <- 0
  
  ### START MCMC
  tryCatch({
    res_BTVTVAR <- BTVT_VAR(iter, H, obs[ii,,], A_marima, burn_in, thin)
    return(list(true_alpha1=true_alpha1, true_alpha2=true_alpha2, true_alpha3=true_alpha3,
                true_gamma=true_gamma, true_A_short=true_A_short, err=err, obs=obs, A_marima=A_marima,
                alpha1=res_BTVTVAR$alpha1, alpha2=res_BTVTVAR$alpha2, alpha3=res_BTVTVAR$alpha3,
                gamma=res_BTVTVAR$gamma, est_gamma=res_BTVTVAR$est_gamma, A_BTVTVAR=res_BTVTVAR$A_BTVTVAR,
                A_BTVTVAR_cmp=res_BTVTVAR$A_BTVTVAR_cmp, A_tvp=A_tvp, A_mgm=A_mgm))
    }, error = function(cond) {
      cat(sprintf('error: seed %i', seeds_list[ii]))
      return(list(true_alpha1=true_alpha1, true_alpha2=true_alpha2, true_alpha3=true_alpha3,
                  true_gamma=true_gamma, true_A_short=true_A_short, err=err, obs=obs, A_marima=A_marima,
                  alpha1=NULL, alpha2=NULL, alpha3=NULL, gamma=NULL, est_gamma=NULL, A_BTVTVAR=NULL,
                  A_BTVTVAR_cmp=NULL, A_tvp=A_tvp, A_mgm=A_mgm))
      }
    )
}
stopCluster(cl)
save.image(file=sprintf('SIMULATION1_SNR2p%ih%i_PARALLEL.RData', P, H))
```

# Metropolis coupled Markov chain Monte Carlo algorithm

MC3.R contains an example of applying BTVT-VAR model with Metropolis coupled Markov chain Monte Carlo (MC^3) algorithm.

```
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
```

For questions and suggestions, please reach out to Wei Zhang(zhngw129) or Michele Guindani (mguindani).
