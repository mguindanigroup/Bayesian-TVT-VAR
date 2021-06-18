library(gtools)
library(marima)
library(tvReg)
library(doSNOW)
library(parallel)
library(foreach)

source('BTVT-VAR_fun.R')
source('SIMULATION1_fun.R')

trials <- 100
seeds_list <- 1:100
N <- 10
true_H <- 3
true_P <- 3
T <- 100

### TRUTH
true_alpha1 <- matrix(NA, true_H, N)
true_alpha2 <- matrix(NA, true_H, N)
true_alpha3 <- matrix(NA, true_H, true_P)
alpha_array <- array(NA, c(true_H, N, N*true_P))
true_gamma <- matrix(NA, true_H, T-true_P)
A <- array(0, c(T-true_P, N, N*true_P))
obs <- matrix(NA, T, N)
max_rg <- rep(NA, true_H)

### SET UP PARALLEL BACKEND
cl <- makePSOCKcluster(100, outfile="", master=nsl(Sys.info()['nodename']))
clusterExport(cl, c('combinations', 'define.model','marima', 'tvVAR', 'rdirichlet', 'adrop', 'BTVT_VAR',
                    'exact_sampling', 'sample_alpha1', 'sample_alpha2', 'sample_alpha3', 'sample_gamma',
                    'psudo_likelihood', 'sample_theta', 'dynamic_coefficients', 'as.tensor', 'cp',
                    'khatri_rao', 'rgig'))
registerDoSNOW(cl)

### SIMULATE COEFFICIENTS WITH STATIONARITY CHECK
output <- foreach (tl=1:100, .multicombine=TRUE) %dopar% {
  set.seed(seeds_list[tl])
  
  ## SIMULATE COEFFICIENTS WITH STATIONARITY CHECK
  for (h in 1:true_H) {
    G <- array(NA, dim = c(2^(h-1), N*true_P, N*true_P))
    while (TRUE) {
      true_alpha1[h,] <- random_coefficient(N, 0.5, 1)
      true_alpha2[h,] <- random_coefficient(N, 0.5, 1)
      for (t in 1:true_P) {
        true_alpha3[h,t] <- random_coefficient(1, 1-1/t, 1/t)
      }
      alpha_array[h,,] <- kronecker(t(true_alpha3[h,]), tcrossprod(true_alpha1[h,], true_alpha2[h,]))
      max_rg[h] <- max(abs(range(alpha_array[h,,])))
      G[1,,] <- rbind(alpha_array[h,,], cbind(diag(1,N*(true_P-1),N*(true_P-1)),
                                              matrix(0,N*(true_P-1),N)))
      if (h>=2) {
        for (hh in 2:2^(h-1)) {
          G[hh,,] <- rbind(alpha_array[h,,]+apply(alpha_array[cbnt(h)[[hh-1]],,,drop=FALSE], 2:3, sum),
                           cbind(diag(1,N*(true_P-1),N*(true_P-1)), matrix(0,N*(true_P-1),N)))
        }
      }
      if (all(sapply(1:2^(h-1), function(x) Mod(eigen(G[x,,], only.values = TRUE)$values[1])) < 0.9) &
          max_rg[h] > 0.5*max_rg[1] & max_rg[h] < 1.5*max_rg[1]) break
    }
  }
  
  ## SIMULATE BINARY INDICATORS
  for (h in 1:true_H) {
    while (TRUE) {
      true_gamma[h,] <- NDARMA(T-true_P, runif(1), runif(1))
      if (mean(true_gamma[h,]) > 0.3) break
    }
  }
  
  ## SIMULATE DATA
  for (t in 1:(T-true_P)) {
    for (h in 1:true_H) {
      A[t,,] <- A[t,,] + alpha_array[h,,]*true_gamma[h,t]
    }
  }
  error <- matrix(rnorm((T-true_P)*N, 0, sd = (1:10)/5), N, T-true_P)
  obs[1:true_P,] <- matrix(rnorm(true_P*N, 0, 1), true_P, N)
  for (t in (true_P+1):T) {
    obs[t,] <- t(A[t-true_P,,] %*% as.vector(t(obs[(t-1):(t-true_P),])) + error[,t-true_P])
  }
  
  ### MARIMA
  P <- 4
  Mod <- define.model(kvar=N, ar=1:P)
  A_marima <- matrix(NA, N, N*P)
  obs.fit <- marima(obs, Mod$ar.pattern)
  A_marima <- -matrix(obs.fit$ar.estimates[,,-1], N, N*P)
  
  ## TVREG PACKAGE
  tvp_var <- tvVAR(obs, P)
  A_tvp <- array(NA, c(T-P, N, N*P))
  for (i in 1:length(tvp_var$coefficients)) {
    A_tvp[,i,] <- tvp_var$coefficients[[i]][,-(N*P+1)]
  }
  
  ### START MCMC
  iter <- 5000
  H <- 4
  res_BTVTVAR <- BTVT_VAR(iter, H, obs, A_marima, 1666, 3)
  
  return(list(true_alpha1=true_alpha1, true_alpha2=true_alpha2, true_alpha3=true_alpha3,
              true_gamma=true_gamma, alpha_array=alpha_array, A=A, error=error, obs=obs, A_marima=A_marima,
              alpha1=res_BTVTVAR$alpha1, alpha2=res_BTVTVAR$alpha2, alpha3=res_BTVTVAR$alpha3,
              gamma=res_BTVTVAR$gamma, est_gamma=res_BTVTVAR$est_gamma, A_BTVTVAR=res_BTVTVAR$A_BTVTVAR,
              A_BTVTVAR_cmp=res_BTVTVAR$A_BTVTVAR_cmp, A_tvp=A_tvp))
}

stopCluster(cl)
save(output, file='SIMULATION1_PARALLEL.RData')