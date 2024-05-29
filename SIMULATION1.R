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
