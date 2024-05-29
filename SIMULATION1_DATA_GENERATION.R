library(tidyverse)
library(marima)
library(rTensor)
library(rlang)
library(gtools)
library(abind)

## import data from txt
Movie_LR_extracted <- read.table("Movie_LR_extracted.txt", sep = "", header = TRUE)
Movie_LR_extracted$run <- as.factor(Movie_LR_extracted$run)

## the 1290 time points corresponding to second subject
s <- 2
ts_Movie_LR <- Movie_LR_extracted[(1290*(s-1)+1):(1290*s),]

## plot the whole time series
ts_Movie_LR_long <- ts_Movie_LR %>%
  gather(key = "variable", value = "value", -c(run,time))
ggplot(ts_Movie_LR_long, aes(x = time, y = value, group = run)) +
  geom_line(aes(color = variable), size = 1) +
  theme_minimal()

## run 1  1-324; run 2 325-661; run 3 662-925; run 4 926-1290
## run 1 11-334; run 2 350-686; run 3 702-965; run 4 981-1345
## plot time series of certain run
run_time <- cbind(c(1,325,662,926), c(324,661,925,1290))
run <- 1
ts_Movie_LR_long <- ts_Movie_LR[run_time[run,1]:run_time[run,2],] %>%
  gather(key = "variable", value = "value", -c(run,time))
ggplot(ts_Movie_LR_long, aes(x = time, y = value, group = run)) +
  geom_line(aes(color = variable), size = 1) +
  theme_minimal()

## scaled observations to estimate VAR coefficients
obs <- scale(ts_Movie_LR[run_time[run,1]:run_time[run,2],-(1:2)], scale = FALSE)
Time <- dim(obs)[1]
N <- dim(obs)[2]

## MARIMA estimates
true_P <- 3
Mod <- define.model(kvar=N, ar=1:true_P)
obs.fit <- marima(obs, Mod$ar.pattern)
A_marima <- -matrix(obs.fit$ar.estimates[,,-1], N, N*true_P)
res_marima <- t(obs.fit$residuals)

## PARAFAC decomposition of the MARIMA coefficient tensor
true_H <- 2
A_marima_tensor <- as.tensor(array(A_marima, dim = c(N, N, true_P)))
cp_A_marima <- cp(A_marima_tensor, num_components = true_H, max_iter = 500, tol = 1e-05)
true_alpha1 <- t(khatri_rao(cp_A_marima$U[[1]], t(cp_A_marima$lambdas^(1/3))))
true_alpha2 <- t(khatri_rao(cp_A_marima$U[[2]], t(cp_A_marima$lambdas^(1/3))))
true_alpha3 <- t(khatri_rao(cp_A_marima$U[[3]], t(cp_A_marima$lambdas^(1/3))))
true_A <- t(true_alpha1) %*% t(khatri_rao(t(true_alpha3), t(true_alpha2)))
est_marima <- array(NA, dim = c(Time-true_P, N))
for (tt in (true_P+1):Time) {
  est_marima[tt-true_P,] <- t(true_A %*% as.vector(t(obs[(tt-1):(tt-true_P),])))
}
err_marima <- obs[(true_P+1):Time,] - est_marima

## true alpha obtained from estimated PARAFAC components of the MARIAM of two runs
env_s2r1 <- env()
load('Movie_LR_extracted_s2r1.RData', envir = env_s2r1)
env_s2r2 <- env()
load('Movie_LR_extracted_s2r2.RData', envir = env_s2r2)

set.seed(2021)
true_H <- 4
cbns <- list()
cbns_count <- 0
for (hh in 1:true_H) {
  for (ii in 1:choose(true_H, hh)) {
    cbns[[cbns_count+ii]] <- combinations(true_H, hh)[ii,]
  }
  cbns_count <- cbns_count + choose(true_H, hh)
}
true_alpha1_tmp <- rbind(env_s2r1$true_alpha1, env_s2r2$true_alpha1)
true_alpha2_tmp <- rbind(env_s2r1$true_alpha2, env_s2r2$true_alpha2)
true_alpha3_tmp <- rbind(env_s2r1$true_alpha3, env_s2r2$true_alpha3)

## define signal to noise ratio as 
## the maximum parameter size of the time-varying parameters divided by the standard deviation of the noise
## SNR=1.8518 in s2r1, SNR=2.7972 in s2r2
## use three SNR values, 0.5, 2.5, 10, denoted by SNR1, SNR2, SNR3 in the saved file name
s <- 2
Rep <- 100
SNR <- 2.5
N <- 27
true_P <- 3

true_alpha1 <- array(NA, c(Rep, true_H, N))
true_alpha2 <- array(NA, c(Rep, true_H, N))
true_alpha3 <- array(NA, c(Rep, true_H, true_P))
for (ii in 1:Rep) {
  G <- array(NA, dim = c(2^true_H-1, N*true_P, N*true_P))
  while (TRUE) {
    true_alpha1[ii,,] <- sapply(rnorm(N*true_H, true_alpha1_tmp, 0.1), function(x) if (abs(x)>0.3) x else 0)  
    true_alpha2[ii,,] <- sapply(rnorm(N*true_H, true_alpha2_tmp, 0.1), function(x) if (abs(x)>0.3) x else 0)
    true_alpha3[ii,,] <- sapply(rnorm(true_P*true_H, true_alpha3_tmp, 0.1), function(x) if (abs(x)>0.3) x else 0)
    for (jj in 1:(2^true_H-1)){
      G[jj,,] <- rbind(t(adrop(true_alpha1[ii,cbns[[jj]],,drop=F], 1)) %*%
                         t(khatri_rao(t(adrop(true_alpha3[ii,cbns[[jj]],,drop=F], 1)),
                                      t(adrop(true_alpha2[ii,cbns[[jj]],,drop=F], 1)))),
                       cbind(diag(1, N*(true_P-1), N*(true_P-1)), matrix(0, N*(true_P-1), N)))
    }
    if (all(sapply(1:(2^true_H-1), function(hh) Mod(eigen(G[hh,,], only.values = TRUE)$values[1]) < 0.9))) break
  }
}

## 400 points segmented into 5 intervals
Time <- 400
time_intervals <- cbind(c(true_P+1, 81, 161, 241, 321), c(80, 160, 240, 320, Time))-true_P
index <- matrix(NA, Rep, 5)
true_A <- array(NA, c(Rep, N, N*true_P, Time-true_P))
obs <- array(NA, c(Rep, Time, N))
err <- array(NA, c(Rep, Time-true_P, N))
for (ii in 1:Rep) {
  for (jj in 1:5) {
    index[ii,jj] <- sample(1:(2^true_H-1), 1, replace = TRUE, prob = rep(1/(2^true_H-1), 2^true_H-1))
    true_A[ii,,,time_intervals[jj,1]:time_intervals[jj,2]] <-
      t(adrop(true_alpha1[ii,cbns[[index[ii,jj]]],,drop=F], 1)) %*%
      t(khatri_rao(t(adrop(true_alpha3[ii,cbns[[index[ii,jj]]],,drop=F], 1)),
                   t(adrop(true_alpha2[ii,cbns[[index[ii,jj]]],,drop=F], 1))))
  }
  sd <- max(abs(true_A[ii,,,]))/SNR
  obs[ii,1:true_P,] <- rnorm(true_P*N, 0, sd)
  err[ii,,] <- rnorm((Time-true_P)*N, 0, sd)
  for (tt in (true_P+1):Time) {
    obs[ii,tt,] <- true_A[ii,,,tt-true_P] %*% c(t(obs[ii,(tt-1):(tt-true_P),])) + err[ii,tt-true_P,]
  }
}
true_A_short <- array(NA, c(Rep, N, N*true_P, 5))
for (ii in 1:Rep) {
  for (jj in 1:5) {
    true_A_short[ii,,,jj] <- true_A[ii,,,jj*80-true_P]
  }
}
true_gamma <- array(0, c(Rep, true_H, Time-true_P))
for (ii in 1:Rep) {
  for (jj in 1:5) {
    true_gamma[ii,cbns[[index[ii,jj]]],time_intervals[jj,1]:time_intervals[jj,2]] <- 1
  }
}

save(cbns, index, time_intervals, err, N, obs, Rep, s, sd, SNR, Time, true_A_short, true_gamma, true_alpha1,
     true_alpha2, true_alpha3, true_H, true_P, file = 'SIMULATION1_SNR2.RData')
