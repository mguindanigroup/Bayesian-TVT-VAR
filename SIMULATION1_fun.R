library(gtools)

### FUNCTION TO GENERATE RANDOM COEFFICIENTS FOR EACH COMPONENT
random_coefficient <- function(n, sparcity, sd) {
  x <- sample(c(0,1), n, prob = c(sparcity, 1-sparcity), replace = TRUE)
  return(x*rnorm(n,0,sd))
}

### FUNCTION TO RETURN ALL COMBINATORIALS OF SET {1,\dots,x-1}
cbnt <- function(x) {
  cbnts <- list()
  for (y in 1:(x-1)) {
    cb <- combinations(x-1,y)
    cb <- lapply(seq_len(nrow(cb)), function(z) cb[z,])
    cbnts <- append(cbnts, cb)
  }
  return(cbnts)
}

### FUNCTION TO SIMULATE BINARY INDICATORS
NDARMA <- function(n, p1, p2) {
  bts <- rep(sample(c(0,1), 1, prob = c(1-p2,p2)), n)
  bts_switch <- sample(c(0,1), n-1, prob = c(1-p1,p1), replace = TRUE)
  for (i in 2:n) {
    bts[i] <- bts_switch[i-1]*bts[i-1]+(1-bts_switch[i-1])*sample(c(0,1), 1, prob = c(1-p2,p2))
  }
  return(bts)
}