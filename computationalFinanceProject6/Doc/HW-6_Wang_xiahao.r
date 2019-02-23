GeneratePricePath <- function(nPath, nSims, S, r, sigma, delta){
  # generate two sets of standard normals that are negatively correlated
  stdNor1 <- matrix(rnorm(nPath*nSims/2), nrow = nPath/2, byrow = TRUE)
  stdNor2 <- -stdNor1
  stdNor_total <- rbind(stdNor1,stdNor2)
  
  priceProcess <- matrix(0, nPath, nSims+1)
  priceProcess[,1] <- S
  
  for(i in 1:nSims){
    priceProcess[, i+1] <- priceProcess[, i] *exp((r - 0.5 * sigma * sigma) * delta + sigma * stdNor_total[,i] * sqrt(delta))
  }
  return(priceProcess)
}

callOptionPrice <- function(nPath, nSims, S, r, sigma){
  paths <- GeneratePricePath(nPath, nSims, S, r, sigma, delta)
  result <- apply(paths, 1, max)
  payoff <- result-strike
  payoff[payoff<0]<-0
  return(mean(payoff))
}

putOptionPrice <- function(nPath, nSims, S, r, sigma){
  paths <- GeneratePricePath(nPath, nSims, S, r, sigma, delta)
  result <- apply(paths, 1, min)
  payoff <- strike - result
  payoff[payoff<0] <- 0
  return(mean(payoff))
}

T <- 1
nPath <- 1000 # row
nSims <- 100 # col
delta <- T/nSims
S <- 98
strike <- 100
r <- 0.03
sigmas <- seq(0.12,0.48, by=0.04)
callPrice <-rep(0.0, length(sigmas))
putPrice <-rep(0.0, length(sigmas))

for(i in 1:length(sigmas)){
  callPrice[i] = callOptionPrice(nPath,nSims, S, r, sigmas[i])
  putPrice[i] = putOptionPrice(nPath,nSims, S, r, sigmas[i])
}

plot(y=callPrice, x= sigmas, type="l", col="red",xlab="volatility")
lines(x= sigmas,y= putPrice, type = "l", col="blue")
abline(h= 3.72244, col="black")
legend("topleft", c("call", "put"), col=c("red", "blue"), lwd=5)




