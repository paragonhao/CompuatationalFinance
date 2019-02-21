Monomial <- function (x, k){
  if (k == 1)
    return(x^0)
  if (k == 2)
    return(x)
  if (k == 3)
    return(x^2)
  if (k == 4)
    return(x^3)
}

Hermit <- function (x, k){
  if (k == 1)
    return(x^0)
  if (k == 2)
    return(2*x)
  if (k == 3)
    return(4*x^2-2)
  if (k == 4)
    return(8*x^3-12*x)
}

Lagurre <- function (x, k){
  if (k == 1)
    return(exp(-x/2))
  if (k == 2)
    return(exp(-x/2)*(1-x))
  if (k == 3)
    return(exp(-x/2)*(1- 2*x+x^2/2))
  if (k == 4)
    return(exp(-x/2)*(1 - 3*x + 3*x^2/2 -x^3/6))
}

InitializeCashFlow <- function(nPath, strike, nSims, priceProcess){
  cashflow <- matrix(0, nPath, nSims)
  isBelowStrike <- strike - priceProcess[,nSims+1]
  isBelowStrike[isBelowStrike<0] <- 0
  cashflow[,nSims] <- cashflow[,nSims] + isBelowStrike
  return(cashflow)
}

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

GetX <- function(priceProcess, col, isBelowStrike){
  x <- priceProcess[, col] * isBelowStrike
  return(x)
}

GetY <- function(currCol, cashflow, nSims, discount_factor , isBelowStrike){
  if (currCol == nSims-1){
    y = cashflow[, (currCol+1):nSims] * discount_factor * isBelowStrike
  }else{
    y = cashflow[, (currCol+1):nSims] %*% discount_factor * isBelowStrike
  }
}

RemoveZeros <- function(isBelowStrike, nPath){
  
}


LSMC <- function (delta, nSims, t, strike, S, r, sigma, nPath, k, METHOD){
  
  # Initialize price process with antithetic variance reduction technique
  priceProcess <- GeneratePricePath(nPath, nSims, S, r, sigma, delta)
  
  # intialize cashflow matrix
  cashflow <- InitializeCashFlow(nPath, strike, nSims, priceProcess)
  
  priceProcess[, 1] = S
  
  zeros = rep(0 , len = nPath)
  zeros_row = rep(0, len = nSims)
  
  for (i in (nSims-1) : 1){
    isBelowStrike = (strike - priceProcess[,i+1] > 0)
    
    discount_factor = exp(seq(-r*delta,-(nSims-i)*r*delta, -r*delta))
    x <- GetX(priceProcess, i+1 , isBelowStrike)
    y <- GetY(i, cashflow, nSims, discount_factor , isBelowStrike)
    
    pos = seq(1, nPath, 1) * isBelowStrike
    
    
    y_clean = y[pos]
    x_clean = x[pos]
    A = matrix(0, k,k)
    B = matrix(0, k)
    X = c()
    for (a in 1:k){
      for (b in 1:k){
        A[a,b] = sum(METHOD(x_clean, a) * METHOD(x_clean, b))
      }
      B[a,1] = sum(METHOD(x_clean,a) * y_clean)
      X = cbind(X, METHOD(x, a))
    }
    reg = solve(A,B)
    continuation = X %*% reg * isBelowStrike
    newCF = apply(cbind(strike - priceProcess[,i+1], zeros), 1, max)
    
    for (j in 1: nPath){
      if (newCF[j] > continuation[j]){
        cashflow[j, ] = zeros_row
        cashflow[j,i] = newCF[j]
      }
    }
  }
  D = exp(seq(-r*delta,-nSims*r*delta, -r*delta))
  return(mean(cashflow%*% D))
}




nSims = 250
S = 40
r = 0.06
sigma = 0.2
strike = 40
t = 0.5
nPath = 1000
k = 3
delta <- t/nSims
out = LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, Monomial)
out
#out2 = LSMC(nSims, t, strike, S, r, sigma, nPath, k, Hermit)
#out3 = LSMC(nSims, t, strike, S, r, sigma, nPath, k, Lagur)

