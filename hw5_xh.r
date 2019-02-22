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
  if ( (nSims - currCol) == 1){
    y <- cashflow[, (currCol+1):nSims] * discount_factor * isBelowStrike
  }else{
    y <- cashflow[, (currCol+1):nSims] %*% discount_factor * isBelowStrike
  }
}

CalculateMatrixA <- function(k, METHOD, x_nonzeros){
  matA <- matrix(0, k,k)
  for (i in 1:k){
    for (j in 1:k){
      matA[i,j] <- sum(METHOD(x_nonzeros, i) * METHOD(x_nonzeros, j))
    }
  }
  return(matA)
}

calculateMatrixB <- function(k, METHOD,x_nonzeros, y_nonzeros){
  vecB <- matrix(0, k,1)
  for (i in 1:k){
    vecB[i,1] <- sum(METHOD(x_nonzeros,i) * y_nonzeros)
  }
  return(vecB)
}

RemoveZeros <- function(isBelowStrike, nPath, variable){
  key <- seq(1, nPath, 1) * isBelowStrike
  return(variable[key])
}

GetRegressionX <- function(x, METHOD){
  LfuncValues <- c()
  for (i in 1:k){
    LfuncValues <- cbind(LfuncValues, METHOD(x, i))
  }
  return(LfuncValues)
}

GetFinalPayOff <- function(currCol, priceProcess, currentPayoff, cashflow, ECV){
  payoff <- strike - priceProcess[,currCol+1]
  payoff[payoff<0] <- 0
  currentPayoff <- cashflow[,currCol] + payoff
  
  for (j in 1: nPath){
    if (currentPayoff[j] > ECV[j]){
      cashflow[j, (currCol+1):nSims] <- 0
      cashflow[j, (currCol-1):1] <- 0
      cashflow[j, currCol] <- currentPayoff[j]
    }
  }
}

################################# LSMC General function for this project #################################
LSMC <- function (delta, nSims, t, strikePrice, S, r, sigma, nPath, k, METHOD, isForwardStart, stopping_time){
  
  # Initialize price process with antithetic variance reduction technique
  priceProcess <- GeneratePricePath(nPath, nSims, S, r, sigma, delta)
  
  # Take the mean of the strike price at small t 
  if(isForwardStart){
    strike <- priceProcess[, stopping_time + 1]
  }
  
  # intialize cashflow matrix
  cashflow <- InitializeCashFlow(nPath, strike, nSims, priceProcess)
  
  priceProcess[, 1] <- S
  
  for (i in (nSims-1) : stopping_time){
    # find out positions where it is below strike
    
    # TODO: trying replacing the strike with if else to differentiate between european
    # can american
    isBelowStrike <- (strike - priceProcess[,i+1] > 0)
    
    discount_factor <- exp(seq(-r*delta,-(nSims-i)*r*delta, -r*delta))
    
    #Get X and discounted Y accordingly
    x <- GetX(priceProcess, i+1 , isBelowStrike)
    y <- GetY(i, cashflow, nSims, discount_factor , isBelowStrike)
    
    # Get rid of zeros in y and x for linear regression analysis
    y_nonzeros <- RemoveZeros(isBelowStrike, nPath, y)
    x_nonzeros <- RemoveZeros(isBelowStrike, nPath, x)
    
    if(length(x_nonzeros) < 2) break;
    
    matA <- CalculateMatrixA(k, METHOD, x_nonzeros)
    vecB <- calculateMatrixB(k, METHOD, x_nonzeros, y_nonzeros)
    
    LfuncValues <- GetRegressionX(x, METHOD)
    
    if(det(matA) > 0.00000000001 && length(matA) >=2){
      # use cholesky decomposition to find coef 
      coef <- chol2inv(chol(matA))%*%vecB
      
      # to find the expected continuation value on path that are not zeros
      ECV <- LfuncValues %*% coef * isBelowStrike
      
      # initialize payoff to current pay on the column/path
      payoff <- strike - priceProcess[,i+1]
      payoff[payoff<0] <- 0
      currentPayoff <- cashflow[,i] + payoff
      # Update cashflow matrix with current payoff
      for (j in 1: nPath){
        if (currentPayoff[j] > ECV[j]){
          cashflow[j, (i+1):nSims] <- 0
          cashflow[j, (i-1):1] <- 0
          cashflow[j,i] <- currentPayoff[j]
        }
      }
    }
  }
  discountFactorForAllPeriods <- exp(seq(-r*delta,-nSims*r*delta, -r*delta))
  finalPayOff <-mean(cashflow%*% discountFactorForAllPeriods)
  return(finalPayOff)
}

# (a)
r <- 0.06
sigma <- 0.2
strike <- 40
nPath <- 100000
stopping_time <- 1
isForwardStart <- FALSE
technique <- Monomial # three methods: Lagurre, Hermit and Monomial, replace with method accordingly
#################### Please uncomment to run the entire scripts #################### 
################################# t= 0.05 ##########################################
# k=2
nSims <- 100
t <- 0.5
k <- 2
delta <- t/nSims
S <-36
payoffK2S36t05 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)
S <-40
payoffK2S40t05 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)
S <-44
payoffK2S44t05 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)

k <- 3
S <-36
payoffK3S36t05 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)
S <-40
payoffK3S40t05 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)
S <-44
payoffK3S44t05 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)

k <- 4
S <-36
payoffK4S36t05 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)
S <-40
payoffK4S40t05 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)
S <-44
payoffK4S44t05 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)

################################# t= 1 ##########################################
nSims <- 200
t <- 1
k <- 2
delta <- t/nSims
S <- 36
payoffK2S36t1 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)
S <-40
payoffK2S40t1 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)
S <-44
payoffK2S44t1 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)

k <- 3
S <-36
payoffK3S36t1 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)
S <-40
payoffK3S40t1 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)
S <-44
payoffK3S44t1 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)

k <- 4
S <-36
payoffK4S36t1 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)
S <-40
payoffK4S40t1 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)
S <-44
payoffK4S44t1 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)

################################# t= 2 ##########################################
nSims <- 400
t <- 2
k <- 2
delta <- t/nSims
S <- 36
payoffK2S36t2 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)
S <-40
payoffK2S40t2 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)
S <-44
payoffK2S44t2 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)

k <- 3
S <-36
payoffK3S36t2 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)
S <-40
payoffK3S40t2 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)
S <-44
payoffK3S44t2 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)

k <- 4
S <-36
payoffK4S36t2 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)
S <-40
payoffK4S40t2 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)
S <-44
payoffK4S44t2 <- LSMC(delta, nSims, t, strike, S, r, sigma, nPath, k, technique, isForwardStart, stopping_time)

# (b)

################################# End of Qn 1  #################################

################################# Qn2  #################################

# 2(a)

ForwardStartOpt <- function(Terminal, t, S, strike, sigma, r){
  delta <- (Terminal / nSims)
  priceProcess <- GeneratePricePath(nPath, nSims, S, r, sigma, delta)
  # intialize cashflow matrix
  strikeCol <- nSims * t + 1
  strikeAtt <- priceProcess[,strikeCol]
  cashflow <- strikeAtt  - priceProcess[, nSims + 1]
  cashflow[cashflow<0] <- 0 
  payoff <- mean(exp(-r*Terminal) * cashflow)
  return(payoff)
}

nPath <- 100000
nSims <- 100
Terminal <- 1
t <- 0.2
S <- 65
strike <- 60
sigma <- 0.2
r <- 0.06
payoff2a <- ForwardStartOpt(Terminal, t, S, strike, sigma, r)

#2(b)


stopping_time <- nSims * t 
isForwardStart <- TRUE
delta <- Terminal/ nSims
k <- 3
payoff2b <- LSMC(delta, nSims, Terminal, strike, S, r, sigma, nPath, k, Hermit, isForwardStart, stopping_time)

