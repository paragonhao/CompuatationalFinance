
GetU <- function(r, h){
  return(exp(r * h + 0.2 * sqrt(h)))
}
GetD <- function(r, h){
  return(exp(r * h - 0.2 * sqrt(h)))
}

GetPStar <- function(r, h, u, d ){
   return((exp(r*h) - d)/(u-d))
}

GenerateBinomialTree <- function(s, N, t, u, d){
  tree <- matrix(0, nrow = N+1, ncol = N + 1)
  
  for (i in 1:(N+1)) {
    for (j in 1:i) {
      tree[i,j] = s * u^(j-1) * d^((i-1)-(j-1))
    }
  }
  
  return(tree)
}

value_binomial_option_euro = function(tree, r, h, k, type) {
  p <- GetPStar(r, h, u, d)
  
  optionTree = matrix(0, nrow=nrow(tree), ncol=ncol(tree))
  
  if(type == 'put') {
    optionTree[nrow(optionTree),] = pmax(k - tree[nrow(tree),], 0)
  } else {
    optionTree[nrow(optionTree),] = pmax(tree[nrow(tree),] - k, 0)
  }
  
  for (i in (nrow(tree)-1):1) {
    for(j in 1:i) {
      optionTree[i, j] = ((1-p)*optionTree[i+1,j] + p*optionTree[i+1,j+1])/exp(r*h)
    }
  }
  
  
  return(optionTree)
}

value_binary_option_euro <- function(tree, r, h, k,u,d){
  p <- GetPStar(r, h, u, d)
  
  optionTree = matrix(0, nrow=nrow(tree), ncol=ncol(tree))
  
  for (i in 1:(N+1)) {
    for (j in 1:i) {
      if( tree[i, j] > k){
        optionTree[i, j] = tree[i, j] 
      }else{
        optionTree[i, j] = 0
      }
    }
  }
  
  for (i in (nrow(tree)-1):1) {
    for(j in 1:i) {
      optionTree[i, j] = ((1-p)*optionTree[i+1,j] + p*optionTree[i+1,j+1])/exp(r*h)
    }
  }
  
  return(optionTree)
}

price_replicating_portfolio <- function(priceTree, payoffTree, N, r, h){
  shares <- matrix(0, nrow = N, ncol = N)
  bond <- matrix(0, nrow = N, ncol = N)
  
  for (i in (nrow(priceTree)-1):1) {
    for(j in 1:i) {
      shares[i,j] = (payoffTree[i+1,j+1] - payoffTree[i+1,j])/(priceTree[i+1,j+1] - priceTree[i+1,j])
      bond[i,j] = (payoffTree[i+1,j+1] - priceTree[i+1,j+1] * shares[i,j]) / exp(r*h)
    }
  }
  return(list(shares,bond))
}

#For 1 a.
t <- 4
r <- 0.02
h <- 0.25
s0 <- 100
k <- 90
u <- GetU(r, h)
d <- GetD(r, h)
N <- t/h

binoTree <- GenerateBinomialTree(s0 , N, t, u, d)
callTree <- value_binomial_option_euro(binoTree, r, h, k, "call")
putTree <- value_binomial_option_euro(binoTree, r, h, k, "put")
premium_straddle1a <- callTree[1,1] + putTree[1,1] # for 1 a. 
replicatingPortfolioCall1a <- price_replicating_portfolio(binoTree, callTree, N, r, h)
replicatingPortfolioPut1a <- price_replicating_portfolio(binoTree, putTree, N, r, h)
cat("Qn 1a Straddle Payoff is:", premium_straddle1a)
cat("Replicating Portfolio for call can be found at replicatingPortfolioCall1a (too big to show)")
cat("Replicating Portfolio for put can be found at replicatingPortfolioPut1a (too big to show)")

#For 1 b
t <- 40
r <- 0.02
h <- 0.025
s0 <- 100
k <- 90
u <- GetU(r, h)
d <- GetD(r, h)
N <- t/h
binoTree2 <- GenerateBinomialTree(s0 , N, t, u, d)
callTree2 <- value_binomial_option_euro(binoTree2, r, h, k, "call")
putTree2 <- value_binomial_option_euro(binoTree2, r, h, k, "put")
premium_straddle1b <- callTree2[1,1] + putTree2[1,1]
replicatingPortfolioCall1b <- price_replicating_portfolio(binoTree2, callTree2, N, r, h)
replicatingPortfolioPut1b <- price_replicating_portfolio(binoTree2, putTree2, N, r, h)
cat("Qn 1b Straddle Payoff is:", premium_straddle1b)
cat("Replicating Portfolio for call can be found at replicatingPortfolioCall1b (too big to show)")
cat("Replicating Portfolio for put can be found at replicatingPortfolioPut1b (too big to show)")

#For 1 c
#Assuming Asset or nothing binary options

t <- 4
r <- 0.02
h <- 0.25
s0 <- 100
k <- 90
u <- GetU(r, h)
d <- GetD(r, h)
N <- t/h
binoTree1c <- GenerateBinomialTree(s0 , N, t, u, d)
binary1c <- value_binary_option_euro(binoTree1c, r, h, k, u,d)
replicatingPortfolioCall1c <- price_replicating_portfolio(binoTree1c, binary1c, N, r, h)
cat("Qn 1c Straddle Payoff is:", binary1c[1, 1])
cat("Replicating Portfolio for call can be found at replicatingPortfolioCall1c (too big to show)")
cat("Replicating Portfolio for put can be found at replicatingPortfolioCall1c (too big to show)")

