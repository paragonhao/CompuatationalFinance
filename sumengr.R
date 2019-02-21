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

LSMC <- function (nSims, t, K, S, r, sigma, nPath, k, FUN){
  priceProcess = matrix(rep(0, len = nPath *(nSims+1)), nrow = nPath)
  rand = matrix(rnorm(nPath*nSims), nrow = nPath)
  CF = matrix(rep(0, len = nPath *(nSims)), nrow = nPath)
  dt = t/nSims
  priceProcess[, 1] = S
  zeros = rep(0 , len = nPath)
  zeros_row = rep(0, len = nSims)
  for (i in 1:nSims){
    priceProcess[, i+1] = priceProcess[, i] *exp((r-sigma^2/2)*dt + sigma*rand[,i]*sqrt(dt))
  }
  CF[, nSims] = apply(cbind(K - priceProcess[,nSims+1], zeros), 1, max)
  for (i in (nSims-1) : 1){
    temp = (K - priceProcess[,i+1] > 0)
    x = priceProcess[, i+1] * temp
    dis_factor = exp(seq(-r*dt,-(nSims-i)*r*dt, -r*dt))
    if (i == nSims-1){
      y = CF[, (i+1):nSims] * dis_factor * temp
    }else{
      y = CF[, (i+1):nSims] %*% dis_factor * temp
    }
    pos = seq(1, nPath, 1) * temp
    y_clean = y[pos]
    x_clean = x[pos]
    A = matrix(rep(0, len = k*k), nrow = k)
    B = matrix(rep(0, len =k), nrow = k)
    X = c()
    for (a in 1:k){
      for (b in 1:k){
        A[a,b] = sum(FUN(x_clean, a) * FUN(x_clean, b))
      }
      B[a,1] = sum(FUN(x_clean,a) * y_clean)
      X = cbind(X, FUN(x, a))
    }
    reg = solve(A,B)
    continuation = X %*% reg * temp
    newCF = apply(cbind(K - priceProcess[,i+1], zeros), 1, max)
    
    for (j in 1: nPath){
      if (newCF[j] > continuation[j]){
        CF[j, ] = zeros_row
        CF[j,i] = newCF[j]
      }
    }
  }
  D = exp(seq(-r*dt,-nSims*r*dt, -r*dt))
  return(mean(CF%*% D))
}




nSims = 250
S = 40
r = 0.06
sigma = 0.2
K = 40
t = 0.5
nPath = 1000
k = 3
out = LSMC(nSims, t, K, S, r, sigma, nPath, k, Monomial)
#out2 = LSMC(nSims, t, K, S, r, sigma, nPath, k, Hermit)
#out3 = LSMC(nSims, t, K, S, r, sigma, nPath, k, Lagur)

