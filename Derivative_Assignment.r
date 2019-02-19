## This project is done collectively by Group 6:
## Group Member: Students: Xiahao Wang, Juan Manuel Ferreyra Maspero, Xinyue Zhu, Yichu Li, Mu Lin

############################ Question 1 ##################################
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
N <- t
binoTree1c <- GenerateBinomialTree(s0 , N, t, u, d)
binary1c <- value_binary_option_euro(binoTree1c, r, h, k, u,d)
replicatingPortfolioCall1c <- price_replicating_portfolio(binoTree1c, binary1c, N, r, h)
cat("Qn 1c Straddle Payoff is:", binary1c[1, 1])
cat("Replicating Portfolio for call can be found at replicatingPortfolioCall1c (too big to show)")
cat("Replicating Portfolio for put can be found at replicatingPortfolioCall1c (too big to show)")
############################ end of question 1 ##################################


############################## QN2 ##############################

# To price an American call/put option
# function to price american put
price_amerput<-function(S_0,K,r,Time,n,return_tree=FALSE){
  dt <- Time/n
  p <- 1/2
  u <- exp(r*dt+0.15*sqrt(dt))
  d <- exp(r*dt-0.15*sqrt(dt))
  
  S_t <- matrix(ncol=n+1,nrow=n+1)
  S_t[1,1] <- S_0
  S_t[n+1,] <- S_0*d^n * cumprod(c(1,rep(u^2,n))) # stock price at t
  
  V_t <- matrix(ncol=n+1,nrow=n+1) # payoff at t
  V_t[n+1,] <- pmax(K - S_t[n+1,],0)
  
  for(i in n:1){
    S_t[i,1:i] <- S_t[i+1,1:i]*u
    V_t[i,1:i] <- pmax(exp(-r*dt)*(p*V_t[i+1,2:(i+1)]+(1-p)*V_t[i+1,-((i+1):(n+1))]),K-S_t[i,1:i])
  }
  
  if(return_tree){
    return(V_t)
  }else{
    return(V_t[1,1])
  }
}

# function to price american call
price_amercall<-function(S_0,K,r,Time,n,return_tree=FALSE){
  dt <- Time/n
  p <- 1/2
  u <- exp(r*dt+0.15*sqrt(dt))
  d <- exp(r*dt-0.15*sqrt(dt))
  
  S_t <- matrix(ncol=n+1,nrow=n+1)
  S_t[1,1] <- S_0
  S_t[n+1,] <- S_0*d^n * cumprod(c(1,rep(u^2,n))) # stock price at t
  
  V_t <- matrix(ncol=n+1,nrow=n+1) # payoff at t
  V_t[n+1,] <- pmax(S_t[n+1,]-K,0)
  
  for(i in n:1){
    S_t[i,1:i] <- S_t[i+1,1:i]*u
    V_t[i,1:i] <- pmax(exp(-r*dt)*(p*V_t[i+1,2:(i+1)]+(1-p)*V_t[i+1,-((i+1):(n+1))]),K-S_t[i,1:i])
  }
  
  if(return_tree){
    return(V_t)
  }else{
    return(V_t[1,1])
  }
}

# To run for an American Option with K=10
put_american=price_amerput(10,10,0.01,1,100)
call_american=price_amercall(10,10,0.01,1,100)
############################ end of question 2 ##################################

############################## QN3 ##############################
r <- 0.02
sig <- 0.2
h <- 1/365
s0 <- 10
t <- 200 # periods
delta <- 0.05
u <- exp(sig*h)
d <- 1/u
p <- (exp(r*h)-d)/(u-d)
D <- c(50,100,150) # div dates

eu_option <- function(s0,r,div,divdates,u,d,h,t,K,type){
  stock   <- matrix(0L, nrow = t+1, ncol = t+1)
  s[1,1]  <- s0
  ndates  <- length(divdates)
  for (j in 2:t+1){
    for (i in 1:j){
      if (div == 0 | j < divdates(1)+1){
        s[i,j] <- s0*u^(j-i)*d^(i-1)
      }
      else if (j < divdates(ndates)+1){
        for (g in 2:ndates){
          if (j >= divdates(g-1)+1 && j < divdates(g)+1){
            s[i,j] <- s0*u^(j-i)*d^(i-1)*(1-div)^(g-1)
          }
        }
      }
      else{
        s[i,j] <- s0*u^(j-i)*d^(i-1)*(1-div)^(ndates)
      }
    }
  }
  
  option <- matrix(0L, nrow = t+1, ncol = t+1)
  for (k in 1:t+1){
    if (type == "C"){
      option[k,t+1] <- max(s[k,t+1]-K,0)
    }
    else if (type == "P"){
      option[k,t+1] <- max(K-s[k,t+1],0)
    }
  }
  
  stock_amt <- matrix(0L, nrow = t, ncol = t)
  bond_amt <- matrix(0L, nrow = t, ncol = t)
  for (j in t:1){
    for (i in i:j){
      option[i,j] <- exp(-r*h)*(p*option[i,j+1]+(1-p)*option[i+1,j+1]);
      stock_amt[i,j] <- (option[i,j+1]-option[i+1,j+1])/(s[i,j]*(u-d));
      bond_amt[i,j] <- option[i,j]-stock_amt[i,j]*s[i,j];
    }
  }
}

dividend <- function(s0,r,div,divdates,u,d,h,t,K,type,nature){
  if (nature=="A"){
    return(eu_option(s0,r,div,divdates,u,d,h,t,K,type))
  }
  else if (nature=="E"){
    return(eu_option(s0,r,div,divdates,u,d,h,t,K,type))
  }
}
############################ end of question 3 ##################################

############################## QN4 ##############################
call_price <- function(N,t,r_f,s0,sigma,n_path,k){
  h=t/N
  ds <- matrix(0,nrow=N,ncol=n_path)
  s <- matrix(s0,nrow=N+1,ncol=n_path)
  c <- vector()
  for (j in 1:n_path){
    .Random.seed
    z=rnorm(n=N)
    for (i in 1:N){
      ds[i,j]=r_f*s[i,j]*h+sigma*s[i,j]*sqrt(h)*z[i]
      s[i+1,j]=s[i,j]+ds[i,j]
    }
    c[j] <- max((mean(s[,j])-k),0)
  }
  
  return(c)
}

s0 <- 200
r <- 0.02
sigma <- 0.2
k <- 220
t=1
n <- 365
n_path <- 10000

call_output <- call_price(N=n,t=t,r_f=r,s0=s0,sigma=sigma,n_path=n_path,k=k)

c_price <- sum(call_output)/n_path
ci_lower <- mean(call_output)-1.96*sd(call_output)
ci_higher <- mean(call_output)+1.96*sd(call_output)
c_price
ci_lower
ci_higher
############################ end of question 4 ##################################

############################ Question 5 ##################################
###FUNCTIONS FOR THIS PROJECT
herm <- function(x,k){
  
  if (k == 1){
    return(rep(1,length(x)))
  }
  if (k == 2){
    return(2*x)
  }
  if (k == 3){
    return(4*x**2-2)
  }
  if (k == 4){
    return(8*x**3-12*x)
  }
}
laguerre <- function(x,k){
  
  if (k == 1){
    return(exp(-x/2))
  }
  if (k == 2){
    return(exp(-x/2)*(1-x))
  }
  if (k == 3){
    return(exp(-x/2)*(1-x*2+x**2/2))
  }
  if (k == 4){
    return(exp(-x/2)*(1-3*x+3*x**2/2-x**3/6))
  }
}
mono <- function(x,k){
  
  if (k == 1){
    return(rep(1,length(x)))
  }
  if (k == 2){
    return(x)
  }
  if (k == 3){
    return(x**2)
  }
  if (k == 4){
    return(x**3)
  }
}

create_m_price <- function(m,n,So,r,sigma,t_path){
  wt <- matrix(rnorm(n*m/2),nrow=m/2,byrow = T)*sqrt(t_path[1])
  wt_path <- t(apply(wt,1,cumsum))
  zt_path <- -wt_path
  matrix1 <- t(apply(wt_path,1,function(x){(So*exp((r-sigma**2/2)*t_path + sigma*x))}))
  matrix2 <- t(apply(zt_path,1,function(x){(So*exp((r-sigma**2/2)*t_path + sigma*x))}))
  s_matrix <- rbind(matrix1,matrix2)
  return(s_matrix)
}

L_selector <- function(method = 'laguerre',x,k){
  if (method == 'laguerre'){
    return(laguerre(x,k))
  }
  if (method == 'hermite'){
    return(herm(x,k))
  }
  if (method == 'monomial'){
    return(mono(x,k))
  }
}
#MAIN FUNCTION TO GENERATE THE LEAST SQUARE APPROXIMATION
LSMC <- function(m,n){
  #The m has to be even in order to continue because I need half rows plus, and half antithetic
  X <- 220
  r <- 0.1
  sigma <- 0.3
  So <- 200
  t <- 1
  k <- 3
  method <- 'laguerre'
  start = 0
  X2 <- X
  
  So <- So/X
  X <- 1
  
  if(m%%2==1){
    stop('m has to be even, please try the function again with an even m')
  }
  #preparing matrices and information for the loop
  delta <- t/n
  t_path <- seq(delta,t,delta)
  set.seed(round(runif(1,1,100000),0))
  s_matrix <- create_m_price(m,n,So,r,sigma,t_path)
  ind_matrix <- matrix(rep(0,m*n),nrow=m,ncol=n)
  ind_matrix[,n] <- as.numeric(s_matrix[,n]<X)
  discount <- exp(-r*delta*seq(1,n))
  a_matrix <- matrix(ncol=k,nrow=k)
  b <- matrix(ncol=1,nrow=k)
  
  if (start != '0' & start > 0){
    start <- (start/t)*n
  } else {
    start <- 1 # it will be executable from 1st column
  }
  
  #start the loop by col from the botton to the beginning
  for (col in (n-1):start){
    
    row_vector <- which((X-s_matrix[,col])>=0)#Case there is no meaning to execute, just skip the row
    
    if (n-col==1){
      y <- ((X-s_matrix[row_vector,(col+1):n])*ind_matrix[row_vector,(col+1):n])*discount[1:(n-col)]
      y <- matrix(y,ncol=1)
    } else {
      y <- ((X-s_matrix[row_vector,(col+1):n])*ind_matrix[row_vector,(col+1):n])%*%discount[1:(n-col)]
    }
    y_max <- apply(y,1,max)
    
    for (i in 1:k){
      for (j in 1:k){
        f <- L_selector(method,s_matrix[row_vector,col],i)*L_selector(method,s_matrix[row_vector,col],j)
        a_matrix[i,j] <- sum(f)
      }
      
      b[i,1] <- sum(y_max*L_selector(method,s_matrix[row_vector,col],i))
      
    }
    a_vector <- chol2inv(chol(a_matrix))%*%b #inverting the matrix with chelosky roots and finding constants from regression
    
    #Start row iteration in order to check continous value vs executable value
    for (row in row_vector){
      
      exp_c_value <- 0
      
      if (k == 2){
        exp_c_value <- sum(a_vector*
                             c(L_selector(method,s_matrix[row,col],1),L_selector(method,s_matrix[row,col],2)))
      }
      if (k == 3){
        exp_c_value <- sum(a_vector*
                             c(L_selector(method,s_matrix[row,col],1),L_selector(method,s_matrix[row,col],2),L_selector(method,s_matrix[row,col],3)))
      }
      if (k == 4){
        exp_c_value <- sum(a_vector*
                             c(L_selector(method,s_matrix[row,col],1),L_selector(method,s_matrix[row,col],2),L_selector(method,s_matrix[row,col],3),L_selector(method,s_matrix[row,col],4)))
      }
      
      if ((X-s_matrix[row,col])>exp_c_value){ #Case the exp value of executing is more than the exp value of continuation
        ind_matrix[row,(col+1):n] <- 0
        ind_matrix[row,col] <- 1
      }
    }
    print(paste('Progress...',round((n-col)/n*100,0),'%'))
  }
  est_lsmc_pp <- sum(((X-s_matrix)*ind_matrix)%*%matrix(discount,ncol=1))/m
  return(est_lsmc_pp*X2)
}

#Initialization
m <- 100000
n <- 250


#a)
option_price <- LSMC(m = m,n = n)

#b)
vector_price <- c(sapply(c(10,100,1000,10000),function(x){LSMC(x,n)}),option_price)

x_ax <- c(10,100,1000,10000,100000)
plot(y=vector_price, x=c(1:5), pch=19,xaxt="n")+
  title('Put price in function of row Qty')+
  axis(1,at= 1:5,labels = x_ax)

#c)
x2_ax <- c(3,10,100,250,1000)
vector_price_2 <- sapply(x2_ax,function(x){LSMC(m,x)})

plot(y=vector_price_2, x=c(1:5), pch=19,xaxt="n")+
  title('Put price in function of row Qty')+
  axis(1,at= 1:5,labels = x2_ax)

############################ end of question 5 ##################################



