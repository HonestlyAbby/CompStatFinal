#Abby Hatch
#CompStat Finals

set.seed(1234567890)

library(SuppDists)
library(mvtnorm)
library(testthat)

#Question 1

time_mat <- matrix(data = NA, nrow=5, ncol = 20)
row.names(time_mat) <- c("r_squared_fun_1", "r_squared_fun_2", "r_squared_fun_3", "r_squared_fun_4", "r_squared_fun_5")

i <- 1
while(TRUE){
  x <- matrix(rnorm(5e5), nrow = 1000, ncol = 500)
  y <- c(rnorm(1000))

  Rprof(filename = "final.out", interval= .001)
  r_squared_fun_1(x, y)
  r_squared_fun_2(x, y)
  r_squared_fun_3(x, y)
  r_squared_fun_4(x, y)
  r_squared_fun_5(x, y)
  Rprof(NULL)
  
  times <- summaryRprof(filename = "final.out")$by.total
  time_mat[1, i] <- times["\"r_squared_fun_1\"",]$total.time
  time_mat[2, i] <- times["\"r_squared_fun_2\"",]$total.time
  time_mat[3, i] <- times["\"r_squared_fun_3\"",]$total.time
  time_mat[4, i] <- times["\"r_squared_fun_4\"",]$total.time
  time_mat[5, i] <- times["\"r_squared_fun_5\"",]$total.time
  
  i <- i+1
  if(i > 20){
    break
  }
}

summary_mat <- matrix(data = NA, nrow = 5, ncol = 6)
row.names(summary_mat) <- c("r_squared_fun_1", "r_squared_fun_2", "r_squared_fun_3", "r_squared_fun_4", "r_squared_fun_5")
colnames(summary_mat) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max")

for(i in 1:5){
  for(j in 1:6){
    summary_mat[i, j] <- summary(time_mat[i,])[j]
  }
}

summary(time_mat[[1, ]])
summary_mat

#The second function is the fastest. My guess is that this is due to the inversion. 
#Inverting such a large matrix is probably computationally difficult, and the second function requires no inversion

#Problem 2

draw_tau = function(data,state){
  shape = 1 + length(data$y)
  rate = 1 + sum(abs(data$y-data$x%*%state$beta))
  return(rgamma(1,shape,rate))
}

draw_beta <- function(data, state){
  lambda <- 1
  nu <- 1/(state$tau*abs(data$y-data$x%*%state$beta))
  
  zvec <- rinvGauss(n = nrow(data$x), nu = nu, lambda = lambda)
  z <- diag(zvec)
  I <- diag(1, nrow = ncol(data$x))
  
  mean <- solve((state$tau^2*t(data$x)%*%z%*%data$x+I)) %*% (state$tau* t(data$x)%*%z%*%as.matrix(data$y))
  sigma <- solve((state$tau^2*t(data$x)%*%z%*%data$x+I))
  
  beta <- rmvnorm(1, mean = mean, sigma=sigma)
  return(t(beta))
}

update_current_state = function(data, state){
  state$beta = draw_beta(data,state)
  state$tau = draw_tau(data,state)
  return(state)
}

gibbs_laplace_regression = function(N_draws,x,y,beta_0,tau_0){
  # data manipulation
  if(is.null(colnames(x))){
    colnames(x) = paste("x",1:ncol(x),sep="_")
  }
  if(!all(abs(x[,1]-1)<.Machine$double.eps)){
    x = cbind(c(1),x)
    colnames(x)[1] = "intercept"
  }
  data = list(y=y,x=x)
  
  # checking for missing initial values and making them if necessary
  if(missing(beta_0)) beta_0 = c(median(data$y),rep(0,ncol(data$x)-1))
  if(missing(tau_0)) tau_0 = 1/mean(abs(data$y-data$x%*%beta_0))
  
  # initial state creation
  current_state = list(beta=beta_0,tau=tau_0)
  
  # output creation
  out = list(beta=matrix(NA,ncol(x),N_draws),tau=rep(NA,N_draws))
  rownames(out$beta) = colnames(x)
  
  # loop for sampling
  for(i in 1:N_draws){
    current_state = update_current_state(data,current_state)
    out$beta[,i] = current_state$beta
    out$tau[i] = current_state$tau
  }
  
  # returning
  return(out)
}

posteriors <- gibbs_laplace_regression(10000, x, y)
betas <- posteriors$beta

#The mean of the betas is 68.87
#The mean of the taus is .0062

apply(betas, 1, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

#Betas 3,4,7,8,9, and 10 are all definitely non-zero

#Problem 3

eigen_fun = function(X,tol=.Machine$double.eps^0.75){
  if(!isSymmetric(X)) stop("x is not symmetric")
  d = ncol(X)
  U = diag(1,d,d,FALSE)
  while(TRUE){
    did_we_update = FALSE
    for(i in 1:(d-1)){
      for(j in (i+1):d){
        # test whether abs(X[i,j]) is too large
        if(abs(X[i,j])>tol){
          # if so, get the Givens Rotation matrix R
          # and update X and U
          
          if((abs(X[i,i] - X[j,j])) < .Machine$double.eps){
            theta <- pi/4
          } else {
            body <- (2*X[j,i])/(X[i,i] - X[j,j])
            theta <- .5*atan(body)
          }
          
          c <- cos(theta)
          s <- sin(theta)
          
          browser()
          
          R <- diag(1, d,d)
          R[i,i] <- c
          R[i,j] <- (s)
          R[j,i] <- (-s)
          R[j,j] <- c
          
          X <- R %*% X %*% t(R)
          
          #X[i,i] <- x_new[i,i]
          #X[i,j] <- x_new[i,j]
          #X[j,i] <- x_new[j,i]
          #X[j,j] <- x_new[j,j]
          
          U <- U %*% R
          # track that in this sweep through we have done at least one update
          did_we_update = TRUE
        }
        
      }
    }
    if(!did_we_update) break
  }
  values=diag(X)
  sorted_values = sort(values,index.return=TRUE,decreasing=TRUE)
  values = sorted_values$x
  vectors = U[,sorted_values$ix]
  return(list(values=values, vectors=vectors))
}

d <- 3
z = matrix(rnorm(d^d),d,d);
X = z + t(z)

eigen(X)
eigen_fun(X)

test_that('same eigen', {
  set.seed(1234567890)
  Z_2 <- matrix(rnorm(2^2),2,2)
  X_2 <- Z_2 + t(Z_2)
  Z_3 <- matrix(rnorm(3^3),3,3)
  X_3 <- Z_3 + t(Z_3)
  Z_4 <- matrix(rnorm(4^4),4,4)
  X_4 <- Z_4 + t(Z_4)
  Z_5 <- matrix(rnorm(5^5),5,5)
  X_5 <- Z_5 + t(Z_5)
  
  expect_equal(eigen(X_2)$values, eigen_fun(X_2)$values)
  expect_equal((as.vector(abs(eigen(X_2)$vectors))), as.vector(abs(eigen_fun(X_2)$vectors)))
  
  expect_equal(eigen(X_3)$values, eigen_fun(X_3)$values)
  expect_equal((as.vector(abs(eigen(X_3)$vectors))), as.vector(abs(eigen_fun(X_3)$vectors)))
  
  expect_equal(eigen(X_4)$values, eigen_fun(X_4)$values)
  expect_equal((as.vector(abs(eigen(X_4)$vectors))), as.vector(abs(eigen_fun(X_4)$vectors)))
  
  expect_equal(eigen(X_5)$values, eigen_fun(X_5)$values)
  expect_equal((as.vector(abs(eigen(X_5)$vectors))), as.vector(abs(eigen_fun(X_5)$vectors)))
})
