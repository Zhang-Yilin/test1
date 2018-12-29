#3
get_l <- function(x, theta, i){
  x1 <- x
  theta1 <- theta
  ret <- exp(x1[i, ]%*%theta1)/(1 + exp(x1[i, ]%*%theta1))
  return(ret)
}




newton <- function(x, y, theta, eps = 0.001, n = 100){
  x1 <- x
  theta1 <- theta
  y1 <- y
  i <- 0
  plottheta <<- list()
  while(i < n){
    Wd <- sapply(1:nrow(x1), function(i){get_l(x1, theta1, i)*(1 - get_l(x1, theta1, i))})
    W <- diag(Wd)
    Pi <- sapply(1:nrow(x1), function(i){get_l(x1, theta1, i)})
    Pi <- matrix(Pi, nrow = length(Pi))
    thetaf <- theta1 + solve(t(x1)%*%W%*%x1) %*% t(x1) %*% (y1 - Pi)
    
    if(sqrt(sum((thetaf - theta1)^2)) < eps){
      print("finished")
      break
    }else{
      i <- i+1
      theta1 <- thetaf
    }
    plottheta[[i]] <<- as.numeric(thetaf)
  }
  print(i)
  if(i == n){
    print("not convergence after n loops")
  }
  return(thetaf)
}







p <- 5
n <- 200 
set.seed(1) 
X <- matrix(rnorm(n*p),ncol=p) 
X <- cbind(rep(1,n),X) 
beta <- matrix(rnorm(p+1),ncol=1) 
eta <- X%*%beta 
Prob <- exp(eta)/(1+exp(eta)) 
y <- as.numeric(runif(n) <= Prob)

theta <- matrix(0, nrow = 6, ncol = 1)

newton(x = X, y = y, theta = theta, eps = 0.01, n = 50)
# [1] "finished"
# [1] 4
# [,1]
# [1,]  1.18613420
# [2,]  0.95005793
# [3,] -0.96895420
# [4,]  0.35054770
# [5,]  0.02557918
# [6,] -1.54630304
plottheta <- as.data.frame(t(as.data.frame(plottheta)))
plottheta$loop <- 1:nrow(plottheta)
df <- data.frame(0,0,0,0,0,0,0)
colnames(df) <- colnames(plottheta)
plottheta <- rbind(df, plottheta)
plot(x = plottheta$loop, y = plottheta[[1]], pch = 1, ylim = c(-1.6, 1.2))
points(x = plottheta$loop, y = plottheta[[2]], pch = 2)
points(x = plottheta$loop, y = plottheta[[3]], pch = 3)     
points(x = plottheta$loop, y = plottheta[[4]], pch = 4)     
points(x = plottheta$loop, y = plottheta[[5]], pch = 5)     
points(x = plottheta$loop, y = plottheta[[6]], pch = 6)


newton(x = X, y = y, theta = theta, eps = 0.001, n = 50)
# [1] "finished"
# [1] 5
# [,1]
# [1,]  1.18613508
# [2,]  0.95005862
# [3,] -0.96895491
# [4,]  0.35054785
# [5,]  0.02557901
# [6,] -1.54630446


plottheta <- as.data.frame(t(as.data.frame(plottheta)))
plottheta$loop <- 1:nrow(plottheta)
df <- data.frame(0,0,0,0,0,0,0)
colnames(df) <- colnames(plottheta)
plottheta <- rbind(df, plottheta)
plot(x = plottheta$loop, y = plottheta[[1]], pch = 1, ylim = c(-1.6, 1.2))
points(x = plottheta$loop, y = plottheta[[2]], pch = 2)
points(x = plottheta$loop, y = plottheta[[3]], pch = 3)     
points(x = plottheta$loop, y = plottheta[[4]], pch = 4)     
points(x = plottheta$loop, y = plottheta[[5]], pch = 5)     
points(x = plottheta$loop, y = plottheta[[6]], pch = 6)

#5
getdata <- function(p){
  n <- 200 
  set.seed(1) 
  X <<- matrix(rnorm(n*p),ncol=p) 
  X <<- cbind(rep(1,n),X) 
  beta <- matrix(rnorm(p+1),ncol=1) 
  eta <- X%*%beta 
  Prob <- exp(eta)/(1+exp(eta)) 
  y <<- as.numeric(runif(n) <= Prob)
  
  theta <<- matrix(0, nrow = p+1, ncol = 1)
}

set.seed(124)

timecost<-list()
result <- list()
lapply(1:6, function(i){
  index <- c(2, 5, 10, 20, 40, 80)
  getdata(index[i])
  timecost[[i]]<<-system.time(newton(x = X, y = y, theta = theta, eps = 0.00001, n = 5))
  result[[i]] <<- newton(x = X, y = y, theta = theta, eps = 0.00001, n = 5)
})
newtontime <- sapply(1:6, function(i){timecost[[i]][3]})
plot(x = c(2, 5, 10, 20, 40, 80), y = newtontime)


#6
newton_ap <- function(x, y, theta, eps = 0.001, n = 100){
  x1 <- x
  theta1 <- theta
  y1 <- y
  i <- 0
  while(i < n){
    Wd <- sapply(1:nrow(x1), function(i){get_l(x1, theta1, i)*(1 - get_l(x1, theta1, i))})
    W <- diag(Wd)
    Pi <- sapply(1:nrow(x1), function(i){get_l(x1, theta1, i)})
    Pi <- matrix(Pi, nrow = length(Pi))
    WW <- diag(as.numeric(diag(t(x1)%*%W%*%x1)))
    thetaf <- theta1 + solve(WW) %*% t(x1) %*% (y1 - Pi)
    if(sqrt(sum((thetaf - theta1)^2)) < eps){
      print("finished")
      break
    }else{
      i <- i+1
      theta1 <- thetaf
    }
  }
  print(i)
  if(i == n){
    print("not convergence after n loops")
  }
  return(thetaf)
}

set.seed(124)
timecost<-list()
result2 <- list()
lapply(1:6, function(i){
  index <- c(2, 5, 10, 20, 40, 80)
  getdata(index[i])
  
  timecost[[i]]<<-system.time(newton_ap(x = X, y = y, theta = theta, eps = 0.00001, n = 5))
  result2[[i]] <<- newton_ap(x = X, y = y, theta = theta, eps = 0.00001, n = 5)
})
newtonaptime <- sapply(1:6, function(i){timecost[[i]][3]})
plot(x = c(2, 5, 10, 20, 40, 80), y = newtontime)
points(x = c(2, 5, 10, 20, 40, 80), y = newtonaptime, pch = 3)

delta <- sapply(1:6, function(i){sqrt(sum((result[[i]] - result2[[i]])^2))})


set.seed(124)
result3 <- list()
lapply(1:6, function(i){
  index <- c(2, 5, 10, 20, 40, 80)
  getdata(index[i])
  df <- as.data.frame(cbind(y,X[, -1]))
  result3[[i]] <<- glm(y ~ ., data = df, family = "binomial")
})

plot (x = c(2, 5, 10, 20, 40, 80), y = delta/c(2, 5, 10, 20, 40, 80), main = "Norm2 of difference between solution of two method")
delta2 <- sapply(1:6, function(i){sqrt(sum((result3[[i]]$coefficients - result2[[i]])^2))})
delta3 <- sapply(1:6, function(i){sqrt(sum((result3[[i]]$coefficients - result[[i]])^2))})
plot (x = c(2, 5, 10, 20, 40, 80), y = delta3, pch = 3, main = "Norm2 of difference between solution and glm")
points(x = c(2, 5, 10, 20, 40, 80), y = delta2)
delta2 > delta3
