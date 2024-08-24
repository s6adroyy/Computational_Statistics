rm(list=ls())
set.seed(100)

library(ROCR)
library(MASS)
library(mvtnorm)
library(ggplot2)
suppressMessages(library(pls))
install.packages("devtools")
install.packages("MASS")
install_github("vqv/ggbiplot")
suppressMessages(library(MASS))
suppressMessages(library(devtools))
suppressMessages(library(ggbiplot))
suppressMessages(library(glmnet))
library(tidyverse)
library(broom)
library(glmnet)
suppressMessages(library(psych))
library (plyr)

n_1 = 124 #Sample size 

eps_sd = sqrt(10)
p <- 9



X_1 <- rnorm (n_1, 6.30, 1.64)
X_2 <- rnorm (n_1, 5.83, 3.84)
X_3 <- rnorm (n_1, 131.32, 45.71)
X_4 <- rnorm (n_1, 86.96, 12.66)
X_5 <- rnorm (n_1, 26023.38, 4548.9)
X_6 <- rnorm (n_1, 61.87, 0.68)
X_7 <- rnorm (n_1, 12.38, 0.43)
X_8 <- rnorm (n_1, 79.96, 5.39)
X_9 <- rnorm (n_1, 23.37, 3.90)


X<- cbind(X_1, X_2, X_3, X_4, X_5, X_6, X_7,X_8, X_9)
eps <- rnorm (n_1 , 0, eps_sd)


beta <- c(1, -1, -1, 1, 0.3, -1, 0, -0.31, -1, 1)

Y<- beta[1]+beta[2]*X_1+beta[3]*X_2+beta[4]*X_3+beta[5]*X_4+beta[6]*X_5+
  beta[7]*X_6+beta[8]*X_7+beta[9]*X_8+beta[10]*X_9 +eps

data<- data.frame(X, Y)

suppressMessages(library(pls))



#### Choosing the number of Components by Cross-Validation ###

index <-sample(1:nrow(data), nrow(data)*0.5)

train.data <- data [index,]
test.data <- data[-index,]

y.train <- train.data$Y
y.test <- test.data$Y

x.train <- data.matrix (train.data[ ,-10]) 
x.test <- data.matrix(test.data[ ,-10])

### Ridge #####

#list of lambda values for the model to try

grid <- 10^seq(2,-2,length.out=100)

## Performing Ridge over the values of lambda ###
cv.ridge <- cv.glmnet(x.train, y.train, alpha = 0, lambda = grid)
summary(cv.ridge)

plot(cv.ridge)

### Finding the optimal Lambda ####

opt.lambda.ridge <- cv.ridge$lambda.min
opt.lambda.ridge

### Finally performing the ridge on the optimal lambda chosen ###

ridge.model <- glmnet(x.train, y.train, lambda = opt.lambda.ridge, alpha = 0)
coef(ridge.model)

ridge.pred <- predict (ridge.model, newx = x.test)
mse.ridge <- mean ((ridge.pred - y.test)^2)
mse.ridge

#### Performing lasso over values of lambda ###

cv.lasso <- cv.glmnet(x.train, y.train, alpha = 1, lambda = grid)
summary(cv.lasso)

plot (cv.lasso)

### Finding the optimal lambda for lasso ###

opt.lambda.lasso <- cv.lasso$lambda.min
opt.lambda.lasso

### Finally performing the lasso on the optimal lambda chosen ###

lasso.model<- glmnet (x.train, y.train, lambda = opt.lambda.lasso, alpha =1)
coef(lasso.model)
lasso.plot <- glmnet(x.train, y.train, family = 'gaussian', alpha = 1)

#plot(lasso.plot)
#coefficients(lasso.model)

lasso.pred <- predict (lasso.model, newx = x.test)
mse.lasso <- mean ((lasso.pred - y.test)^2)
mse.lasso


#### OLS regression ####

lm.fit <- lm(Y ~ ., data = train.data)
summary(lm.fit)
lm.fit.predict <- predict(lm.fit, newdata = test.data)

### OLS MSE #####

mse.ols <- mean((lm.fit.predict-y.test)^2)
mse.ols

all.mse <- cbind(mse.ridge, mse.lasso, mse.ols)
colnames(all.mse) <- c("Ridge", "Lasso", "OLS")

all.mse


#########################################################

###### DGP SIMULATION #####

### Deleting the insignificant covariates X_6 and X_7 ## 


set.seed(122)

n = 124
p1 = 9 



### Using the values of the mean and the coefficients form the base paper####

####Number of repetitions ####


rep <- 50

#### Set values ####
X_1 <- rnorm (n, 6.30, 1.64)
X_2 <- rnorm (n, 5.83, 3.84)
X_3 <- rnorm (n, 131.32, 45.71)
X_4 <- rnorm (n, 86.96, 12.66)
X_5 <- rnorm (n, 26023.38, 4548.9)
X_8 <- rnorm (n, 79.96, 5.39)
X_9 <- rnorm (n, 23.37, 3.90)

X.dgp1<- cbind(X_1, X_2, X_3, X_4, X_5,X_8, X_9)

cor (X.dgp1, method = "pearson")

eps_sd = sqrt(10)

eps <- rnorm (n , 0, eps_sd)



beta.dgp1 <- c(1, -1, -1, 1, 0.3, -1, 0, -0.31, -1, 1)

Y.dgp1<- beta.dgp1[1]*X_1+beta.dgp1[2]*X_2+beta.dgp1[3]*X_3+beta.dgp1[4]*X_4+beta.dgp1[5]*X_5+beta.dgp1[6]*X_8+
          beta.dgp1[7]*X_9 +eps

data.dgp1 <- data.frame (X.dgp1 , Y.dgp1)

index <-sample(1:nrow(data.dgp1), nrow(data.dgp1)*0.5)

train.data.dgp1 <- data.dgp1 [index,]
test.data.dgp1 <- data.dgp1[-index,]

y.train.dgp1 <- train.data.dgp1$Y.dgp1
y.test.dgp1 <- test.data.dgp1$Y.dgp1

x.train.dgp1 <- data.matrix (train.data.dgp1[ ,-10]) 
x.test.dgp1 <- data.matrix (test.data.dgp1[, -10]) 

### Set grid ##
grid <- 10^seq(2,-2,length.out=100)

MSE.sim <- matrix(NA,rep,3)
colnames(MSE.sim) <- c("OLS", "Ridge", "Lasso")

####################################################################3

for (i in 1:rep){
  
  ## training sample
  
  X_1 <- rnorm (n, 6.30, 1.64)
  X_2 <- rnorm (n, 5.83, 3.84)
  X_3 <- rnorm (n, 131.32, 45.71)
  X_4 <- rnorm (n, 86.96, 12.66)
  X_5 <- rnorm (n, 26023.38, 4548.9)
  
  X_8 <- rnorm (n, 79.96, 5.39)
  X_9 <- rnorm (n, 23.37, 3.90)
  
  X.dgp1<- cbind(X_1, X_2, X_3, X_4, X_5,X_8, X_9)
  
  eps_sd = sqrt(10)
 
  beta.dgp1 <- c(-1, -1, 1, 1, 0.3, -0.31, 1)

  eps <- rnorm (n , 0, eps_sd)
  
 
  Y.dgp1<- beta.dgp1[1]*X_1+beta.dgp1[2]*X_2+beta.dgp1[3]*X_3+beta.dgp1[4]*X_4+beta.dgp1[5]*X_5+beta.dgp1[6]*X_8+
    beta.dgp1[7]*X_9 +eps
  
  
  data.dgp1 <- data.frame (X.dgp1 , Y.dgp1)
  
  index <-sample(1:nrow(data.dgp1), nrow(data.dgp1)*0.5)
  
  train.data.dgp1 <- data.dgp1 [index,]
  test.data.dgp1 <- data.dgp1[-index,]
  
  y.train.dgp1 <- train.data.dgp1$Y.dgp1
  y.test.dgp1 <- test.data.dgp1$Y.dgp1
  
  x.train.dgp1 <- data.matrix (train.data.dgp1[ ,-10]) #glmnet requires it to be a matrix
  x.test.dgp1 <- data.matrix (test.data.dgp1[, -10]) #glmnet requires it to be a matrix
 
  ## Implement the 3 methods 
  
 ### OLS ####
  
  lm.fit.dgp1 <- lm(Y.dgp1 ~ ., data = train.data.dgp1)
  summary(lm.fit.dgp1)
  lm.fit.predict.dgp1 <- predict(lm.fit.dgp1, newdata = test.data.dgp1)
  
  ### OLS MSE #####
  
  MSE.sim[i,1]<- mean((lm.fit.predict.dgp1-y.test.dgp1)^2)
  MSE.sim[i,1]
  
### Ridge #  
  
  cv_ridge.sim <- cv.glmnet(x.train.dgp1,y.train.dgp1 , alpha = 0, lambda = grid)
  
  opt_lambda.ridge.sim <- cv_ridge.sim$lambda.min
  
  ridge.model.dgp1 <- glmnet(x.train.dgp1, y.train.dgp1, lambda = opt_lambda.ridge.sim, alpha = 0)
  
  ridge.model.dgp1.graph <- glmnet(x.train.dgp1, y.train.dgp1,family = "gaussian", alpha = 0)
 
   plot(  ridge.model.dgp1.graph)
  
  coef(ridge.model.dgp1)
  
  ridge.pred.dgp1 <- predict (ridge.model.dgp1, newx = x.test.dgp1)
  
  ### Ridge MSE ###
  
  MSE.sim[i,2] <- mean ((ridge.pred.dgp1 - y.test.dgp1)^2)
  MSE.sim[i,2]
  
  #### Lasso ####
  
  cv_lasso.sim <- cv.glmnet(x.train.dgp1, y.train.dgp1, alpha = 1, lambda = grid)
  
  opt_lambda.lasso.sim <- cv_lasso.sim$lambda.min
  
  lasso.model.dgp1<- glmnet (x.train.dgp1, y.train.dgp1, lambda = opt_lambda.lasso.sim, alpha =1)
  
  
  lasso.model.dgp1.graph<- glmnet (x.train.dgp1, y.train.dgp1, family = "gaussian", alpha =1)
  plot(lasso.model.dgp1.graph)
  
  coef(lasso.model.dgp1)
  
  lasso.pred.dgp1 <- predict (lasso.model.dgp1, newx = x.test.dgp1)
  
  ### Lasso MSE ###
  
MSE.sim[i,3] <- mean ((lasso.pred.dgp1 - y.test.dgp1)^2)
MSE.sim[i,3]

}


  
#Compute row means for the 3 methods ######

colMeans(MSE.sim)
  
##############################

###### Including a covariance matrix of my choice #####

## Cov matrix function


set.seed (300)

p=9

offdiag <- function (x, at =0) {
  if (is.matrix(x)) {
    result <- x[row(x) == col(x) - at]
  } else {
    len <- length(x)
    result <- matrix(0.8, nrow = len + abs(at), ncol = len + abs(at))
    result[row(result) == col(result) - at] <- x
  }
  return(result)
}

cov.sim<-offdiag(x= rep(1,p), at=0)
true.means <- c(6.30, 5.83, 131.32, 86.96, 26023.38, 61.87, 12.38, 79.96, 23.37)
beta.dgp1 <- c(1, -1, -1, 1, 0.3, -1, 0, -0.31, -1, 1)
n = 500

eps_sd = sqrt(10)
eps <- rnorm (n , 0, eps_sd)


rep <- 50

grid <- 10^seq(2,3,length.out=100)

MSE.sim.dgp2 <- matrix(NA,rep,3)
colnames(MSE.sim.dgp2) <- c("OLS", "Ridge", "Lasso")

for (i in 1:rep){
  
  ## training sample
  cov.sim<-offdiag(x= rep(1,p), at=0)
  df<- mvrnorm(n=n, mu=true.means, Sigma=cov.sim)
  X_f<- rep(1, nrow(df))
  X.train.dgp2 <- as.data.frame(cbind(X_f,df))
  eps.train <- rnorm (n , 0, eps_sd)
  Y.train.dgp2 <- as.matrix(X.train.dgp2) %*% beta.dgp1 + eps.train
  data.train.dgp2 <- data.frame ("Y.sim"=Y.train.dgp2, "X.sim"=X.train.dgp2)
  
  ## test sample
  df_test<- mvrnorm(n=n, mu=true.means, Sigma=cov.sim)
  X_f_test<-rep(1, nrow(df_test))
  X.test.dgp2 <- as.data.frame(cbind(X_f_test,df_test))
  eps.test <- rnorm (n , 0, eps_sd)
  Y.test.dgp2<- as.matrix(X.test.dgp2) %*% beta.dgp1+ eps.test
  data.test.dgp2 <- data.frame ("Y.sim.test"=Y.test.dgp2, "X.sim.test"=X.test.dgp2)
  

  #### OLS ####
  
  lm.fit.dgp2<-glmnet(as.matrix(X.train.dgp2), as.matrix(Y.train.dgp2), alpha = 0, lambda = 0,family="gaussian")
  lm.fit.predict.dgp2<-predict(lm.fit.dgp2,newx=as.matrix(X.test.dgp2))
  
  
  # OLS MSE 
  MSE.sim.dgp2[i,1] <- mean((lm.fit.predict.dgp2-Y.test.dgp2)^2)

  #
  
  ########## Ridge
  
  cv_ridge.dgp2 <- cv.glmnet(as.matrix(X.train.dgp2), as.matrix(Y.train.dgp2), alpha = 0, lambda = grid)
  opt_ridgelambda.dgp2 <- cv_ridge.dgp2$lambda.min
  
  ridge.sim.dgp2 <- glmnet(as.matrix(X.train.dgp2), as.matrix(Y.train.dgp2), alpha = 0, family = 'gaussian', lambda = opt_ridgelambda.dgp2, intercept = FALSE)
  
  coef(ridge.sim.dgp2)
  
  ridge.predict.dgp2 <- predict(ridge.sim.dgp2, newx = as.matrix(X.test.dgp2))
  
  # Ridge MSE
  MSE.sim.dgp2[i,2] <- mean(( ridge.predict.dgp2 - Y.test.dgp2)^2)
  
  ########## LASSO
  
  cv_lasso.dgp2 <- cv.glmnet(as.matrix(X.train.dgp2), as.matrix(Y.train.dgp2), alpha = 1, lambda = grid)
  opt_lassolambda.dgp2 <- cv_lasso.dgp2$lambda.min
  
  lasso.sim.dgp2 = glmnet(as.matrix(X.train.dgp2), as.matrix(Y.train.dgp2), alpha = 1, family = 'gaussian', lambda = opt_lassolambda.dgp2, intercept = FALSE, standardize = TRUE)
  
  coef( lasso.sim.dgp2)
  
  lasso.predict.dgp2 <- predict(lasso.sim.dgp2, newx = as.matrix(X.test.dgp2))
  
  # LASSO MSE
  MSE.sim.dgp2[i,3] <- mean(( lasso.predict.dgp2 - Y.test.dgp2)^2)
  
}


  
######### Compute row means for the three methods for the new dgp with own cov #####

colMeans(MSE.sim.dgp2)


######### DGP with increase in  lambda 

set.seed (300)

p=7

offdiag <- function (x, at =0) {
  if (is.matrix(x)) {
    result <- x[row(x) == col(x) - at]
  } else {
    len <- length(x)
    result <- matrix(0.6, nrow = len + abs(at), ncol = len + abs(at))
    result[row(result) == col(result) - at] <- x
  }
  return(result)
}

cov.sim<-offdiag(x= rep(1,p), at=0)
true.means <- c(6.30, 5.83, 131.32, 86.96, 26023.38, 79.96, 23.37)
beta.dgp1 <-   c(1, -1, -1, 1, 0.3, -1, -1, 1)
n = 124

eps_sd = sqrt(10)
eps <- rnorm (n , 0, eps_sd)



rep <- 50

grid <- 5^seq(3,4,length.out=100)

MSE.sim.dgp2 <- matrix(NA,rep,3)
colnames(MSE.sim.dgp2) <- c("OLS", "Ridge", "Lasso")

for (i in 1:rep){
  
  ## training sample
  cov.sim<-offdiag(x= rep(1,p), at=0)
  df<- mvrnorm(n=n, mu=true.means, Sigma=cov.sim)
  X_f<- rep(1, nrow(df))
  X.train.dgp2 <- as.data.frame(cbind(X_f,df))
  eps.train <- rnorm (n , 0, eps_sd)
  Y.train.dgp2 <- as.matrix(X.train.dgp2) %*% beta.dgp1 + eps.train
  data.train.dgp2 <- data.frame ("Y.sim"=Y.train.dgp2, "X.sim"=X.train.dgp2)
  
  ## test sample
  df_test<- mvrnorm(n=n, mu=true.means, Sigma=cov.sim)
  X_f_test<-rep(1, nrow(df_test))
  X.test.dgp2 <- as.data.frame(cbind(X_f_test,df_test))
  eps.test <- rnorm (n , 0, eps_sd)
  Y.test.dgp2<- as.matrix(X.test.dgp2) %*% beta.dgp1+ eps.test
  data.test.dgp2 <- data.frame ("Y.sim.test"=Y.test.dgp2, "X.sim.test"=X.test.dgp2)
  
  
  #### OLS ####
  
  lm.fit.dgp2<-glmnet(as.matrix(X.train.dgp2), as.matrix(Y.train.dgp2), alpha = 0, lambda = 0,family="gaussian")
  lm.fit.predict.dgp2<-predict(lm.fit.dgp2,newx=as.matrix(X.test.dgp2))
  
  
  # OLS MSE 
  MSE.sim.dgp2[i,1] <- mean((lm.fit.predict.dgp2-Y.test.dgp2)^2)
  
  
  ########## Ridge
  
  cv_ridge.dgp2 <- cv.glmnet(as.matrix(X.train.dgp2), as.matrix(Y.train.dgp2), alpha = 0, lambda = grid)
  opt_ridgelambda.dgp2 <- cv_ridge.dgp2$lambda.min
  
  ridge.sim.dgp2 <- glmnet(as.matrix(X.train.dgp2), as.matrix(Y.train.dgp2), alpha = 0, family = 'gaussian', lambda = opt_ridgelambda.dgp2, intercept = FALSE)
  
  ridge.sim.dgp2plot <- glmnet(as.matrix(X.train.dgp2), as.matrix(Y.train.dgp2), alpha = 0, family = 'gaussian')
  plot (ridge.sim.dgp2plot)
  
  coef(ridge.sim.dgp2)
  
  ridge.predict.dgp2 <- predict(ridge.sim.dgp2, newx = as.matrix(X.test.dgp2))
  
  # Ridge MSE
  MSE.sim.dgp2[i,2] <- mean(( ridge.predict.dgp2 - Y.test.dgp2)^2)
  
  ########## LASSO
 
  cv_lasso.dgp2 <- cv.glmnet(as.matrix(X.train.dgp2), as.matrix(Y.train.dgp2), alpha = 1, lambda = grid)
  opt_lassolambda.dgp2 <- cv_lasso.dgp2$lambda.min
  
  lasso.sim.dgp2 = glmnet(as.matrix(X.train.dgp2), as.matrix(Y.train.dgp2), alpha = 1, family = 'gaussian', lambda = opt_lassolambda.dgp2, intercept = FALSE, standardize = TRUE)
  
  lasso.sim.dgp2graph <- glmnet(as.matrix(X.train.dgp2), as.matrix(Y.train.dgp2), alpha = 1, family = 'gaussian')
  plot (lasso.sim.dgp2graph)
  
  coef( lasso.sim.dgp2)
  
  lasso.predict.dgp2 <- predict(lasso.sim.dgp2, newx = as.matrix(X.test.dgp2))
  
  # LASSO MSE
  MSE.sim.dgp2[i,3] <- mean(( lasso.predict.dgp2 - Y.test.dgp2)^2)
  
}



######### Compute row means for the three indicators for the new dgp with own cov #####

colMeans(MSE.sim.dgp2)


  ###### DGP 3 ##### Evaluating the performance of ridge , lasso and ols based on the true parameters 
  
  set.seed(200)
  
  n = 200
  p = 20
  
  beta.vec.dgp3 <- c(rep(1, 5), rep(0.3, 5), rep(0, 10))
  
  gen_datasets.dgp3 = function(p, sample_size){
    eps <- rnorm(n = sample_size, mean = 0, sd = sqrt(10))
    beta.vec.dgp3 <- c(rep(1, 5), rep(0.3, 5), rep(0, 10))
    X3 <- rmvnorm(n= sample_size, mean=rep(0, p), sigma=diag(sample(1:2, p, replace = T), nrow = p, ncol = p))
    Y3 <- X3 %*% beta.vec.dgp3 + eps
    return(data.frame(Y3, X3))  
  }
  
  test_set.3 = gen_datasets.dgp3(20, sample_size = 200)
  test_setX3 = test_set.3[, -1]
  test_setY3 = test_set.3[,  1]
  
  # Expected Prediction Error OLS
  epe.OLS3 = function(p, sample_size){
    # Creating a training set 
    train_set.3 = gen_datasets.dgp3(p = p, sample_size = sample_size)
    train_set.3$X3 = train_set.3[, -1]
    beta.ols3 = solve(t(as.matrix(train_set.3$X3)) %*% as.matrix(train_set.3$X3)) %*% t(as.matrix(train_set.3$X3)) %*% as.matrix(train_set.3$Y3)
    # Sum of errors for OLS
    ols.error = mean((as.matrix(test_setY3) - (as.matrix(test_setX3) %*% beta.ols3))**2)
    # Number of Correct zero parameters
    nzc = sum(beta.ols3[11:20] == beta.vec.dgp3[11:20])
    
    frame = cbind(nzc, ols.error)
    colnames(frame) <- c("Number of correct zero coefficients", "Prediction Error OLS")
    return(frame)
  }
  
  library (plyr)
  replication.OLS = replicate(100, epe.OLS3(p = 20, sample_size = 200), simplify = FALSE)
  replication.OLS.df <- ldply(replication.OLS, data.frame)
  sim.epe = apply(replication.OLS.df, 2, mean)
  simulation.result.OLS = data.frame(sim.epe)
  rownames(simulation.result.OLS) = c("Number of Correct Zero Coefficients", "Expected Prediction Error OLS")
  colnames(simulation.result.OLS) = c("Simulation Results: OLS")
  simulation.result.OLS
  
  epe.REG = function(p, sample_size){
    train_set.3 = gen_datasets.dgp3(p = p, sample_size = sample_size)
    train_set.3$X3 = train_set.3[, -1]
    grid = seq(0, 1, by = 0.02)
    beta.vec.dgp3 <- c(rep(1, 5), rep(0.3, 5), rep(0, 10))
    ridge.model = glmnet(as.matrix(train_set.3$X3), as.matrix(train_set.3$Y3), alpha = 0, lambda = grid, intercept = FALSE, 
                         standardize = F)
    lasso.model = glmnet(as.matrix(train_set.3$X3), as.matrix(train_set.3$Y3), alpha = 1, lambda = grid, intercept = FALSE, 
                         standardize = F)
    cv.ridge.model = cv.glmnet(as.matrix(train_set.3$X3), as.matrix(train_set.3$Y3), alpha = 0, lambda = grid)
    cv.lasso.model = cv.glmnet(as.matrix(train_set.3$X3), as.matrix(train_set.3$Y3), alpha = 1, lambda = grid)
    best.lambda.ridge = cv.ridge.model$lambda.min
    best.lambda.lasso = cv.lasso.model$lambda.min
    ridge.pred.opt = predict(ridge.model, s = best.lambda.ridge, newx = as.matrix(test_setX3))
    lasso.pred.opt = predict(lasso.model, s = best.lambda.lasso, newx = as.matrix(test_setX3))
    ridge.coef.opt = predict(ridge.model, type = "coefficients", s = best.lambda.ridge)[1:21, ]
    lasso.coef.opt = predict(lasso.model, type = "coefficients", s = best.lambda.lasso)[1:21, ]
    nzc.ridge = sum(ridge.coef.opt[12:21] == beta.vec.dgp3[11:20])
    nzc.lasso = sum(lasso.coef.opt[12:21] == beta.vec.dgp3[11:20])
    pe.ridge = mean((test_setY3 - ridge.pred.opt)**2)
    pe.lasso = mean((test_setY3 - lasso.pred.opt)**2)
    frame.REG = cbind(nzc.ridge, pe.ridge, nzc.lasso, pe.lasso)
    colnames(frame.REG) = c("Number of Correct Zero Coefficients Ridge", "Prediction Error Ridge",
                            "Number of Correct Zero Coefficients LASSO", "Prediction Error LASSO") 
    return(frame.REG)
  }
  simulation.REG = apply(replicate(100, epe.REG(20, 200)), 2, mean)
  ridge.sim = data.frame(simulation.REG[1:2])
  lasso.sim = data.frame(simulation.REG[3:4])
  sim.res.reg = cbind(ridge.sim, lasso.sim)
  colnames(sim.res.reg) = c("Simulation Results: Ridge", "Simulation Results: LASSO")
  rownames(sim.res.reg) = c("Number of Correct Zero Coefficients", "Expected Prediction Error")
  sim.res.reg
  
  
  ################################
  
  set.seed (200)
  beta.vec <- c(rep(1, 5), rep(0.03, 5), rep(0, 10))
  
  gen_datasets = function(p, sample_size){
    eps <- rnorm(n = sample_size, mean = 0, sd = 1)
    beta.vec <- c(rep(1, 5), rep(0.03, 5), rep(0, 10))
    X <- rmvnorm(n= sample_size, mean=rep(0, p), sigma=diag(sample(1:2, p, replace = T), nrow = p, ncol = p))
    Y <- X %*% beta.vec + eps
    return(data.frame(Y, X))  
  }
  
  
  test_set = gen_datasets(20, sample_size = 100)
  test_setX = test_set[, -1]
  test_setY = test_set[,  1]
  
  # Expected Prediction Error OLS
  epe.OLS = function(p, sample_size){
    # Creating a training set 
    train_set = gen_datasets(p = p, sample_size = sample_size)
    train_set$X = train_set[, -1]
    beta.ols = solve(t(as.matrix(train_set$X)) %*% as.matrix(train_set$X)) %*% t(as.matrix(train_set$X)) %*% as.matrix(train_set$Y)
    # Sum of errors for OLS
    ols.error = mean((as.matrix(test_setY) - (as.matrix(test_setX) %*% beta.ols))**2)
    # Number of Correct zero parameters
    nzc = sum(beta.ols[11:20] == beta.vec[11:20])
    
    frame = cbind(nzc, ols.error)
    colnames(frame) <- c("Number of correct zero coefficients", "Prediction Error OLS")
    return(frame)
  }
  library (plyr)
  replication.OLS = replicate(100, epe.OLS(p = 20, sample_size = 100), simplify = FALSE)
  replication.OLS.df <- ldply(replication.OLS, data.frame)
  sim.epe = apply(replication.OLS.df, 2, mean)
  simulation.result.OLS = data.frame(sim.epe)
  rownames(simulation.result.OLS) = c("Number of Correct Zero Coefficients", "Expected Prediction Error OLS")
  colnames(simulation.result.OLS) = c("Simulation Results: OLS")
  simulation.result.OLS
  
  library(glmnet)
  epe.REG = function(p, sample_size){
    train_set = gen_datasets(p = p, sample_size = sample_size)
    train_set$X = train_set[, -1]
    grid = seq(0, 1, by = 0.02)
    beta.vec <- c(rep(1, 5), rep(0.03, 5), rep(0, 10))
    ridge.model = glmnet(as.matrix(train_set$X), as.matrix(train_set$Y), alpha = 0, lambda = grid, intercept = FALSE, 
                         standardize = F)
    lasso.model = glmnet(as.matrix(train_set$X), as.matrix(train_set$Y), alpha = 1, lambda = grid, intercept = FALSE, 
                         standardize = F)
    cv.ridge.model = cv.glmnet(as.matrix(train_set$X), as.matrix(train_set$Y), alpha = 0, lambda = grid)
    cv.lasso.model = cv.glmnet(as.matrix(train_set$X), as.matrix(train_set$Y), alpha = 1, lambda = grid)
    best.lambda.ridge = cv.ridge.model$lambda.min
    best.lambda.lasso = cv.lasso.model$lambda.min
    ridge.pred.opt = predict(ridge.model, s = best.lambda.ridge, newx = as.matrix(test_setX))
    lasso.pred.opt = predict(lasso.model, s = best.lambda.lasso, newx = as.matrix(test_setX))
    ridge.coef.opt = predict(ridge.model, type = "coefficients", s = best.lambda.ridge)[1:21, ]
    lasso.coef.opt = predict(lasso.model, type = "coefficients", s = best.lambda.lasso)[1:21, ]
    nzc.ridge = sum(ridge.coef.opt[12:21] == beta.vec[11:20])
    nzc.lasso = sum(lasso.coef.opt[12:21] == beta.vec[11:20])
    pe.ridge = mean((test_setY - ridge.pred.opt)**2)
    pe.lasso = mean((test_setY - lasso.pred.opt)**2)
    frame.REG = cbind(nzc.ridge, pe.ridge, nzc.lasso, pe.lasso)
    colnames(frame.REG) = c("Number of Correct Zero Coefficients Ridge", "Prediction Error Ridge",
                            "Number of Correct Zero Coefficients LASSO", "Prediction Error LASSO") 
    return(frame.REG)
  }
  simulation.REG = apply(replicate(100, epe.REG(20, 100)), 2, mean)
  ridge.sim = data.frame(simulation.REG[1:2])
  lasso.sim = data.frame(simulation.REG[3:4])
  sim.res.reg = cbind(ridge.sim, lasso.sim)
  colnames(sim.res.reg) = c("Simulation Results: Ridge", "Simulation Results: LASSO")
  rownames(sim.res.reg) = c("Number of Correct Zero Coefficients", "Expected Prediction Error")
  sim.res.reg