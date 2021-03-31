# NONPARAMETRIC STATISTICS FINAL ASSIGNMENT
# DANYU ZHANG & DANIEL ALONSO
# MASTER IN STATISTICS FOR DATA SCIENCE


###################################################################################
###################################################################################
## CATEGORY A, PROBLEM 6: Exercise 5.11 ##
########################### PART A
# Reading the challenger.txt data
challenger <- read.table('https://raw.githubusercontent.com/egarpor/handy/master/datasets/challenger.txt', header = TRUE)

X = challenger$temp # Assign temperature column to X
Y = challenger$fail.field # Assign fail.field column to X
plot(X, Y) # plotting a scatterplot of temp vs fail.field

logistic <- function(x) 1 / (1 + exp(-x)) # define logistic function

# Bandwidth that undersmooths
h <- 0.2 # bandwidth
x <- seq(10, 28, l = 5000) # x grid of points to cover

suppressWarnings(
  fit_glm <- sapply(x, function(x) {
    # setting the weights for each data point using Gaussian kernels
    K <- dnorm(x = x, mean = X, sd = h)
    # fitting the local logistic regression and taking the coefficient
    glm.fit(x = cbind(1, X - x), y = Y, weights = K,
            family = binomial())$coefficients[1]
  })
)

# plotting the estimated probabilities of having an incident in function of temperature when overfitting
plot(x, logistic(fit_glm)) 

# Bandwidth that oversmooths
h <- 10
x <- seq(10, 28, l = 5000)

suppressWarnings(
  fit_glm <- sapply(x, function(x) {
    # setting the weights for each data point using Gaussian kernels
    K <- dnorm(x = x, mean = X, sd = h) 
    # fitting local logistic regression and taking the coefficient
    glm.fit(x = cbind(1, X - x), y = Y, weights = K, 
            family = binomial())$coefficients[1]
  })
)

# plotting the estimated probabilities of having an incident in function of temperature when underfitting
plot(x, logistic(fit_glm)) 

# Adequate Bandwidth
h <- 2
x <- seq(10, 28, l = 5000)

suppressWarnings(
  fit_glm <- sapply(x, function(x) {
    # setting the weights for each data point using Gaussian kernels
    K <- dnorm(x = x, mean = X, sd = h)
    # fitting local logistic regression and taking the coefficient
    glm.fit(x = cbind(1, X - x), y = Y, weights = K,
            family = binomial())$coefficients[1]
  })
)

# plotting the estimated probabilities of having an incident in function of temperature when fitting with an adequate bandwidth
plot(x, logistic(fit_glm)) 

########################### PART B
n <- length(Y) # number of data points 
h <- seq(1.5, 2.3, by=0.1) # grid of bandwidths' values that we want to test in order to have a maximized value of likelihood

suppressWarnings(
  LCV <- sapply(h, function(h) {
    sum(sapply(1:n, function(i) {
      # kernel
      K <- dnorm(x = X[i], mean = X[-i], sd = h)
      # optimizing function
      nlm(f = function(beta) { 
        # Maximum likelihood function
        -sum(K * (Y[-i] * (beta[1] + beta[2] * (X[-i] - X[i])) -
                    log(1 + exp(beta[1] + beta[2] * (X[-i] - X[i])))))
      }, p = c(0, 0))$minimum
    }))
  })
)

# plotting the likelihoods in function of bandwidths' values
plot(h, LCV, type = "o")
abline(v = h[which.max(LCV)], col = 2)


########################### PART C
h <- 2
# setting the grid with the values that we want to predict in order to visualize it from the plot temperature vs probability
x <- seq(-0.6, 11.67, l = 500) 

suppressWarnings(
  fit_glm <- sapply(x, function(x) {
    K <- dnorm(x = x, mean = X, sd = h)
    glm.fit(x = cbind(1, X - x), y = Y, weights = K,
            family = binomial())$coefficients[1]
  })
)

# plotting the estimated probabilities of having an incident in function of temperature when fitting with the optimized bandwidth value
plot(x, logistic(fit_glm))

# Find probability at x=-0.6 and x=11.67
pred0.6 = logistic(fit_glm)[1]
pred11.67 = logistic(fit_glm)[length(fit_glm)]

pred0.6
pred11.67

########################### PART D
h <- 2
x1 <- -0.6
x2 <- 11.67

suppressWarnings(
  fit_glm1 <- sapply(x1, function(x1) {
    K <- dnorm(x = x1, mean = X, sd = h) # fitting at the point x1 = -0.6
    glm.fit(x = cbind(1, X - x1), y = Y, weights = K,
            family = binomial())$coefficients[1]
  })
)

suppressWarnings(
  fit_glm2 <- sapply(x2, function(x2) {
    K <- dnorm(x = x2, mean = X, sd = h) # fitting at the point x2 = 11.67
    glm.fit(x = cbind(1, X - x2), y = Y, weights = K,
            family = binomial())$coefficients[1]
  })
)

fit_glm1 # -0.6 fit
fit_glm2 # 11.67 fit

###################################################################################
###################################################################################
## CATEGORY B, PROBLEM 4: Exercise 4.9 ##
########################### PART A
# Local cubic estimator implementation
lce <- function(x, data, h, p=3) {
    ## INITIALIZING VARIABLES
    # Resulting vector initilization
    result <- c()

    # Predictors initialization (copied from data col 1)
    predictors <- data[,1]

    # Response initialization (copied from data col 2)
    Y <- data[,2]

    # e_1 (vector of size p+1 x 1 of zeros with the exception of the first entry)
    e_1 <-  matrix(c(c(1), rep(0,p)),nrow=p+1,ncol=1)

    # X matrix
    X <- list()

    ## MAIN LOOP
    for (k in 1:length(x)){
        # k-th matrix of the list of X, whose size will correspond with the length
        # of the values the user wants to predict (length of x)
        X[[k]] <- matrix(,nrow=length(predictors),ncol=p+1)

        # Filling up the k-th X matrix entry by entry
        for (i in 1:dim(X[[k]])[1]) {
            for (j in 1:dim(X[[k]])[2]) {
                # the k-th X matrix's i,j-th positions are filled
                # subtracting the i-th predictor's value with the k-th
                # element in the values the user wants to evaluate
                # all to the power j-1, where the first value obtained
                # is a 0 in the exponent, rendering the first column of each
                # X matrix a column of ones, and the last one to the power
                # of the p defined by the user (default p=3)
                X[[k]][i,j] <- (predictors[i] - x[k])^(j-1)
            }
        }

        # Weights, these change in every loop, therefore must be initialized
        # in every loop. 
        # Each weight corresponds to one element of the x vector provided
        # by the user
        weights <- c()
        for (i in 1:length(predictors)) {
            # We use the K_h normal kernel of the difference between
            # the i-th predictor and the k-th element of x
            weights[i] <- pnorm((predictors[i] - x[k])/h)/h
        }
        # We create a matrix where the entries of the weights vector
        # are in the diagonal of a matrix (where all other values are 0)
        weights <- diag(weights)

        # We calculate W_i^P (x[i]), that is the matrix product of:
        # the transpose of e_1 times the inverse of the matrix product
        # of the transpose of the k-th X times the weights times the k-th X
        # times the transpose of the k-th X times the weights times the response
        # We then assign this to the k-th element of the result vector initialized before
        result[k] <- t(e_1) %*% solve(t(X[[k]]) %*% weights %*% X[[k]]) %*% t(X[[k]]) %*% weights %*% Y
    }
    
    ## RETURNING VALUES
    # We return the estimator evaluated at the vector x
    return(result)
}

########################### PART B
# function to estimate
m = function(x) ((x-1)^2)

# bandwidth
h = 0.5

# data simulation
pred <- rnorm(500, mean=1, sd=1) # simulating 500 normally distributed values with mean = 1 and std = 1
resp <- c() # creating an empty vector to be appended values for the response
for (i in 1:length(pred)) {
    # generating the response adding the epsilon value
    # the epsilon error value will be a randomly generated
    # normally distributed value with mean 0 and std = 0.5
    resp <- c(resp, m(pred[i])) + rnorm(1, mean=0, sd=0.5)
}

# appending to list object
simulated_dataset <- matrix(,nrow=length(pred),ncol=2) # generating empty matrix with nrow = lenght of predictor vector
simulated_dataset[,1] <- pred # adding the predictors to column 1 of the matrix
simulated_dataset[,2] <- resp # adding the response to column 2 of the matrix

# Running the custom implementation
test_sample <- rnorm(500, mean=1, sd=1) # we generate 500 random normally distributed values with mean = 1 and sd = 1 to test
imp <- lce(x=test_sample, data=simulated_dataset, h=h) # we run the implementation using the dataset and bandwidth supplied

# plotting the points to trace the line
plot(test_sample, imp)

# generate grid to plot and test the accuracy of the implementation
x_grid <- seq(min(pred),max(pred), l=500)
plot(pred,resp) # plot the predictors vs the response

# tracing the regresion and custom local cubic estimator lines
lines(x_grid, m(x_grid), col=1)
lines(x_grid, lce(x=x_grid, data=simulated_dataset, h=h), col=2)
legend("top", legend = c("True reg","LocPoly (3rd deg)"), lwd=2, col=1:2) 

###################################################################################
###################################################################################
## CATEGORY B, PROBLEM 4: Exercise 4.9 ##
########################### PART A
# Importing necessary libraries
library(ks)
library(dplyr)
library(ggplot2)

# loading the data
load(url('https://raw.githubusercontent.com/egarpor/handy/master/datasets/ovals.RData'))

ovals$labels = as.factor(ovals$labels) # converting the labels column to factor
train <- ovals[1:2000,] # splitting the dataset into train 
test <- ovals[2001:3000,] # and test set

# plotting the grid
ggplot(data = train, aes(x.1, x.2, color = labels)) +
  geom_point() +   
  scale_color_manual(values = c("1" = "red", "2" = "blue", "3" = "yellow"))


########################### PART B
# Filtering x1, x2 and x3 by their respective grouping (determined by the labels column)
x1 <- ovals %>% filter(labels=="1") %>%
  select(x.1,x.2)
x2 <- ovals %>% filter(labels=="2") %>%
  select(x.1,x.2)
x3 <- ovals %>% filter(labels=="3") %>%
  select(x.1,x.2)

# Computing the plug in bandwidth for each subset (x1, x2 and x3)
h1 <- ks::Hpi(x1)
h2 <- ks::Hpi(x2)
h3 <- ks::Hpi(x3)

cont <- seq(0, 0.05, l = 20) # Defining a grid of 20 elements between 0 and 0.05
col <- viridis::viridis

par(mfrow = c(1, 3)) # creating a figure grid of 1x3
# plotting the KDEs for each filtered dataset utilizing the plug-in bandwidths
plot(ks::kde(x = x1, H = h1), display = "filled.contour2",
     abs.cont = cont, col.fun = col, main = "Plug in bandwidth for group1")
plot(ks::kde(x = x2, H = h2), display = "filled.contour2",
     abs.cont = cont, col.fun = col, main = "Plug in bandwidth for group2")
plot(ks::kde(x = x3, H = h3), display = "filled.contour2",
     abs.cont = cont, col.fun = col, main = "Plug in bandwidth for group3")

# Obtaining diagonal plug in bandwidths
h1_diag <- ks::Hpi.diag(x1)
h2_diag <- ks::Hpi.diag(x2)
h3_diag <- ks::Hpi.diag(x3)

par(mfrow = c(1, 3)) # creating a figure grid of 1x3
# plotting the KDEs for each filtered dataset utilizing the diagonal plug-in bandwidths
plot(ks::kde(x = x1, H = h1_diag), display = "filled.contour2",
     abs.cont = cont, col.fun = col, main = "Diagonal Plug in bandwidth for group1")
plot(ks::kde(x = x2, H = h2_diag), display = "filled.contour2",
     abs.cont = cont, col.fun = col, main = "Diagonal Plug in bandwidth for group2")
plot(ks::kde(x = x3, H = h3_diag), display = "filled.contour2",
     abs.cont = cont, col.fun = col, main = "Diagonal Plug in bandwidth for group3")

########################### PART C
x <- train[,1:2]
groups<-train$labels

# Hpi bandwidths are computed
kda <- ks::kda(x=x, x.group = groups) 

########################### PART D
par(mfrow = c(1, 1)) # Fixing the grid for individual plots
# We plot the KDA with the parameters as follows
plot(kda, 
     col = rainbow(3), 
     lwd = 2, 
     col.pt = 1, 
     cont = seq(5, 85, by = 20), 
     col.part = rainbow(3, alpha = 0.25), 
     drawpoints = TRUE)

########################### PART E
new_data=test[1:2] # dividing the test set to predict
pred_kda <- predict(kda, x = new_data) # predicting over the new_data section of test
ks::compare(x.group = test$labels, est.group = pred_kda)$cross # confusion matrix for the KDA model

########################### PART F
lda_model = MASS::lda(x, grouping = groups, data=train) # running LDA model 
pred_lda = predict(lda_model, newdata = new_data)$class # predicting using LDA
ks::compare(test$labels, pred_lda) # confusion matrix for LDA

########################### PART G
qda_model = MASS::qda(x, grouping = groups, data=train) # running QDA model
pred_qda = predict(qda_model, newdata =new_data)$class # predicting using QDA
ks::compare(test$labels, pred_qda) # confusion matrix for QDA