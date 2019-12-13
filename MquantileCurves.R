# Required Packages 
library(fda)

# Preliminarily required variables
domain = seq(from=0, to=1, by=0.001) # domain that function will be defined 
obs.domain = seq(from=0, to=1, by=0.01) # observed time points (t_{ij}) 
num.knots <- 20 # The number of knots 
knots <- seq(from=0, to=1, by=1/num.knots) # vector of knots (including endpoints) 
b.basis <- bsplineS(obs.domain, knots) # Design matrix of B-spline bases (with knots)
b.expand <- bsplineS(domain, knots) # B-spline basis matrix at domain 

rho <- function(x, tau, c=1.345){ # Asymmetric Huber loss function
	# tau : skewness paramaeter
	# c : cutoff value 
  if(x>c){
    return(tau*c*(2*x-c))
  }else if(x>0){
    return(tau*x^2)
  }else if(x>(-c)){
    return((1-tau)*x^2)
  }else{
    return(-(1-tau)*c*(2*x+c))
  }
}

psi <- function(x, tau, c=1.345){ # Derivative of rho 
	# tau : skewness paramaeter
	# c : cutoff value 
  if(x>c){
    return(2*tau*c)
  }else if(x>0){
    return(2*tau*x)
  }else if(x>(-c)){
    return(2*(1-tau)*x)
  }else{
    return(-2*(1-tau)*c)
  }
}

mquant.curve <- function(y, tau, c=1.345, pen=FALSE, lambda=0){
	# y : observed data points (y_{ij}) (on time points (t_{ij})) 
	# tau : expectile level 
	# c : cutoff value
	# pen : a logical argument whether to use penalized spline. If TRUE, ridge penalized P-spline fit is used 
	# lambda : regularization parameter 
  coeff <- rep(0, num.knots+3)
  n <- length(y)
  if(!pen){
    weight.new <- rep(tau, n); weight <- vector(length = n)
    sigma <- sqrt(mean(lm(y ~ b.basis -1)$residuals^2))
    iter <- 1
    while(sum((weight-weight.new)^2)>1e-6){ # iteratively re-weighted least squares: until the weight converges 
      weight <- weight.new
      #print(paste(iter, "th iteration"))
      lmfit <- lm(y ~ b.basis - 1, weights = weight)
      sigma <- sqrt(mean(lmfit$residuals^2))
      weight.new <- sapply(lmfit$residuals/sigma, psi, tau=tau, c=c)/(lmfit$residuals)
      iter <- iter + 1 
      if(iter == 500){ # maximum iteration number : to be customized
        print(paste("procedure does not converge")); 
        break; 
      }
      if(!is.na(coeff[1])){
        coeff <- lmfit$coefficients
      }
    }
  }else{
    weight.new <- rep(tau, n); weight <- vector(length = n)
    sigma <- sqrt(mean(lm(y ~ b.basis -1)$residuals^2))
    iter <- 1
     # D is a matrix of difference of B-spline bases
    D <- matrix(0, nrow=num.knots+1, ncol=num.knots+3)
    index <- c(1:(num.knots+1))
    D[cbind(index, index)] <- 1; D[cbind(index, index+1)] <- -2; D[cbind(index, index+2)] <- 1
    while(sum((weight-weight.new)^2)>1e-6){
      weight <- weight.new
      #print(paste(iter, "th iteration"))
      W <- diag(weight)
      coeff <- solve(t(b.basis) %*% W %*% b.basis+2*sigma*lambda*t(D)%*%D, t(b.basis) %*% W %*% y)
      iter <- iter + 1 
      y.fit <- b.basis %*% coeff
      sigma <- sqrt(mean((y-y.fit)^2))
      weight.new <- as.vector(sapply((y-y.fit)/sigma, psi, tau=tau, c=c)/(y-y.fit))
      if(iter == 500){ # maximum iteration number : to be customized
        print(paste("procedure does not converge")); 
        break; 
      }
    }
    coeff <- list(coeff=coeff, weight=weight)
  }
  return(coeff)
}
D <- matrix(0, nrow=num.knots+1, ncol=num.knots+3)
index <- c(1:(num.knots+1))
D[cbind(index, index)] <- 1; D[cbind(index, index+1)] <- -2; D[cbind(index, index+2)] <- 1

gcv.select.m <- function(y, tau, c, lambda=exp(seq(-3, 5, length=100))){
	# gcv.select.m function returns the best regularization parameter value based on GCV. 
	# y : observed data points (y_{ij}) (on time points (t_{ij})) 
	# tau : expectile level 
	# c : cutoff value
	# lambda : vector of candidates of lambda values 
  gcv.result <- vector(length=length(lambda))
  for(i in 1:length(lambda)){
    lamb <- lambda[i]
    mquant.fit <- mquant.curve(y, tau, c=c, pen=T, lambda=lamb)
    y.fit <- b.basis %*% mquant.fit$coeff
    sigma <- sqrt(mean((y-y.fit)^2))
    W <- diag(mquant.fit$weight)
    S.mat <- b.basis %*% solve(t(b.basis) %*% W %*% b.basis+2*sigma*lamb*t(D)%*%D) %*% t(b.basis) %*% W 
    gcv.result[i] <- sum(mquant.fit$weight*(y-y.fit)^2)/(length(y)-sum(diag(S.mat)))^2
  }
  #print(gcv.result)
  return(lambda[which.min(gcv.result)])
}

cv.select.m <- function(y, tau, c,  lambda=exp(seq(-3, 5, length=100))){
	# cv.select.m function returns the best regularization parameter value based on CV. 
	# y : observed data points (y_{ij}) (on time points (t_{ij})) 
	# tau : expectile level 
	# c : cutoff value
	# lambda : vector of candidates of lambda values 
  cv.result <- vector(length=length(lambda))
  for(i in 1:length(lambda)){
    lamb <- lambda[i]
    mquant.fit <- mquant.curve(y, tau, c,  pen=T, lambda=lamb)
    y.fit <- b.basis %*% mquant.fit$coeff
    sigma <- sqrt(mean((y-y.fit)^2))
    W <- diag(mquant.fit$weight)
    S.mat <- diag(sqrt(diag(W))) %*% b.basis %*% solve(t(b.basis) %*% W %*% b.basis+2*sigma*lamb*t(D)%*%D) %*% t(b.basis) %*% diag(sqrt(diag(W))) 
    cv.result[i] <- sum(mquant.fit$weight*(y-y.fit)^2/(1-diag(S.mat))^2)
  }
  #print(cv.result)
  return(lambda[which.min(cv.result)])
}

mquant.curve.pspline <- function(y, tau, c=1.345, gcv=TRUE){
	# mquant.curve.pspline function returns P-spline-based expectile fit. 
	# y : observed data points (y_{ij}) (on time points (t_{ij})) 
	# tau : expectile level 
	# gcv : a logical argument. If TRUE, GCV is used to select lambda; otherwise, CV is used. 
  if(!gcv){
    lambda <- cv.select.m(y, tau, c)
    coeff <- mquant.curve(y, tau, c, pen=T, lambda=lambda)$coeff
    return(coeff)
  }else{
    lambda <- gcv.select.m(y, tau, c)
    coeff <- mquant.curve(y, tau, c, pen=T, lambda=lambda)$coeff
    return(coeff)
  }
}

# Example 
set.seed(1213)
library(lqmm)
N <- 200
# simul.norm is a matrix of N discretely observed functional data (y_{ij}) on time points (t_{ij}) 
simul.norm <- matrix(nrow=length(obs.domain), ncol=N)
true1 <- 2-cos(pi*obs.domain) # true mean function 
true.cluster <- sample(c(1,2), size=N, rep=T, prob=c(0.6, 0.4)) # true cluster alignment
poo <- function(s){return(ral(1, sigma=sig, tau=s)-sig*(1/s-1/(1-s)))} # Asymmetric Laplace-distributed error generator 
for(i in 1:N){
	if(true.cluster[i]==1){
		simul.norm[,i] <- true1+sapply(rep(0.1, length=101), poo)
	}else{
		simul.norm[,i] <- true1+sapply(rep(0.9, length=101), poo)
	}
}

coeff.tau0.5 <- apply(simul.norm, 2, mquant.curve.pspline, tau=0.5)
coeff.tau0.1 <- apply(simul.norm, 2, mquant.curve.pspline, tau=0.1)
coeff.tau0.9 <- apply(simul.norm, 2, mquant.curve.pspline, tau=0.9)
fit.mquant0.5 <- b.expand%*%coeff.tau0.5
fit.mquant0.1 <- b.expand%*%coeff.tau0.1
fit.mquant0.9 <- b.expand%*%coeff.tau0.9

	i1 <- 4; i2 <- 8
	par(mfrow=c(1,2))
	plot(obs.domain, simul.norm[,i1], type="o", ylim=range(simul.norm[,c(i1,i2)]), xlab="domain", ylab="data level", main="Data with Positive Skewness")
	lines(obs.domain, true1, col="blue", lwd=1.5)
	lines(domain, fit.mquant0.5[,i1], col="red")
	lines(domain, fit.mquant0.1[,i1], col="violet")
	lines(domain, fit.mquant0.9[,i1], col="orange")
	legend("bottomright", c("0.1th M-quantile", "0.5th M-quantile", "0.9th M-quantile"), lty=1, col=c("violet", "red", "orange"))
	plot(obs.domain, simul.norm[,i2], type="o", ylim=range(simul.norm[,c(i1,i2)]), xlab="domain", ylab="data level", main="Data with Negative Skewness")
	lines(obs.domain, true1, col="blue", lwd=1.5)
	lines(domain, fit.mquant0.5[,i2], col="red")
	lines(domain, fit.mquant0.1[,i2], col="violet")
	lines(domain, fit.mquant0.9[,i2], col="orange")    
	legend("bottomright", c("0.1th M-quantile", "0.5th M-quantile", "0.9th M-quantile"), lty=1, col=c("violet", "red", "orange"))
