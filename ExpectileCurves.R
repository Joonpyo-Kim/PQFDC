# Required Packages 
library(fda)

# Preliminarily required variables
domain = seq(from=0, to=1, by=0.001) # domain that function will be defined 
obs.domain = seq(from=0, to=1, by=0.01) # observed time points (t_{ij}) 
num.knots <- 20 # The number of knots 
knots <- seq(from=0, to=1, by=1/num.knots) # vector of knots (including endpoints) 
b.basis <- bsplineS(obs.domain, knots) # Design matrix of B-spline bases (with knots)
b.expand <- bsplineS(domain, knots) # B-spline basis matrix at domain 

expec.curve <- function(y, tau, pen=FALSE, lambda=0){
	# y : observed data points (y_{ij}) (on time points (t_{ij})) 
	# tau : expectile level 
	# pen : a logical argument whether to use penalized spline. If TRUE, ridge penalized P-spline fit is used 
	# lambda : regularization parameter 
  coeff <- rep(0, num.knots+3)
  n <- length(y)
  if(!pen){
    weight.new <- rep(tau, n); weight <- vector(length = n)
    iter <- 1
    while(sum((weight-weight.new)^2)!=0){ # iteratively re-weighted least squares: until the weight converges 
      weight <- weight.new
      #print(paste(iter, "th iteration"))
      lmfit <- lm(y ~ b.basis - 1, weights = weight)
      weight.new <- tau * (y > lmfit$fitted.values) + (1-tau) * (y <= lmfit$fitted.values)
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
    iter <- 1
    # D is a matrix of difference of B-spline bases
    D <- matrix(0, nrow=num.knots+1, ncol=num.knots+3)
    index <- c(1:(num.knots+1))
    D[cbind(index, index)] <- 1; D[cbind(index, index+1)] <- -2; D[cbind(index, index+2)] <- 1
    while(sum((weight-weight.new)^2)!=0){
      weight <- weight.new
      #print(paste(iter, "th iteration"))
	  W <- diag(weight)
	  coeff <- solve(t(b.basis) %*% W %*% b.basis+lambda*t(D)%*%D, t(b.basis) %*% W %*% y)
      iter <- iter + 1 
      y.fit <- b.basis %*% coeff
      weight.new <- as.vector(tau * (y > y.fit) + (1-tau) * (y <= y.fit))
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

gcv.select <- function(y, tau, lambda=exp(seq(-2, 5, length=100))){
	# gcv.select function returns the best regularization parameter value based on GCV. 
	# y : observed data points (y_{ij}) (on time points (t_{ij})) 
	# tau : expectile level 
	# lambda : vector of candidates of lambda values 
	gcv.result <- vector(length=length(lambda))
	for(i in 1:length(lambda)){
		lamb <- lambda[i]
		exp.fit <- expec.curve(y, tau, pen=T, lambda=lamb)
		y.fit <- b.basis %*% exp.fit$coeff
		W <- diag(exp.fit$weight)
		S.mat <- b.basis %*% solve(t(b.basis) %*% W %*% b.basis+lamb*t(D)%*%D) %*% t(b.basis) %*% W 
		gcv.result[i] <- sum(exp.fit$weight*(y-y.fit)^2)/(length(y)-sum(diag(S.mat)))^2
	}
#		print(gcv.result)
	return(lambda[which.min(gcv.result)])
}

expec.curve.pspline <- function(y, tau){
	# expec.curve.pspline function returns P-spline-based expectile fit. 
	# y : observed data points (y_{ij}) (on time points (t_{ij})) 
	# tau : expectile level 
	lambda <- gcv.select(y, tau)
	coeff <- expec.curve(y, tau, pen=T, lambda=lambda)$coeff
	return(coeff)
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

coeff.tau0.5 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.5)
coeff.tau0.1 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.1)
coeff.tau0.9 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.9)
fit.expec0.5 <- b.expand%*%coeff.tau0.5
fit.expec0.1 <- b.expand%*%coeff.tau0.1
fit.expec0.9 <- b.expand%*%coeff.tau0.9

	i1 <- 4; i2 <- 8
	par(mfrow=c(1,2))
	plot(obs.domain, simul.norm[,i1], type="o", ylim=range(simul.norm[,c(i1,i2)]), xlab="domain", ylab="data level", main="Data with Positive Skewness")
	lines(obs.domain, true1, col="blue", lwd=1.5)
	lines(domain, fit.expec0.5[,i1], col="red")
	lines(domain, fit.expec0.1[,i1], col="violet")
	lines(domain, fit.expec0.9[,i1], col="orange")
	legend("bottomright", c("0.1th expectile", "0.5th expectile", "0.9th expectile"), lty=1, col=c("violet", "red", "orange"))
	plot(obs.domain, simul.norm[,i2], type="o", ylim=range(simul.norm[,c(i1,i2)]), xlab="domain", ylab="data level", main="Data with Negative Skewness")
	lines(obs.domain, true1, col="blue", lwd=1.5)
	lines(domain, fit.expec0.5[,i2], col="red")
	lines(domain, fit.expec0.1[,i2], col="violet")
	lines(domain, fit.expec0.9[,i2], col="orange")    
	legend("bottomright", c("0.1th expectile", "0.5th expectile", "0.9th expectile"), lty=1, col=c("violet", "red", "orange"))
