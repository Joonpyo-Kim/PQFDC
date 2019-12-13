library(fda)
library(lqmm)
library(Funclustering)
library(mclust)
library(glmnet)
library(gam)

###1. skewness is constant of time
N <- 200
num.knots <- 20
set.seed(1221)
result.vec <- matrix(nrow=7, ncol=5)
result.vec.rand <- matrix(nrow=7, ncol=5)
for(sn in c(1,3,5)){
	sig <- sn*0.1
	domain = seq(from=0, to=1, by=0.001)
	obs.domain = seq(from=0, to=1, by=0.01)
	knots <- seq(from=0, to=1, by=1/num.knots)
	b.basis <- bsplineS(obs.domain, knots)
	b.expand <- bsplineS(domain, knots)
	simul.norm <- matrix(nrow=length(obs.domain), ncol=N)
	true1 <- 2-cos(pi*obs.domain)
	#true1 <- (2*obs.domain-1)^{3} 
	true.cluster <- sample(c(1,2), size=N, rep=T, prob=c(0.6, 0.4))
	poo <- function(s){return(ral(1, sigma=sig, tau=s)-sig*(1/s-1/(1-s)))}
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
	coeff.tau0.01 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.01)
	coeff.tau0.99 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.99)
	fit.expec0.5 <- b.expand%*%coeff.tau0.5
	fit.expec0.1 <- b.expand%*%coeff.tau0.1
	fit.expec0.9 <- b.expand%*%coeff.tau0.9
	fit.expec0.01 <- b.expand%*%coeff.tau0.01
	fit.expec0.99 <- b.expand%*%coeff.tau0.99
	
	
	#kmeans 
	#1-1) mean based (0.5 expectile)
	clust.obj <- kmeans(t(coeff.tau0.5), centers=2)
	tpr.clust <- vector(length=2)
	per <- permutations(2,2)
	for(i in 1:2){
		poo <- (per[i,])[clust.obj$cluster]
		tpr.clust[i] <- mean(poo == true.cluster)
	}
	result.vec[1,sn] <- max(tpr.clust) 
	result.vec.rand[1,sn] <- adjustedRandIndex(clust.obj$cluster, true.cluster)
		
	#2-6) weighted average   
    J <- 5
    w1 <- w2 <- matrix(nrow=N, ncol=J)
    weighted.expec <- matrix(nrow=length(domain), ncol=N)
    for(i2 in 1:N){
      for(i in 1:J){
        search.vec <- seq(from=floor(1000/J*(i-1))+1, to=floor(1000/J*i))
        w1[i2,i] <- abs(0.001*sum(fit.expec0.9[search.vec, i2]-fit.expec0.5[search.vec, i2]))
        w2[i2,i] <- abs(0.001*sum(fit.expec0.5[search.vec, i2]-fit.expec0.1[search.vec, i2]))
      }
      #w1[i2,] <- 1/w1[i2,]; w2[i2,] <- 1/w2[i2,]
      w11 <- w1[i2,]/(w1[i2,]+w2[i2,]); w22 <- w2[i2,]/(w1[i2,]+w2[i2,]); w1[i2,] <- w11; w2[i2,] <- w22
      poo1 <- w1[i2,1]; poo2 <- w2[i2,1]
      for(i in 1:J){
        search.vec <- seq(from=floor(1000/J*(i-1))+1, to=floor(1000/J*i))
        poo1 <- c(poo1, rep(w1[i2,i], length(search.vec)))
        poo2 <- c(poo2, rep(w2[i2,i], length(search.vec)))
      }
      weighted.expec[,i2] <-poo1*fit.expec0.9[,i2]+poo2*fit.expec0.1[,i2]
    }
    coeff.tau.al <- matrix(nrow=N, ncol=num.knots+3)
    for(i in 1:N){
      coeff.tau.al[i,] <- solve(t(b.expand)%*%b.expand, t(b.expand)%*%weighted.expec[,i])
    }    
	clust.obj <- kmeans(coeff.tau.al, centers=2)
    per <- permutations(2,2)
	for(i in 1:2){
		poo <- (per[i,])[clust.obj$cluster]
		tpr.clust[i] <- mean(poo == true.cluster)
	}
	result.vec[2,sn] <- max(tpr.clust) 
	result.vec.rand[2,sn] <- adjustedRandIndex(clust.obj$cluster, true.cluster)

    
	#2-6-2) weighted average: 0.01, 0.99   
   J <- 5
    w1 <- w2 <- matrix(nrow=N, ncol=J)
    weighted.expec <- matrix(nrow=length(domain), ncol=N)
    for(i2 in 1:N){
      for(i in 1:J){
        search.vec <- seq(from=floor(1000/J*(i-1))+1, to=floor(1000/J*i))
        w1[i2,i] <- abs(0.001*sum(fit.expec0.99[search.vec, i2]-fit.expec0.5[search.vec, i2]))
        w2[i2,i] <- abs(0.001*sum(fit.expec0.5[search.vec, i2]-fit.expec0.01[search.vec, i2]))
      }
      #w1[i2,] <- 1/w1[i2,]; w2[i2,] <- 1/w2[i2,]
      w11 <- w1[i2,]/(w1[i2,]+w2[i2,]); w22 <- w2[i2,]/(w1[i2,]+w2[i2,]); w1[i2,] <- w11; w2[i2,] <- w22
      poo1 <- w1[i2,1]; poo2 <- w2[i2,1]
      for(i in 1:J){
        search.vec <- seq(from=floor(1000/J*(i-1))+1, to=floor(1000/J*i))
        poo1 <- c(poo1, rep(w1[i2,i], length(search.vec)))
        poo2 <- c(poo2, rep(w2[i2,i], length(search.vec)))
      }
      weighted.expec[,i2] <-poo1*fit.expec0.99[,i2]+poo2*fit.expec0.01[,i2]
    }
    coeff.tau.al <- matrix(nrow=N, ncol=num.knots+3)
    for(i in 1:N){
      coeff.tau.al[i,] <- solve(t(b.expand)%*%b.expand, t(b.expand)%*%weighted.expec[,i])
    }    
	clust.obj <- kmeans(coeff.tau.al, centers=2)
    per <- permutations(2,2)
	for(i in 1:2){
		poo <- (per[i,])[clust.obj$cluster]
		tpr.clust[i] <- mean(poo == true.cluster)
	}
	result.vec[3,sn] <- max(tpr.clust) 
	result.vec.rand[3,sn] <- adjustedRandIndex(clust.obj$cluster, true.cluster)
	
	
	#2-6-5) weighted average: error distribution
	w <- matrix(nrow=2, ncol=N)
	coeff.taus <- array(dim=c(2, num.knots+3, N))
	coeff.taus[1,,] <- coeff.tau0.1
	coeff.taus[2,,] <- coeff.tau0.9
	weighted.expec <- matrix(nrow=length(obs.domain), ncol=N)
	residual <- matrix(nrow=length(obs.domain), ncol=N)
	for(i in 1:N){
		residual[,i] <- simul.norm[,i]-(b.basis %*% coeff.tau0.5)[,i]
		den <- density(residual[,i])
		w[1,i] <-  approx(den$x, den$y, xout=quantile(residual[,i], probs=0.1))$y
		w[2,i] <-  approx(den$x, den$y, xout=quantile(residual[,i], probs=0.9))$y
	}
	coeff.tau.al <- matrix(nrow=N, ncol=num.knots+3)
	for(i in 1:N){
	w.avg <- function(v){
		return(sum(w[,i]*v)/sum(w[,i]))
	}		
		coeff.tau.al[i,] <- apply(coeff.taus[,,i], 2, w.avg)
	}     
	clust.obj <- kmeans(coeff.tau.al, centers=2)
	tpr.clust <- vector(length=2)
	per <- permutations(2,2)
	for(i in 1:2){
		poo <- (per[i,])[clust.obj$cluster]
		tpr.clust[i] <- mean(poo == true.cluster)
	}
	result.vec[4,sn] <- max(tpr.clust) 
	result.vec.rand[4,sn] <- adjustedRandIndex(clust.obj$cluster, true.cluster)
	
	
	#2-6-6) weighted average: error distribution
	w <- matrix(nrow=2, ncol=N)
	coeff.taus <- array(dim=c(2, num.knots+3, N))
	coeff.taus[1,,] <- coeff.tau0.01
	coeff.taus[2,,] <- coeff.tau0.99
	weighted.expec <- matrix(nrow=length(obs.domain), ncol=N)
	residual <- matrix(nrow=length(obs.domain), ncol=N)
	for(i in 1:N){
		residual[,i] <- simul.norm[,i]-(b.basis %*% coeff.tau0.5)[,i]
		den <- density(residual[,i])
		w[1,i] <-  approx(den$x, den$y, xout=quantile(residual[,i], probs=0.01))$y
		w[2,i] <-  approx(den$x, den$y, xout=quantile(residual[,i], probs=0.99))$y
	}
	coeff.tau.al <- matrix(nrow=N, ncol=num.knots+3)
	for(i in 1:N){
	w.avg <- function(v){
		return(sum(w[,i]*v)/sum(w[,i]))
	}		
		coeff.tau.al[i,] <- apply(coeff.taus[,,i], 2, w.avg)
	}     
	clust.obj <- kmeans(coeff.tau.al, centers=2)
	tpr.clust <- vector(length=2)
	per <- permutations(2,2)
	for(i in 1:2){
		poo <- (per[i,])[clust.obj$cluster]
		tpr.clust[i] <- mean(poo == true.cluster)
	}
	result.vec[5,sn] <- max(tpr.clust) 
	result.vec.rand[5,sn] <- adjustedRandIndex(clust.obj$cluster, true.cluster)
	
	
	#2-8)funclust 
	bsp <- create.bspline.basis(breaks=seq(0, 1, by=1/num.knots))
	fd.obj <- fd(coef=coeff.tau0.5, basisobj=bsp)
	tpr <- vector(length=100)
	tpr2 <- vector(length=100)
		clust.obj <- funclust(fd.obj, K=2)
		tpr.clust <- vector(length=2)
		per <- permutations(2,2)
		for(i in 1:2){
			poo <- (per[i,])[clust.obj$cls]
			tpr.clust[i] <- mean(poo == true.cluster)
		}
	result.vec[6,sn] <- max(tpr.clust)
	result.vec.rand[6,sn] <- adjustedRandIndex(clust.obj$cls, true.cluster)
	
	#2-9)curveclust
	simul.norm.curv <- matrix(nrow=128, ncol=200)
	simul.norm.list <- list()
	simul.norm.curv[14:114,] <- simul.norm
	for(i in 1:13){
		simul.norm.curv[i,] <- simul.norm.curv[28-i,]
	}
	for(i in 115:128){
		simul.norm.curv[i,] <- simul.norm.curv[228-i,]
	}
	for(i in 1:200){
		simul.norm.list[[i]] <- simul.norm.curv[,i]
	}
	CCD <- apply(simul.norm.curv, 2, FUN = function(y) {
            wd(y, filter.number = 1, family = "DaubExPhase")
        })
   CCD <- new("CClustData", Y=simul.norm.list, filter.number=1)
   CCDred <- getUnionCoef(CCD)
   CCO            = new("CClustO")
   CCO["nbclust"] = 2
   CCR <- getFCMM(CCDred, CCO)
   summary(CCR)
   groups <- apply(CCR["Tau"], 1, which.max)
   	for(i in 1:2){
			tpr.clust[i] <- mean(groups == true.cluster)
		}
	max(tpr.clust)
	result.vec[7,sn] <- max(tpr.clust)
	result.vec.rand[7,sn] <- adjustedRandIndex(groups, true.cluster)

}

###2. skewness is not extreme & constant of time
N <- 200
num.knots <- 20
set.seed(1221)
result.vec <- matrix(nrow=7, ncol=5)
result.vec.rand <- matrix(nrow=7, ncol=5)
for(sn in c(1,3,5)){
	sig <- sn*0.1
	domain = seq(from=0, to=1, by=0.001)
	obs.domain = seq(from=0, to=1, by=0.01)
	knots <- seq(from=0, to=1, by=1/num.knots)
	b.basis <- bsplineS(obs.domain, knots)
	b.expand <- bsplineS(domain, knots)
	simul.norm <- matrix(nrow=length(obs.domain), ncol=N)
	true1 <- 2-cos(pi*obs.domain)
	#true1 <- (2*obs.domain-1)^{3} 
	true.cluster <- sample(c(1,2), size=N, rep=T, prob=c(0.6, 0.4))
	poo <- function(s){return(ral(1, sigma=sig, tau=s)-0.2*(1/s-1/(1-s)))}
	for(i in 1:N){
		if(true.cluster[i]==1){
			simul.norm[,i] <- true1+sapply(rep(0.3, length=101), poo)
		}else{
			simul.norm[,i] <- true1+sapply(rep(0.7, length=101), poo)
		}
	}
	
	coeff.tau0.5 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.5)
	coeff.tau0.1 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.1)
	coeff.tau0.9 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.9)
	coeff.tau0.01 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.01)
	coeff.tau0.99 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.99)
	fit.expec0.5 <- b.expand%*%coeff.tau0.5
	fit.expec0.1 <- b.expand%*%coeff.tau0.1
	fit.expec0.9 <- b.expand%*%coeff.tau0.9
	fit.expec0.01 <- b.expand%*%coeff.tau0.01
	fit.expec0.99 <- b.expand%*%coeff.tau0.99
	
	
	#kmeans 
	#1-1) mean based (0.5 expectile)
	clust.obj <- kmeans(t(coeff.tau0.5), centers=2)
	tpr.clust <- vector(length=2)
	per <- permutations(2,2)
	for(i in 1:2){
		poo <- (per[i,])[clust.obj$cluster]
		tpr.clust[i] <- mean(poo == true.cluster)
	}
	result.vec[1,sn] <- max(tpr.clust) 
	result.vec.rand[1,sn] <- adjustedRandIndex(clust.obj$cluster, true.cluster)
		
	#2-6) weighted average   
    J <- 5
    w1 <- w2 <- matrix(nrow=N, ncol=J)
    weighted.expec <- matrix(nrow=length(domain), ncol=N)
    for(i2 in 1:N){
      for(i in 1:J){
        search.vec <- seq(from=floor(1000/J*(i-1))+1, to=floor(1000/J*i))
        w1[i2,i] <- abs(0.001*sum(fit.expec0.9[search.vec, i2]-fit.expec0.5[search.vec, i2]))
        w2[i2,i] <- abs(0.001*sum(fit.expec0.5[search.vec, i2]-fit.expec0.1[search.vec, i2]))
      }
      #w1[i2,] <- 1/w1[i2,]; w2[i2,] <- 1/w2[i2,]
      w11 <- w1[i2,]/(w1[i2,]+w2[i2,]); w22 <- w2[i2,]/(w1[i2,]+w2[i2,]); w1[i2,] <- w11; w2[i2,] <- w22
      poo1 <- w1[i2,1]; poo2 <- w2[i2,1]
      for(i in 1:J){
        search.vec <- seq(from=floor(1000/J*(i-1))+1, to=floor(1000/J*i))
        poo1 <- c(poo1, rep(w1[i2,i], length(search.vec)))
        poo2 <- c(poo2, rep(w2[i2,i], length(search.vec)))
      }
      weighted.expec[,i2] <-poo1*fit.expec0.9[,i2]+poo2*fit.expec0.1[,i2]
    }
    coeff.tau.al <- matrix(nrow=N, ncol=num.knots+3)
    for(i in 1:N){
      coeff.tau.al[i,] <- solve(t(b.expand)%*%b.expand, t(b.expand)%*%weighted.expec[,i])
    }    
	clust.obj <- kmeans(coeff.tau.al, centers=2)
    per <- permutations(2,2)
	for(i in 1:2){
		poo <- (per[i,])[clust.obj$cluster]
		tpr.clust[i] <- mean(poo == true.cluster)
	}
	result.vec[2,sn] <- max(tpr.clust) 
	result.vec.rand[2,sn] <- adjustedRandIndex(clust.obj$cluster, true.cluster)

    
	#2-6-2) weighted average: 0.01, 0.99   
   J <- 5
    w1 <- w2 <- matrix(nrow=N, ncol=J)
    weighted.expec <- matrix(nrow=length(domain), ncol=N)
    for(i2 in 1:N){
      for(i in 1:J){
        search.vec <- seq(from=floor(1000/J*(i-1))+1, to=floor(1000/J*i))
        w1[i2,i] <- abs(0.001*sum(fit.expec0.99[search.vec, i2]-fit.expec0.5[search.vec, i2]))
        w2[i2,i] <- abs(0.001*sum(fit.expec0.5[search.vec, i2]-fit.expec0.01[search.vec, i2]))
      }
      #w1[i2,] <- 1/w1[i2,]; w2[i2,] <- 1/w2[i2,]
      w11 <- w1[i2,]/(w1[i2,]+w2[i2,]); w22 <- w2[i2,]/(w1[i2,]+w2[i2,]); w1[i2,] <- w11; w2[i2,] <- w22
      poo1 <- w1[i2,1]; poo2 <- w2[i2,1]
      for(i in 1:J){
        search.vec <- seq(from=floor(1000/J*(i-1))+1, to=floor(1000/J*i))
        poo1 <- c(poo1, rep(w1[i2,i], length(search.vec)))
        poo2 <- c(poo2, rep(w2[i2,i], length(search.vec)))
      }
      weighted.expec[,i2] <-poo1*fit.expec0.99[,i2]+poo2*fit.expec0.01[,i2]
    }
    coeff.tau.al <- matrix(nrow=N, ncol=num.knots+3)
    for(i in 1:N){
      coeff.tau.al[i,] <- solve(t(b.expand)%*%b.expand, t(b.expand)%*%weighted.expec[,i])
    }    
	clust.obj <- kmeans(coeff.tau.al, centers=2)
    per <- permutations(2,2)
	for(i in 1:2){
		poo <- (per[i,])[clust.obj$cluster]
		tpr.clust[i] <- mean(poo == true.cluster)
	}
	result.vec[3,sn] <- max(tpr.clust) 
	result.vec.rand[3,sn] <- adjustedRandIndex(clust.obj$cluster, true.cluster)
	
	
	#2-6-5) weighted average: error distribution
	w <- matrix(nrow=2, ncol=N)
	coeff.taus <- array(dim=c(2, num.knots+3, N))
	coeff.taus[1,,] <- coeff.tau0.1
	coeff.taus[2,,] <- coeff.tau0.9
	weighted.expec <- matrix(nrow=length(obs.domain), ncol=N)
	residual <- matrix(nrow=length(obs.domain), ncol=N)
	for(i in 1:N){
		residual[,i] <- simul.norm[,i]-(b.basis %*% coeff.tau0.5)[,i]
		den <- density(residual[,i])
		w[1,i] <-  approx(den$x, den$y, xout=quantile(residual[,i], probs=0.1))$y
		w[2,i] <-  approx(den$x, den$y, xout=quantile(residual[,i], probs=0.9))$y
	}
	coeff.tau.al <- matrix(nrow=N, ncol=num.knots+3)
	for(i in 1:N){
	w.avg <- function(v){
		return(sum(w[,i]*v)/sum(w[,i]))
	}		
		coeff.tau.al[i,] <- apply(coeff.taus[,,i], 2, w.avg)
	}     
	clust.obj <- kmeans(coeff.tau.al, centers=2)
	tpr.clust <- vector(length=2)
	per <- permutations(2,2)
	for(i in 1:2){
		poo <- (per[i,])[clust.obj$cluster]
		tpr.clust[i] <- mean(poo == true.cluster)
	}
	result.vec[4,sn] <- max(tpr.clust) 
	result.vec.rand[4,sn] <- adjustedRandIndex(clust.obj$cluster, true.cluster)
	
	
	#2-6-6) weighted average: error distribution
	w <- matrix(nrow=2, ncol=N)
	coeff.taus <- array(dim=c(2, num.knots+3, N))
	coeff.taus[1,,] <- coeff.tau0.01
	coeff.taus[2,,] <- coeff.tau0.99
	weighted.expec <- matrix(nrow=length(obs.domain), ncol=N)
	residual <- matrix(nrow=length(obs.domain), ncol=N)
	for(i in 1:N){
		residual[,i] <- simul.norm[,i]-(b.basis %*% coeff.tau0.5)[,i]
		den <- density(residual[,i])
		w[1,i] <-  approx(den$x, den$y, xout=quantile(residual[,i], probs=0.01))$y
		w[2,i] <-  approx(den$x, den$y, xout=quantile(residual[,i], probs=0.99))$y
	}
	coeff.tau.al <- matrix(nrow=N, ncol=num.knots+3)
	for(i in 1:N){
	w.avg <- function(v){
		return(sum(w[,i]*v)/sum(w[,i]))
	}		
		coeff.tau.al[i,] <- apply(coeff.taus[,,i], 2, w.avg)
	}     
	clust.obj <- kmeans(coeff.tau.al, centers=2)
	tpr.clust <- vector(length=2)
	per <- permutations(2,2)
	for(i in 1:2){
		poo <- (per[i,])[clust.obj$cluster]
		tpr.clust[i] <- mean(poo == true.cluster)
	}
	result.vec[5,sn] <- max(tpr.clust) 
	result.vec.rand[5,sn] <- adjustedRandIndex(clust.obj$cluster, true.cluster)
	
	
	#2-8)funclust 
	bsp <- create.bspline.basis(breaks=seq(0, 1, by=1/num.knots))
	fd.obj <- fd(coef=coeff.tau0.5, basisobj=bsp)
	tpr <- vector(length=100)
	tpr2 <- vector(length=100)
		clust.obj <- funclust(fd.obj, K=2)
		tpr.clust <- vector(length=2)
		per <- permutations(2,2)
		for(i in 1:2){
			poo <- (per[i,])[clust.obj$cls]
			tpr.clust[i] <- mean(poo == true.cluster)
		}
	result.vec[6,sn] <- max(tpr.clust)
	result.vec.rand[6,sn] <- adjustedRandIndex(clust.obj$cls, true.cluster)
	
	#2-9)curveclust
	simul.norm.curv <- matrix(nrow=128, ncol=200)
	simul.norm.list <- list()
	simul.norm.curv[14:114,] <- simul.norm
	for(i in 1:13){
		simul.norm.curv[i,] <- simul.norm.curv[28-i,]
	}
	for(i in 115:128){
		simul.norm.curv[i,] <- simul.norm.curv[228-i,]
	}
	for(i in 1:200){
		simul.norm.list[[i]] <- simul.norm.curv[,i]
	}
	CCD <- apply(simul.norm.curv, 2, FUN = function(y) {
            wd(y, filter.number = 1, family = "DaubExPhase")
        })
   CCD <- new("CClustData", Y=simul.norm.list, filter.number=1)
   CCDred <- getUnionCoef(CCD)
   CCO            = new("CClustO")
   CCO["nbclust"] = 2
   CCR <- getFCMM(CCDred, CCO)
   summary(CCR)
   groups <- apply(CCR["Tau"], 1, which.max)
   	for(i in 1:2){
			tpr.clust[i] <- mean(groups == true.cluster)
		}
	max(tpr.clust)
	result.vec[7,sn] <- max(tpr.clust)
	result.vec.rand[7,sn] <- adjustedRandIndex(groups, true.cluster)
	
}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               



### 3. skewness is varying for time
N <- 200
num.knots <- 20
set.seed(1221)
result.vec <- matrix(nrow=7, ncol=5)
result.vec.rand <- matrix(nrow=7, ncol=5)
for(sn in c(1,3,5)){
	sig <- sn*0.1
	domain = seq(from=0, to=1, by=0.001)
	obs.domain = seq(from=0, to=1, by=0.01)
	knots <- seq(from=0, to=1, by=1/num.knots)
	b.basis <- bsplineS(obs.domain, knots)
	b.expand <- bsplineS(domain, knots)
	simul.norm <- matrix(nrow=length(obs.domain), ncol=N)
	true1 <- 2-cos(pi*obs.domain)
	true.cluster <- sample(c(1,2), size=N, rep=T, prob=c(0.6, 0.4))
	poo <- function(s){return(ral(1, sigma=sig, tau=s)-0.2*(1/s-1/(1-s)))}
	tau.vec <- exp(seq(from=log(1/9), to=log(9), length=101))/(1+exp(seq(from=log(1/9), to=log(9), length=101)))
	for(i in 1:N){
		if(true.cluster[i]==1){
			simul.norm[,i] <- true1+sapply(tau.vec, poo)
		}else{
			simul.norm[,i] <- true1+sapply(rev(tau.vec), poo)
		}
	}

	
	coeff.tau0.5 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.5)
	coeff.tau0.1 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.1)
	coeff.tau0.9 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.9)
	coeff.tau0.01 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.01)
	coeff.tau0.99 <- apply(simul.norm, 2, expec.curve.pspline, tau=0.99)
	fit.expec0.5 <- b.expand%*%coeff.tau0.5
	fit.expec0.1 <- b.expand%*%coeff.tau0.1
	fit.expec0.9 <- b.expand%*%coeff.tau0.9
	fit.expec0.01 <- b.expand%*%coeff.tau0.01
	fit.expec0.99 <- b.expand%*%coeff.tau0.99
	
	
	#kmeans 
	#1-1) mean based (0.5 expectile)
	clust.obj <- kmeans(t(coeff.tau0.5), centers=2)
	tpr.clust <- vector(length=2)
	per <- permutations(2,2)
	for(i in 1:2){
		poo <- (per[i,])[clust.obj$cluster]
		tpr.clust[i] <- mean(poo == true.cluster)
	}
	result.vec[1,sn] <- max(tpr.clust) 
	result.vec.rand[1,sn] <- adjustedRandIndex(clust.obj$cluster, true.cluster)
		
	#2-6) weighted average   
    J <- 5
    w1 <- w2 <- matrix(nrow=N, ncol=J)
    weighted.expec <- matrix(nrow=length(domain), ncol=N)
    for(i2 in 1:N){
      for(i in 1:J){
        search.vec <- seq(from=floor(1000/J*(i-1))+1, to=floor(1000/J*i))
        w1[i2,i] <- abs(0.001*sum(fit.expec0.9[search.vec, i2]-fit.expec0.5[search.vec, i2]))
        w2[i2,i] <- abs(0.001*sum(fit.expec0.5[search.vec, i2]-fit.expec0.1[search.vec, i2]))
      }
      #w1[i2,] <- 1/w1[i2,]; w2[i2,] <- 1/w2[i2,]
      w11 <- w1[i2,]/(w1[i2,]+w2[i2,]); w22 <- w2[i2,]/(w1[i2,]+w2[i2,]); w1[i2,] <- w11; w2[i2,] <- w22
      poo1 <- w1[i2,1]; poo2 <- w2[i2,1]
      for(i in 1:J){
        search.vec <- seq(from=floor(1000/J*(i-1))+1, to=floor(1000/J*i))
        poo1 <- c(poo1, rep(w1[i2,i], length(search.vec)))
        poo2 <- c(poo2, rep(w2[i2,i], length(search.vec)))
      }
      weighted.expec[,i2] <-poo1*fit.expec0.9[,i2]+poo2*fit.expec0.1[,i2]
    }
    coeff.tau.al <- matrix(nrow=N, ncol=num.knots+3)
    for(i in 1:N){
      coeff.tau.al[i,] <- solve(t(b.expand)%*%b.expand, t(b.expand)%*%weighted.expec[,i])
    }    
	clust.obj <- kmeans(coeff.tau.al, centers=2)
    per <- permutations(2,2)
	for(i in 1:2){
		poo <- (per[i,])[clust.obj$cluster]
		tpr.clust[i] <- mean(poo == true.cluster)
	}
	result.vec[2,sn] <- max(tpr.clust) 
	result.vec.rand[2,sn] <- adjustedRandIndex(clust.obj$cluster, true.cluster)

    
	#2-6-2) weighted average: 0.01, 0.99   
   J <- 5
    w1 <- w2 <- matrix(nrow=N, ncol=J)
    weighted.expec <- matrix(nrow=length(domain), ncol=N)
    for(i2 in 1:N){
      for(i in 1:J){
        search.vec <- seq(from=floor(1000/J*(i-1))+1, to=floor(1000/J*i))
        w1[i2,i] <- abs(0.001*sum(fit.expec0.99[search.vec, i2]-fit.expec0.5[search.vec, i2]))
        w2[i2,i] <- abs(0.001*sum(fit.expec0.5[search.vec, i2]-fit.expec0.01[search.vec, i2]))
      }
      #w1[i2,] <- 1/w1[i2,]; w2[i2,] <- 1/w2[i2,]
      w11 <- w1[i2,]/(w1[i2,]+w2[i2,]); w22 <- w2[i2,]/(w1[i2,]+w2[i2,]); w1[i2,] <- w11; w2[i2,] <- w22
      poo1 <- w1[i2,1]; poo2 <- w2[i2,1]
      for(i in 1:J){
        search.vec <- seq(from=floor(1000/J*(i-1))+1, to=floor(1000/J*i))
        poo1 <- c(poo1, rep(w1[i2,i], length(search.vec)))
        poo2 <- c(poo2, rep(w2[i2,i], length(search.vec)))
      }
      weighted.expec[,i2] <-poo1*fit.expec0.99[,i2]+poo2*fit.expec0.01[,i2]
    }
    coeff.tau.al <- matrix(nrow=N, ncol=num.knots+3)
    for(i in 1:N){
      coeff.tau.al[i,] <- solve(t(b.expand)%*%b.expand, t(b.expand)%*%weighted.expec[,i])
    }    
	clust.obj <- kmeans(coeff.tau.al, centers=2)
    per <- permutations(2,2)
	for(i in 1:2){
		poo <- (per[i,])[clust.obj$cluster]
		tpr.clust[i] <- mean(poo == true.cluster)
	}
	result.vec[3,sn] <- max(tpr.clust) 
	result.vec.rand[3,sn] <- adjustedRandIndex(clust.obj$cluster, true.cluster)
	
	
	#2-6-5) weighted average: error distribution
	w <- matrix(nrow=2, ncol=N)
	coeff.taus <- array(dim=c(2, num.knots+3, N))
	coeff.taus[1,,] <- coeff.tau0.1
	coeff.taus[2,,] <- coeff.tau0.9
	weighted.expec <- matrix(nrow=length(obs.domain), ncol=N)
	residual <- matrix(nrow=length(obs.domain), ncol=N)
	for(i in 1:N){
		residual[,i] <- simul.norm[,i]-(b.basis %*% coeff.tau0.5)[,i]
		den <- density(residual[,i])
		w[1,i] <-  approx(den$x, den$y, xout=quantile(residual[,i], probs=0.1))$y
		w[2,i] <-  approx(den$x, den$y, xout=quantile(residual[,i], probs=0.9))$y
	}
	coeff.tau.al <- matrix(nrow=N, ncol=num.knots+3)
	for(i in 1:N){
	w.avg <- function(v){
		return(sum(w[,i]*v)/sum(w[,i]))
	}		
		coeff.tau.al[i,] <- apply(coeff.taus[,,i], 2, w.avg)
	}     
	clust.obj <- kmeans(coeff.tau.al, centers=2)
	tpr.clust <- vector(length=2)
	per <- permutations(2,2)
	for(i in 1:2){
		poo <- (per[i,])[clust.obj$cluster]
		tpr.clust[i] <- mean(poo == true.cluster)
	}
	result.vec[4,sn] <- max(tpr.clust) 
	result.vec.rand[4,sn] <- adjustedRandIndex(clust.obj$cluster, true.cluster)
	
	
	#2-6-6) weighted average: error distribution
	w <- matrix(nrow=2, ncol=N)
	coeff.taus <- array(dim=c(2, num.knots+3, N))
	coeff.taus[1,,] <- coeff.tau0.01
	coeff.taus[2,,] <- coeff.tau0.99
	weighted.expec <- matrix(nrow=length(obs.domain), ncol=N)
	residual <- matrix(nrow=length(obs.domain), ncol=N)
	for(i in 1:N){
		residual[,i] <- simul.norm[,i]-(b.basis %*% coeff.tau0.5)[,i]
		den <- density(residual[,i])
		w[1,i] <-  approx(den$x, den$y, xout=quantile(residual[,i], probs=0.01))$y
		w[2,i] <-  approx(den$x, den$y, xout=quantile(residual[,i], probs=0.99))$y
	}
	coeff.tau.al <- matrix(nrow=N, ncol=num.knots+3)
	for(i in 1:N){
	w.avg <- function(v){
		return(sum(w[,i]*v)/sum(w[,i]))
	}		
		coeff.tau.al[i,] <- apply(coeff.taus[,,i], 2, w.avg)
	}     
	clust.obj <- kmeans(coeff.tau.al, centers=2)
	tpr.clust <- vector(length=2)
	per <- permutations(2,2)
	for(i in 1:2){
		poo <- (per[i,])[clust.obj$cluster]
		tpr.clust[i] <- mean(poo == true.cluster)
	}
	result.vec[5,sn] <- max(tpr.clust) 
	result.vec.rand[5,sn] <- adjustedRandIndex(clust.obj$cluster, true.cluster)
	
	
	#2-8)funclust 
	bsp <- create.bspline.basis(breaks=seq(0, 1, by=1/num.knots))
	fd.obj <- fd(coef=coeff.tau0.5, basisobj=bsp)
	tpr <- vector(length=100)
	tpr2 <- vector(length=100)
		clust.obj <- funclust(fd.obj, K=2)
		tpr.clust <- vector(length=2)
		per <- permutations(2,2)
		for(i in 1:2){
			poo <- (per[i,])[clust.obj$cls]
			tpr.clust[i] <- mean(poo == true.cluster)
		}
	result.vec[6,sn] <- max(tpr.clust)
	result.vec.rand[6,sn] <- adjustedRandIndex(clust.obj$cls, true.cluster)
	
	#2-9)curveclust
	simul.norm.curv <- matrix(nrow=128, ncol=200)
	simul.norm.list <- list()
	simul.norm.curv[14:114,] <- simul.norm
	for(i in 1:13){
		simul.norm.curv[i,] <- simul.norm.curv[28-i,]
	}
	for(i in 115:128){
		simul.norm.curv[i,] <- simul.norm.curv[228-i,]
	}
	for(i in 1:200){
		simul.norm.list[[i]] <- simul.norm.curv[,i]
	}
	CCD <- apply(simul.norm.curv, 2, FUN = function(y) {
            wd(y, filter.number = 1, family = "DaubExPhase")
        })
   CCD <- new("CClustData", Y=simul.norm.list, filter.number=1)
   CCDred <- getUnionCoef(CCD)
   CCO            = new("CClustO")
   CCO["nbclust"] = 2
   CCR <- getFCMM(CCDred, CCO)
   summary(CCR)
   groups <- apply(CCR["Tau"], 1, which.max)
   	for(i in 1:2){
			tpr.clust[i] <- mean(groups == true.cluster)
		}
	max(tpr.clust)
	result.vec[7,sn] <- max(tpr.clust)
	result.vec.rand[7,sn] <- adjustedRandIndex(groups, true.cluster)

	
}
