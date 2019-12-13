# Use the same dataset as in Example 2

N <- 300
num.knots <- 20
domain = seq(from=0, to=1, by=0.001)
obs.domain = seq(from=0, to=1, by=0.01)
knots <- seq(from=0, to=1, by=1/num.knots)
b.basis <- bsplineS(obs.domain, knots)
b.expand <- bsplineS(domain, knots)

result.vec <- matrix(nrow=5 ncol=5)
result.vec.rand <- matrix(nrow=5, ncol=5)

c.value <- 0.5 # cutoff value

    coeff.tau0.5 <- apply(simul.norm, 2, mquant.curve.pspline, tau=0.5, c=c.value)
    coeff.tau0.1 <- apply(simul.norm, 2, mquant.curve.pspline, tau=0.1, c=c.value)
    coeff.tau0.9 <- apply(simul.norm, 2, mquant.curve.pspline, tau=0.9, c=c.value)
    coeff.tau0.01 <- apply(simul.norm, 2, mquant.curve.pspline, tau=0.01, c=c.value)
    coeff.tau0.99 <- apply(simul.norm, 2, mquant.curve.pspline, tau=0.99, c=c.value)
    fit.mquant0.5 <- b.expand%*%coeff.tau0.5
    fit.mquant0.1 <- b.expand%*%coeff.tau0.1
    fit.mquant0.9 <- b.expand%*%coeff.tau0.9
    fit.mquant0.01 <- b.expand%*%coeff.tau0.01
    fit.mquant0.99 <- b.expand%*%coeff.tau0.99
    
    #kmeans 
    #1-1) mean based (0.5 mquanttile)
	clust.obj <- kmeans(t(coeff.tau0.5), centers=3)
	tpr.clust <- vector(length=6)
	per <- permutations(3,3)
	for(i in 1:6){
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
        w1[i2,i] <- abs(0.001*sum(fit.mquant0.9[search.vec, i2]-fit.mquant0.5[search.vec, i2]))
        w2[i2,i] <- abs(0.001*sum(fit.mquant0.5[search.vec, i2]-fit.mquant0.1[search.vec, i2]))
      }
      #w1[i2,] <- 1/w1[i2,]; w2[i2,] <- 1/w2[i2,]
      w11 <- w1[i2,]/(w1[i2,]+w2[i2,]); w22 <- w2[i2,]/(w1[i2,]+w2[i2,]); w1[i2,] <- w11; w2[i2,] <- w22
      poo1 <- w1[i2,1]; poo2 <- w2[i2,1]
      for(i in 1:J){
        search.vec <- seq(from=floor(1000/J*(i-1))+1, to=floor(1000/J*i))
        poo1 <- c(poo1, rep(w1[i2,i], length(search.vec)))
        poo2 <- c(poo2, rep(w2[i2,i], length(search.vec)))
      }
      weighted.expec[,i2] <-poo1*fit.mquant0.9[,i2]+poo2*fit.mquant0.1[,i2]
    }
    coeff.tau.al <- matrix(nrow=N, ncol=num.knots+3)
    for(i in 1:N){
      coeff.tau.al[i,] <- solve(t(b.expand)%*%b.expand, t(b.expand)%*%weighted.expec[,i])
    }     
	clust.obj <- kmeans(coeff.tau.al, centers=3)
	tpr.clust <- vector(length=6)
	per <- permutations(3,3)
	for(i in 1:6){
		poo <- (per[i,])[clust.obj$cluster]
		tpr.clust[i] <- mean(poo == true.cluster)
	}
	result.vec[2,sn] <- max(tpr.clust) 
	result.vec.rand[2,sn] <- adjustedRandIndex(clust.obj$cluster, true.cluster)

    
    #2-6-2) weighted average: 0.01, 0.99   
    w1 <- w2 <- matrix(nrow=N, ncol=J)
    weighted.expec <- matrix(nrow=length(domain), ncol=N)
    for(i2 in 1:N){
      for(i in 1:J){
        search.vec <- seq(from=floor(1000/J*(i-1))+1, to=floor(1000/J*i))
        w1[i2,i] <- abs(0.001*sum(fit.mquant0.99[search.vec, i2]-fit.mquant0.5[search.vec, i2]))
        w2[i2,i] <- abs(0.001*sum(fit.mquant0.5[search.vec, i2]-fit.mquant0.01[search.vec, i2]))
      }
      #w1[i2,] <- 1/w1[i2,]; w2[i2,] <- 1/w2[i2,]
      w11 <- w1[i2,]/(w1[i2,]+w2[i2,]); w22 <- w2[i2,]/(w1[i2,]+w2[i2,]); w1[i2,] <- w11; w2[i2,] <- w22
      poo1 <- w1[i2,1]; poo2 <- w2[i2,1]
      for(i in 1:J){
        search.vec <- seq(from=floor(1000/J*(i-1))+1, to=floor(1000/J*i))
        poo1 <- c(poo1, rep(w1[i2,i], length(search.vec)))
        poo2 <- c(poo2, rep(w2[i2,i], length(search.vec)))
      }
      weighted.expec[,i2] <-poo1*fit.mquant0.99[,i2]+poo2*fit.mquant0.01[,i2]
    }
    coeff.tau.al <- matrix(nrow=N, ncol=num.knots+3)
    for(i in 1:N){
      coeff.tau.al[i,] <- solve(t(b.expand)%*%b.expand, t(b.expand)%*%weighted.expec[,i])
    }     
	clust.obj <- kmeans(coeff.tau.al, centers=3)
	tpr.clust <- vector(length=6)
	per <- permutations(3,3)
	for(i in 1:6){
		poo <- (per[i,])[clust.obj$cluster]
		tpr.clust[i] <- mean(poo == true.cluster)
	}
	result.vec[2,sn] <- max(tpr.clust) 
	result.vec.rand[2,sn] <- adjustedRandIndex(clust.obj$cluster, true.cluster)

    
    #2-6-5) weighted average: error distribution
    w <- matrix(nrow=2, ncol=N)
    coeff.taus <- array(dim=c(2, num.knots+3, N))
    # coeff.tau0.5 <- apply(simul.norm, 2, mquant.curve.pspline, tau=0.5, c=c.value)
    coeff.taus[1,,] <- coeff.tau0.1 #apply(simul.norm, 2, mquant.curve.pspline, tau=0.1, c=c.value)
    coeff.taus[2,,] <- coeff.tau0.9 #apply(simul.norm, 2, mquant.curve.pspline, tau=0.9, c=c.value)
    weighted.mquant <- matrix(nrow=length(obs.domain), ncol=N)
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
	clust.obj <- kmeans(coeff.tau.al, centers=3)
	tpr.clust <- vector(length=6)
	per <- permutations(3,3)
	for(i in 1:6){
		poo <- (per[i,])[clust.obj$cluster]
		tpr.clust[i] <- mean(poo == true.cluster)
	}
	result.vec[2,sn] <- max(tpr.clust) 
	result.vec.rand[2,sn] <- adjustedRandIndex(clust.obj$cluster, true.cluster)


    #2-6-6) weighted average: error distribution
    w <- matrix(nrow=2, ncol=N)
    coeff.taus <- array(dim=c(2, num.knots+3, N))
    # coeff.tau0.5 <- apply(simul.norm, 2, mquant.curve.pspline, tau=0.5, c=c.value)
    coeff.taus[1,,] <- coeff.tau0.01 #apply(simul.norm, 2, mquant.curve.pspline, tau=0.01, c=c.value)
    coeff.taus[2,,] <- coeff.tau0.99 #apply(simul.norm, 2, mquant.curve.pspline, tau=0.99, c=c.value)
    weighted.mquant <- matrix(nrow=length(obs.domain), ncol=N)
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
	clust.obj <- kmeans(coeff.tau.al, centers=3)
	tpr.clust <- vector(length=6)
	per <- permutations(3,3)
	for(i in 1:6){
		poo <- (per[i,])[clust.obj$cluster]
		tpr.clust[i] <- mean(poo == true.cluster)
	}
	result.vec[2,sn] <- max(tpr.clust) 
	result.vec.rand[2,sn] <- adjustedRandIndex(clust.obj$cluster, true.cluster)