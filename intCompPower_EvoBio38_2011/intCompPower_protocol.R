require(mvtnorm); require(MASS);  require(corpcor)

# The first section below provides a list of functions that are used repeatedly.
# The second section below provides the protocols for each analysis and figure.

########### Functions
# This function simulates a correlation matrix with a given integration level
# int is the desired integration level of the matrix
# p is number of variables 
# int.type determines whether int is given as mean absolute correlation 
# (mean.cor) or as the relative variance of eigenvalues (rel.ve)
# cv is the heterogeneity of the off-diagonal elements
# if int.type is rel.ve then both q and cv are determined by the function by 
# maximizing cv maximizing cv follows the equations given in Pavlicev et al. 2009 

simulateR <- function(int, p, cv=min((1-q)/q, 0.9), int.type=c("mean.cor", "rel.ve")) {
	R <- matrix(1,p,p)
	if (int.type=="rel.ve") {
		func <- function(v, p, rve) {rve-v/(p-1)}
		op <- optimize(f=func, interval=c(0,1), p=p, rve=int, maximum=TRUE)
		q <- sqrt(op$object)
		cv <- sqrt(op$maximum)/op$objective
		} else {q <- int; cv <- cv}
	qz <- 0.5*log((1+q)/(1-q))
	cv <- ifelse((1-q)/q<cv, (1-q)/q, cv)
	z <- rnorm(p*(p-1)/2, mean=qz, sd=cv*qz)	
	R[lower.tri(R)] <- (exp(2*z)-1)/(exp(2*z)+1)
	rt <- t(R)
	R[upper.tri(R)] <- rt[upper.tri(rt)]
	R}
	

### This function is for shrinking (bending) of singular matrices following Jorjani et al. 2003

shrinkR <- function(R, tol=10^-8) {
    ei <- eigen(R, symmetric=TRUE)
    d <- ei$values
    if (min(d) < tol) {
       di <- d
       di[d<tol] <- 2*tol
       di <- di*(sum(d)/sum(di))
       M <- ei$vectors %*% (di * t(ei$vectors))
       M <- cov2cor(M)
       dimnames(M) <- dimnames(R)
       return(M)} else {return(R)}
   }
   
### This function is for calculating coefficient of variation of the absolute values of the off-diagonal elements of a correlation matrix R (Pavlicev et al. 2009)
cv.r <- function(R) {
	r <- abs(R[lower.tri(R)])
	sqrt(var(r))/mean(r)}

### A function for calculating the mean coeeficient of determination of a correlation matrix R

IIr2 <- function(R) {
		r <- R[lower.tri(R)]
		mean(r^2)}

### A function for caluclating the mean of the absolute values of the off-diagonal elements of a correlation matrix R (Cane 1993)

IIr <- function(R) {
		r <- abs(R[lower.tri(R)])
		mean(r)}

### A function for calculating the relative standard deviation of eigenvalues of a correlation matrix R (Pavlicev et al. 2009)

IIsde <- function(R) {
	d <- eigen(R)$values
	p <- length(d)
	sqrt(sum((d-1)^2)/(p*(p-1)))
	}

### A function for calculating Van Valen's (1974) redundancy index

IIred <- function(R) {
	p <- ncol(R)
	R2 <- numeric(p)
	for (j in 1:p) {R2[j] <- t(R[-j,j])%*%ginv(R[-j,-j])%*%R[-j,j]}
	mean(R2)*(p-1)/p}

### A function for calculating the first approximation for the Hansen and Houle (2008) integration index for correlation matrix R including the correction from Hansen and Houle 2009

IIhh1 <- function(R) {
	d <- eigen(R)$values
	k <- length(d)
	hd <- (1/mean(1/d))/mean(d)
	i <- (var(d)*(k-1)/k)/mean(d)^2
	ir <- (var(1/d)*(k-1)/k)/mean(1/d)^2
	1-hd*(1+2*(i+ir-1+hd+2*i*ir/(k+2))/(k+2))}

### A function for calculating the second approximation for the Hansen and Houle (2008) integration index for correlation matrix R

IIhh2 <- function(R) {
	1-mean(1/diag(ginv(R)))}


############# Protocols
###### Protocol for generating matrices

ifuns <- list("IIr"=IIr, "IIr2"=IIr2, "IIsde"=IIsde,  "IIred"=IIred, "IIhh1"=IIhh1, "IIhh2"=IIhh2)
ni <- length(ifuns)
nm <- 2500
q <- runif(nm,0.02,0.98)
p <- round(runif(nm, 5, 35))
cv.intl <- runif(nm,0.02,2)

RR <- RR3 <- RR1 <- RF <- list()
for (i in 1:nm) {
	RF[[i]] <- Ri <- simulateR(int=q[i], p=p[i], cv=cv.intl[i], int.type="mean.cor")
	RR[[i]] <- shrinkR(Ri, tol=10^-6)
	RR3[[i]] <- shrinkR(Ri, tol=10^-3)
	RR1[[i]] <- shrinkR(Ri, tol=10^-1)
	}

si <- which(!sapply(RF, is.positive.definite))

### calculating stats

cv <- sapply(RF, cv.r)
cvs <- sapply(RR, cv.r)

It <- It6 <- It3 <- It1 <- c()
for (i in 1:ni) {
	ifun <- ifuns[[i]]
	It <- cbind(It, sapply(RF, ifun))
	It6 <- cbind(It6, sapply(RR, ifun))
	It3 <- cbind(It3, sapply(RR3, ifun))
	It1 <- cbind(It1, sapply(RR1, ifun))
	}
colnames(It1) <- colnames(It3) <- colnames(It6) <- colnames(It)  <- names(ifuns)

save.image("ws-properties mixed-p.Rdata")
  
######## Figure 1
## cv vs cvs (without the empirical data)
quartz(width=5, height=4.5)
par(mar=c(4,4,0.5,0.5))
plot(It6[,"IIr"], cvs,  xlab=NA, ylab=NA, cex=1.2, cex.axis=1.2, ylim=c(0,1), xlim=c(0,1), pch=19)
mtext("Mean absolute correlation",cex=1.2, adj=NA, line=2.8, side=1)
mtext("CV absolute correlation",cex=1.2, adj=NA, line=2.8, side=2)

## Effect of shrinking by p
plot(It6[,"IIr"], cvs)
points(It6[which(p>=5 & p<=15),"IIr"], cvs[which(p>=5 & p<=15)], pch=19, col="red")
points(It6[which(p>=16 & p<=25),"IIr"], cvs[which(p>=16 & p<=25)], pch=19, col="green")
points(It6[which(p>=26 & p<=35),"IIr"], cvs[which(p>=26 & p<=35)], pch=19, col="blue")
legend("topright", legend=c("5-15", "16-25", "26-35"), text.col=c("red", "green", "blue"), cex=1.5)

#### Figure 2
## IIrve and IIr2 are exactly the same
quartz(width=10, height=5)
layout(matrix(1:2,1,2))
par(mar=c(4,4,4,0.5))
plot(It6[,"IIr"],It6[,"IIsde"], xlim=c(0,1), ylim=c(0,1), ylab=NA, xlab=NA, cex=1.2, cex.axis=1.2)
segments(0,0,1,1,col="grey60", lwd=4, lty=6)
mtext("Mean absolute correlation",cex=1.1, adj=NA, line=2.3, side=1)
mtext("Relative variance/SD of eigenvalues",cex=1, adj=NA, line=2.7, side=2)
points(It6[,"IIr"],It6[,"IIsde"]^2, cex=1.2, pch=2, col="blue")
points(It6[,"IIr"],It6[,"IIr"]^2, col="grey60", pch=19, cex=0.5)
legend("topleft", legend=c("rSDE", "rVE"), pch=c(1,2), col=c("black", "blue"), pt.cex=1.7, pt.lwd=1.5, cex=1.1, bty="n")
mtext("A.",cex=1.2, adj=-0.1, line=1.8, side=3, font=2)
plot(It6[,"IIr2"],It6[,"IIsde"], xlim=c(0,1), ylim=c(0,1), ylab=NA, xlab=NA, cex=1.2, cex.axis=1.2)
mtext("Mean squared correlation",cex=1.1, adj=NA, line=2.3, side=1)
mtext("Relative variance/SD of eigenvalues",cex=1, adj=NA, line=2.7, side=2)
points(It6[,"IIr2"],It6[,"IIsde"]^2, cex=1.2, pch=2, col="blue")
segments(0,0,1,1,col="grey60", lwd=4, lty=6)
legend("topleft", legend=c("rSDE", "rVE"), pch=c(1,2), col=c("black", "blue"), pt.cex=1.7, pt.lwd=1.5, cex=1.1, bty="n")
mtext("B.",cex=1.2, adj=-0.1, line=1.8, side=3, font=2)


##### Figure 3
# IIhh2 and IIhh1 are above and below IIred, respectively
quartz(width=9, height=6)
layout(matrix(1:9, 3, 3, byrow=TRUE))
par(mar=c(4,4.5,0,0.5), oma=c(0,0,3,0))
plot(It6[,"IIred"], It6[,"IIhh1"], xlim=c(0,1), ylim=c(0,1), cex.axis=1.4, xlab=NA, ylab=NA, col="black",lwd=1.6)
mtext("HH1",cex=1, adj=NA, line=2.7, side=2)
segments(0,0,1,1, lwd=4, col="grey60")
plot(It3[,"IIred"], It3[,"IIhh1"], xlim=c(0,1), ylim=c(0,1), cex.axis=1.4, xlab=NA, ylab=NA, col="black", lwd=1.6)
mtext("HH1",cex=1, adj=NA, line=2.7, side=2)
segments(0,0,1,1, lwd=4, col="grey60")
plot(It1[,"IIred"], It1[,"IIhh1"], xlim=c(0,1), ylim=c(0,1), cex.axis=1.4, xlab=NA, ylab=NA, col="black", lwd=1.6)
mtext("HH1",cex=1, adj=NA, line=2.7, side=2)
segments(0,0,1,1, lwd=4, col="grey60")
mtext(expression("Bending level: 10"^"-6"), adj=0.04, line=1, side=3, outer=TRUE)
mtext(expression("Bending level: 10"^"-3"), adj=0.45, line=1, side=3, outer=TRUE)
mtext(expression("Bending level: 10"^"-1"), adj=0.86, line=1, side=3, outer=TRUE)

par(mar=c(4,4.5,0,0.5))
plot(It6[,"IIred"], It6[,"IIhh2"], xlim=c(0,1), ylim=c(0,1), cex.axis=1.4, xlab=NA, ylab=NA, col="black",lwd=2)
mtext("HH2",cex=1, adj=NA, line=2.7, side=2)
segments(0,0,1,1, lwd=4, col="grey60")
plot(It3[,"IIred"], It3[,"IIhh2"], xlim=c(0,1), ylim=c(0,1), cex.axis=1.4, xlab=NA, ylab=NA, col="black", lwd=2)
mtext("HH2",cex=1, adj=NA, line=2.7, side=2)
segments(0,0,1,1, lwd=4, col="grey60")
plot(It1[,"IIred"], It1[,"IIhh2"], xlim=c(0,1), ylim=c(0,1), cex.axis=1.4, xlab=NA, ylab=NA, col="black", lwd=2)
mtext("HH2",cex=1, adj=NA, line=2.7, side=2)
segments(0,0,1,1, lwd=4, col="grey60")

# IIhh2 is exactly the same as IIred when corrected for the effective number of units of information:
pp <- (p-1)/p
par(mar=c(4,4.5,0,0.5))
plot(It6[,"IIred"], It6[,"IIhh2"]*pp, xlim=c(0,1), ylim=c(0,1), cex.axis=1.4, xlab=NA, ylab=NA, col="black",cex=1.5, lwd=2)
mtext("HH2 scaled",cex=1, adj=NA, line=2.7, side=2)
mtext("VVred",cex=1, adj=NA, line=2.7, side=1)
segments(0,0,1,1, lwd=4, col="grey60")
plot(It3[,"IIred"], It3[,"IIhh2"]*pp, xlim=c(0,1), ylim=c(0,1), cex.axis=1.4, xlab=NA, ylab=NA, col="black", cex=1.5, lwd=2)
mtext("HH2 scaled",cex=1, adj=NA, line=2.7, side=2)
mtext("VVred",cex=1, adj=NA, line=2.7, side=1)
segments(0,0,1,1, lwd=4, col="grey60")
plot(It1[,"IIred"], It1[,"IIhh2"]*pp, xlim=c(0,1), ylim=c(0,1), cex.axis=1.4, xlab=NA, ylab=NA, col="black", cex=1.5, lwd=2)
mtext("VVred",cex=1, adj=NA, line=2.7, side=1)
mtext("HH2 scaled",cex=1, adj=NA, line=2.7, side=2)
segments(0,0,1,1, lwd=4, col="grey60")

##### Figure 4
## IIred (hence IIhh1/IIhh2) vs. IIsde (without the empirical data)
quartz(width=10, height=5)
layout(matrix(1:2,1,2))
par(mar=c(4,4.5,4,0.5))
plot(It6[,"IIsde"], It6[,"IIred"], xlim=c(0,1), ylim=c(0,1), cex.axis=1, xlab=NA, ylab=NA)
mtext("Relative SD of eigenvalues",cex=1.2, adj=NA, line=2.3, side=1)
mtext("Van Valen's redundancy index",cex=1.2, adj=NA, line=2.7, side=2)
mtext("A.",cex=1.2, adj=-0.2, line=1.8, side=3, font=2)
plot(It1[,"IIsde"], It1[,"IIred"], xlim=c(0,1), ylim=c(0,1), cex.axis=1, xlab=NA, ylab=NA)
mtext("Relative SD of eigenvalues",cex=1.2, adj=NA, line=2.3, side=1)
mtext("Van Valen's redundancy index",cex=1.2, adj=NA, line=2.7, side=2)
mtext("B.",cex=1.2, adj=-0.2, line=1.8, side=3, font=2)

##### Figure 5
## IIred (hence IIhh1/IIhh2) vs. IIsde
cvsets <- list(which(cvs<0.1), which(cvs>0.2 & cvs<0.3), which(cvs>0.4 & cvs<0.5), which(cvs>0.6))
lcv <- length(cvsets)
cvtxt <- list("cv: 0 - 0.1", "cv: 0.2 - 0.3", "cv: 0.4 - 0.5", "cv: 0.6 - 1")
psets <- list(which(p<15), which(p>=15 & p<25), which(p>=25))
lp <- length(psets)
ptxt <- list("p: 5 - 14", "p: 15 - 24", "p: 25 - 35")
quartz(width=12, height=6)
par(mar=c(2,2,1,1), oma=c(2,2,1.5,0))
layout(matrix(1:(lcv*lp), nr=lp, nc=lcv))
for(i in 1:lcv){
	for (j in 1:lp) {
		ic <- cvsets[[i]]
		jp <- psets[[j]]
		h <- ic[ic%in%jp]
		plot(It6[,"IIsde"], It6[,"IIred"], xlim=c(0,1), ylim=c(0,1), col="grey60", cex.axis=1.4, xlab=NA, ylab=NA)
		points(It6[h,"IIsde"], It6[h,"IIred"], cex.axis=1.4, pch=19)
		text(x=0.7, y=0.2, labels=cvtxt[i], cex=1.2)
		text(x=0.7, y=0.09, labels=ptxt[j], cex=1.2)
		}
	mtext("rSDE", side=1, line=0.5, outer=TRUE, cex=0.8, adj=0.12+(i-1)*0.257)
	}

for (i in 1:lp) {
	mtext("VVred", side=2, line=0.5, outer=TRUE, cex=1, adj=0.86-(i-1)*0.355)
	}

#### stronger bending
cvsets <- list(which(cvs<0.1), which(cvs>0.2 & cvs<0.3), which(cvs>0.4 & cvs<0.5), which(cvs>0.6))
lcv <- length(cvsets)
cvtxt <- list("cv: 0 - 0.1", "cv: 0.2 - 0.3", "cv: 0.4 - 0.5", "cv: 0.6 - 1")
psets <- list(which(p<15), which(p>=15 & p<25), which(p>=25))
lp <- length(psets)
ptxt <- list("p: 5 - 14", "p: 15 - 24", "p: 25 - 35")
quartz(width=12, height=6)
par(mar=c(2,2,1,1), oma=c(2,2,1.5,0))
layout(matrix(1:(lcv*lp), nr=lp, nc=lcv))
for(i in 1:lcv){
	for (j in 1:lp) {
		ic <- cvsets[[i]]
		jp <- psets[[j]]
		h <- ic[ic%in%jp]
		plot(It1[,"IIsde"], It1[,"IIred"], xlim=c(0,1), ylim=c(0,1), col="grey60", cex.axis=1.4, xlab=NA, ylab=NA)
		points(It1[h,"IIsde"], It1[h,"IIred"], cex.axis=1.4, pch=19)
		text(x=0.7, y=0.2, labels=cvtxt[i], cex=1.2)
		text(x=0.7, y=0.09, labels=ptxt[j], cex=1.2)
		}
	mtext("rSDE", side=1, line=0.5, outer=TRUE, cex=1, adj=0.12+(i-1)*0.257)
	}

for (i in 1:lp) {
	mtext("VVred", side=2, line=0.5, outer=TRUE, cex=1, adj=0.86-(i-1)*0.355)
	}

save.image("ws-properties.Rdata")
ws=ls()
rm(list=ws[-which(ws%in%c("IIhh1","IIhh2","IIr","IIr2","IIred","IIsde","shrinkR",        "simulateR"))])

########### sampling simulations (Figs 6-7)

## one nonparametric bootstrap round for the bootNPcor function
sampnp <- function(N, X1, func, shtol){
		Xb1 <- X1[sample(1:nrow(X1), N[1], replace=TRUE),]
		func(shrinkR(cor(Xb1), tol=shtol))}
	
## generating distribution using parametric bootstrap
# X1 is a sample matrix of n observations and p variables; func is the function to operate on the matrix (e.g., an integration index)
# shtol is the shrinking level of the matrix
bootNPcor <- function(X1, n, func, n.boot=500, shtol) {
		N <- rep(n, n.boot)
		dist <- sort(sapply(N, sampnp, X1=X1, func=func, shtol=shtol))
		dist}

## simulations for p=35; repeat for p=10
ifuns <- list("rSDE"=IIsde, "VVred"=IIred)
ni <- length(ifuns)

p <- 35
figs <- c("6a","6b")
# p=10
# figs <- c("7a","7b")
qq <- c(0.10, 0.22, 0.4, 0.65, 0.85)
nq <- length(qq)
nx <- 100 #sample size of the initial matrix
sht <- c(10^-6, 10^-1) # matrix bending level
NN <- seq(30, 94, 8) 
ln <- length(NN)

n.boot <- 500

RF <- array(sapply(qq, simulateR, p=p, cv=0.3, int.type="mean.cor"), dim=c(p, p, nq))
inds <- expand.grid(1:ln, 1:nq, 1:ni)

XL <- distINL <- ItrL <- list()
distIn <- matrix(NA, n.boot, ln*nq*ni)
for (s in 1:2) {
	XL[[s]] <- XX <- array(apply(RF, 3, function(R){rmvnorm(nx, sigma=shrinkR(R, tol=sht[s]))}), dim=c(nx, p, nq))
	RR <- array(apply(XX, 3, function(X) {shrinkR(cor(X),tol=sht[s])}), dim=c(p, p, nq))
	itr <- c()
	for (j in 1:nrow(inds)) {
			nj <- inds[j,1]
			xj <- inds[j,2]
			ifj <- inds[j,3]
			distIn[,j] <- bootNPcor(X1=XX[,,xj], n=NN[nj], func=ifuns[[ifj]], shtol=sht[s], n.boot=n.boot)
			itr <- c(itr, ifuns[[ifj]](RR[,,xj]))
			print(c(s,j))}
	ItrL[[s]] <- itr <- c(matrix(itr, ln, nq*ni)[1,], apply(RR, 3, IIr))
	distINL[[s]] <- distIN <- array(distIn, dim=c(n.boot, ln, nq*ni), dimnames=list(NULL,as.character(NN),NULL))
		}

for (s in 1:2) {
	distIN <- distINL[[s]]
	itr <- ItrL[[s]]
	quartz(width=12, height=3.5, file=paste("Fig", figs[s], ".pdf"), type="pdf")
	#quartz(width=12, height=3.5)
	par(mar=c(2,2,1,1), oma=c(2,2.8,2,0))
	layout(matrix(1:(nq*ni),nr=ni,nc=nq, byrow=TRUE))
	for (h in 1:(nq*ni)) {
		m <- colMeans(distIN[,,h])
		cil <- apply(distIN[,,h], 2, quantile, probs=0.025)
		cih <- apply(distIN[,,h], 2, quantile, probs=0.975)
		plot(NN, m, pch="-", cex=2, lwd=1.5, xaxt="n", ylim=c(0,1))
		axis(1, at=NN[seq(1,length(NN),2)], labels=NN[seq(1,length(NN),2)], cex.axis=1.2)
		segments(NN, cil, NN, cih, lwd=2)
		segments(0,itr[h],100,itr[h], col="red", lwd=1.5)}
		for (i in 1:ni) {mtext(names(ifuns)[i], side=2, line=0.8, outer=TRUE, cex=1, adj=0.8-(i-1)*0.58)}
		for (i in 1:nq) {
			mtext("Sample size", side=1, line=0.5, outer=TRUE, cex=1, adj=0.07+(i-1)*0.217)
			mtext(paste("Mean |r| = ",itr[nq*ni+i]), side=3, line=-0.5, outer=TRUE, cex=1, adj=0.06+(i-1)*0.222)}
			mtext(c("A.","B.")[s], side=3, outer=TRUE, adj=-0.02, line=0.5, font=2)
		dev.off()
	}
save.image("ws-dist 35p.Rdata")
# save.image("ws-dist 10p.Rdata")

ws=ls()
rm(list=ws[-which(ws%in%c("IIhh1","IIhh2","IIr","IIr2","IIred","IIsde","shrinkR",        "simulateR"))])

######## Power analysis Figures 8 and 9

## one parametric bootstrap round for the bootPARcor function
samppa <- function(N, R, func, shtol) {
		X1 <- rmvnorm(N, sigma=R)
		func(shrinkR(cor(X1), tol=shtol))}
		
## generating distribution using parametric bootstrap
# R is a correlation matrix; n is sample size, func is the function to operate on the matrix (e.g., an integration index)
# shtol is the shrinking level of the matrix
bootPARcor <- function(R, n, n.boot=500, shtol=10^-6, func) {
		N <- rep(n, n.boot)
		dist <- sort(sapply(N, samppa, S1=S1, func=func, shtol=shtol))
		dist}


ifun <- IIsde
p=35
# p=10 # repeat for 10 variables
# p=60 # repeat for 60 variables
shtol=10^-6

n.boot=500

sde <- seq(0.15, 0.85, 0.02)
rve <- sde^2
rve <- rep(rve, each=10)
nm <- length(rve)

NN <- seq(10, 82, 6)
ln <- length(NN)

RNc <-  expand.grid(1:nm, 1:ln)# each row is a vector of indices, the first element indicates the R matrix and the second indicates the N element to be included in one combination of d and N

############################################
RR <- array(dim=c(p,p,nm))
for (i in 1:nm) {
	RR[,,i] <- shrinkR(simulateR(int=rve[i], int.type="rel.ve", p=p), tol=shtol)
	}

II <- apply(RR, 3, ifun)
o <- order(II)
RR <- RR[,,o]
II <- II[o]
#Ir <- apply(RR, 3, IIr)

R0l <- shrinkR(simulateR(int=min(rve)-0.1*min(rve), int.type="rel.ve", p=p), tol=shtol) # weakyly integrated reference matrix
ifun(R0l) 
# 0.1516977
min(II) 
# 0.1554438
dI1 <- colMeans(matrix(abs(II-ifun(R0l)), 10, nm/10))

R0h <- shrinkR(simulateR(int=max(rve)+0.01*max(rve), int.type="rel.ve", p=p), tol=shtol) # strongly integrated reference matrix
ifun(R0h)
# 0.8540006
 max(II)
# 0.850202
dI2 <- colMeans(matrix(abs(II-ifun(R0h)), 10, nm/10))
orderdI2 <- order(dI2) # check that the order is correct and consistent across collumns
dI2i <- dI2[orderdI2] # switches direction to match that of the first set

quartz(width=10, height=5)
layout(matrix(1:2,1,2))
par(mar=c(4,4,4,0.5))
hist(dI1)
hist(dI1[dI1<0.1], add=TRUE)
hist(dI2)
hist(dI2[dI2<0.1], add=TRUE)

###################

getd0d1 <- function(N, R0, R1) {
	Rp <- (R0+R1)/2
	R1b <- shrinkR(cor(rmvnorm(N, sigma=R1)), tol=shtol)
	Rpb1 <- shrinkR(cor(rmvnorm(N, sigma=Rp)), tol=shtol)
	Rpb2 <- shrinkR(cor(rmvnorm(N, sigma=Rp)), tol=shtol)
	R0b <- shrinkR(cor(rmvnorm(N, sigma=R0)), tol=shtol)
	d0 <- abs(ifun(Rpb1)-ifun(Rpb2))
	d1 <- abs(ifun(R1b)-ifun(R0b))
	c(d0,d1)
	} # one sample with one combination of d and N

######### reference matrix is the most weakly integrated
pw <- c()
for (i in 1:nrow(RNc)) {
	ir <- RNc[i,1]; inn <- RNc[i,2]
	R <- RR[,,ir]
	nn <- rep(NN[inn], n.boot)
	d0d1 <- matrix(sapply(nn, getd0d1, R1=R, R0=R0l), n.boot, 2, byrow=TRUE)
	ci0 <- quantile(d0d1[,1], 0.95)
	d1 <- d0d1[,2]
	pw <- c(pw, length(d1[d1>=ci0])/(n.boot+1))
	print(i)}

PW1 <- pw
save.image("ws-power analysis1.R")

PWm <- array(pw, dim=c(10, nm/10, ln))
SP1 <- matrix(colMeans(PWm), nr=nm/10, nc=ln)
order(dI1) # check that the order is correct

DNP1p35 <- DNP1 <- list(dI1, NN, SP1); save(DNP1p35, file=paste("DNP1 p", p, ".Rdata", sep=""))
# DNP1p10 <- DNP1 <- list(dI1, NN, SP1); save(DNP1p10, file=paste("DNP1 p", p, ".Rdata", sep=""))
# DNP1p60 <- DNP1 <- list(dI1, NN, SP1); save(DNP1p60, file=paste("DNP1 p", p, ".Rdata", sep=""))

######### reference matrix is the most strongly integrated

pw <- c()
for (i in 1:nrow(RNc)) {
	ir <- RNc[i,1]; inn <- RNc[i,2]
	R <- RR[,,ir]
	nn <- rep(NN[inn], n.boot)
	d0d1 <- matrix(sapply(nn, getd0d1, R1=R, R0=R0h), n.boot, 2, byrow=TRUE)
	ci0 <- quantile(d0d1[,1], 0.95)
	d1 <- d0d1[,2]
	pw <- c(pw, length(d1[d1>=ci0])/(n.boot+1))
	print(i)}

PW2 <- pw
save.image("ws-power analysis2.R")

PWm <- array(pw, dim=c(10, nm/10, ln))
SP2 <- matrix(colMeans(PWm), nr=nm/10, nc=ln)
SP2i<-SP2[orderdI2,]
DNP2p35 <- DNP2 <- list(dI2i, NN, SP2i); save(DNP2p35, file=paste("DNP2 p", p, ".Rdata", sep=""))
# DNP2p10 <- DNP2 <- list(dI2i, NN, SP2i); save(DNP2p10, file=paste("DNP2 p", p, ".Rdata", sep=""))
# DNP2p60 <- DNP2 <- list(dI2i, NN, SP2i); save(DNP2p60, file=paste("DNP2 p", p, ".Rdata", sep=""))

save.image("ws-power analysis.R")

## Figure 8
quartz(width=10, height=5)
layout(matrix(1:2,1,2))
par(mar=c(4,4,4,0.5))
contour(dI1, NN, SP1, lwd=1.5, labcex=0.8, vfont=c("serif","bold"), cex.axis=1.2, col=c(rep("black",8),"red","black"))
axis(1,at=seq(0.12,0.28,0.02), labels=TRUE, cex.axis=0.55, tcl=-0.4, padj=-2.3)
segments(0,30,1,30, lty=3)
segments(0,40,1,40, lty=3)
segments(0.11,0,0.11,90, lty=3)
segments(0.18,0,0.18,90, lty=3)
segments(0.22,0,0.22,90, lty=3)
mtext("Integration difference",cex=1.2, adj=NA, line=2.3, side=1)
mtext("Sample size",cex=1.2, adj=NA, line=2.7, side=2)
mtext("A.",cex=1.2, adj=-0.1, line=1.8, side=3, font=2)
contour(dI2i, NN, SP2i, lwd=1.5, labcex=0.8, vfont=c("serif","bold"), cex.axis=1.2, col=c(rep("black",8),"red","black"))
axis(1,at=seq(0.12,0.28,0.02), labels=TRUE, cex.axis=0.55, tcl=-0.4, padj=-2.3)
segments(0,30,1,30, lty=3)
segments(0,40,1,40, lty=3)
segments(0.10,0,0.10,90, lty=3)
segments(0.17,0,0.17,90, lty=3)
segments(0.21,0,0.21,90, lty=3)
mtext("Integration difference",cex=1.2, adj=NA, line=2.3, side=1)
mtext("Sample size",cex=1.2, adj=NA, line=2.7, side=2)
mtext("B.",cex=1.2, adj=-0.1, line=1.8, side=3, font=2)


## Figure 9
quartz(width=10, height=5)
layout(matrix(1:2,1,2))
par(mar=c(4,4,4,0.5))
contour(DNP1p10[[1]], DNP1p10[[2]], DNP1p10[[3]], xlim=c(0,0.6), lwd=2.5, labcex=0.8, vfont=c("serif","bold"), cex.axis=1.2, levels=0.8, lty=6)
contour(DNP1p35[[1]], DNP1p35[[2]], DNP1p35[[3]], xlim=c(0,0.6), lwd=2.5, labcex=0.8, vfont=c("serif","bold"), cex.axis=1.2, levels=0.8, add=TRUE, lty=1)
contour(DNP1p60[[1]], DNP1p60[[2]], DNP1p60[[3]], xlim=c(0,0.6), lwd=2.5, labcex=0.8, vfont=c("serif","bold"), cex.axis=1.2, levels=0.8, add=TRUE, lty=3)
legend("topright", legend=c("p = 10", "p = 35", "p = 60"), lty=c(6,1,3), lwd=2)

mtext("Integration difference",cex=1.2, adj=NA, line=2.3, side=1)
mtext("Sample size",cex=1.2, adj=NA, line=2.7, side=2)
mtext("A.",cex=1.2, adj=-0.1, line=1.8, side=3, font=2)

segments(0,30,1,30, lty=3)
segments(0,40,1,40, lty=3)

contour(DNP2p10[[1]], DNP2p10[[2]], DNP2p10[[3]], xlim=c(0,0.6), lwd=2.5, labcex=0.8, vfont=c("serif","bold"), cex.axis=1.2, levels=0.8, lty=6)
contour(DNP2p35[[1]], DNP2p35[[2]], DNP2p35[[3]], xlim=c(0,0.6), lwd=2.5, labcex=0.8, vfont=c("serif","bold"), cex.axis=1.2, levels=0.8, add=TRUE, lty=1)
contour(DNP2p60[[1]], DNP2p60[[2]], DNP2p60[[3]], xlim=c(0,0.6), lwd=2.5, labcex=0.8, vfont=c("serif","bold"), cex.axis=1.2, levels=0.8, add=TRUE, lty=3)
legend("topright", legend=c("p = 10", "p = 35", "p = 60"), lty=c(6,1,3), lwd=2)
mtext("Integration difference",cex=1.2, adj=NA, line=2.3, side=1)
mtext("Sample size",cex=1.2, adj=NA, line=2.7, side=2)
mtext("B.",cex=1.2, adj=-0.1, line=1.8, side=3, font=2)

segments(0,30,1,30, lty=3)
segments(0,40,1,40, lty=3)















