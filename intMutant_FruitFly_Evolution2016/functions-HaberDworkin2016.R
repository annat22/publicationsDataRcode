#### functions to go along with protocol for Haber & Dworkin 2016

#### kmN2Nkm
# Converts an array of k x m x N to a matrix of N x k*m.
# Input (XX) is an array of k x m x N (e.g., N=number of specimens, k=number of landmarks, m=number of dimensions). 
# Output is a matrix of N specimens, where each specimen is a row vector arranged as {x1, y1, z1, x2, y2, z2...}
kmN2Nkm <- function(XX) {
	m <- ncol(XX)
	X <- t(apply(XX,3,t))
	if (!is.null(dimnames(XX)[[1]])) {dimnames(X) <- list(dimnames(XX)[[3]], paste(rep(rownames(XX), each=m), rep(colnames(XX),m), sep="_"))}
	X}

##### Nkm2kmN
# Converts a matrix of N x k*m to an array of k x m x N. 
# Input (X) is a matrix of N specimens, where each specimen is a row vector arranged as {x1, y1, z1, x2, y2, z2...}
Nkm2kmN <- function(X, m=3) {
		k <- ncol(X)/m
		N <- nrow(X)
		AA <- array(apply(array(t(X), dim=c(m,k,N)), 3, t), dim=c(k,m,N))
		cn <- colnames(X)[seq(1,ncol(X),m)]
		if (m==2) mm <- c("x","y") else mm <- c("x","y","z")
		#dimnames(AA) <- list(substr(cn,1,nchar(cn)-2), mm, rownames(X))
		dimnames(AA) <- list(paste("LM", 1:k, sep=""), mm, rownames(X))
		AA}
		
####### centsize
# calculates centroid size as the mean euclidean distance of b/w every landmark and the centroid
# X is one specimen matrix of k landmarks by m dimensions
centsize <- function(X) {
	sqrt(sum(apply(X, 2, var)*(nrow(X)-1)))}

######### shrink
# shrinks the matrix to avoid singularity
# M is a covariance matrix; tol is shrinking tollerance
# returns the shrunk matrix
# for both correlation and covariance matrices
# Based on Jorjani, H., L. Klei, and U. Emanuelson. 2003. A Simple Method for Weighted Bending of Genetic (Co)variance Matrices. Journal of Dairy Science 86:677-679. 

shrink <- function(M, tol=10^-8) {
    ei <- eigen(M, symmetric=TRUE)
    d <- ei$values
    rtol <- tol * mean(d)
    if (min(d) < rtol) {
    	if (sum(diag(M))==ncol(M)) {
    		di <- d
    		di[d<rtol] <- 2*rtol
    		di <- di*(sum(d)/sum(di))
    		Ms <- ei$vectors %*% (di * t(ei$vectors))
    		Ms <- cov2cor(Ms)
    		} else {
    			d[d<rtol] <- rtol
    			Ms <- ei$vectors %*% (d * t(ei$vectors))
    			}
    	dimnames(Ms) <- dimnames(M)
       	return(Ms)} else {return(M)}
   }

######## bootNP.BCa
# Calculates the observed theta, bootstraps the sample, and calculates BCa CI's
# x is a vector of either univariate data or 1:nrow(X) for a multivariate dataset X 
# (see example below), in which case X would be passed on as an additional argument of theta.fun
# For theta that includes comparison of two multivariate datasets X would 
# be rbind(X1,X2) and theta.fun would have an additional argument specifying grouping factor. 
# If replace=FALSE then data is permuted instead of bootstrapped

# Output is a list including: 
# the observed estimate (thetahat), 
# a vector of pseudovalues (thetastar), 
# the confidence interval for the specified alpha, including BCa if BCa=TRUE,
# the BCa parameters, z0, acc, and u, if BCa=TRUE,

### modified from http://www-rohan.sdsu.edu/~babailey/stat672/bootstrap.r
### and http://www.rohan.sdsu.edu/~babailey/stat672/bcanon.r

# To bootstrap bi- and multivariate datasets,
# write theta.fun so that its argument x
# is a set of observation indices
# and simply pass as data to bootstrap the vector 1,2,..n.
# For example, to bootstrap
# the correlation coefficient from a set of 15 data pairs:
#       xdata <- matrix(rnorm(30),ncol=2)
#       n <- 15
#       theta.funCor <- function(x,xdata){cor(xdata[x,1],xdata[x,2])} 
#       results <- bootNP.BCa(x=1:n, nboot=20, theta.fun=theta.funCor, xdata)


bootNP.BCa <- function(x, nboot, alpha=c(0.025,0.975), replace=TRUE, BCa=TRUE, theta.fun, ...){
	call <- match.call()
	thetahat <- theta.fun(x,...) # observed value
	n <- length(x)
	
	# Checking if a grouping factor is supplied for theta.fun (for comparing two multivariate datasets) 
	# and resample accordingly 
	args <- list(...)
	if (!"gf"%in%names(args)) { # no grouping factor, so all one group
		bootsam <- apply(matrix(x, nrow=n, nc=nboot), MARGIN=2, sample, size=n, replace=replace)
		} else {
			gf <- args["gf"] # grouping factor supplied for theta.fun
			ngf <- table(gf) # number of groups in the grouping factor
			bootsam <- c()
			i=1
			for (j in 1:length(ngf)){
				bj <- apply(matrix(i:(i+ngf[j]-1), nrow=ngf[j], nc=nboot), MARGIN=2, sample, size=ngf[j], replace=replace)
				i <- i+ngf[j]
				bootsam <- rbind(bootsam, bj)
				} # ensuring that original sample size is maintained within each group
			}
		
	thetastar <- apply(bootsam, MARGIN=2, theta.fun,...) # pseudovalues
 	confpoints <- NULL; z0 <- NULL; acc <- NULL; u=NULL
 	if (BCa==TRUE) {
 		z0 <- qnorm(sum(thetastar<thetahat, na.rm=TRUE)/(nboot+1))
	   	u <- rep(0,n)
   		for(i in 1:n) {u[i] <- theta.fun(x[-i],...)}  # Jackknife pseudovalues
   		bias <- mean(u, na.rm=TRUE)-u
   		acc <- sum(bias^3)/(6*(sum(bias^2))^1.5) # acceleration
   		zalpha <- qnorm(alpha)
   		cf <- pnorm(z0+(z0+zalpha)/(1-acc*(z0+zalpha))) # correction factors for both alphas
   		confpoints <- quantile(thetastar, probs=cf, na.rm=TRUE)
   		names(confpoints) <- paste("BCaCI", alpha, sep="")
 	} else confpoints <- quantile(thetastar, probs=alpha, na.rm=TRUE)
 	return(list(thetahat=thetahat, thetastar=thetastar, confpoints=confpoints, z0=z0, acc=acc, u=u, call=call, theta.args=names(list(...)), groupingFactor=ifelse(exists("gf"), gf, NA)))}


##### rSDE
# calculates the relative standard deviation of eigenvalues based on Van Valen (1974)
# for either covariances or correlation matrices (M)
# Covariance matrices are standardized by the total variance
# For correlation matrices this reduces to rSDE of Pavlicev et al. (2009) 

rSDE <- function(M) {
	d <- eigen(M, symmetric=TRUE)$values
	p <- length(d)
	sqrt(sum((d-mean(d))^2)*p/(sum(d)*sum(d)*(p-1)))}


##### Eccentricity 
# calculated as the inverse of effective number of dimensions from Kirkpatrick 2009 
eccentNd <- function(V) {
	ev <- eigen(V)$values
	ev[1]/sum(ev)} # inverse of effective number of dimensions; Kirkpatrick 2009


####### jackknife
# generates a jackknifed distribution of the statistic given in theta.fun
# theta.fun is the function applied to each pseudosample
# x is either a vector of univariate data or a vector of row indices for multivariate data
# xdata is the matrix of multivariate data that matches x if x is row indices
# jack.gf is grouping factor in case the jackknifed units are not the individual
# (e.g., family, locality, etc.). Default is to jackknife by individual
# ... additional arguments transfered to theta.fun
# CI calculation is from 
# http://www.math.ntu.edu.tw/~hchen/teaching/LargeSample/references/R-bootstrap.pdf 
# and http://darwin.phyloviz.net/ComparingPartitions/index.php?link=Tut9

jackknife <- function(x, theta.fun, jack.gf=as.character(1:length(x)), alpha=c(0.025,0.975), ...) {
	call <- match.call()
	thetahat <- theta.fun(x, ...)
	n <- length(unique(jack.gf))
	jackdist0 <- c()
	for (gfl in unique(jack.gf)) {
		xl <- x[-which(jack.gf==gfl)]
		jackdist0[gfl] <- theta.fun(xl, ...)
		}
	jackdist <- n*thetahat-(n-1)*jackdist0 # corrected for estimation bias
	
	est <- mean(jackdist)
	ci <- est+qt(alpha,n-1)*sqrt(var(jackdist)/n)
	return(list(thetahat=thetahat, est=est, confpoints=ci, jackdist=jackdist, jackdist.or=jackdist0, call=call))}
	
######### randsk
# calculates similarity (here distance) between two covariance matrices (V1 and V2) 
# using the random skewer method,following Cheverud 1996 J.Evol.Biol 9:5-42 
# as explained in Cheverud and Marroig 2007 Genet.Mol.Biol. 30(2):461-469 
# and Marroig et al 2009 Evol.Biol 36:136-148. 
# Skewers are drawn from a normal distribution following Marroig et al. 2012 Evo.Bio. 38:225-241
# The similarity values usually provided by random skewers are here transformed to distances

randsk <- function(V1, V2, n.it=5000, B=NULL) {	
	if (is.null(B)) {B <- matrix(rnorm (ncol(V1)*n.it, mean = 0, sd = 1), ncol(V1), n.it)} # generating a sample of selection vectors
	B <- t(t(B)/sqrt(colSums(B^2))) # scaling them to unit length
	Z1 <- V1%*%B # response vectors of first VCV matrix
	Z2 <- V2%*%B # response vectors of second VCV matrix
	Z1 <- t(t(Z1)/sqrt(colSums(Z1^2))) # scaling response vectors of V1
	Z2 <- t(t(Z2)/sqrt(colSums(Z2^2))) # scaling response vectors of V2
	r <- diag(t(Z1)%*%Z2) # caluculating their dot-product
	z <- mean(0.5*log((1+r)/(1-r))) # taking their mean using fisher's transformation
	s <- (exp(z/0.5)-1)/(exp(z/0.5)+1) # un-transforming the mean
	d <- sqrt(1-s^2) # converting to a distance
	0.5*log((1+d)/(1-d)) # fisher transformation to normalize the distances
	}


####### comsubsp
# Common subspace, or Krzanowski's method for comparing two covariance matrices
# following Zelditch et al., 2006, Evolution & Development, 8, 46-60 and green book p.308
# Returns a distance metric, whereas Krzanowski's original metric 
# is a similarity one (Blows et al., 2004, American Naturalist, 163, 329-340)
# This is the fastest and most stable calculation
# V1, V2 are the covaraince matrices, q is number of dimensions to retain

comsubsp <- function(V1,V2, q=floor(ncol(V1)/2), ...) {
 	E <- eigen(V1, symmetric=TRUE)$vectors
 	Q <- E[,1:q]%*%solve(E)[1:q,]
 	E <- eigen(V2, symmetric=TRUE)$vectors
 	R <- E[,1:q]%*%solve(E)[1:q,]
	j <- try(eigen(Q-R)$values, silent=TRUE)
	if (is.character(j)) NA else sqrt(sum(asin(as.numeric(j[1:q]))^2))
 	} # distance; Zelditch version 
		

###### MBreleig
# A measure of the relative eigenvalues 
# following Mitteroecker and Bookstein, 2009, Evolution, 63, 727-737
# also {Bookstein and Mitteroecker, 2014, Evol Biol, 41, 336-350}

MBreleig <- function(V1, V2, ...) {
	# to avoid halting the loop when the matrices are un-solvable the whole thing is in a "try" function:
	ev <- try(eigen(solve(V2,V1))$values, silent=TRUE) 
	if (is.character(ev)) NA else sqrt(sum(log(as.numeric(ev))^2))
	}

