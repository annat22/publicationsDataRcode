### FUNCTIONS needed for prtcl-MdivII.R

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
		dimnames(AA) <- list(substr(cn,1,nchar(cn)-2), mm, rownames(X))
		#dimnames(AA) <- list(paste("LM", 1:k, sep=""), mm, rownames(X))
		AA}

######### shrink
# shrinks the matrix to avoid singularity
# M is a covariance matrix; tol is shrinking tollerance
# for both correlation and covariance matrices

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

########## ILMD
# Computes interlandmark distances
# X is a specimen matrix of k landmarks by m dimensions; works for both 2D and 3D
# LM is a matrix of kd x 2 where kd is the number of desired distances; each row in LM specifies the 
# row indices (or landmark names if applicable) in X for the pair of landmarks between which the
# distance is to be calculated. by default computes all pairwise distances
# last updated July 20 2010

ILMD <- function(X, LM=NULL) {
		k <- nrow(X)
		if(is.null(rownames(X))){lnames <- paste("l",1:k,sep=".")} else {(lnames <- rownames(X))}
		D <- as.matrix(dist(X))
		l1 <- lnames[combn(k,2)[1,]]
		l2 <- lnames[combn(k,2)[2,]]		
		Dnm <- paste(l1, l2, sep="-")
		d <- D[lower.tri(D)]; names(d) <- Dnm
		if (is.null(LM)) {d} else {
			a <- b <- cbind(match(LM[,1],lnames),match(LM[,2],lnames))
			o <- t(apply(a,1,order))
			a[which(o[,1]==2),1] <- b[which(o[,1]==2),2]
			a[which(o[,1]==2),2] <- b[which(o[,1]==2),1]
			d[paste(lnames[a[,1]],lnames[a[,2]], sep="-")]
			}}
			
####### ILMDtes
# generates a definition matrix for interlandmark distances based on a tesselation representation of a 2D or 3D configuration.
# X is k x m matrix of landmark coordinates of one configuration; k landmarks and m dimensions. 
# (presumably mean configuration of a sample)
# if X is symmetric then only one side of the configuration should be included 
# (e.g., only midline and right bilateral landmarks); otherwise symmetry is ignored.
# output is a matrix with two columns. each row is a pair of landmarks that define one interlandmark distance
# use ILMD.R function in order to calculate the interlandmark distances themselves from the definition matrix generated here.

require(geometry)
ILMDtes <- function(X) {
	if (is.null(rownames(X))) rownames(X) <- paste("l",1:nrow(X),sep=".")
	tes <- delaunayn(X)
	ild0 <- rbind(tes[,1:2], tes[,2:3], tes[,3:4], tes[,c(1,4)])
	ir <- which(ild0[,1]>ild0[,2])
	ild0[ir,] <- ild0[ir,2:1]
	ild.nm <- rownames(ild0) <- apply(ild0,1,paste, collapse="-")
	uil <- unique(ild.nm)
	ild <- ild0[uil,]
	ild <- ild[order(ild[,1]),]
	cbind(rownames(X)[ild[,1]], rownames(X)[ild[,2]])}
	

##### evorateM
# Calculates the evolutionary rate matrix following Revell (2009) and Revell and Collar (2010)
# Input is the data matrix (X) and the tree in ape's phylo format
# rownames(X) includes tip labels of tree but order does not need to be the same
# both tree and X can contain taxa that do not appear in the other object, these taxa will be ignored)
# Output is a symmetric matrix with nrow and ncol equals to ncol(X)

require(ape)
evorateM <- function(X, phy) {
			tx <-  phy$tip.label
			Nt <- length(tx)
			Xi <- X[tx,]
			Se <- vcv.phylo(phy)[tx,tx] # C in Revell (2009)
			invS <- solve(Se)
			one <- matrix(1,Nt,1)
			a <- t(t(one)%*%invS%*%Xi)*sum(invS)^-1
			Re <- t(Xi-one%*%t(a))%*%invS%*%(Xi-one%*%t(a))/(Nt-1)
			Re}

##### A reduced version of ace from package ape
ace <- function(phy, X) {
		Nt <- length(phy$tip.label)
		Nn <- phy$Nnode
		dis <- dist.nodes(phy)
		MRCA <- mrca(phy,full=TRUE)
		M <- dis[Nt+1, MRCA]
		dim(M) <- rep(sqrt(length(M)), 2)
		Sho<-M[-(1:Nt), 1:Nt] 	
		S <- vcv.phylo(phy)
		invS <- solve(S)
		a0 <- colSums(invS)%*%X*1/sum(invS)
		A <- rbind(a0, Sho[-1,]%*%invS%*%(X-(rep(1,Nt)%*%a0))+rep(1,(Nn-1))%*%a0)
		rownames(A) <- phy$node.label
		A}

######## bootNP.BCa
# Calculates the observed theta, bootstraps the sample, and calculates BCa CI's
# x is a vector of either univariate data or 1:nrow(X) for a multivariate dataset X 
# (see example below), in which case X would be passed on as an additional argument of theta.fun
# For theta that includes comparison of two multivariate datasets X would 
# be rbind(X1,X2) and theta.fun would have an additional argument specifying grouping factor. 
### modified from http://www-rohan.sdsu.edu/~babailey/stat672/bootstrap.r
### and http://www.rohan.sdsu.edu/~babailey/stat672/bcanon.r

# To bootstrap bi- and multivariate datasets,
# write theta.fun so that its argument x
# is the set of observation indices
# and simply pass as data to bootstrap the vector 1,2,..n.
# For example, to bootstrap
# the correlation coefficient from a set of 15 data pairs:
#       xdata <- matrix(rnorm(30),ncol=2)
#       n <- 15
#       theta.funCor <- function(x,xdata){cor(xdata[x,1],xdata[x,2])} # x is a vector specifying row indices
#       results <- bootNP.BCa(x=1:n, nboot=20, theta.fun=theta.funCor, xdata)


bootNP.BCa <- function(x, nboot, alpha=c(0.025,0.975), replace=TRUE, BCa=TRUE, theta.fun, gf=as.factor(rep("gf1", length(x))), ...){
	call <- match.call()
	n <- length(x)
	thetahat <- theta.fun(x,...) # observed value
	#bootsam <- matrix(sample(x,size=n*nboot,replace=TRUE), nrow=n, ncol=n.boot)
	ngf <- table(gf) # number of groups in the grouping factor, for comparing two multivariate datasets
	if (length(ngf)==1) {bootsam <- apply(matrix(x, nrow=n, nc=nboot), MARGIN=2, sample, size=n, replace=replace)
		} else {
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
 	return(list(thetahat=thetahat, thetastar=thetastar, confpoints=confpoints, z0=z0, acc=acc, u=u, call=call))}


##### rSDE
# calculates the relative standard deviation of eigenvalues based on Van Valen (1974)
# for either covariances or correlation matrices (M)
# Covariance matrices are standardized by the total variance
# For correlation matrices this reduces to rSDE of Pavlicev et al. (2009) 


rSDE <- function(M, ...) {
	d <- eigen(M, symmetric=TRUE)$values
	p <- length(d)
	sqrt(sum((d-mean(d))^2)*p/(sum(d)*sum(d)*(p-1)))}

##### AvFlex
# Calculates average flexibility 
# following Rolian 2009 and Marroig et al. 2009
# M is a Variance-covariance matrix
# B is a set of random vectors (skewers)
# B can be either provided or generated within (defult)

Avflex <- function(M, n.it=1000, B=NULL) {	
	if (is.null(B)) {
		B <- matrix(rnorm (ncol(M)*n.it, mean = 0, sd = 1), ncol(M), n.it)}
	B <- t(t(B)/sqrt(colSums(B^2)))
	Z <- M%*%B
	Z <- t(t(Z)/sqrt(colSums(Z^2)))
	r <- diag(t(Z)%*%B)
	z <- mean(0.5*log((1+r)/(1-r)))
	(exp(z/0.5)-1)/(exp(z/0.5)+1)}
	
#### AvCondE
# Calculated average conditional evolvability following Hansen and Haoule 2008 using simulations
# V is a VCV matrix; the set of random vectors (B) can be either determined in advance or generated within (defult)
require(MASS)
AvCondE <- function(V, n.it=1000, B=NULL, ...) {	
	if (is.null(B)) {
		B <- matrix(rnorm (ncol(V)*n.it, mean = 0, sd = 1), ncol(V), n.it)
		B <- t(t(B)/sqrt(colSums(B^2))) # standardized to unit length
		}
		mean(1/diag(t(B)%*%ginv(V)%*%B))
		}
	
	
#### AvE
# Calculated average evolvability following Hansen and Haoule 2008 using simulations
# V is a VCV matrix; the set of random vectors (B) can be either determined in advance or generated within (defult)
AvE <- function(V, n.it=1000, B=NULL, ...) {	
	if (is.null(B)) {
		B <- matrix(rnorm (ncol(V)*n.it, mean = 0, sd = 1), ncol(V), n.it)
		B <- t(t(B)/sqrt(colSums(B^2))) # standardized to unit length
		}
		mean(diag(t(B)%*%V%*%B))
		}

#### aceRoot
# calclating root sate assuming Brownian motion (i.e., phylogenetically-weighted mean)
aceRoot <- function(phy, Xt) {
		tx <- phy$tip.label
		Nt <- length(tx)
		Xt <- Xt[tx,]
		Se <- vcv.phylo(phy)
		invS <- solve(Se)
		one <- matrix(1,Nt,1)
		c(t(one)%*%invS%*%Xt*sum(invS)^-1)
		}
