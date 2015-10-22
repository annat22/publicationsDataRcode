### FUNCTIONS needed for prtcl-EvoOfInt_EvoBio.R

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
		
####### centsize
# calculates centroid size as the mean euclidean distance of b/w every landmark and the centroid
# X is one specimen matrix of k landmarks by m dimensions
centsize <- function(X) {
	sqrt(sum(apply(X, 2, var)*(nrow(X)-1)))}

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
	
########## randsk
#calculates simmilarity between two covariance matrices (V1 and V2) using the random skewer method,
# following Cheverud 1996 J.Evol.Biol 9:5-42 
# explained in Cheverud and Marroig 2007 Genet.Mol.Biol. 30(2):461-469 
# and Marroig et al 2009 Evol.Biol 36:136-148 
# skewers are drawn from a normal distribution following Marroig et al. 2012 Evo.Bio. 38:225-241

randsk <- function(V1, V2, n.it=1000, B=NULL) {	
	if (is.null(B)) {
		B <- matrix(rnorm (ncol(V1)*n.it, mean = 0, sd = 1), ncol(V1), n.it) # generating a sample of random vectors
		B <- t(t(B)/sqrt(colSums(B^2))) # scaling them to unit length
		}
	Z1 <- V1%*%B # response vectors of first VCV matrix
	Z2 <- V2%*%B # response vectors of second VCV matrix
	Z1 <- t(t(Z1)/sqrt(colSums(Z1^2))) # scaling response vectors of V1
	Z2 <- t(t(Z2)/sqrt(colSums(Z2^2))) # scaling response vectors of V2
	r <- diag(t(Z1)%*%Z2) # caluculating their dot-product
	z <- mean(0.5*log((1+r)/(1-r))) # taking their mean using fisher's transformation
	(exp(z/0.5)-1)/(exp(z/0.5)+1) # re-transforming the mean to vary between -1 and 1
	}

#### matcor
# matrix correlation
# V1 and V2 are either covariances or correlation matrices; covariances are converted to correlations
matcor <- function(V1, V2, ...) {
	if (any(diag(V1)!=1)) {V1 <- cov2cor(V1)}
	if (any(diag(V2)!=1)) {V2 <- cov2cor(V2)}
	cor(V1[lower.tri(V1)], V2[lower.tri(V2)])}

##### tip.disp
# calculates disparity of tips descending from a given node
tip.disp <- function (node, data, phy, disp=c("trace","avg.sq")) {
    Ntip <- length(phy$tip.label)
    nb.node <- which(phy$node.label==node)
    Xt <- data[tips(phy, Ntip+nb.node), ]
    if (disp=="trace") {d <- sum(diag(var(Xt)))}
    if (disp=="avg.sq") {d <- mean(dist(Xt, method = "euclidean")^2)}   
    d}

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

###############################################################
#### decdiv 
# All of the following functions are taken from 
# Pavoine et al. (2010) Ecological Monographs 2010, 80:485-507
# http://www.esajournals.org/doi/abs/10.1890/09-1290.1
# http://www.esapubs.org/archive/mono/M080/017/supfileA.R

decdiv <- function(phy, df, dis = NULL, tol = 1e-08){
	
	if(is.vector(df)){
        df <- cbind.data.frame(df)
    }
    if(!is.data.frame(df)) stop("df should be a data frame")
    if(any(apply(df, 2, sum)<tol)) stop("null column in df")
    if(any(df < -tol)) stop("negative values in df")
    df[df < tol] <- 0
    df <- as.data.frame(apply(df, 2, function(x) x/sum(x)))
    
    disc2 <- function(samples, dis = NULL, structures = NULL, tol = 1e-08) 
    {
        if (!inherits(samples, "data.frame")) 
            stop("Non convenient samples")
        if (any(samples < 0)) 
            stop("Negative value in samples")
        if (any(apply(samples, 2, sum) < 1e-16)) 
            stop("Empty samples")
        if (!is.null(dis)) {
            if (!inherits(dis, "dist")) 
                stop("Object of class 'dist' expected for distance")
            if (!is.euclid(dis)) 
                warning("Euclidean property is expected for distance")
            dis <- as.matrix(dis)
            if (nrow(samples) != nrow(dis)) 
                stop("Non convenient samples")
        }
        if (is.null(dis)) 
            dis <- (matrix(1, nrow(samples), nrow(samples)) - diag(rep(1, 
                nrow(samples)))) * sqrt(2)
        if (!is.null(structures)) {
            if (!inherits(structures, "data.frame")) 
                stop("Non convenient structures")
            m <- match(apply(structures, 2, function(x) length(x)), 
                ncol(samples), 0)
            if (length(m[m == 1]) != ncol(structures)) 
                stop("Non convenient structures")
            m <- match(tapply(1:ncol(structures), as.factor(1:ncol(structures)), 
                function(x) is.factor(structures[, x])), TRUE, 0)
            if (length(m[m == 1]) != ncol(structures)) 
                stop("Non convenient structures")
        }
        Structutil <- function(dp2, Np, unit) {
            if (!is.null(unit)) {
                modunit <- model.matrix(~-1 + unit)
                sumcol <- apply(Np, 2, sum)
                Ng <- modunit * sumcol
                lesnoms <- levels(unit)
            }
            else {
                Ng <- as.matrix(Np)
                lesnoms <- colnames(Np)
            }
            sumcol <- apply(Ng, 2, sum)
            Lg <- t(t(Ng)/sumcol)
            colnames(Lg) <- lesnoms
            Pg <- as.matrix(apply(Ng, 2, sum)/nbhaplotypes)
            rownames(Pg) <- lesnoms
            deltag <- as.matrix(apply(Lg, 2, function(x) t(x) %*% 
                dp2 %*% x))
            ug <- matrix(1, ncol(Lg), 1)
            dg2 <- t(Lg) %*% dp2 %*% Lg - 1/2 * (deltag %*% t(ug) + 
                ug %*% t(deltag))
            colnames(dg2) <- lesnoms
            rownames(dg2) <- lesnoms
            return(list(dg2 = dg2, Ng = Ng, Pg = Pg))
        }
        Diss <- function(dis, nbhaplotypes, samples, structures) {
            structutil <- list(0)
            structutil[[1]] <- Structutil(dp2 = dis, Np = samples, 
                NULL)

            ###
            diss <- list(as.dist(structutil[[1]]$dg2))
            fun1 <- function(x){
                y <- x
                y[y<tol] <- 0
                return(y)
            }
            diss <- lapply(diss, fun1)
            diss <- lapply(diss, function(x) sqrt(2*x))
            ###

            if (!is.null(structures)) {
                for (i in 1:length(structures)) {
                    structutil[[i + 1]] <- Structutil(structutil[[1]]$dg2, 
                    structutil[[1]]$Ng, structures[, i])
                }
                ###
                diss <- c(diss, tapply(1:length(structures), factor(1:length(structures)),
                    function(x) (as.dist(structutil[[x + 1]]$dg2))))
                diss <- lapply(diss, fun1)
                diss <- lapply(diss, function(x) sqrt(2*x))
                ###
            }
            return(diss)
        }
        nbhaplotypes <- sum(samples)
        diss <- Diss(dis^2, nbhaplotypes, samples, structures)
        names(diss) <- c("samples", names(structures))
        if (!is.null(structures)) {
            return(diss)
        }
        return(diss$samples)
    }
   
    decdivV <- function(freq){
       nsp <- length(phy$leaves)
	   nnodes <- length(phy$nodes)
	   matno <- as.data.frame(matrix(0, nnodes, nsp))
	   rownames(matno) <- names(phy$nodes)
	   names(matno) <- names(phy$leaves)
	   for(i in 1:nsp){
		  matno[phy$path[[i]][-length(phy$path[[i]])], i] <- 1
	   }
	   matfr <- as.matrix(matno) %*% diag(freq)
	   matfr2 <- as.data.frame(t(matfr))
    	divno <- divc(matfr2, dis)
        matfr3 <- cbind.data.frame(matfr2, diag(freq))
        names(matfr3) <- c(names(matfr2), names(phy$leaves))
        matfr4 <- matfr3[, apply(matfr3, 2, sum)!=0]
        if(ncol(matfr4)==0) stop("only one species considered")
        discno <- disc2(matfr4, dis, tol = tol)
        lambdano <- apply(matfr4, 2, sum)
        prdist <- diag(lambdano)%*%as.matrix(discno^2/2)%*%diag(lambdano)
        colnames(prdist) <- rownames(prdist) <- names(matfr4)
        fun1 <- function(x){
            x <- x[apply(matfr3[x], 2, sum)!=0]
            if(length(x) == 1) return(0)
            else return(sum(prdist[x, x])/2)
        }
        res <- unlist(lapply(phy$parts, fun1))
        lambdano <- apply(matfr3, 2, sum)
        lambdano[lambdano < tol] <- 1
        res <- res * 1/as.vector(lambdano)[1:nnodes]
	   return(res)
    }
    return(apply(df, 2, decdivV))
}

plot.decdiv <- 
function (x, y = NULL, f.phylog = 0.5, cleaves = 1, cnodes = 0, vnodes = NULL, vcolor = "white", adjust = 1, vmax = max(vnodes),
    labels.leaves = names(x$leaves), clabel.leaves = 1, labels.nodes = names(x$nodes), 
    clabel.nodes = 0, sub = "", csub = 1.25, possub = "topleft", 
    draw.box = FALSE, clegend = 1, ...) 
{
    if (!inherits(x, "phylog")) 
        stop("Non convenient data")
    leaves.number <- length(x$leaves)
    leaves.names <- names(x$leaves)
    nodes.number <- length(x$nodes)
    nodes.names <- names(x$nodes)
    if (length(labels.leaves) != leaves.number) 
        labels.leaves <- names(x$leaves)
    if (length(labels.nodes) != nodes.number) 
        labels.nodes <- names(x$nodes)
    leaves.car <- gsub("[_]", " ", labels.leaves)
    nodes.car <- gsub("[_]", " ", labels.nodes)
    mar.old <- par("mar")
    on.exit(par(mar = mar.old))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    if (f.phylog < 0.05) 
        f.phylog <- 0.05
    if (f.phylog > 0.95) 
        f.phylog <- 0.95
    maxx <- max(x$droot)
    plot.default(0, 0, type = "n", xlab = "", ylab = "", xaxt = "n", 
        yaxt = "n", xlim = c(-maxx * 0.15, maxx/f.phylog), ylim = c(-0.05, 
            1), xaxs = "i", yaxs = "i", frame.plot = FALSE)
    x.leaves <- x$droot[leaves.names]
    x.nodes <- x$droot[nodes.names]

    reglage <- maxx * 0.10 * adjust
    cex0 <- par("cex") * clegend
    ###
        legender <- function(br0, sq0, clegend) {
        sq0 <- round(sq0, dig = 3)
        cha <- as.character(sq0[1])
        for (i in (2:(length(sq0)))) cha <- paste(cha, sq0[i], 
            sep = " ")
        yh <- max(c(strheight(cha, cex = cex0), br0))
        h <- strheight(cha, cex = cex0)
        y0 <- par("usr")[3] + yh / (par("usr")[2] - par("usr")[1]) * (par("usr")[4] - par("usr")[3])
        x0 <- par("usr")[1] + h/2
        for (i in (1:(length(br0)))) {
            cha <- sq0[i]
            cha <- paste(" ", cha, sep = "")
            xh <- strwidth(cha, cex = cex0)
            text(x0 + xh/2, y0, cha, cex = cex0)
            z0 <- br0[i]
            x0 <- x0 + xh + z0
            symbols(x0, y0, circles = z0, bg = vcolor, 
                    fg = "black", add = TRUE, inch = FALSE)
            x0 <- x0 + z0
        }
        invisible()
    }
    if (clegend > 0) {
            sq0 <- pretty(c(vnodes, vmax), 4)
            l0 <- length(sq0)
            sq0 <- (sq0[1:(l0 - 1)] + sq0[2:l0])/2
            br0 <- (sq0 / vmax * reglage)
            legender(br0, sq0, clegend = clegend)
    }
    ###
    
    
    ###
    if(clegend > 0){
        cha <- as.character(sq0[1])
        yh <- max(c(strheight(cha, cex = cex0), br0*2))
        decalage <- yh / (par("usr")[2]-par("usr")[1]) * (par("usr")[4]-par("usr")[3]) - 0.05
    }
    else decalage <- 0
    ypos <- seq(decalage, 1, length = leaves.number + 2)[(leaves.number:1)+1]
    if(!is.null(y)) ypos <- ypos[y] 
    ###
    
    names(ypos) <- leaves.names
    xcar <- maxx * 1.05
    xx <- c(x.leaves, x.nodes)
    if (clabel.leaves > 0) {
        for (i in 1:leaves.number) {
            text(xcar, ypos[i], leaves.car[i], adj = 0, cex = par("cex") * 
                clabel.leaves)
            segments(xcar, ypos[i], xx[i], ypos[i], col = grey(0.7))
        }
    }
    yleaves <- ypos[1:leaves.number]
    xleaves <- xx[1:leaves.number]

    if (cleaves > 0) {
        for (i in 1:leaves.number) {
            points(xx[i], ypos[i], pch = 21, bg = 1, cex = par("cex") * 
                cleaves)
        }
    }
    yn <- rep(0, nodes.number)
    names(yn) <- nodes.names
    y <- c(ypos, yn)
    for (i in 1:length(x$parts)) {
        w <- x$parts[[i]]
        but <- names(x$parts)[i]
        ypos[but] <- mean(ypos[w])
        b <- range(ypos[w])
        segments(xx[but], b[1], xx[but], b[2])
        x1 <- xx[w]
        y1 <- ypos[w]
        x2 <- rep(xx[but], length(w))
        segments(x1, y1, x2, y1)
    }
    if (cnodes > 0 | !is.null(vnodes)) {
        sq <- vnodes /vmax
        names(sq) <- nodes.names
        for (i in nodes.names) {
            if(length(vnodes)!= length(x$nodes))
                points(xx[i], ypos[i], pch = 21, bg = "white", cex = cnodes)
            else{
                ###
                radi <- sq[i] * reglage
                if(sq[i] > 0)
                symbols(xx[i], ypos[i], circles = radi, 
                   bg = vcolor, fg = "black", add = TRUE, inch = FALSE)
                ###
                }
        }
    }

    if (clabel.nodes > 0) {
        scatterutil.eti(xx[names(x.nodes)], ypos[names(x.nodes)], 
            labels.nodes, clabel.nodes)
    }
    x <- (x.leaves - par("usr")[1])/(par("usr")[2] - par("usr")[1])
    ypos <- ypos[leaves.names]
    xbase <- (xcar - par("usr")[1])/(par("usr")[2] - par("usr")[1])
    if (csub > 0) 
        scatterutil.sub(sub, csub = csub, possub = possub)
    if (draw.box) 
        box()
    if (cleaves > 0) 
        points(xleaves, yleaves, pch = 21, bg = 1, cex = par("cex") * 
            cleaves)

    return(invisible(list(xy = data.frame(x = x, y = ypos), xbase = xbase, 
        cleaves = cleaves)))
}

rtest.decdiv <- function(phy, freq, dis = NULL, nrep = 99, vranking = "complexity", ties.method = "average", option = 1:3, optiontest = NULL, tol = 1e-08)
{

    #*******************************************************************************#
    #                         Checking of the parameters                            #
    #*******************************************************************************#
    
    if(!is.vector(freq)) stop("freq must be a unique vector")
    if (!is.numeric(nrep) | nrep <= 1) 
        stop("Non convenient nrep")
    if(sum(freq) < tol) stop("empty sample")
    if(any(freq < -tol)) stop("negative values in df")
    
    #*******************************************************************************#
    #                               Basic notations                                 #
    #*******************************************************************************#
    
    freq[freq < tol] <- 0
    freq <- freq / sum(freq)
    
    nsp <- length(phy$leaves)
  	nnodes <- length(phy$nodes)
    if(is.null(dis))
    dis <- as.dist(sqrt(2*matrix(1, nsp, nsp) - diag(rep(1, nsp))))
    
    #*******************************************************************************#
    #                               Node ranking                                    #
    #*******************************************************************************#

    complexity <- function(phy){   

	    matno <- as.data.frame(matrix(0, nnodes, nnodes))
	    rownames(matno) <- names(phy$nodes)
	    names(matno) <- names(phy$nodes)
        pathnodes <- phy$path[-(1:nsp)]
	    for(i in 1:nnodes){
	        matno[pathnodes[[i]], i] <- 1
	    }
        listno <- lapply(1:nnodes, function(i) names(matno)[matno[i, ] > 0])
        names(listno) <- names(phy$nodes)
        nbdes <- cbind.data.frame(lapply(phy$parts, function(x) prod(1:length(x))))
        compl <- lapply(listno, function(x) prod(nbdes[x]))
        compltab <- cbind.data.frame(compl)
        compltab <- cbind.data.frame(t(compltab))
        names(compltab) <- "complexity"
        return(compltab)
        
    }

    droot <- function(phy){
        roottab <- cbind.data.frame(phy$droot[-(1:nsp)])
        names(roottab) <- "droot"
        return(roottab)
    }

    if(is.numeric(vranking)){
        vrank <- as.data.frame(rank(vranking, ties.method = ties.method))
        names(vrank) <- "free"
    }
    else
    vrank <- sapply(vranking, function(x) rank(get(x)(phy), ties.method = ties.method))
   
    if(!any(option == 3))
        r1 <- length(option)
    else
        r1 <- length(option) + length(vranking) - 1

    #*******************************************************************************#
    #                       Field observations                                      #
    #*******************************************************************************#
    
    vobs <- decdiv(phy, freq, dis, tol = 1e-08)

    #*******************************************************************************#
    #                       Statistics for the four tests                           #
    #*******************************************************************************#    
    
    namnodes <- rownames(vobs)
    stat1 <- function(v){
        v <- v/sum(v)
        return(max(v))
    }
    stat2 <- function(v){
        v <- v/sum(v)
        fun1 <- function(m){
            return(abs(sum(sort(v)[1:m]) - m/nnodes))
        }
        return(max(unlist(lapply(1:nnodes, fun1))))
    }
    stat3 <- function(v){
        # this statistics has been sightly changed because the consideration of ties in Ollier et al.
        # was not explained, althought ties always happen with such a methodology. 
        funstat3 <- function(vrank1){
            v <- v/sum(v)
            return(sum(rank(vrank1, ties.method = ties.method)*v)/nnodes)
        }
        return(apply(vrank, 2, funstat3))
    }

    methods <- c("stat1", "stat2", "stat3")[option]
    
    #*******************************************************************************#
    #                      Statistics on field observations                         #
    #*******************************************************************************#
    
    statobs <- unlist(sapply(methods, function(x) get(x)(vobs[, 1])))

    #*******************************************************************************#
    #                              Permutation scheme                               #
    #*******************************************************************************#     
    
    funperm <- function(i){
        e <- sample(1:nsp)
        vtheo <- decdiv(phy, freq[e], as.dist(as.matrix(dis)[e, e]), tol = tol)
        stattheo <- unlist(sapply(methods, function(x) get(x)(vtheo[, 1])))
        return(stattheo)
    }
    
    tabsimu <- as.data.frame(t(cbind.data.frame(lapply(1:nrep, funperm))))
    rownames(tabsimu) <- paste("t", 1:nrep, sep="")
    if(r1 == 2 & methods[1] == "stat3")
        names(tabsimu) <- paste("stat3", names(tabsimu), sep=".")
    
    #*******************************************************************************#
    #                                     End                                       #
    #*******************************************************************************# 
    
    optiondefault <- c("greater", "greater", "two-sided", "two-sided", "two-sided")
    names(optiondefault) <- c("stat1", "stat2", "stat3.complexity", "stat3.droot", "stat3.free")
    
    if(r1 == 1)
    {
    if(!is.null(optiontest))
    return(as.randtest(obs = statobs, sim = tabsimu[, 1], alter = optiontest, call = "rtest.decdiv"))
    else
    return(as.randtest(tabsimu[, 1], statobs, alter = optiondefault[names(tabsimu)], call = "rtest.decdiv"))
    }
    
    if(!is.null(optiontest))
    return(as.krandtest(obs = statobs, sim = tabsimu, alter = optiontest, call = "rtest.decdiv"))
    else
    return(as.krandtest(obs = statobs, sim = tabsimu, alter = optiondefault[names(tabsimu)], 
        call = "rtest.decdiv"))

}

################ screeplot.ppca
# function modified from package adephylo
# plot the eigenvalues from a summary object of pPCA
# with a decomposition of their Moran's I value plotted against the variance associated with each eigenvalue

screeplot.ppca <- function (sumry, ...) {
    opar <- par("las")
    on.exit(par(las = opar))
    labels <- lapply(1:nrow(sumry$ppca), function(i) bquote(lambda[.(i)]))
    par(las = 1)
    xmax <- max(sumry$ppca[, "var"]) ## original code had $pca
    I0 <- unlist(sumry$Istat[1])
    Imin <- unlist(sumry$Istat[2])
    Imax <- unlist(sumry$Istat[3])
    cexl <- rep(1, length=nrow(sumry$ppca))
    cexl[sumry$ppca[, "var"]<10^-06] <- 0.7
    plot(x = sumry$ppca[, 2], y = sumry$ppca[, 3], type = "n", 
        xlab = "Variance", ylab = "Moran's I", 
        xlim = c(0, xmax*1.1), ylim = c(Imin * 1.1, Imax * 1.1), 
        yaxt = "n", ...)
    text(x = sumry$ppca[, 2], y = sumry$ppca[, 3], do.call(expression, 
        labels), cex=cexl)
    ytick <- c(I0, round(seq(Imin, Imax, le = 5), 1))
    ytlab <- as.character(round(seq(Imin, Imax, le = 5), 1))
    ytlab <- c(as.character(round(I0, 1)), as.character(round(Imin, 
        1)), ytlab[2:4], as.character(round(Imax, 1)))
    axis(side = 2, at = ytick, labels = ytlab)
    rect(0, Imin, xmax*1.1, Imax, lty = 2)
    segments(0, I0, xmax*1.1, I0, lty = 2)
    abline(v = 0)
    title("Eigenvalues decomposition")
    return(invisible(match.call()))
}

