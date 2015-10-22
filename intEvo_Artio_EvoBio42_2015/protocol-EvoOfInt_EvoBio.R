### R Protocol for reproducing analyses from 
### Haber (2014) "The Evolution of Morphological Integration in the Ruminant Skull"
### Evolutionary Biology
####### PART 1: GENERATING ILMD VCV MATRICES from the A_NkmSym dataset
require(geomorph)

### Requires the following functions, which can be downloaded from DRYAD or gitthub

source("functions-EvoOfInt_EvoBio.R")

### Requires the following objects, which can be downloaded from DRYAD or gitthub
load("A_NkmSymm.Rdata") # data; Procrustes coordinates
load("specimenInfo.Rdata") # information about specimens
load("VerL.Rdata") # error vcv matrices
load("ildefL.Rdata") # definition tables for the interlandmark distances
load("LtI.Rdata") # phylogenetic trees

XX <- Nkm2kmN(A) # converting the data matrix of N specimens by km variables to an array of k landmarks by m dimensions by N specimens

txnm <- LtI[["FV"]]$tip.label # taxa names, same for all trees


# calculating interlandmark distances for each of the definition sets ("IL32", "ILtes")
 

ll <- which(substr(ildefL[["IL32"]][,2],start=nchar(ildefL[["IL32"]][,2])-1,stop=nchar(ildefL[["IL32"]][,2]))=="_L") # ilmd's that span both sides need to be cut in the middle because only half of the symmetric configuration is included

tab.txE <- expand.grid(txnm, names(ildefL)) # all relevant combinations of taxon and ILMD definition; saves a loop
nn <- EL <- list()
for (i in 1:nrow(tab.txE)) {
	txi <- tab.txE[i,1]
	def <- tab.txE[i,2]
	ii <- which(specimenInfo[,"species.FV"]==txi)
	mf.txi <- specimenInfo[ii,"mf"] # sex of each specimen
	spp.txi <- specimenInfo[ii,"tx"] # lowest taxonomic id of each specimen
	EDc <- ED <- t(apply(XX[,,ii], 3, ILMD, LM=ildefL[[def]]))
	if (def=="IL32") {ED[,ll] <- 0.5*ED[,ll]}
	gm <- apply(ED,2,mean) # grand mean of the sample
	for (tj in unique(spp.txi)) {
		# centering each sex within each spp around its own mean 
		# and recentering around grand mean:
		ft <- which(mf.txi=="F" & spp.txi==tj)
		mt <- which(mf.txi=="M" & spp.txi==tj)
		ED[ft,] <- t(t(scale(ED[ft,], scale=FALSE))+gm)
		ED[mt,] <- t(t(scale(ED[mt,], scale=FALSE))+gm)
		}
	EL[[as.character(def)]][[as.character(txi)]] <- ED
	nn[[as.character(txi)]] <- nrow(ED) # sample size
	print(paste(def, txi))}

# Generating one set of P matrices based on mean-standardized data (covariances;"EVms") 
# and one based on variance-standardized data (correlations; "ER") 
# for each ILMD definition (IL32, ILtes)

V2R2 <- list()
tab.EM <- expand.grid(names(EL), c("EVms","ER"))
for (i in 1:nrow(tab.EM)) {
	Ei <- tab.EM[i,1]
	Mi <- tab.EM[i,2]
	for (tx in txnm) {
		X <- EL[[Ei]][[tx]]
		mcol <- colMeans(X)
		V <- var(X)-VerL[[Ei]][[tx]]
		if (Mi=="EVms") M <- shrink(V/mcol%*%t(mcol)) else M <- cov2cor(shrink(V))
		V2R2[[paste(Ei, Mi, sep="")]][[tx]] <- M
	}
}

save(V2R2, file="V2R2.Rdata")
save(nn, file="nn.Rdata")

######################################################################################
####### PART 2: PAIRWISE SIMILARITIES


require(mvtnorm); require(vegan)

source("functions-EvoOfInt_EvoBio.R")
load("V2R2.Rdata")
load("nn.Rdata")

txnm <- names(V2R2[[1]])
Nt <- length(txnm)

simfun <- list(randsk, randsk, matcor, matcor) # random skewers for vcv and matrix correltion for cor matrices
names(simfun) <- names(V2R2)

n.it=1000
n.boot=999

## Calculating repeatability for each P matrix (one for each taxon)

reptL <- list()
for(dt in names(V2R2)) {
	for (tx in txnm) {
		M <- V2R2[[dt]][[tx]]
		p <- ncol(M)
		N <- rep(nn[[tx]], n.boot)		
		VB <- array(sapply(N, function(n){X <- rmvnorm(n=n, sigma=M); shrink(var(X))}), dim=c(p,p,n.boot))
		tB <- apply(VB, 3, simfun[[dt]], V2=M, n.it=n.it)
		zb <- 0.5*log((1+tB)/(1-tB))
		t <- mean(zb)
		reptL[[dt]][[tx]] <- (exp(t/0.5)-1)/(exp(t/0.5)+1) # repeatability; between the observed matrix and every bootstrapped matrix
		rm(VB, M, N)
		#print(c(date(), tx,dt))
		write.table(t(c(date(), "rept", dt, tx)), file="progressP1t1.tab", append=TRUE, sep="\t", col.names=FALSE, row.names=FALSE)
		}
	}

save(reptL, file="reptL.Rdata")

# Generating expected distribution under H0 of no divergence, based on each P matrix
B0L <- list()
for(dt in names(V2R2)) {
	for (tx in txnm) {
		M <- V2R2[[dt]][[tx]]
		p <- ncol(M)
		N <- rep(nn[[tx]], n.boot)		
		VB1 <- array(sapply(N, function(n){X <- rmvnorm(n=n, sigma=M); shrink(var(X))}), dim=c(p,p,n.boot))
		VB2 <- array(sapply(N, function(n){X <- rmvnorm(n=n, sigma=M); shrink(var(X))}), dim=c(p,p,n.boot))
		b0 <- c(); for (b in 1:n.boot) {b0[b] <- simfun[[dt]](VB1[,,b], VB2[,,b])}
		B0L[[dt]][[tx]] <- b0 # between pairs of bootstraped matrices yielding distribution for H0 of no divergence
		rm(VB1, VB2, N)
		#print(c(date(), tx,dt))
		write.table(t(c(date(), "B0", dt, tx)), file="progressP1t1.tab", append=TRUE, sep="\t", col.names=FALSE, row.names=FALSE)
		}
	}

save(B0L, file="B0L.Rdata")

# Generating expected distribution for H0 of no similarity for Random Skewers
rb0L <- list()
for (dt in names(V2R2)[1:2]) {
	# one set of random vectors
	p <- ncol(V2R2[[dt]][[1]])
	B1 <- matrix(rnorm(p*n.it*2, mean = 0, sd = 1), p, (n.it-1)]*2)
	B1 <- t(t(B1)/sqrt(colSums(B1^2))) # scaled to unit length
	# another set 
	B2 <- matrix(rnorm(p*n.it*2, mean = 0, sd = 1), p, (n.it-1)*2)
	B2 <- t(t(B2)/sqrt(colSums(B2^2)))
	# expected distribution between two random vectors
	rb0L[[dt]] <- abs(diag(t(B1)%*%B2)) 
	rm(B1, B2, N)
	}

#### calculating observed similarities and testing them against both H0
IJ <- t(combn(txnm,2))
pwP1P0L <- list()


for (dt in names(V2R2)) {
	pwP1P0 <- rnall <- c()
	p <- ncol(V2R2[[dt]][[1]])
	for (h in 1:nrow(IJ)) {
		tx1 <- IJ[h,1]; tx2 <- IJ[h,2]
		r.obs <- simfun[[dt]](V2R2[[dt]][[tx1]], V2R2[[dt]][[tx2]]) # observed similarity
		rept1 <- reptL[[dt]][[tx1]]; rept2 <- reptL[[dt]][[tx2]]
		r.adj <- r.obs/sqrt(rept1*rept2)
		# testing H0 of no divergence
		bt1 <- B0L[[dt]][[tx1]]; bt2 <- B0L[[dt]][[tx2]]
		p1t1 <- length(bt1[bt1<=r.obs])/(length(bt1)+1)
		p1t2 <- length(bt2[bt2<=r.obs])/(length(bt2)+1)
		# testing H0 of no similarity
		if (dt=="ILtesEVms" | dt=="IL32EVms") {
			#  random vectors for vcv 
			rb0 <- rb0L[[dt]]
			p0 <- length(rb0[rb0>=r.obs])/(n.boot+1)			
			} else {
				# mantel for cor
				p0 <- mantel(V2R2[[dt]][[tx1]], V2R2[[dt]][[tx2]])$signif}		
		pwP1P0 <- rbind(pwP1P0, c(r.obs=r.obs,p1t1=p1t1,p1t2=p1t2, p0=p0, r.adj=r.adj))
		rnall <- c(rnall, paste(tx1,tx2,sep="-"))
		#print(c(date(),h,dt))
		write.table(t(c(date(),"pw", dt, tx1, tx2, h)), file="progressP1t1.tab", append=TRUE, sep="\t", col.names=FALSE, row.names=FALSE)		
		}
	rownames(pwP1P0) <- rnall
	pwP1P0L[[dt]] <- pwP1P0
	}

save(pwP1P0L, file="pwP1P0L.Rdata")
save.image("ws-pwP1P0.R")

###  percentages of significant pw comparisons; Table 1
splt <- sapply(rownames(pwP1P0L[[1]]), strsplit, "-")
pdf <- numeric(length=length(splt)) # identifies the family for pairing
for (l in 1:length(splt)) {
	if (splt[[l]][1]%in%txnm[3:19] & splt[[l]][2]%in%txnm[3:19]) pdf[l] <- 1
	else if (splt[[l]][1]%in%txnm[20:51] & splt[[l]][2]%in%txnm[20:51]) pdf[l] <- 2}

a = 0.001
pwtab <- c()
for (dt in names(pwP1P0L)) {
	M <- pwP1P0L[[dt]]
	p0 <- M[,"p0"]
	psig.all0 <- round(length(p0[p0<=a])*100/length(p0))
	psig.C0 <- round(length(p0[p0<=a & pdf==1])*100/length(p0[pdf==1]))
	psig.B0 <- round(length(p0[p0<=a & pdf==2])*100/length(p0[pdf==2]))
	#p1 <- M[,"p1"]
	p1 <- ifelse(M[,"p1t1"]>M[,"p1t2"], M[,"p1t1"], M[,"p1t1"])
	psig.all1 <- round(length(p1[p1<=a])*100/length(p1))
	psig.C1 <- round(length(p1[p1<=a & pdf==1])*100/length(p1[pdf==1]))
	psig.B1 <- round(length(p1[p1<=a & pdf==2])*100/length(p1[pdf==2]))
	pwtab <- rbind(pwtab, c(paste(psig.all0,psig.all1,sep="/"), paste(psig.C0,psig.C1,sep="/"), paste(psig.B0,psig.B1,sep="/")))
	}
colnames(pwtab) <- c(paste("psig.all (",length(p0),")",sep=""), paste("psig.C (",length(p0[pdf==1]),")",sep=""), paste("psig.B (",length(p0[pdf==2]),")",sep=""))
rownames(pwtab) <- names(pwP1P0L)
write.table(pwtab, file="Table1_Psig.tab", sep="\t", quote=FALSE) # Table 1

save.image("ws-pwP1P0.R")

### rept boxplots Fig.1
quartz(width=6, height=3.5, file="Fig1_reptBXP.jpg", dpi=150, type="jpg", bg="white")
# postscript(width=6, height=3.5, file="Fig1_reptBXP.eps", bg="white")
par(mar=c(4,4,1,1))
boxplot(reptL, names=c("IL32", "ILtes", "IL32", "ILtes"), pars=list(cex.axis=1.1))
mtext(c("Covariances", "Covariance", "Correlations", "Correlations"), side=1, outer=FALSE, line=2, cex=1, at=1:4, adj=0.5)
mtext("Repeatabilities", side=2, outer=FALSE, line=2.5, cex=1.3)
dev.off()

### adj pw boxplots BC Fig. 2
quartz(width=6, height=3.5, file="Fig2_pw-adjBXP-BC.jpg", dpi=150, type="jpg")
#postscript(width=6, height=3.5, file="Fig2_pw-adjBXP-BC.eps", horizontal=TRUE, bg="white")
par(mar=c(4,4,1,1))
plot(x=0:16, xlim=c(0,16), ylim=c(0,1), pch=NA, xaxt="n", yaxt="n", xlab=NA, ylab=NA)
i=0
for (dt in 1:4) {
	boxplot(list(pwP1P0L[[dt]][pdf==0,"r.adj"],pwP1P0L[[dt]][pdf==1,"r.adj"],pwP1P0L[[dt]][pdf==2,"r.adj"]), names=c("A", "C", "B"), pars=list(cex.axis=1), add=TRUE, at=(dt+i):(dt+i+2))
	i=i+3}

mtext(paste(c("IL32", "ILtes"), c("cov", "cov", "cor", "cor"), sep=", "), side=1, outer=FALSE, line=2, cex=1, at=seq(2,16,4), adj=0.5)
mtext("Pairwise similarities", side=2, outer=FALSE, line=2.5, cex=1.3)
dev.off()


# adj vs obs Fig S2
quartz(width=6.5,height=6.5, file="FigS2_pwP0P1obs-adj.jpg", dpi=150, type="jpg")
layout(matrix(1:4,2,2))
par(mar=c(2,2,0.5,0.5), oma=c(2.5,2.5,2.5,2))
for (l in names(pwP1P0L)) {
		M <- pwP1P0L[[l]]
		plot(M[,"r.obs"],M[,"r.adj"], xlab="raw r", ylab="adj r")
		segments(0,0,1,1,col="red")
		}

mtext("Observed pairwise similarities", side=1, line=0.8, outer=TRUE, cex=1.2, adj=0.5)
mtext("Adjusted pairwise similarities", side=2, line=0.8, outer=TRUE, cex=1.2, adj=0.5)
for (i in 1:2) {mtext(c("Covariances","Correlations")[i], side=3, outer=TRUE, cex=1.2, adj=c(0.07, 0.68)[i])}
for (i in 1:2) {mtext(c("IL32","ILtes")[i], side=4,outer=TRUE, cex=1.2, adj=c(0.95, 0.40)[i])}

dev.off()

# Integration similarity matrix for each of the four datasets

SL <- list()
for (dt in names(pwP1P0L)) {
	S <- diag(1, Nt,Nt)
	#S <- diag(sapply(reptL[[dt]], mean), Nt,Nt)
	S[lower.tri(S)] <- pwP1P0L[[dt]][,"r.adj"]; S <- t(S)
	S[lower.tri(S)] <- pwP1P0L[[dt]][,"r.adj"]
	dimnames(S) <- list(txnm,txnm)
	SL[[dt]] <- round(S, 2)
	}

# Convert similarities to distances
DL <-Dz <- Df <- list()
for (dt in names(SL)) {
	S <- SL[[dt]]
	Df[[dt]] <- as.dist(sqrt(2*(1-S))) # Foote's notes + web
	Dz[[dt]] <- as.dist(sqrt(1-S^2)) # Zar's from Jamniczky et al 2009
	}

### check that distances are euclidean/metric
require(ade4)
sapply(Df, is.euclid)
#  IL32EVms ILtesEVms    IL32ER   ILtesER 
#   FALSE     FALSE      TRUE      TRUE 

sapply(Dz, is.euclid)
# IL32EVms ILtesEVms    IL32ER   ILtesER 
#    TRUE     TRUE      TRUE      TRUE 


# Random Skewers isn't euclidean sometimes probably because of the random element
## Shrinking the Random Skewers matrices makes them euclidean
for (l in 1:2) {
	SL[[l]] <- shrink(SL[[l]])}

DL <- Dz <- Df <- list()
for (dt in names(SL)) {
	S <- SL[[dt]]
	Df[[dt]] <- as.dist(sqrt(2*(1-S))) # Foote's notes + web
	Dz[[dt]] <- as.dist(sqrt(1-S^2)) # Zar's alienation (Jamniczky et al 2009)
	DL[[dt]] <- sqrt(2*(1-S)) # Foote's notes + web
	}

sapply(Df, is.euclid)
#  IL32EVms ILtesEVms    IL32ER   ILtesER 
#   TRUE      TRUE       TRUE      TRUE 

save(SL, file="SL.Rdata")
save(DL, file="DL.Rdata")
save.image("ws-pattern_intsim.Rdata")

#############################################################################
########### PART 3: TESTING FOR OVERALL PHYLOGENETIC SIGNAL in D
#### using decdiv from Pavoine et al (2010); Table 2 and Figure S3

require(ade4); require(ape)

source("functions-EvoOfInt_EvoBio.R")
load("LtI.Rdata")
load("DL.Rdata")
load("comptreesL.Rdata") # phylogenetic trees
load("specimenInfo.Rdata") # phylogenetic trees

## Pruning trees to include only taxa with inegration information
require(ape)
LtI <- list()
for (tr in names(comptreesL)) {
	tree <- comptreesL[[tr]]
	txi <- unique(specimenInfo[specimenInfo[,"IntM"]=="IM",paste("species",tr,sep=".")])
	LtI[[tr]] <- drop.tip(tree, tip=which(!tree$tip.label%in%txi)) 
	}

save(LtI, file="LtI.Rdata")


decdivResL <- list()
stats <- c()
for (tr in names(LtI)) {
	tree <- LtI[[tr]]
	#tree$tip.label <- paste("t", 1:51, sep="")
	phy <- newick2phylog(write.tree(tree), add.tools = TRUE)
	for (dt in names(DL)) {
		D <- DL[[dt]]
		for (bl in c("complexity", "droot")) {
			print(c(tr, dt, bl))
			quartz(width=6, height=8, file=paste("FigS3_decdivplot_",dt,tr,".jpg",sep=""), type="jpg", dpi=150)
			plot.decdiv(phy, vnodes=decdiv(phy, rep(1,nrow(D)), as.dist(D)))
			dev.off()
			res <- rtest.decdiv(phy, freq=rep(1,nrow(D)), dis=as.dist(D), nrep = 999, vranking=bl)
			decdivResL[[tr]][[dt]][[bl]] <- res 
			stats <- rbind(stats, cbind(paste(tr, dt, bl, sep="-"), res$obs, res$expvar, res$alter, res$pvalue))}
			}
	write.table(stats, "Table2_decdiv.stats.tab", sep="\t", append=TRUE) # Table 2
	}

save.image("ws-pattern_decdiv.Rdata")

#############################################################################
########### PART 4:  pPCA of INTEGRATION SPACE
#### following Jombart et al (2010)  Bioinformatics Bioinformatics 2010, 26:1907-1909
#### and Jombart et al (2010) Journal of Theoretical Biology 2010, 264:693-701

require(adephylo); require(phytools); require(phylobase)

source("functions-EvoOfInt_EvoBio.R")
load("LtI.Rdata")
load("DL.Rdata")

# generating integration space using Principal Coordinate Analysis
pk=50 # keeping all 50 dimensions initially
eigL <- CL <- list()
for (dt in names(DL)) {
	mds <- cmdscale(DL[[dt]], k=pk, eig=TRUE)
	CL[[dt]] <- mds$points # scores in integration space
	eigL[[dt]] <- mds$eig # eigenvalues of integration space
	}

# screeplot Figure S4
quartz(width=6.5,height=6.5, file="FigS4_screeplotMDS.jpg", dpi=150, type="jpg")
layout(matrix(1:4,2,2))
par(mar=c(2,2,0.5,0.5), oma=c(2.5,2.5,2.5,2))

for (dt in names(eigL)) {
	rev <- 100*eigL[[dt]]/sum(eigL[[dt]])
	plot(rev, lwd=1.5, cex.lab=1.2)
	th <- which(cumsum(rev)>95)[1]
	segments(0,th,51,th, lty=2, lwd=2)}

mtext("Eigenvalue index", side=1, line=0.8, outer=TRUE, cex=1.2, adj=0.5)
mtext("Precent of variance", side=2, line=0.8, outer=TRUE, cex=1.2, adj=0.5)
for (i in 1:2) {mtext(c("Covariances","Correlations")[i], side=3, outer=TRUE, cex=1.2, adj=c(0.07, 0.68)[i])}
for (i in 1:2) {mtext(c("IL32","ILtes")[i], side=4,outer=TRUE, cex=1.2, adj=c(0.95, 0.40)[i])}

dev.off()

# flipping some axes for comparability among datasets
# CL0 <- CL
CL[["ILtesEVms"]][,1] <- -CL[["ILtesEVms"]][,1]
#CL[["IL32EVms"]][,2] <- -CL[["IL32EVms"]][,2]
CL[["IL32EVms"]][,1] <- -CL[["IL32EVms"]][,1]
CL[["IL32ER"]][,1] <- -CL[["IL32ER"]][,1]

save(CL, file="CL.Rdata")

# Figure S5; non-phylogenetic integration space 

x=1;y=2
quartz(width=6.5,height=6.5, file=paste("FigS5_IntSpace_PCo",x,y,"metric.jpg", sep=""), dpi=150, type="jpg")
layout(matrix(1:4,2,2))
par(mar=c(2,2,0.5,0.5), oma=c(2.5,2.5,2.5,2))

for (l in 1:4) {
	C <- CL[[l]]
	plot(C[,x],C[,y], pch=NA, cex.axis=1, cex.lab=1.2)#, xlim=c(-0.3, c(0.4,0.4,0.2,0.2)[l]))
	points(C[1:2,x],C[1:2,y], lwd=1.5,cex=2, pch=8, font=2)
	points(C[3:19,x], C[3:19,y], pch=24, bg="grey50", cex=1.5, lwd=1.5, font=2)
	points(C[20:51,x], C[20:51,y], pch=1, cex=1.5, lwd=1.6, font=2)

mtext("PCo1", side=1, line=0.8, outer=TRUE, cex=1.2, adj=0.5)
mtext("PCo2", side=2, line=0.8, outer=TRUE, cex=1.2, adj=0.5)
for (i in 1:2) {mtext(c("Covariances","Correlations")[i], side=3, outer=TRUE, cex=1.2, adj=c(0.07, 0.68)[i])}
for (i in 1:2) {mtext(c("IL32","ILtes")[i], side=4,outer=TRUE, cex=1.2, adj=c(0.95, 0.40)[i])}
}
dev.off()

#plot(C[,x], C[,y], pch=NA)
#legend("topright", legend=c("Bovids", "Cervids", "Tragulids"), pch=c(1,24,8), pt.bg="grey50", cex=2, pt.lwd=2.5)

###  pPCA 

# scale tree branches by height of tree
for (tr in names(LtI)) {
	tree <- LtI[[tr]]
	tree$edge.length <- tree$edge.length/max(nodeHeights(tree))	
	LtI[[tr]] <- tree
	}


ppcaSumL <- ppcaResL <- list() # summary and full results of pPCA
for (tr in names(LtI)) {
	tree <- LtI[[tr]]
	for (dt in names(CL)) {
		X <- CL[[dt]]
		tr4d <- phylo4d(tree, X)
		for (bl in c("patristic", "Abouheif")) {
			print(c(bl, dt, tr))
			ppcaResL[[tr]][[dt]][[bl]] <- res <- ppca(tr4d, scale=FALSE, scannf=FALSE, nfposi=2, nfnega=1, method=bl)
			ppcaSumL[[tr]][[dt]][[bl]] <- summary(res, printres=TRUE)
			}	
		}
	}


save(ppcaResL, file="ppcaResL.Rdata")
save.image("ws-pattern_cmds-ppca.Rdata")

# Figures 3 and 4; Figure S6, S8, S10, S12, S14

for (tr in names(LtI)) {
	for (bl in c("patristic", "Abouheif")) {
		for (i in 1:2) {
			#quartz(width=7, height=6, file=paste("ppcaEig",c("IL32","ILtes")[i],bl,tr,".jpg",sep="-"), type="jpg", dpi=150)
			postscript(width=7, height=6, file=paste("ppcaEig",c("IL32","ILtes")[i],bl,tr,".eps",sep="-"), bg="white")
			par(oma=c(0,0,3,0))
			layout(matrix(1:4, 2, 2))

			# Covariances
			if (bl=="patristic") par(mar=c(5,4.5,2.8,1.5)) else par(mar=c(3.5,4.5,2.8,1.5))
			res <- ppcaSumL[[tr]][[names(CL)[i]]][[bl]]
			eig <- res$ppca$eig		
			bp <- barplot(eig, beside=TRUE, ylim=c(0,max(eig)*1.2), main="Eigenvalues", cex.axis=1)
			#title(main="Eigenvalues", outer=TRUE, line=1)
			arrows(bp[1],1.25*max(eig), bp[1], 1.05*max(eig), lwd=1.5, length=0.08)
			mtext("G1", side=3, at=bp[1], font=2, line=0, cex=0.9)
			arrows(bp[length(eig)], 0.22*max(eig), bp[length(eig)], 0.05*max(eig), lwd=1.5, length=0.08)
			text(x=bp[length(eig)], y=0.3*max(eig), labels="L1", font=2)
			box("figure")
			par(mar=c(4.5,4.5,2.8,1.5))
			screeplot.ppca(res, cex.lab=1.2) # modified function provided in the accompanied functions file
			box("figure")
			mtext("A. Covariances", side=3, line=1, outer=TRUE, adj=0.02)
			
			#Correlations
			if (bl=="patristic") par(mar=c(5,4.5,2.8,1.5)) else par(mar=c(3.5,4.5,2.8,1.5))
			res <- ppcaSumL[[tr]][[names(CL)[i+2]]][[bl]]		
			eig <- res$ppca$eig		
			bp <- barplot(eig, beside=TRUE, ylim=c(0,max(eig)*1.2), main="Eigenvalues", cex.axis=1)
			#title(main="Eigenvalues", outer=TRUE, line=1)
			arrows(bp[1],1.25*max(eig), bp[1], 1.05*max(eig), lwd=1.5, length=0.08)
			mtext("G1", side=3, at=bp[1], font=2, line=0, cex=0.9)
			arrows(bp[length(eig)], 0.22*max(eig), bp[length(eig)], 0.05*max(eig), lwd=1.5, length=0.08)
			text(x=bp[length(eig)], y=0.3*max(eig), labels="L1", font=2)
			box("figure")
			par(mar=c(4.5,4.5,2.8,1.5))
			screeplot.ppca(res, cex.lab=1.2)
			box("figure")
			mtext("B. Correlations", side=3, line=1, outer=TRUE, adj=0.6)
			dev.off()
			}
		}
	}
		


##### pPCA Scores Figure 5 and Figures S7, S9, S11, S13, S15

load("nn.Rdata")
txnm <- LtI[["FV"]]$tip.label
Nt <- length(txnm)
txnm.full <- names(LtI[["FV"]]$tip.label)

for (tr in names(LtI)) {
	tree <- LtI[[tr]]
	tree$edge.length <- tree$edge.length/max(nodeHeights(tree))
	for(dt in 1:2) {
		#quartz(width=9.4, height=8, file=paste("ppca-",c("IL32","ILtes")[dt],tr,".jpg", sep=""), type="jpg", dpi=300)
		postscript(width=9.4, height=8, file=paste("ppca-",c("IL32","ILtes")[dt],tr,".eps", sep=""), bg="white")
		par(mar=c(2,0,3,1))
		plot.phylo(tree, show.tip.label=FALSE, no.margin=FALSE, x.lim=3.2)
		segments(rep(1.01, Nt), 1:Nt, rep(2.295), 1:Nt, col="grey")
		sc <- cbind(ppcaResL[[tr]][[dt]][["Abouheif"]][["li"]], ppcaResL[[tr]][[dt+2]][["Abouheif"]][["li"]], ppcaResL[[tr]][[dt]][["patristic"]][["li"]], ppcaResL[[tr]][[dt+2]][["patristic"]][["li"]])
		sc <- scale(sc[,-seq(2,12,3)])
		#xi <- c(1.05, 1.2, 1.45, 1.55, 1.75, 1.9, 2.1, 2.25)
		xplus <- c(0, rep(c(0.1, 0.24),3), 0.1)
		xi <- 1.1+cumsum(xplus)
		for (i in 1:ncol(sc)) {
			p <- which(sc[,i]>0)
			n <- which(sc[,i]<0)
			cols <- bgs <- character(length=Nt)
			cols[p] <- bgs[n] <- "white"
			cols[n] <- bgs[p] <- "black"
			points(rep(xi[i], Nt), 1:Nt, pch=21, col=cols, cex=abs(sc[,i])+0.3/1.2, bg=bgs)
		}
		#text(2.3, 1:Nt, labels=paste(txnm.full," (",nn,")", sep=""), cex=0.85, adj=0, vfont=c("sans serif", "bold"))
		text(2.3, 1:Nt, labels=paste(txnm.full," (",nn,")", sep=""), cex=0.85, adj=0)
		segments(xi[c(2,4,6)]+0.12, 1, xi[c(2,4,6)]+0.12, Nt, lty=2)
		mtext(rep(c("Covariances", "Correlations"), 2), at=c(1.07,1.41,1.74,2.08), adj=0.2, cex=0.8, line=-0.5)
		mtext(c("A", "B"), at=c(1,1.68), adj=0, cex=1, line=0.8, font=2)
		mtext(c("Abouheif", "Patristic"), at=c(1.1,1.78), adj=0, cex=1, line=0.8)
		mtext(rep(c("G1", "L1"),4), side=1, at=xi, adj=0.5, cex=0.8, line=-0.7)
		text(0.24, 28.8, labels="Bovidae", cex=1, font=2, adj=0)
		text(0.24, 6.8, labels="Cervidae", cex=1, font=2, adj=0)
		text(0.06, 2.2, labels="Tragulidae", cex=1, font=2, adj=0)
		axis(side=1, pos=0,cex=0.1, at=seq(0, 0.8, 0.2), labels=50-seq(0,40,10), cex.axis=0.7, font=2, tck=-0.01, padj=-2)
		mtext("Time (Mya)", adj=0, at=c(-0.01,-1), side=1, line=0.5, cex=0.8, font=2)
		dev.off()
		}
	}

#### correlation between cophenetic distance and scores on G1

nodes <- c("Bovidae"="15FV", "Cervidae"="25FV")
dts <- c("IL32EVms"="Covariances", "IL32ER"="Correlations")

for (tr in names(LtI)) {
	quartz(width=7, height=6, file=paste("G1~pd_withinFamily",tr ,".jpg",sep=""), type="jpg", dpi=150)
	layout(matrix(1:4,2,2))
	par(mar=c(4.5,4,1,1))
	for (nd in names(nodes)){
		for (dt in c("IL32EVms", "IL32ER")) {
			ndtree <- extract.clade(LtI[[tr]], node=nodes[nd])
			sc <- ppcaResL[[tr]][[dt]][["Abouheif"]][["li"]][ndtree$tip.label,1]
			#sc <- scale(sc)
			pd <- cophenetic.phylo(ndtree)
			pd <- pd[lower.tri(pd)]
			ds <- as.matrix(dist(sc))
			ds <- ds[lower.tri(ds)]
			plot(pd, ds, xlab="Phylogenetic distance", ylab="G1 scores", cex.lab=1.3, cex.axis=1.1)
			legend("topleft", legend=paste(nd, dts[dt]), bty="n", adj=0.1, cex=1.2)	
		}
	}
	dev.off()
}

######## LOADINGS; Figure 6

load("V2R2.Rdata")
load("ildefL.Rdata")
load("ppcaResL.Rdata")

VL <- V2R2[["IL32ER"]] # using only correlations of the IL32 dataset
#VL <- V2R2[["IL32EVms"]]
p <- ncol(VL[[1]])
pname <- colnames(VL[[1]])
p1 <- pname[combn(p,2)[1,]]
p2 <- pname[combn(p,2)[2,]]		

def <- ildefL[["IL32"]]

di <- which(paste(def[,2],def[,1],sep="-")%in%pname)
if (length(di)>0) {	deff <- def
					def[di,1] <- deff[di,2]
					def[di,2] <- deff[di,1]
					} # making sure order of LM in ILMD's in def is consistent with pname

nv <- p*(p-1)/2
Vv <- sapply(VL, function(V){ V[lower.tri(V)]})
Vv[Vv<0] <- -Vv[Vv<0] #taking absolute values of the correlations
a1 <- apply(Vv,1,cor,y=ppcaResL[["FV"]][["IL32ER"]][["Abouheif"]][["li"]][,1]) # loadings
a1 <- a1[i1 <- order(a1,decreasing=TRUE)]
p1 <- p1[i1]
p2 <- p2[i1]
a1 <- round(a1, 2)

m1 <- def[match(p1, paste(def[,1],def[,2],sep="-")),3]
m2 <- def[match(p2, paste(def[,1],def[,2],sep="-")),3]
mm <- paste(m1, m2, sep=" = ")
Rnm <- paste(p1, p2, sep=" = ")
L1 <- cbind(Rnm[1:40], mm[1:40], a1[1:40], Rnm[nv:(nv-39)], mm[nv:(nv-39)], a1[nv:(nv-39)]) # the 40 highest negative and highest positive loadings
write.table(L1, sep="\t", file="TableS6_loadings-abs-ILtes32ER.tab", quote=FALSE)

hsta <- hist(a1, plot=FALSE)
hc <- hsta$counts
hb <- hsta$breaks
tbm1 <- c()
for (i in 2:length(hb)) {
	tbi <- numeric(length(unique(m1)))
	names(tbi) <- unique(m1)
	tb <- table(m1[a1>hb[i-1] & a1<=hb[i]])
	tbi[names(tb)] <- tb
	tbm1 <- rbind(tbm1, tbi)}

tbm2 <- c()
for (i in 2:length(hb)) {
	tbi <- numeric(length(unique(m1)))
	names(tbi) <- unique(m1)
	tb <- table(m2[a1>hb[i-1] & a1<=hb[i]])
	tbi[names(tb)] <- tb
	tbm2 <- rbind(tbm2, tbi)}

tbm <- tbm1+tbm2
tbm <- tbm[,c(4,2,1,6,3,5)]
tbm <- tbm/rowSums(tbm)
#barplot(t(tbm), names.arg=hsta$mids, beside=FALSE, legend=TRUE, col=c("black", "gray25", "gray35", "gray55", "gray80", "gray92"), border=NA)
#barplot(t(tbm), names.arg=hsta$mids, beside=FALSE, legend=TRUE, col=c("black", "gray30", "gray50", "gray65", "gray85", "white"))

## plot with boxplot
quartz(width=7.5, height=6.5, file="module-loadingsIL32ER.jpg", type="jpg", dpi=150)
layout(matrix(c(1,2,2,2),4,1))
par(mar=c(2,7,2,2.6))
boxplot(a1, frame.plot=FALSE, horizontal=TRUE, xaxt="n")
axis(1, at=round(seq(-0.7, 0.5, 0.2),2), line=-1.5, cex.axis=1.4)
mtext("A.", side=3, adj=-0.1, font=2)
par(mar=c(5,7,4,0))
barplot(t(tbm), names.arg=hsta$mids, beside=FALSE, legend=TRUE, col=c("black", "gray30", "gray50", "gray65", "gray85", "white"), args.legend=list(x=15.2, y=0.42, bg="white", cex=1.6), cex.axis=1.6, cex.names=1.6, xlab="Loadings", cex.lab=2)
mtext("B.", side=3, adj=-0.1, line=2, font=2)
title(ylab="Module representation", adj=0.5, line=3.5, cex.lab=2)
dev.off()

# plot without boxplot
quartz(width=10.5, height=6.5, file="Fig6_module-loadingsIL32ER.jpg", type="jpg", dpi=150)
#postscript(width=10.5, height=6.5, file="Fig6_module-loadingsIL32ER.eps", bg="white")
par(mar=c(5,7,0,0), oma=c(0,0,3,0))
barplot(t(tbm), names.arg=hsta$mids, beside=FALSE, legend=TRUE, col=c("black", "gray30", "gray50", "gray65", "gray85", "white"), args.legend=list(x=15.3, y=0.44, bg="white", cex=1.6), cex.axis=1.7, cex.names=1.6, xlab="Loadings", cex.lab=2)
mtext(hc, side=3, at=(hsta$mids+0.96)/1.6, line=0.5, font=2, cex=1.4, outer=TRUE)
title(ylab="Module representation", adj=0.5, line=3.5, cex.lab=2)
dev.off()


save.image("ws-pattern_ppcaLoadings.Rdata")

#################################################################
##### PART 5: Disparity of integration space; Figure 7

require(geiger); require(phytools)

source("functions-EvoOfInt_EvoBio.R")
load("LtI.Rdata")
load("CL.Rdata")

nsims=1000
nodes <- c("Bovidae"="15FV", "Cervidae"="25FV", "Caprinae"="139FV")
DispL <- DispExpL <- list()

for (tr in names(LtI)) {	
	tree <- LtI[[tr]]
	tree$edge.length <- tree$edge.length/max(nodeHeights(tree)) # it doesn't make a difference
	Disp <- array(dim=c(3,length(nodes),length(CL)), dimnames=list(c("obs","expL","expU"),names(nodes), names(CL)))
	DispExp <- array(dim=c(nsims,length(nodes),length(CL)), dimnames=list(1:nsims,names(nodes), names(CL)))
	for (dt in names(CL)) {
		X <- CL[[dt]]
		XS <- sim.char(tree, evorateM(X, tree), nsims, model="BM")
		for (nd in names(nodes)) {
			Xi <- X[tree$tip.label,]
			#if (nd==2 & tr!="FV") {nodes[i]="45FV"}
			obs <- tip.disp(node=nodes[nd], data=Xi, phy=tree, disp="trace")
			DispExp[,nd,dt] <- distexp <- apply(XS, 3, function(Xs){Xx <- Xs[tree$tip.label,]
			tip.disp(node=nodes[nd], data=Xx, phy=tree, disp="trace")})
			expci <- quantile(distexp, probs=c(0.025, 0.975))
			Disp[,nd,dt] <- c(obs, expci)
			print(c(tr,dt,nd))}	
		}
	DispL[[tr]] <- Disp
	DispExpL[[tr]] <- DispExp
	}


save.image("ws-IntDisparity_incsub_scaled.Rdata")

legs <- c("IL32EVms"="IL32, cov", "ILtesEVms"="ILtes, cov", "IL32ER"="IL32, cor", "ILtesER"="ILtes, cor")
for (tr in names(LtI)) {
	quartz(width=10, height=6, file=paste("dispInt.Scaled",tr,".jpg"), type="jpg", dpi=150)
	layout(matrix(1:4, 2, 2))
	par(mar=c(3,4.5,1,1))
	x <- 1:length(nodes)
	for (dt in names(CL)) {
		d.exp <- DispExpL[[tr]][,,dt]
		d.obs <- DispL[[tr]][1,,dt]
		plot(x+0.2, d.obs, pch="*", cex=2.5, lwd=2, xaxt="n", xlab=NA, xlim=c(0.5,length(nodes)+0.7), font=2, ylim=c(0, max(d.exp)), cex.axis=1, ylab="Integration disparity", cex.lab=1.3)
		for (i in 1:length(x)) {boxplot(d.exp[,i], at=x[i], boxwex=0.4, add=TRUE)}
		#points(x+0.1, d.obs, pch="*", cex=1, lwd=2)
		legend("topright", legend=legs[dt], cex=1.2, bty="n")
		axis(1, at=x, tick=FALSE, labels=names(nodes), cex.axis=1.2, font=2, line=-0.5)
		}
	dev.off()}
	

#quartz(width=5.5, height=3.5, file="Fig7_IntDisp_incsub_FV32ER.jpg", type="jpg", dpi=150)
postscript(width=5.5, height=3.5, file="Fig7_IntDisp_incsub_FV32ER.eps", bg="white")
par(mar=c(3,4.5,2.5,1))
tr <- "FV"
dt <- "IL32ER"
x <- 1:length(nodes)
		d.exp <- DispExpL[[tr]][,,dt]
		d.obs <- DispL[[tr]][1,,dt]
		plot(x+0.2, d.obs, pch="*", cex=2.5, lwd=2, xaxt="n", xlab=NA, xlim=c(0.5,length(nodes)+0.7), font=2, ylim=c(0, max(d.exp)), cex.axis=1, ylab="Integration disparity", cex.lab=1.3)
		for (i in 1:length(x)) {boxplot(d.exp[,i], at=x[i], boxwex=0.4, add=TRUE)}
		axis(1, at=x, tick=FALSE, labels=names(nodes), cex.axis=1.2, line=-0.5)
dev.off()

