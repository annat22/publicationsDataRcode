### R Protocol reproducing analyses from 
### Haber (2015) "The role of intrinsic constraints in the diversification of ruminants"

############################################################################################ 
####### PART 1: GENERATING ILMD VCV MATRICES AND SPECIES MEANS from the A_NkmSym dataset
require(ape)

### Requires the following functions, which can be downloaded from DRYAD or gitthub

source("functions-MdivII.R")

### Requires the following objects, which can be downloaded from DRYAD or gitthub
load("A_NkmSymm.Rdata") # data; Procrustes coordinates
load("specimenInfo.Rdata") # information about specimens
load("ildefL.Rdata") # definition tables for the interlandmark distances
load("comptreesL.Rdata") # phylogenetic trees, complete

### Pruning trees to include only species with shape and vcv information
# one set for species with vcv as well as mean shape (for the within-pop variation analysis)
# and one set with all species with mean shape info with or without vcv 
# (for the among-pop morphological divergence analysis)
# the only difference between the trees in terms of terminal taxa 
# is that Cervus elaphus is split into two species in CT and MR trees  

LtI <- LtM <- list()

for (tr in names(comptreesL)) {
	tree <- comptreesL[[tr]]
	sptr <- paste("species", tr, sep=".")
	txi <- droplevels(unique(specimenInfo[specimenInfo$IntM=="IM",sptr])) # species with covariance matrix as well as mean shape
	LtI[[tr]] <- drop.tip(tree, tip=which(!tree$tip.label%in%txi))
	txm <- droplevels(unique(specimenInfo[,sptr])) # all species with mean shape, with or without vcv
	LtM[[tr]] <- drop.tip(tree, tip=which(!tree$tip.label%in%txm))
	}
rm(tr)

save(LtM, file="LtM.Rdata")
save(LtI, file="LtI.Rdata")

### calculating interlandmark distances for each of the definition sets ("IL32", "ILtes")

AA <- Nkm2kmN(A) # converting the data matrix of N specimens by km variables to an array of k landmarks by m dimensions by N specimens

EDL <- list("IL32"=t(apply(AA, 3, ILMD, LM=ildefL[["IL32"]])), "ILtes"=t(apply(AA, 3, ILMD, LM=ildefL[["ILtes"]])))

ll <- which(substr(ildefL[["IL32"]][,2],start=nchar(ildefL[["IL32"]][,2])-1,stop=nchar(ildefL[["IL32"]][,2]))=="_L") # ilmd's that span both sides need to be cut in the middle because only half of the symmetric configuration is included; see below

nnL <- XmL <- EL <- list()
tab <- expand.grid(names(LtM), names(EDL)) # all combinations of trees and ILMD definition; saves a loop
for (i in 1:nrow(tab)) {
	tr <- as.character(tab[i,1])
	txnm <- LtM[[tr]]$tip.label # taxa names
	def <- as.character(tab[i,2])
	ED <- EDL[[def]]
	if (def=="IL32") {ED[,ll] <- 0.5*ED[,ll]} # cutting in half the ilmd's that span both sides
	Xm <- c()
	for (tx in txnm) {
		sptr <- paste("species", tr, sep=".")
		ii <- which(specimenInfo[,sptr]==tx) # all specimens in this species
		mf.txi <- specimenInfo[ii,"mf"] # sex of each specimen
		spp.txi <- specimenInfo[ii,"tx"] # lowest taxonomic id of each specimen
		E <- as.matrix(ED[ii,])
		if (length(ii)==1) E <- t(E)
		gm <- apply(E,2,mean) # grand mean of the species
		for (tj in unique(spp.txi)) {
			# centering each sex within each subspecies around its own mean 
			# and recentering around grand mean:
			ft <- which(mf.txi=="F" & spp.txi==tj)
			mt <- which(mf.txi=="M" & spp.txi==tj)
			if (length(ft)>1) E[ft,] <- t(t(scale(E[ft,], scale=FALSE))+gm)
			if (length(mt)>1) E[mt,] <- t(t(scale(E[mt,], scale=FALSE))+gm)
			}
		EL[[tr]][[def]][[tx]] <- E
		nnL[[tr]][tx] <- nrow(E) # sample size
		Xm <- rbind(Xm, colMeans(E))
		print(paste(tr, def, txi))
		}
	rownames(Xm) <- txnm
	XmL[[tr]][[def]] <- Xm
	}

save(XmL, file="XmL.Rdata")
#### calculating vcv matrices for species with sample size > 26
# taxa are the same for all three trees.
# however, the specimens included in Cervus elaphus are different for CT and MR because it is split

nn <- EpL <- PL <- list()
for (tr in names(LtI)) {
	for (def in names(EDL)) {
		txnm <- LtI[[tr]]$tip.label # taxa with >26 specimens
		for (tx in txnm) {
			E <- EL[[tr]][[def]][[tx]]
			ii <- which(specimenInfo[rownames(E),"IntM"]=="IM")
			Ei <- scale(E[ii,], center=TRUE, scale=colMeans(E[ii,])) # mean scaled
			EpL[[tr]][[def]][[tx]] <-	Ei
			PL[[tr]][[def]][[tx]] <- shrink(var(Ei))
			nn[[tr]][[tx]] <- length(ii)
		}
	}
}

save(EpL, file="EpL.Rdata") # mean scaled data for each taxon (that has P) in each tree
save(PL, file="PL.Rdata") # within-pop vcv for each taxon in each tree
save(nn, file="nn.Rdata") # sample size for each taxon in each tree
save.image("ws-MdivII_EVL.Rdata")

##################################################################################
###### PART 2; Within-pop variation properties (e.bar, c.bar, f.bar, rSDE(P))

source("functions-MdivII.R")
load("EpL.Rdata")

theta.ii <- function(x, Xdata, statfun, permcol) {
	Xth <- Xdata[x,]
	if (permcol==TRUE) {
		Xth <- apply(Xth, MARGIN=2, sample, size=nrow(Xth), replace=FALSE)
		Xth <- apply(Xth, MARGIN=2, sample, size=nrow(Xth), replace=FALSE)} # one round doesn't seem to randomize it properly
	statfun(shrink(var(Xth)))}

ifunsL <- list("rSDE"=rSDE, "AvFlex"=AvFlex, "AvE"=AvE, "AvCondE"=AvCondE)

n.perm=1000
n.boot=999

BresL <- IIL <- list()

for (tr in names(EpL)){
	for (dt in names(EpL[[tr]])) {
		II <- array(NA, dim=c(length(EpL[[tr]][[dt]]), 10, length(ifunsL)), dimnames=list(names(EpL[[tr]][[dt]]), c("obs","BCa.L","BCa.U", "bL", "bU", "pm","pL","pU","N","Iadj"), names(ifunsL)))
		for (iif in names(ifunsL))	{	
			for (tx in names(EpL[[tr]][[dt]])) {
				print(paste(c(date(), dt, iif, tx))
				#write.table(t(c(date(), dt, iif, tx)), file="progress_II.tab", sep="\t", append=TRUE, row.name=FALSE, col.name=FALSE)
				E <- EpL[[tr]][[dt]][[tx]]
				if (iif=="rSDE") {
					permsam <- matrix(1:nrow(E), nrow(E), n.perm)
					permdist <- apply(X=permsam, MARGIN=2, FUN=theta.ii, Xdata=E, statfun=ifunsL[[iif]], permcol=TRUE)
					pm <- mean(permdist)
					pI <- quantile(permdist, probs=c(0.025,0.975))
					Bres <- bootNP.BCa(1:nrow(E), nboot=n.boot, theta.fun=theta.ii, Xdata=E, statfun=ifunsL[[iif]], permcol=FALSE)
					bI <- quantile(Bres$thetastar, probs=c(0.025,0.975)) # regular confidence intervals
					io <- Bres$thetahat
					bc <- Bres$confpoints # BCa corrected CI's
					i.adj <- io-pm
					BresL[[dt]][[iif]][[tx]] <- Bres
					rm(Bres)
				} else {
					Bres <- bootNP.BCa(1:nrow(E), nboot=n.boot, theta.fun=theta.ii, Xdata=E, statfun=ifunsL[[iif]], permcol=FALSE)				
					bI <- quantile(Bres$thetastar, probs=c(0.025,0.975))
					io <- Bres$thetahat
					bc <- Bres$confpoints
					i.adj <- pm <- NA; pI <- c(NA, NA)
					BresL[[dt]][[iif]][[tx]] <- Bres
					rm(Bres)}
				II[tx,,iif] <- c(io, bc, bI, pm, pI, nrow(E), i.adj)
				}
			}	
	IIL[[tr]][[dt]] <- II
	}
}

save(IIL, file="IIL.Rdata")
save.image(file="ws-II_bootBCA+perm.R")

##### Comparing Il32 to ILtes; FigureA3
 
#load("IIL.Rdata")

#axlb <- c("rSDE"="rSDE(P)", "AvFlex"=expression(bar("f")), "AvE"=expression(bar("e")), "AvCondE"=expression(bar("c")))

iflb <- c("rSDE"="rSDE(P)", "AvFlex"="Average flexibility", "AvE"="Average evolvability", "AvCondE"="Average cond. evolvability")

quartz(width=4, height=10, file="FigureA3_ILtes~IL32_IECF.jpg", type="jpg", dpi=150)
layout(matrix(1:4, 4, 1))
par(mar=c(4.5,5,1,1))

for (iif in names(iflb)) {
	X <- IIL[["FV"]][[1]][,1:3,iif]
	Y <- IIL[["FV"]][[2]][,1:3,iif]
	plot(X[,"obs"], Y[,"obs"], xlab="IL32", ylab="ILtes", pch=NA, ylim=range(Y, na.rm=TRUE), xlim=range(X, na.rm=TRUE), cex.axis=1.5, cex.lab=1.5)
	abline(lm(Y[,"obs"] ~ X[,"obs"]))
	if (iif!="AvCondE"){
		segments(X[,"obs"], Y[,"BCa.L"], X[,"obs"], Y[,"BCa.U"], col="grey60", lwd=1)
		segments(X[,"BCa.L"], Y[,"obs"], X[,"BCa.U"], Y[,"obs"], col="grey60", lwd=1)
		}	
points(X[,"obs"], Y[,"obs"], pch=21, bg="grey80", cex=1.6)
	legend("bottomright", legend=iflb[iif], bty="n", cex=1.5)
}

dev.off()


#### observed vs adjusted rSDE(P); Figure A4

quartz(width=4, height=6, file="Figure A4_rSDE_obs~adj.jpg", type="jpg", dpi=150)
layout(matrix(1:2, 2, 2))
par(mar=c(5, 5, 1, 1))
for (i in 1:2) {
	II <- IIL[["FV"]][[i]][,,"rSDE"]
	plot(II[,"obs"],II[,"Iadj"], xlab="rSDE; observed", ylab="rSDE; adjusted", cex.axis=1.2, cex.lab=1.4, cex=1.4, pch=21, bg="grey80")
	abline(lm(II[,"Iadj"]~II[,"obs"]))
	legend("bottomright", legend=c("IL32","ILtes")[i], bty="n", cex=1.5)}
dev.off()

### Comparing the different matrix properties; Figure 1

combs <- cbind(t(expand.grid("rSDE", names(iflb)[-1])), c("AvCondE", "AvE"))

quartz(width=4, height=10, file="Figure2_IECFcomc_32EVms.jpg", type="jpg", dpi=150)
layout(matrix(1:4, 4, 1))
par(mar=c(4,5,1.5,1.5))


for (j in 1:ncol(combs)) {
	xe <- combs[1,j]
	ye <- combs[2,j]
		X <- IIL[["FV"]][["IL32EVms"]][,1:3,xe]
		Y <- IIL[["FV"]][["IL32EVms"]][,1:3,ye]
		plot(X[,"obs"], Y[,"obs"], xlab=iflb[xe], ylab=iflb[ye], pch=NA, ylim=range(Y, na.rm=TRUE), xlim=range(X, na.rm=TRUE), cex.axis=1.5, cex.lab=1.5)
		abline(lm(Y[,"obs"] ~ X[,"obs"]))
		segments(X[,"obs"], Y[,"BCa.L"], X[,"obs"], Y[,"BCa.U"], col="grey60", lwd=1)
		segments(X[,"BCa.L"], Y[,"obs"], X[,"BCa.U"], Y[,"obs"], col="grey60", lwd=1)
		points(X[,"obs"], Y[,"obs"], pch=21, bg="grey80", cex=1.6)
		#mtext(c("A","B","C","D")[j], side=3, line=1, adj=-0.25, font=2)
		}
dev.off()
	
#### Fitting evolutionary models to matrix properties; Tables 2 and 3

require(maticce); require(geiger)

load("IIL.Rdata")
load("LtI.Rdata")
iflb <- c("rSDE"="rSDE(P)", "AvE"="Average evolvability", "AvCondE"="Average cond. evolvability")


tx <- rownames(IIL[[1]][[1]])
nodes <- c("139FV", "15FV", "27FV", "25FV", "8FV") # effectively "regimes"
txndL <- list()
for (ni in nodes) {
	ndn <- which(LtI[[1]]$node.label==ni)+length(tx)
	tipsndn <- tips(LtI[[1]], node=ndn)
	names(tipsndn) <- NULL
	txndL[[ni]] <-tipsndn
	} # which tip taxa belong to each node/regime

for (tl in names(LtI)) {LtI[[tl]]$edge.length <- LtI[[tl]]$edge.length/50 } # scaled to tree height
trees <- lapply(list("FV"=LtI[[1]], "CT"=LtI[[2]]), ape2ouch) # converted to ouch trees

FitL <- THmL <- parsL <- list()

for (tr in names(trees)) {
	for (dt in names(IIL[[tr]])) {
		fit <- c()
		for (iif in names(iflb)) {
			xi <- IIL[[tr]][[dt]][,"obs",iif] # character to be modeled
			if (iif=="rSDE") {xi <- xi - IIL[[tr]][[dt]][,"pm",iif]} # adjusting rSDE by sample size

			# fitting OU
			resOU <- runBatchHansen(ouchTrees=trees[[tr]], characterStates=log(xi), cladeMembersList=txndL)
			ou.ic <- informationCriterion.hansenBatch(resOU)[[1]]
			Mod <- resOU$regMatrix$overall # matrix specifying regimes for each model
			nm <- nrow(Mod)
	
			# fitting BM
			resBM <- fitContinuous(multi2di(LtI[[tr]]), log(xi), SE = 0, model = "BM")	
	
	 		# combining OU and BM results for calculating weights
			aicc.all <- c(ou.ic$AICc, resBM$opt$aicc)
			rel.lk <- exp(-0.5*(aicc.all-min(aicc.all)))
			weight.all <- round(rel.lk/sum(rel.lk), 4) 

			# organizing results table
			ou.fit <- cbind(round(ou.ic$AICc,2), weight.all[1:nrow(Mod)])
			bm.fit <- c(round(resBM$opt$aicc,2), weight.all[nrow(Mod)+1])
			fit <- cbind(fit, rbind(ou.fit, "BM"=bm.fit))
			
			Th <- resOU$theta[[1]]
			trdf <- as(trees[[tr]],"data.frame")
			colnames(Th) <-  as.vector(trdf[["labels"]])
			thm <- colSums(Th[,tx]*weight.all[1:nrow(Mod)]) # weighted averaged theta values for all nodes and tips
			THmL[[tr]][[dt]][[iif]] <- thm[tx]
			
			alpha <- resOU$hansens[[1]][,"theta / alpha"]
			sigma.sq <- resOU$hansens[[1]][,"sigma.squared"]
			dec1 <- c(); for (l in nodes) {dec1[l] <- txndL[[l]][[1]]}
			dec1["8FV"] <- "25FV"
			pars <- cbind(exp(Th[,dec1]), alpha, sigma.sq) # parameters for all models for each element in <nodes>
			colnames(pars) <- c(nodes, "alpha", "sigma2")
			pars <- rbind(pars, "Average"=c(exp(thm[dec1]), sum(alpha*weight.all[1:nrow(Mod)]), sum(sigma.sq*weight.all[1:nrow(Mod)])))
			pars <- rbind(pars, "BM"=c(rep(exp(resBM$opt$z0),6), round(resBM$opt$sigsq,2)))
			parsL[[tr]][[dt]][[iif]] <- pars
			write.table(pars[,c(5:1,6:7)], file=paste("pars",tr,dt,iif,".tab",sep=""), sep="\t") # table 3
			}
		Fit <- cbind(c(ou.ic$K, resBM$opt$k), rbind(Mod, rep(0,length(nodes))), fit)
		dimnames(Fit) <- list(c(1:nm,"BM"), c("K", nodes, rep(c("AICc", "AICcwi"), 3)))
		write.table(Fit[,c(1,6:2,7:12)], file=paste("fit",tr,dt,".tab",sep="_"), sep="\t") # table 2
		FitL[[tr]][[dt]] <- Fit
		}
	}

save(THmL, file="THmL.Rdata")
save(FitL, file="FitL.Rdata")
save.image("ws-maticce_IILadj.Rdata")
	
##### FIGURE 2
## Plotting rSDE(P) against the FV tree
## rSDE(P) values are adjusted for sample size by subtracting the value of a permuted matrix
require(ape)
load("LtI.Rdata")
load("IIL.Rdata")
load("THmL.Rdata")

nn <- IIL[["FV"]][["IL32EVms"]][,,"rSDE"][,"N"]
II <- IIL[["FV"]][["IL32EVms"]][,,"rSDE"]+1 # for ploting against the tree, needs to start at 1 instead of 0
th <- exp(THmL[["FV"]][["IL32EVms"]][["rSDE"]])+1

tree <- LtI[[1]]
tree$edge.length <- tree$edge.length/50 # scaled by tree height

tx <- names(tree$tip.label)
tx[10:15] <- c("Odocoileus h. californicus", "Odocoileus h. columbianus", "Odocoileus h. hemionus", "Odocoileus v. borealis", "Odocoileus v. couesi", "Odocoileus v. leucurus")
Nt <- length(tx)

quartz(width=7.7, height=8.5, file="Figure2_32EVms_adj.jpg", type="jpg", dpi=150)
par(mar=c(2.3,0,0,0))
plot.phylo(tree, show.node.label=FALSE, show.tip.label=FALSE, use.edge.length=TRUE, edge.width=2.5, x.lim=2.4)
segments(rep(1,Nt), 1:Nt, rep(1.5, Nt), lwd=0.5, col="grey30")
points(II[,"Iadj"], 1:Nt, pch=19, font=2)
segments(II[,"BCa.L"]-II[,"pm"]+1, 1:Nt, II[,"BCa.U"]-II[,"pm"]+1, 1:Nt, lwd=2.2, font=2)
segments(th, 1:51, th, c(2:51,51), lty=3, col="grey5", lwd=2, font=2)

axis(side=1, pos=0, cex=0.1, at=seq(1, 1.5, 0.1), labels=seq(0,0.5,0.1), cex.axis=0.9, tck=-0.01, padj=-1.2)
mtext("rSDE(P)", adj=0, at=c(1.2,-1), side=1, line=0.8, cex=1)

text(1.53, 1:Nt, labels=paste(tx," (",nn,")", sep=""), cex=0.9, adj=0)

axis(side=1, pos=0, cex=0.1, at=seq(0, 0.8, 0.2), labels=50-seq(0,40,10), cex.axis=0.9, tck=-0.01, padj=-1.2)
mtext("Time (Mya)", adj=0, at=c(0,-1), side=1, line=0.8, cex=1)

text(0.20, 28.5, labels="Bovidae", cex=1, adj=0)
text(0.29, 6.8, labels="Cervidae", cex=1, adj=0)
text(0.06, 2.2, labels="Tragulidae", cex=1, adj=0)
text(0.21, 18.5, labels="8FV", cex=1, adj=0)
text(0.44, 18.2, labels="Cervinae", cex=0.85, adj=0)
text(0.39, 49.6, labels="Caprinae", cex=1, adj=0)

dev.off()

###############################################################################################################
#########  PART 3: Morphological diversification

require(geiger); require(ape); require(mvtnorm)

source("functions-MdivII.R")
load("XmL.Rdata")
load("LtM.Rdata")
load("THmL.Rdata")


# a parametric sampling function
# N is sample size, M is the covariance matrix underlying the sampling distribution
# shtol is shrinking tolerance for the shrink function
parsamp <- function(N, M, shtol=10^-8) {
		Xs <- rmvnorm(N, sigma=M)
		shrink(var(Xs), tol=shtol)
		}

# tip disparity calculated as the mean squared deviation from the root
tip.disp <- function (Xt, phy) {
 	a0 <- aceRoot(phy, Xt)  # root state 
 	Xto <-t(t(Xt)-a0)
   sum(diag(Xto%*%t(Xto)))/(nrow(Xto)-1)
   }


n.boot=500
nodes <- c("Cervidae_woA"="25FV", "Bovidae_woCp"="15FV", "Bovidae"="15FV", "Cervidae"="25FV", "Cervinae"="27FV", "Caprinae"="139FV") # Cervida is analysed with and without Alce, and Bovida is analysed with and without Caprinae

## calculating the among-pop vcv, and its properties
ApopEL <- DicL <- list()

for (tr in names(LtM)) {
	tree <- LtM[[tr]]
	for (dt in names(XmL[[tr]])) {	
		tx <- tree$tip.label
		X <- XmL[[tr]][[dt]][tx,]
		Disp <- Rates <- Eccent <- matrix(nr=length(nodes),nc=3, dimnames=list(names(nodes), c("obs","ciL","ciU")))		
		for (nd in names(nodes)) {
			print(c(tr,dt,nd))
			ndtree <- extract.clade(tree, node=nodes[nd])
			if (nd=="Cervidae_woA") {ndtree <- drop.tip(ndtree, tip=c("Alc_al2","Alc_am4"))}
			if (nd=="Bovidae_woCp") {ndtree <- drop.tip(ndtree, tip=prop.part(ndtree)[[which(ndtree$node.label=="139FV")]])}
			txnd <- ndtree$tip.label
			Xnd <- X[txnd,]
			a0 <- aceRoot(ndtree, Xnd) # clade's phylogenetically-weighted mean; inferred root state
			Xnd <- t(t(Xnd)/a0) # scale matrix of species mean by the clade's mean
			Dic <- shrink(evorateM(Xnd, ndtree), tol=10^-6) # among-pop matrix (Dic in Table 1)
			XS <- sim.char(ndtree, Dic, nsim=n.boot, model="BM")
			disp.expdist <- apply(XS, 3, function(Xs){tip.disp(Xs, ndtree)})
			Disp[nd,] <- c(tip.disp(Xnd, ndtree), quantile(disp.expdist, probs=c(0.025, 0.975)))
			BDic <- array(sapply(rep(length(txnd),n.boot), parsamp, M=Dic), dim=c(nrow(Dic),ncol(Dic),n.boot))
			rates.dist <- apply(BDic, 3, function(V){sum(diag(V))})
			Rates[nd,] <- c(sum(diag(Dic)), quantile(rates.dist, probs=c(0.025, 0.975)))
			Eccent.dist <- apply(BDic, 3, rSDE)
			Eccent[nd,] <- c(rSDE(Dic), quantile(Eccent.dist, probs=c(0.025, 0.975)))
			DicL[[tr]][[dt]][[nd]] <- Dic
			print(c(tr, dt, nd))
			rm(ndtree, BDic, disp.expdist, rates.dist, Eccent.dist)
			}
		ApopEL[[tr]][[dt]][["Disp"]] <- Disp # average distance from clade's root (same as matrix size of the uncorrected among-pop vcv, D)
		ApopEL[[tr]][[dt]][["Rates"]] <- Rates # matrix size of Dic
		ApopEL[[tr]][[dt]][["Eccent"]] <- Eccent # matrix shape of Dic
		}
	}

save(DicL, file="DicL.RData")
save.image("ws-Mdiv.all.RData")

### FIGURE 3; disparity, rates, and eccentrcity of among-pop matrix (Dic)

#tr="FV"; dt="IL32"

dec1 <- c("Cervidae"="Hyd_in9", "Cervinae"="Rus_un1", "Bovidae"="Bis_bi9", "Caprinae"="Cap_ib8")

ylb <- c("Disparity", "Rate of dispersion", "Eccentricity of dispersion"); names(ylb) <- names(ApopEL[[1]][[1]])
#names(ylb) <- names(ApopEL[[1]])

x <- c(1, 1.25, 2, 2.25, 3:4)


for (tr in names(LtM)[1:2]) {
	for (dt in names(XmL[[tr]])) {	
		th.ii <- round(exp(THmL[[tr]][[dt]][["rSDE"]][dec1]),2)
		th.e <- round(exp(THmL[[tr]][[dt]][["AvE"]][dec1]),4)
		names(th.ii) <- names(th.e) <- names(dec1)
		th.ii <- sort(th.ii)
		th.e <- th.e[names(th.ii)]

		ci <- which(names(th.ii)=="Cervidae")
		xord <- c(names(th.ii)[1:ci], "Cervidae_woA", names(th.ii)[-(1:ci)])
		bi <- which(xord=="Bovidae")
		xord <- c(xord[1:bi], "Bovidae_woCp", xord[-(1:bi)])
		quartz(width=5.2, height=7, file=paste("Figure3_ApopE_",dt,tr,".jpg", sep=""), type="jpg", dpi=150)
		layout(matrix(1:3, 3, 1))
		par(mar=c(0,5,1,1), oma=c(6,0,0,0))
		for (ef in names(ApopEL[[tr]][[dt]])) {
			X <- ApopEL[[tr]][[dt]][[ef]]
			plot(x, ylim=c(min(X)*0.75,max(X)*1.05), pch=NA, xaxt="n", xlab=NA, xlim=c(0.5,4.5), cex.axis=1.2, ylab=ylb[ef], cex.lab=1.6)
			points(x, X[xord,"obs"], pch=c(4,1,4,1,4,4), cex=1.8, lwd=2)
			if (ef=="Disp") segx <- x-0.1 else segx <- x
			segments(segx, X[xord,"ciL"], segx, X[xord,"ciU"], lwd=2)
			}
			axis(1, at=1:4, tick=FALSE, labels=names(th.ii), cex.axis=1.9, line=-0.4, outer=TRUE)
			axis(1, at=1:4, labels=th.ii, tick=FALSE, cex.axis=1.5, line=1.3, lty=0, outer=TRUE)
			mtext("rSDE(P)", 1, line=2.3, at=0.16, outer=TRUE, cex=0.8, adj=1)
			axis(1, at=1:4, labels=th.e, tick=FALSE, cex.axis=1.5, line=2.9, lty=0, outer=TRUE)
			mtext(expression(bar("e")), 1, line=3.8, at=0.16, outer=TRUE, cex=0.9, adj=1)
		dev.off()
		}
	}

###############################################################################################################
#########  PART 4: Evolvability along specific directions of divergence, accounting for orientation of Dic

require(ape)
require(phytools)
 
source("functions-MdivII.R")

load("LtM.Rdata") # three alternative trees
load("PL.Rdata") # P
load("XmL.Rdata")
load("nn.Rdata")
load("DicL.RData")

nodes <- c("Cervidae_woA"="25FV", "Bovidae_woCp"="15FV", "Bovidae"="15FV", "Cervidae"="25FV", "Caprinae"="139FV")


######## FIGURE 4
### e(d.cl) in Table 1
###  Dic ~ P; evolvabilities along eigenvectors of among-pop Dic; P is the species within-pop P


EvDicPspL <- list()

for (tr in names(LtM)) {
	tree <- LtM[[tr]]
	for (dt in names(DicL[[tr]])) {	
		PLdt <- PL[[tr]][[dt]] # using covariances only for each dataset
		for (nd in names(nodes)){
			ndtree <- extract.clade(tree, node=nodes[nd]) # specific clade only; one clade at a time
			if (nd=="Cervidae_woA") {ndtree <- drop.tip(ndtree, tip=c("Alc_al2","Alc_am4"))}
			if (nd=="Bovidae_woCp") {ndtree <- drop.tip(ndtree, tip=prop.part(ndtree)[[which(ndtree$node.label=="139FV")]])}
			txnd <- ndtree$tip.label # species within the clade who have morphological information
			txnd <- txnd[which(txnd%in%names(PLdt))]
			Dic <- DicL[[tr]][[dt]][[nd]]
			eig <- eigen(Dic)
			ev <- 100*eig$values/sum(eig$values)
			q <- which(cumsum(ev)<96) # keeping 95% of variation in Dic
			tab <- expand.grid(txnd, q)
			E <- matrix(nr=length(txnd), nc=length(q))
			colnames(E) <- as.character(round(ev[q]))
			rownames(E) <- txnd
			for (i in 1:nrow(tab)) {
				P <- PLdt[[tab[i,1]]]
				x <- eig$vectors[,tab[i,2]]
				E[tab[i,1], tab[i,2]] <- t(x)%*%P%*%x
				rm(P)
				}
			EvDicPspL[[tr]][[dt]][[nd]] <- E
			}
		}
	}

save.image("ws-EvDicPsp.RData")
	
## Plots

for (tr in names(LtM)) {
	for (dt in names(EvDicPspL[[tr]])) {	
		quartz(width=15, height=3.8, file=paste("EvVePsp", tr, dt, "Fig.jpg", sep="_"), type="jpg", dpi=150)
		layout(matrix(1:3,1,3))
		par(mar=c(5.5,0,1,1), oma=c(0,7,0,0))
		ymx <- max(sapply(EvDicPspL[[tr]][[dt]], max))
		for (nd in names(nodes)[-(1:2)]) {
			E <- EvDicPspL[[tr]][[dt]][[nd]]
			boxplot(as.data.frame(E), yaxt="n", cex.axis=2.2, cex.lab=2.4, ylim=c(0,ymx*1.05), ylab=NA, xlab="Eigenvalue (%)")
			legend("topright", legend=nd, cex=2.5, bty="n")
			}
		axis(2, outer=TRUE, line=0, cex.axis=2.2)
		mtext("Evolvability", 2, line=4, outer=TRUE, cex=1.8)
		dev.off()
		}
	}

# again excluding Alces doesn't make any difference
tr="FV"
	for (dt in names(XmL[[tr]])) {
		quartz(width=11, height=4, file=paste("EvDicPsp", tr, dt, "_woA.jpg", sep="_"), type="jpg", dpi=150)
		layout(matrix(1:2,1,2))
		par(mar=c(5,5,1,1))
		ymx <- max(sapply(EvDicPspL[[tr]][[dt]], max))*100
		for (nd in names(nodes)[c(4,1)]) {
			E <- EvDicPspL[[tr]][[dt]][[nd]]*100
			boxplot(as.data.frame(E), cex.axis=1.2, cex.lab=1.4, ylim=c(0,ymx*1.05), ylab="Evolvability (%)", xlab="Eigenvalue (%)")
			if (nd=="Cervidae") {
					legend("topright", nd, cex=1.2, bty="n")
			} else {legend("topright", "Cervidae without Alces", cex=1.2, bty="n")} 
		}
		dev.off()}

# Excluding caprines from Bovidae doesn't make a difference either; median is a little higher for ILtes only
tr="FV"
	for (dt in names(XmL[[tr]])) {
		quartz(width=11, height=4, file=paste("EvDicPspL", tr, dt, "_woCp.jpg", sep="_"), type="jpg", dpi=150)
		layout(matrix(1:2,1,2))
		par(mar=c(5,5,1,1))
		ymx <- max(sapply(EvDicPspL[[tr]][[dt]], max))*100
		for (nd in names(nodes)[c(3,2)]) {
			E <- EvDicPspL[[tr]][[dt]][[nd]]*100
			boxplot(as.data.frame(E), cex.axis=1.2, cex.lab=1.4, ylim=c(0,ymx*1.05), ylab="Evolvability (%)", xlab="Eigenvalue (%)")
			if (nd=="Bovidae") {
					legend("topright", nd, cex=1.2, bty="n")
			} else {legend("topright", "Bovidae without Caprinae", cex=1.2, bty="n")} 
		}
		dev.off()}


######## FIGURE 5
### e(d.sp) in Table 1
###  d.sp ~ Pav; evolvabilities along divergence vector between each species and the clade root
# Pav is the pooled within-pop matrix of the clade (averaged across species; clades are given in 'nodes' above)

EvPavL <- list()

for (tr in names(LtM)) {
	tree <- LtM[[tr]]
	for (dt in names(XmL[[tr]])) {	
		PLdt <- PL[[tr]][[dt]] # using covariances only for each dataset
		X <- XmL[[tr]][[dt]][tree$tip.label,] # multivariate mean of each species in the original trait space (not PC scores), ordered as in the tree 
		for (nd in names(nodes)){
			ndtree <- extract.clade(tree, node=nodes[nd]) # specific clade only; one clade at a time
			if (nd=="Cervidae_woA") {ndtree <- drop.tip(ndtree, tip=c("Alc_al2","Alc_am4"))}
			if (nd=="Bovidae_woCp") {ndtree <- drop.tip(ndtree, tip=prop.part(ndtree)[[which(ndtree$node.label=="139FV")]])}
			txnd <- ndtree$tip.label # species within the clade who have morphological information
			Xnd <- X[txnd,]
			a0 <- aceRoot(ndtree, Xnd)
			Xnd <- t(t(Xnd)/a0) # mean-standardized within the clade using phy-weighted mean
			txndi <- txnd[txnd%in%names(PLdt)] # species that have both morphological and integration information
			# calculating Pav
			Pav.ss <- matrix(0, nrow(PLdt[[1]]), ncol(PLdt[[1]]))
			for (tx in txndi) {
				Ptx.ss <- PLdt[[tx]]*(nn[[tr]][[tx]]-1) # sum of squares
				Pav.ss <- Pav.ss + Ptx.ss
				}
			Pav <- Pav.ss/(sum(unlist(nn))-length(txndi)) # pooled ms
			# evolvabilities alond d.sp
			E <- matrix(nr=length(txndi),nc=7, dimnames=list(txndi, c("e.min", "e.av", "e.obs", "e.max", "c.av", "c.obs", "Mdiv")))
			ev <- eigen(Pav)$values
			e.min <- min(ev)
			e.av <- mean(ev)
			c.av <- 1/mean(1/ev)*(1+(2*var(1/ev)/mean(1/ev)^2)/(length(ev)+2))
			e.max <- max(ev)
			for (tx in txndi) {
				x <- Xnd[tx,]-(a0/a0)
				mdiv <- sqrt(sum(x^2)) # morphological divergence (Euclidean distance in morphospace) between tx and the clade's root
				x <- x/norm(as.matrix(x),"F") # standardized to unit length
				e.obs <- t(x)%*%Pav%*%x
				c.obs <- 1/(t(x)%*%solve(Pav)%*%x)
				E[tx,] <- c(e.min, e.av, e.obs, e.max, c.av, c.obs, mdiv)
				}
		EvPavL[[tr]][[dt]][[nd]] <- E
		print(c(tr, dt, nd))
		}
	}
}

save.image("ws-HHevPav.RData")

## Plots
for (tr in names(LtM)) {
	for (dt in names(XmL[[tr]])) {
			quartz(width=5, height=9, file=paste("EvHH",dt,tr,"PavFig.jpg",sep="_"), type="jpg", dpi=150)
			layout(matrix(1:3,3,1))
			par(mar=c(5,5,1,1), oma=c(0,0,3,0))
			for (nd in names(nodes)[-(1:2)]) {
				E <- EvPavL[[tr]][[dt]][[nd]]*100
				E <- E[o <- order(E[,"Mdiv"]),]
				E[,"Mdiv"] <- log(E[,"Mdiv"]/100)
				plot(E[,"Mdiv"], E[,"e.obs"], pch=21, xlab="Morphological divergence", ylab="Evolvability (%)", ylim=c(0,max(E[,"e.max"])*1.2), cex.lab=1.7, cex.axis=1.4, cex=1.5, lwd=1.5, bg="grey60")
				segments(-1, E[,"e.min"], max(E[,"Mdiv"])*1.1, E[,"e.min"], lwd=1)
				segments(-1, E[,"e.max"], max(E[,"Mdiv"])*1.1, E[,"e.max"], lwd=1)
				segments(-1, E[,"e.av"], max(E[,"Mdiv"])*1.1, E[,"e.av"], lty=2, lwd=1.5)
				points(E[,"Mdiv"], E[,"c.obs"], pch=22, cex=1.5, lwd=1.5, bg="grey90")
				segments(-1, E[,"c.av"], max(E[,"Mdiv"])*1.1, E[,"c.av"], lty=3)
				legend("topright", nd, cex=1.8, bty="n")
				legend("topleft", c("e", "c"), pch=c(21,22), cex=1.8, pt.bg=c("grey60", "grey90"), bty="n", ncol=2)
				#mtext(c("e.min", "e.ave", "e.max"), side=4, at=c(E[,"e.min"],E[,"e.av"],E[,"e.max"]), las=1)
				if (nd=="Bovidae") {
					i <- which(rownames(E)%in%rownames(EvPavL[[tr]][[dt]][["Caprinae"]]))
					points(E[i,"Mdiv"], E[i,"e.obs"], pch=19, cex=1.5, col="black")
					}
				}
			dev.off()
			}
		}

## Excluding Alces from Cervidae makes no difference

tr="FV"
dt="IL32"
			quartz(width=5, height=7, file=paste("EvHH",dt,tr,"Pa_woA.jpg",sep="_"), type="jpg", dpi=150)
			layout(matrix(1:2,2,1))
			par(mar=c(4,4.5,1,1), oma=c(0,0,3,0))
			for (nd in names(nodes)[c(4,1)]) {
				E <- EvPavL[[tr]][[dt]][[nd]]*100
				E <- E[o <- order(E[,"Mdiv"]),]
				E[,"Mdiv"] <- log(E[,"Mdiv"]/100)
				plot(E[,"Mdiv"], E[,"e.obs"], pch=21, xlab="Morphological divergence", ylab="Evolvability (%)", ylim=c(0,max(E[,"e.max"])*1.2), cex.lab=1.5, cex.axis=1.4, cex=1.5, lwd=1.5, bg="grey60")
				segments(-1, E[,"e.min"], max(E[,"Mdiv"])*1.1, E[,"e.min"], lwd=1)
				segments(-1, E[,"e.max"], max(E[,"Mdiv"])*1.1, E[,"e.max"], lwd=1)
				segments(-1, E[,"e.av"], max(E[,"Mdiv"])*1.1, E[,"e.av"], lty=2, lwd=1.5)
				points(E[,"Mdiv"], E[,"c.obs"], pch=22, cex=1.5, lwd=1.5, bg="grey90")
				segments(-1, E[,"c.av"], max(E[,"Mdiv"])*1.1, E[,"c.av"], lty=3)
				if (nd=="Cervidae") {
					legend("topright", nd, cex=1.2, bty="n")
					} else {legend("topright", "Cervidae without Alces", cex=1.2, bty="n")} 
				legend("topleft", c("e", "c"), pch=c(21,22), pt.bg=c("grey60", "grey90"), pt.cex=1.6, bty="n")
				}
			dev.off()

## Excluding caprines from Bovidae makes no difference 

tr="FV"
dt="IL32"
			quartz(width=5, height=7, file=paste("EvHH",dt,tr,"Pa_woCap.jpg",sep="_"), type="jpg", dpi=150)
			layout(matrix(1:2,2,1))
			par(mar=c(4,4.5,1,1), oma=c(0,0,3,0))
			for (nd in names(nodes)[c(3,2)]) {
				E <- EvPavL[[tr]][[dt]][[nd]]*100
				E <- E[o <- order(E[,"Mdiv"]),]
				E[,"Mdiv"] <- log(E[,"Mdiv"]/100)
				plot(E[,"Mdiv"], E[,"e.obs"], pch=21, xlab="Morphological divergence", ylab="Evolvability (%)", ylim=c(0,max(E[,"e.max"])*1.2), cex.lab=1.5, cex.axis=1.4, cex=1.5, lwd=1.5, bg="grey60")
				segments(-1, E[,"e.min"], max(E[,"Mdiv"])*1.1, E[,"e.min"], lwd=1)
				segments(-1, E[,"e.max"], max(E[,"Mdiv"])*1.1, E[,"e.max"], lwd=1)
				segments(-1, E[,"e.av"], max(E[,"Mdiv"])*1.1, E[,"e.av"], lty=2, lwd=1.5)
				points(E[,"Mdiv"], E[,"c.obs"], pch=22, cex=1.5, lwd=1.5, bg="grey90")
				segments(-1, E[,"c.av"], max(E[,"Mdiv"])*1.1, E[,"c.av"], lty=3)
				if (nd=="Bovidae") {
					legend("topright", nd, cex=1.2, bty="n")
					} else {legend("topright", "Bovidae without Caprinae", cex=1.2, bty="n")} 
				legend("topleft", c("e", "c"), pch=c(21,22), pt.bg=c("grey60", "grey90"), pt.cex=1.6, bty="n")
				}
			dev.off()




