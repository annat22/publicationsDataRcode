######### Protocol to reproduce all analyses and generate  
######### all figures and tables in Haber and Dworkin 2016

## read in data etc.
## see README file for information on these files
load("inf_spec.Rdata")
load("inf_gt.Rdata")
load("XallData.Rdata") # all genotypes superimposed together
source("functions-HaberDworkin2016.R")

# regression model correcting for replicate effect and allometry within each genotype

Xall.res <- Xall

for (gt in rownames(inf_gt)) {
	ind <- which(inf_spec$genotype==gt)
	Xgt <- Xall[ind,]
	repID <- droplevels(inf_spec$repID[ind]) # replicate factor
	lcs <- log(inf_spec$centsize[ind]) # log centroid size
	reg <- lm(Xgt ~ repID + lcs) # regression model
	Xm <- colMeans(predict(reg)) # # predicted mean configuration to add back to the residuals 
	#Xm <- coef(reg)[1,] # predicted mean configuration to add back to the residuals
	Xall.res[ind,] <- reg$residuals + matrix(Xm,nr=length(ind),nc=length(Xm),byrow=TRUE) # recentered on the genotype mean
	}


save(Xall.res, file="Xall.res.Rdata")
write.table(Xall.res, file="Xall.res.txt", sep="\t", row.names=FALSE, col.names=FALSE) # for LORY


##########################################################################################
#### comparing shape scores with and without the model
eig <- eigen(var(Xall.res))
C1 <- Xall.res%*%eig$vectors[,1:30] # shape scores with model

load("/Users/annat/Dropbox/AnnatHaber/Annat_WingIntegrationV1/dataCode/CsrL.Rdata")
C0 <- c()
for (gt in names(CsrL[["P"]])) {
	C0 <- rbind(C0, CsrL[["P"]][[gt]])} # shape scores previously computed without model

C0 <- C0[rownames(C1),]

quartz(width=6, height=6, file="PCcompared.jpg", type="jpg")
layout(matrix(1:4,2,2))
par(mar=c(4.5,4.5,0.5,0.5))
for (i in 1:4) {
	plot(C1[,i], C0[,i], xlab=paste("PC",i," with model", sep=""), ylab=paste("PC",i," without model", sep=""))}
dev.off()

rm(list=ls())

########################################################################################
########### load LORY output and combine all datasets with reduced dimensionality; Figure S1
require(vegan)
load("inf_spec.Rdata")
load("inf_gt.Rdata")
load("Xall.res.Rdata")
#source("functions-HaberDworkin2016.R")


p <- 30 # number of dimensions to keep

quartz(width=6, height=10, file="FigS1_screeplot.jpg", type="jpg", dpi=150)
layout(matrix(1:3,3,1))
par(mar=c(3,3,0.5,0.5))
CsrL <- XsrL <- list()
for (dt in c("P","Jtps","Jebs")) {
		if (dt=="P") {
			Xdt <- Xall.res
			} else if (dt=="Jtps") {
				Xdt <- as.matrix(read.table("JxcTPSsr.dat", sep="\t", header=TRUE))
				} else {
					Xdt <- as.matrix(read.table("JxcEBSsr.dat", sep="\t", header=TRUE))}
	rownames(Xdt) <- rownames(inf_spec)

	eig <- eigen(var(Xdt))
	Cdt <- Xdt%*%eig$vectors[,1:p]
	ev <- round(eig$values/sum(eig$values),3)[1:58]
	plot(ev, cex=1.2, pch=19,cex.axis=1.4, xlab="PC index", ylab="% eigenvalue")
	be <- bstick(length(ev))
	points(be, col="gray15", cex=2,pch="*", font=2,lwd=2)
	legend("topright", legend=paste(dt,"size-regressed"), cex=2)
	cumev <- round(100*cumsum(ev))
	ei <- c(which(cumev>94)[1], which(cumev>98)[1], which(cumev==100)[1])
	text(ei,c(0.165,0.125,0.075), labels=paste(cumev[ei], "%", sep=""), cex=1.8)
	arrows(ei, c(0.155,0.115,0.065), ei, c(0.07,0.0450,0.02), lwd=1.5, length=0.15, angle=25)

	for (gt in rownames(inf_gt)) {
		ind <- which(inf_spec$genotype==gt)
		XsrL[[dt]][[gt]] <- Xdt[ind,] # Procrustes and LORY original data, size-regressed
		CsrL[[dt]][[gt]] <- Cdt[ind,] # shape scores, Procrustes and LORY, size-regressed
		}
	}
dev.off()
	
inf_specL <- split(inf_spec, inf_spec$genotype)

JKco <- matrix(as.numeric(read.table("JxcCo.dat", sep="\t", header=TRUE)[1,]), nc=2, byrow=TRUE)
colnames(JKco) <- c("x", "y")
save(JKco, file="JKco.RData") # Coordinates of Jacobians

save(CsrL, file="CsrL.Rdata") # shape scores, size-regressed
save(XsrL, file="XsrL.Rdata") # original data, size-regressed
save(inf_specL, file="inf_specL.Rdata")

#save.image("ws-prelim.RData")

rm(list=ls())

###################################################################################################
########## Matrix properties: rSDE, eccentricity (inverse of Nd), total variance; with BCa CI's

load("inf_spec.Rdata")
load("inf_gt.Rdata")
source("functions-HaberDworkin2016.R")
load("CsrL.Rdata")
CL <- CsrL; rm(CsrL) # morphospace scores

funsL <- list("s2"=function(V){sum(diag(V))}, "rSDE"=rSDE, "eccentNd"=eccentNd) # list of functions to calculate matrix properties

## theta function to pass to bootNP.BCa
theta.iis2gv <- function(x, Xdata, statfun) {
	Xth <- Xdata[x,] # data resampled (or original if x=1:nrow(Xdata)); see bootNP.BCa
	V <- shrink(var(Xth), tol=10^-8)
	statfun(V)} 


n.boot=999 # number of bootstrap iterations

IIS2L <- list()
for (dt in c("P","Jtps","Jebs")) {
	for (df in names(funsL)) {
		IIS2 <- c()
		for (gt in rownames(inf_gt)) {
			print(c(dt,gt,df))
			Xdata <- CL[[dt]][[gt]]
			bres <- bootNP.BCa(1:nrow(Xdata), nboot=n.boot, theta=theta.iis2gv, Xdata=Xdata, statfun=funsL[[df]])
			IIS2 <- rbind(IIS2, c("obs"=bres$thetahat, bres$confpoints))
			write.table(t(c(date(),dt,gt,df)), file="progressIIS2.tab", quote=FALSE, sep="\t", append=TRUE, col.names=FALSE, row.names=FALSE)
			#print(c(date(),dt,gt,df))
			}
		rownames(IIS2) <- rownames(inf_gt)
		IIS2L[[dt]][[df]] <- IIS2
		}
	}
	
save(IIS2L, file="IIS2L.RData")
save.image("ws-IIS2.Rdata")

rm(list=ls())

#################################################################################
########## All pairwise shape differences; with BCa CI's

load("inf_spec.Rdata")
load("inf_gt.Rdata")
source("functions-HaberDworkin2016.R")
load("CsrL.Rdata")
CL <- CsrL; rm(CsrL) # morphospace scores

n.boot=999
p <- 30

IJ <- t(combn(rownames(inf_gt),2)) # all pairwise combinations of genoypes (lower triangle)

theta.shapedist <- function(x, Xdata, gf) {
	Xth <- Xdata[x,] # data resampled (or original if x=1:nrow(Xdata); see bootNP.BCa)
	gf <- gf[x] # grouping factor resampled (otherwise gets permuted)
	gfl <- levels(gf)
	#breaking down the concatinated dataset Xdata into the two samples to be compared (two genotypes)
	X1 <- Xth[gf==gfl[1],]
	X2 <- Xth[gf==gfl[2],]
	#calculating mean shapes
	m1 <- colMeans(X1)
	m2 <- colMeans(X2)
	sqrt(sum((m1-m2)^2)) # euclidean distance on tangent space
	}

shapeDL <- resL.shapeD <- list()

for (dt in c("P","Jtps", "Jebs")) {
	D <- c()
	for (i in 1:nrow(IJ)) {
		gt1 <- IJ[i,1] # genotype 1
		gt2 <- IJ[i,2] # genotype 2
		Xgt1 <- CL[[dt]][[gt1]]
		Xgt2 <- CL[[dt]][[gt2]]
		Xdata <- rbind(Xgt1, Xgt2) # data to pass to theta in bootNP.BCa
		gff <- as.factor(c(rep(gt1, nrow(Xgt1)), rep(gt2, nrow(Xgt2)))) #  grouping factor to pass to theta in bootNP.BCa
		bres <- bootNP.BCa(x=1:nrow(Xdata), nboot=n.boot, theta.fun=theta.shapedist, Xdata=Xdata, gf=gff)
		d.obs <- bres$thetahat # observed value
		D <- rbind(D, c("d.obs"=d.obs, bres$confpoints, quantile(bres$thetastar, probs=c(0.025,0.975))))	
		#print(c(dt, i, gt1, gt2))
		write.table(t(c(date(),dt, i, gt1, gt2)), file="progressShapeD.tab", quote=FALSE, sep="\t", append=TRUE, col.names=FALSE, row.names=FALSE)		
		resL.shapeD[[dt]][[paste(gt1,gt2,sep="-")]] <- bres
		}
	rownames(D) <-  apply(IJ, 1, paste, collapse="-")
	shapeDL[[dt]] <- D
	rm(bres)
	}


save(shapeDL, file="shapeDL.RData")
save.image("ws-shapeD_sizeReg.RData")

rm(list=ls())

#################################################################################
######## All pairwise differences in matrix orientation; with jackknife CI's

load("inf_spec.Rdata")
load("inf_gt.Rdata")
source("functions-HaberDworkin2016.R")
load("CsrL.Rdata")
CL <- CsrL; rm(CsrL) # morphospace scores

distfunL <-  list("randsk"=randsk, "comsubsp"= comsubsp, "MBreleig"= MBreleig) # three different metrics

theta.intdist <- function(x, Xdata, gf, distfun) {
	Xth <- Xdata[x,] # data resampled (or original if x=1:nrow(Xdata)); see bootNP.BCa
	gf <- gf[x]
	gfl <- levels(gf)
	#breaking down the concatinated dataset Xdata into the two samples to be compared (e.g, two genotypes)
	X1 <- Xth[gf==gfl[1],]
	X2 <- Xth[gf==gfl[2],]
	distfun(shrink(var(X1), tol=10^-8), shrink(var(X2), tol=10^-8))}

IJ <- t(combn(rownames(inf_gt),2)) # all pairwise combinations of genoypes (lower triangle only so no 

intDL <- resL.intD <- list()

p=ncol(CL[[1]][[1]])

for (dt in c("P","Jtps","Jebs")) {
	D <- array(dim=c(nrow(IJ), 5, length(distfunL)), dimnames=list(apply(IJ, 1, paste, collapse="-"), c("d.obs", "jack.orCI.025", "jack.orCI.975", "jack.bcCI.025", "jack.bcCI.975"), names(distfunL)))
	for (i in 1:nrow(IJ)) {
		gt1 <- IJ[i,1]; gt2 <- IJ[i,2]
		Xgt1 <- CL[[dt]][[gt1]][,1:p]
		Xgt1 <- t(t(Xgt1)-colMeans(Xgt1)) # centering on the mean for b0
		Xgt2 <- CL[[dt]][[gt2]][,1:p]
		Xgt2 <- t(t(Xgt2)-colMeans(Xgt2)) # centering on the mean for b0
		Xdata <- rbind(Xgt1, Xgt2) # concatenated data to pass to theta.fun in jackknife
		gff <- as.factor(c(rep(gt1, nrow(Xgt1)), rep(gt2, nrow(Xgt2)))) #  grouping factor to pass to theta.fun in jackknife
		for (df in names(distfunL)) {
			#print(c(date(),dt, i, gt1, gt2, df))
			jres <- jackknife(1:nrow(Xdata), theta.fun=theta.intdist, gf=gff, Xdata=Xdata, distfun=distfunL[[df]])
			d.obs <- jres$thetahat
			D[i,,df] <- c(d.obs, quantile(jres$jackdist.or, probs=c(0.025,0.975)), jres$confpoints)
			write.table(t(c(date(),dt, i, gt1, gt2, df)), file="progressIntD.tab", quote=FALSE, sep="\t", append=TRUE, col.names=FALSE, row.names=FALSE)		
			resL.intD[[dt]][[paste(gt1,gt2,sep="-")]][[df]] <- jres
			}
		}
	intDL[[dt]] <- D
	rm(jres)
	}

save.image("ws-intD.RData")
save(intDL, file="intDL.RData")

rm(list=ls())

#################################################################################
######### PCoA spaces compared using Procrustes superimposition
require("vegan")

load("inf_gt.Rdata")
load("intDL.Rdata")
load("shapeDL.RData")

Ngt <- nrow(inf_gt)

# pathway grouping factor
pthw.gf <- inf_gt$pathway[inf_gt$allele!="WT"]
pthw.gf[pthw.gf=="Hh"] <- "TGF-b"
pthw.gf <- droplevels(pthw.gf)

procMDSresL <- pathwayMANOVAresL <- list()

k=20

for (dt in names(shapeDL)) {
	
	## Shape space
	d <- shapeDL[[dt]][,1]
	D <- matrix(0,Ngt, Ngt)
	D[lower.tri(D)] <- d; D <- D + t(D)
	rownames(D) <- colnames(D) <- rownames(inf_gt)
	Csh <- cmdscale(D, k=k)
	ptwMsh <- manova(Csh[inf_gt$allele!="WT",] ~ pthw.gf)
	pathwayMANOVAresL[[dt]][["shapeD"]] <- ptwMsh
	
	# Covariance space (orientation)
	for (dfun in dimnames(intDL[[dt]])[[3]]){
		d <- intDL[[dt]][,1,dfun]
		D <- matrix(0,Ngt, Ngt)
		D[lower.tri(D)] <- d; D <- D + t(D)
		rownames(D) <- rownames(inf_gt)
		Cint <- cmdscale(D, k=k)
		ptwMint <- manova(Cint[inf_gt$allele!="WT",] ~ pthw.gf)
		pathwayMANOVAresL[[dt]][[dfun]] <- ptwMint
		
	### Superimposing shape space and orientation spaces
	procSO <- protest(Csh, Cint, scale=TRUE, symmetric=TRUE)
	procMDSresL[[dt]][[dfun]] <- procSO
	
	### Disparity
  	disparity <- function(X) {sum(diag(var(X)))}
 	# disparity <- function(X) {mean(dist(X, method = "euclidean")^2)}
	sh.jc <- procSO$X # shape scores in the joint space
	disp.sh <- disparity(sh.jc)
	int.jc <- procSO$Yrot # covariance scores in the joint space
		}
	}

save.image("ws-procMDS.RData")
save(procMDSresL, file="procMDSresL.RData")

rm(list=ls())

############################################################################################
######## FIGURES and TABLES

## Figure 1: illustrating matrix properties

require(ellipse)

quartz(width=7, height=7, file="Figure1ellipses.jpg", type="jpg", dpi=150)
par(mar=c(0,0,0,0))

e1 <- ellipse(0.7) #ref
e2 <- ellipse(0.92, centre=c(-3.3,3.3)) # higher eccentricity
e3 <- ellipse(0.7, centre=c(-5.5,-5.5), scale=c(0.65,0.65)) # smaller size
e4 <- ellipse(-0.7, centre=c(6.7,0)) # perpendicular orientation
eL <- list(e1,e2,e3,e4[c(75:100,1:74),])


plot(c(-9,9), c(-9,9), asp=1, pch=NA, xaxt="n", yaxt="n", frame=FALSE, ann=FALSE)
for (ei in 1:length(eL)) {
	segments(eL[[ei]][-100,1], eL[[ei]][-100,2], eL[[ei]][-1,1], eL[[ei]][-1,2], lwd=1.8)
	segments(eL[[ei]][100,1], eL[[ei]][100,2], eL[[ei]][1,1], eL[[ei]][1,2], lwd=1.8)
	arrows(eL[[ei]][51,1], eL[[ei]][51,2], eL[[ei]][2,1], eL[[ei]][2,2], length=0.15, angle=20, col="grey40")	
	arrows(eL[[ei]][26,1], eL[[ei]][26,2], eL[[ei]][76,1], eL[[ei]][76,2], length=0.15, angle=20, col="grey40")
	#arrows(arL[[ei]],lwd=4, col="grey30")
}

text(x=c(-5, -7, 8.5), y=c(4.5, -4.5, 1.5), labels=c("C", "A", "B"), font=2, cex=1.2)

dev.off()

########## Figures 4-7 
load("intDL.Rdata")
load("IIS2L.RData")
load("shapeDL.RData")
load("inf_gt.Rdata")
load("procMDSresL.Rdata")

# specifying comparison types between genotypes
IJ <- t(combn(rownames(inf_gt),2)) # all pairwise combinations of genotypes (lower triangle only so no redundancy)
nr.pw <- rownames(shapeDL[[1]])[which(IJ=="Sam", arr.ind=TRUE)[,1]] # comparisons between sam and everything else
nr <- unlist(strsplit(nr.pw,"-")); nr <- nr[nr!="Sam"] # names of all genotypes except Sam ordered as in nr.pw 
nr.pw.mt <- nr.pw[inf_gt[nr,"WT"]=="Sam"] # comparisons between sam and all mutants as in nr.pw
nr.mt <- nr[inf_gt[nr,"WT"]=="Sam"] # names of all mutants ordered as in nr.pw 

# pch by wild type 
pchs <- c()
pchs[inf_gt$WT=="NC"] <- 13
pchs[inf_gt$WT=="MA"] <- 21
pchs[inf_gt$WT=="Sam"] <- 19
pchs[rownames(inf_gt)=="Sam"] <- 8
names(pchs) <- rownames(inf_gt)

# pathway grouping factor
pthw.gf <- inf_gt$pathway[inf_gt$allele!="WT"]
pthw.gf[pthw.gf=="Hh"] <- "TGF-b"
pthw.gf <- droplevels(pthw.gf)

##############
#### correlation between TPS and EBS for shapeD and intD; Figure S2
quartz(width=6,height=8, file="FigS2_EBSvsTPS_shapeDintD_randsk.jpg", type="jpg", dpi=150)
layout(matrix(1:2, 2, 1))
par(mar=c(4.5,4.5,0.8,0.8))

plot(shapeDL[["Jtps"]][,1], shapeDL[["Jebs"]][,1], xlab="TPS", ylab="EBS", pch=NA, cex.lab=1.4, cex.axis=1.2)
segments(shapeDL[["Jtps"]][,2], shapeDL[["Jebs"]][,1], shapeDL[["Jtps"]][,3], shapeDL[["Jebs"]][,1], col="grey40", lwd=1.6)
segments(shapeDL[["Jtps"]][,1], shapeDL[["Jebs"]][,2], shapeDL[["Jtps"]][,1], shapeDL[["Jebs"]][,3], col="grey40", lwd=1.6)
points(shapeDL[["Jtps"]][,1], shapeDL[["Jebs"]][,1], pch=21, bg="grey80", cex=1.2)
segments(0,0,2,2)
legend("topleft", "All pairwise shape distances", bty="n", cex=1.2, adj=0.05)

plot(intDL[["Jtps"]][,1,"randsk"], intDL[["Jebs"]][,1,"randsk"], xlab="TPS", ylab="EBS", pch=NA, cex.lab=1.4, cex.axis=1.2)
segments(intDL[["Jtps"]][,2,"randsk"], intDL[["Jebs"]][,1,"randsk"], intDL[["Jtps"]][,3,"randsk"], intDL[["Jebs"]][,1,"randsk"], col="grey40", lwd=1.6)
segments(intDL[["Jtps"]][,1,"randsk"], intDL[["Jebs"]][,2,"randsk"], intDL[["Jtps"]][,1,"randsk"], intDL[["Jebs"]][,3,"randsk"], col="grey40", lwd=1.6)
points(intDL[["Jtps"]][,1,"randsk"], intDL[["Jebs"]][,1,"randsk"], pch=21, bg="grey80", cex=1.2)
legend("topleft", "All pairwise covariance distances", bty="n", cex=1.2, adj=0.05)
segments(0,0,2,2)

dev.off()

#### correlation between Jacobian-based (EBS) and Proc-based shapeD and intD; Figure S3

quartz(width=6,height=8, file="FigS3_PvsJebs_shapeDintD_randsk.jpg", type="jpg", dpi=150)
layout(matrix(1:2, 2, 1))
par(mar=c(4.5,4.5,0.8,0.8))

plot(shapeDL[["P"]][,1], shapeDL[["Jebs"]][,1], xlab="Procrustes", ylab="Jacobians (EBS)", pch=NA, cex.lab=1.4, cex.axis=1.2)
segments(shapeDL[["P"]][,2], shapeDL[["Jebs"]][,1], shapeDL[["P"]][,3], shapeDL[["Jebs"]][,1], col="grey40", lwd=1.6)
segments(shapeDL[["P"]][,1], shapeDL[["Jebs"]][,2], shapeDL[["P"]][,1], shapeDL[["Jebs"]][,3], col="grey40", lwd=1.6)
points(shapeDL[["P"]][,1], shapeDL[["Jebs"]][,1], pch=21, bg="grey80", cex=1.2)
#abline(lm(shapeDL[["Jebs"]][,1]~shapeDL[["P"]][,1]))
#segments(0,0,2,2)
legend("topleft", "All pairwise shape distances", bty="n", cex=1.2, adj=0.05)

plot(intDL[["P"]][,1,"randsk"], intDL[["Jebs"]][,1,"randsk"], xlab="Procrustes", ylab="Jacobians (EBS)", pch=NA, cex.lab=1.4, cex.axis=1.2)
segments(intDL[["P"]][,2,"randsk"], intDL[["Jebs"]][,1,"randsk"], intDL[["P"]][,3,"randsk"], intDL[["Jebs"]][,1,"randsk"], col="grey40", lwd=1.6)
segments(intDL[["P"]][,1,"randsk"], intDL[["Jebs"]][,2,"randsk"], intDL[["P"]][,1,"randsk"], intDL[["Jebs"]][,3,"randsk"], col="grey40", lwd=1.6)
points(intDL[["P"]][,1,"randsk"], intDL[["Jebs"]][,1,"randsk"], pch=21, bg="grey80", cex=1.2)
legend("topleft", "All pairwise covariance distances", bty="n", cex=1.2, adj=0.05)
segments(0,0,2,2)

dev.off()


#### shpaeD vs total variance; all relative to Sam only; Figure 4 and S8

quartz(width=5,height=11,file="Fig4_shapeD~variance.jpg", type="jpg", dpi=150)
layout(matrix(1:3,3,1))
par(mar=c(5,5,4,1))
for (dt in names(shapeDL)) {
	y <- shapeDL[[dt]][nr.pw,]
	x <- IIS2L[[dt]][["s2"]][nr,]
	# e <- s2/p
	plot(x[,1], y[,1], pch=NA, ylab="Shape distance from Sam", xlab="Total variance", xlim=c(min(x,na.rm=TRUE), max(x,na.rm=TRUE)), ylim=c(0, max(y,na.rm=TRUE)), cex.lab=1.6, cex.axis=1.4)
	segments(x[,2], y[,1], x[,3], y[,1], col="grey30")
	segments(x[,1], y[,2], x[,1], y[,3], col="grey30")
	points(x[,1], y[,1], pch=pchs[nr], bg="white", cex=2)
	segments(IIS2L[[dt]][["s2"]]["Sam",1], 0, IIS2L[[dt]][["s2"]]["Sam",1], 4, lty=2, lwd=2)
	segments(0, 0, 1, 0, lty=2, lwd=2)
	legend("bottomright", legend=c("MA","NC","Mutants"), pch=c(21,13,19), bg="white", cex=1.4)
	mtext(dt, line=2, adj=-0.12, cex=1.4, font=2)
	#text(x=0.18, y=1.17,labels="3045")
	text(x=x[c("3045","10413"),1], y=y[c("Sam-3045","Sam-10413"),1],labels=c("3045","10413"), adj=c(-0.3,-0.5))
	}
dev.off()
### intD vs shapeD; all relative to Sam only; Figure 5
# one figure for each of the similarity measure (randsk, comsubsp, MB)
# three datasets in each figure

for (df in c("randsk", "comsubsp", "MBreleig")) {
	quartz(width=5, height=11, file=paste("Fig5_intD~shapeD_", df, ".jpg", sep=""), type="jpg", dpi=150)
	layout(matrix(1:3,3,1))
	par(mar=c(6.5,6.5,4,1))
	for (dt in names(shapeDL)) {
		y <- shapeDL[[dt]][nr.pw,1:3]
		x <- intDL[[dt]][nr.pw,1:3,df]
		plot(x[,1], y[,1], pch=NA, ylab="Shape distance from Sam", xlab="Covariance distance from Sam", ylim=c(0, max(y,na.rm=TRUE)), xlim=c(min(x,na.rm=TRUE), max(x,na.rm=TRUE)), cex.lab=1.5, cex.axis=1.4)
		#segments(0, 0, 1, 0, lty=2, lwd=2)
		#segments(0, 0, 0, 1, lty=2, lwd=2, col="grey30")
		segments(x[,2], y[,1], x[,3], y[,1], col="grey60")
		segments(x[,1], y[,2], x[,1], y[,3], col="grey60")
		points(x[,1], y[,1], pch=pchs[nr], bg="white", cex=1.8)
		legend("bottomright", legend=c("Mutants", "MA", "NC"), pch=c(19,21,13), pt.bg="white", cex=1.4)
		mtext(dt, line=2, adj=-0.12, cex=1.2, font=2)
		#legend("topleft", legend=c("P"="Procrustes", "Jtps"="TPS", "Jebs"="EBS")[dt], cex=1.5, adj=0.3, bty="n")
		#legend("topright", legend=c("Sam"), lty=2, bty="n", cex=1.6)
		text(x=x[c("Sam-3045","Sam-10413"),1], y=y[c("Sam-3045","Sam-10413"),1],labels=c("3045","10413"), adj=c(-0.3,-0.5))
		}
	dev.off()
	}


### PCoA spaces superimposed; Figure 6 and S9-S12

for (dt in names(shapeDL)) {
		
	for (dfun in dimnames(intDL[[1]])[[3]]){	## Covariance space (orientation)
		sh <- procMDSresL[[dt]][[dfun]]$X[,1:2]
		int <-procMDSresL[[dt]][[dfun]]$Yrot[,1:2]
		ev <- procMDSresL[[dt]][[dfun]]$svd$d; ev <- round(ev*100/sum(ev))
		x=1;y=2
		quartz(width=5, height=8, type="jpg", file=paste("Fig6_mdsShapeInt12_", dt, "_", dfun,".jpg", sep=""), dpi=300) 
		layout(matrix(1:2, 2, 1))
		par(mar=c(4,4,0.5,0.5))
		plot(rbind(sh,int), pch=NA, cex.lab=1, cex.axis=1, xlab=paste0("PCo",x," (",ev[x],"%)"), ylab=paste0("PCo",y," (",ev[2],"%)"))
		segments(0,-0.5,0,0.5, lty=2, col="grey40")
		segments(-0.5,0,0.5,0, lty=2, col="grey40")
		segments(sh[,1], sh[,2], int[,1], int[,2], col="grey60")
		points(sh, pch=ifelse(inf_gt$WT=="Sam", 19, 15))
		points(int, pch=ifelse(inf_gt$WT=="Sam", 21, 22), bg="grey90")
		legend("bottomright", legend=c("shape nWT", "shape mt", "covariance nWT", "covariance mt"), pch=c(15,19,22,21), pt.bg="grey90", cex=1)

		sh.mt <- sh[inf_gt$allele!="WT",]; int.mt <- int[inf_gt$allele!="WT",]
		plot(rbind(sh.mt,int.mt), pch=NA, cex.lab=1, cex.axis=1, xlab="PCo1", ylab="PCo2", xlim=c(min(sh.mt[,1],int.mt[,1]),0))
		segments(0,-0.5,0,0.5, lty=2, col="grey40")
		segments(-0.5,0,0.5,0, lty=2, col="grey40")
		segments(sh.mt[,1], sh.mt[,2], int.mt[,1], int.mt[,2], col="grey60")
		points(sh.mt, pch=19, cex=1.5)
		points(int.mt, pch=21, bg="grey90", cex=1.5)
		points(rbind(sh["Sam",], int["Sam",]), pch=c(19,21), bg="grey90", cex=1.5)
		legend("bottomleft", legend=c("shape mt", "covariance mt"), pch=c(19,21), pt.bg="grey90", cex=1)
		text(sh[c("Sam", "3045", "10413"),], labels=c("Sam", "3045", "10413"), cex=0.8, pos=c(3,2,2), offset=0.5)
		text(sh.mt[pthw.gf=="TGF-b",], labels="t", cex=0.8, offset=0, vfont=c("serif", "bold italic"), col="white")
		text(int.mt[pthw.gf=="TGF-b",], labels="t", cex=0.8, offset=0, vfont=c("serif", "bold italic"))
		dev.off()
		rm(sh,int)
		}
	}



#### Association between eccentricity (inverse of effective number of dimensions) and rSDE; figure S4

	quartz(width=5, height=10, file="FigS4_eccentNd~rSDE.jpg", type="jpg", dpi=150)
	layout(matrix(1:3,3,1))
	par(mar=c(5.2,5.2,1,1))
	for (dt in names(shapeDL)) {
		y <- IIS2L[[dt]][["eccentNd"]][nr,]
		x <- IIS2L[[dt]][["rSDE"]][nr,]
		plot(x[,1], y[,1], pch=NA, ylab="Eccentricity", xlab="rSDE", xlim=c(min(x,na.rm=TRUE), max(x,na.rm=TRUE)), ylim=c(min(y,na.rm=TRUE), max(y,na.rm=TRUE)), cex.lab=2, cex.axis=1.6)
		segments(0, IIS2L[[dt]][["eccentNd"]]["Sam",1], 5, IIS2L[[dt]][["eccentNd"]]["Sam",1], col="grey30", lty=2, lwd=2)
		segments(IIS2L[[dt]][["rSDE"]]["Sam",1], 0, IIS2L[[dt]][["rSDE"]]["Sam",1], 5, col="grey30", lt=2, lwd=2)
		segments(x[,2], y[,1], x[,3], y[,1], col="grey60")
		segments(x[,1], y[,2], x[,1], y[,3], col="grey60")
		points(x[,1], y[,1], pch=pchs[nr], cex=2, bg="white")
		legend("topleft", legend=c("Mutants", "MA", "NC"), pch=c(19,1,13), cex=1.4, box.col="white")
		legend("topright", legend=dt, cex=1.4, text.font=2, bty="n")
		}
dev.off()

#ylab=expression(paste("eccentricity; SD(", lambda, ")"))
#ylab=expression(paste("eccentricity; ", lambda[1]/sum(lambda[i],i=2,n)))


#### Table 2: all matrix measures
# highlighting those that differ siginifacntly from Sam

M <- c()
for (df in cn <- c("s2","rSDE","eccentNd")) {
	I <- round(IIS2L[["Jebs"]][[df]], 2)
	m <- paste(I[,1], " (", I[,2], ",", I[,3], ")", sep="")
	is <- which(I[,3] < I["Sam",1] | I[,2] > I["Sam",1])
	m[is] <- paste("*",  m[is], sep="")
	M <- cbind(M, m)
	}

rownames(M) <- rownames(IIS2L[[1]][[1]])
colnames(M) <- cn
O <- round(intDL[["Jebs"]][nr.pw,1:3,"randsk"],2)
S <- round(shapeDL[["Jebs"]][nr.pw,1:3],2)
O <- rbind("0", O)
S <- rbind("0", S)
M <- cbind(M, "cov.dis"=paste(O[,1], " (", O[,2], ",", O[,3], ")", sep=""), "shape.dis"=paste(S[,1], " (", S[,2], ",", S[,3], ")", sep=""))
M[-1,4:5] <- paste("*", M[-1,4:5], sep="")

write.table(M, file="Table2_matrixMeasures.tab", sep="\t", quote=FALSE)

### eccentricity vs variance/shapeD/intD; figure 7 and S14-S15
# one figure for each dataset
# intD is only random skewers

for (dt in names(IIS2L)) {
	quartz(width=5, height=10, file=paste("Fig7_eccentNd~s2shapeDintD_",dt,".jpg",sep=""), type="jpg", dpi=150)
	layout(matrix(1:3,3,1))
	par(mar=c(5.5,7,2.5,1))
	for (ci in c("s2", "shapeD", "intD")) {
	y <- IIS2L[[dt]][["eccentNd"]][nr,]
	ylb <- "Eccentricity"
		if (ci=="s2") {
			x <- IIS2L[[dt]][["s2"]][nr,]
			xlb <- "Total variance"
			} else if (ci=="shapeD"){
				x <- shapeDL[[dt]][nr.pw,]
				xlb <- "Shape distance from Sam"
				} else {
					x <- intDL[[dt]][nr.pw,,"randsk"]
					xlb <- "Covariance distance from Sam"
					}
		plot(x[,1], y[,1], pch=NA, xlab=xlb, ylab=ylb, xlim=c(min(x[,1:3],na.rm=TRUE), max(x[,1:3],na.rm=TRUE)), ylim=c(min(y,na.rm=TRUE), max(y,na.rm=TRUE)), cex.lab=1.6, cex.axis=1.4)
		segments(0, IIS2L[[dt]][["eccentNd"]]["Sam",1], 4, IIS2L[[dt]][["eccentNd"]]["Sam",1], col="grey30", lty=2, lwd=2)
		if (ci=="s2") {
			segments(IIS2L[[dt]][["s2"]]["Sam",1], 0, IIS2L[[dt]][["s2"]]["Sam",1], 1, col="grey30", lty=2, lwd=2)
			text(x=x["2513",1], y=y["2513",1],labels=c("2513"), adj=c(1.3,-0.5))
			} else {	
				text(x=x[c("Sam-2513","Sam-12772"),1], y=y[c("2513","12772"),1],labels=c("2513","12772"), adj=c(1.3,1))} 
		segments(x[,2], y[,1], x[,3], y[,1], col="grey60")
		segments(x[,1], y[,2], x[,1], y[,3], col="grey60")
		points(x[,1], y[,1], pch=pchs[nr], bg="white", cex=2)
		legend("topleft", legend=c("Mutants", "MA", "NC"), pch=c(19,21,13), pt.bg="white", bg="white", cex=1.4)
		mtext(c(s2="A",shapeD="B", intD="C")[ci], side=2, line=5.5, font=2, las=1, at=max(y)*1.10)
		}
	dev.off()
	}

### figure S13: rSDE (distance from sam) vs variance/shapeD/intD
# only Jebs
# intD is only random skewers

dt="Jebs"
	quartz(width=5, height=10, file="FigS12_rSDE~s2shapeDintD_Jebs.jpg", type="jpg", dpi=150)
	layout(matrix(1:3,3,1))
	par(mar=c(5.5,7,2.5,1))
	for (ci in c("s2", "shapeD", "intD")) {
		if (ci=="s2") {
			x <- IIS2L[[dt]][["s2"]][nr,]
			xlb <- "Total variance"
			y <- IIS2L[[dt]][["rSDE"]][nr,]
			ylb <- "rSDE"
			} else if (ci=="shapeD"){
				x <- shapeDL[[dt]][nr.pw,]
				xlb <- "Shape distance from Sam"
				y <- IIS2L[[dt]][["rSDE"]][nr,] - IIS2L[[dt]][["rSDE"]]["Sam",] # distance from sam
				ylb <- "rSDE distance from Sam"
				} else {
					x <- intDL[[dt]][nr.pw,,"randsk"]
					xlb <- "Covariance distance from Sam"
					y <- IIS2L[[dt]][["rSDE"]][nr,] - IIS2L[[dt]][["rSDE"]]["Sam",] # distance from sam
					ylb <- "rSDE distance from Sam"
					}
		plot(x[,1], y[,1], pch=NA, xlab=xlb, ylab=ylb, xlim=c(min(x[,1:3],na.rm=TRUE), max(x[,1:3],na.rm=TRUE)), ylim=c(min(y,na.rm=TRUE), max(y,na.rm=TRUE)), cex.lab=1.5, cex.axis=1.4)
		if (ci=="s2") {
			segments(0, IIS2L[[dt]][["rSDE"]]["Sam",1], 4, IIS2L[[dt]][["rSDE"]]["Sam",1], col="grey30", lty=2, lwd=2)
			segments(IIS2L[[dt]][["s2"]]["Sam",1], 0, IIS2L[[dt]][["s2"]]["Sam",1], 1, col="grey30", lty=2, lwd=2)
			} else {segments(0, 0, max(x)*2, 0, lty=2, lwd=2, col="grey30")}
		segments(x[,2], y[,1], x[,3], y[,1], col="grey60")
		segments(x[,1], y[,2], x[,1], y[,3], col="grey60")
		points(x[,1], y[,1], pch=pchs[nr], bg="white", cex=2)
		legend("topleft", legend=c("Mutants", "MA", "NC"), pch=c(19,21,13), pt.bg="white", bg="white", cex=1.4)
		mtext(c(s2="A",shapeD="B", intD="C")[ci], side=2, line=5.5, font=2, las=1, at=max(y)*1.10)
		}
	dev.off()

##### regression of shape distance over covariance distance
load("intDL.Rdata")
load("shapeDL.RData")
load("inf_gt.Rdata")

x <- intDL[["Jebs"]][1:24,"d.obs","randsk"]
y <- shapeDL[["Jebs"]][1:24,"d.obs"]
reg <- lm( y ~ x + I(x^2))

# > summary(reg)
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.48427 -0.09175 -0.01381  0.21137  0.45430 
# Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2.8364     0.9226  -3.074  0.00575 **
# x            9.0695     2.8216   3.214  0.00416 ** .  
# I(x^2)      -5.1904     2.1052  -2.466  0.02238 *     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.2532 on 21 degrees of freedom
# Multiple R-squared:  0.6352,	Adjusted R-squared:  0.6005
# F-statistic: 18.29 on 2 and 21 DF,  p-value: 2.519e-05
# > confint(reg)
#                2.5 %     97.5 %
# (Intercept) -4.755009 -0.9177331
# x            3.201710 14.9371946
# I(x^2)      -9.568443 -0.8124214




