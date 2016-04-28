																	# most of the major functions are in this file, _Euc version produces the distance matrix via Wills 98 instead of original foote 93
library(calibrate) 																					# needed for a function for producing alternative strat data
library(zoo)																						# needed for funtion na.approx that is used to estimate values between points in rarefaction curves
library(earth)																					# needed for function staxlab which can be used to rotate labels and stack them)	
library(moments)																					# needed for the skewness function that is calculated for each time bin (skewness of the bootstrap replicate distribution)
library(nlme) 
library(paleoTS) 



### Disparity analysis using SoV and SoR
DtT <- function(MS,PA,boo) {
	skewness.b.sple <- c()			# the skewness of the distribution of the bootstrap resample at everybin
	rar <- min(apply(PA,2,sum))			# the diversity value the Sum of Ranges is rarefied to
	if (rar==1) rar <- 2
	Variance <- c() ; Range <- c()
	
	MeanBootCentroid<- matrix(NA,nrow=ncol(PA),ncol=ncol(MS))		###### ADDED BY MARTIN HUGHES 13/04/2015
	SDBootCentroid<- c()											###### ADDED BY MARTIN HUGHES 13/04/2015
	
	for (tps in 1:ncol(PA)) {			# a loop for computing disparity for each time interval
		occ <- which(PA[,tps]==1)
		if (length(occ)==1) {occ <- c(occ,occ) ; rar <- 2 }
		vari <- c() ; ran <- c()
		centri<- matrix(NA,ncol=ncol(MS),nrow=boo)					###### ADDED BY MARTIN HUGHES 13/04/2015
		for (i in 1:boo) {
			b.sple <- sample(occ,length(occ),replace=T)
			r.sple <- sample(occ,rar,replace=F)
			vari <- c(vari,SoV(MS[b.sple,]))		# records disparity values of boostrap replicates
			ran <- c(ran,SoR(MS[r.sple,]))		# same but for range-based estimates
			centri[i,]<- apply(MS[b.sple,],2,mean)					###### ADDED BY MARTIN HUGHES 13/04/2015
			}
		outvar <- which(vari==0) ; outran <- which(ran==0)
		ifelse(length(outvar!=0),Variance <- cbind(Variance,rbind(mean(vari[-c(outvar)]),sd(vari[-c(outvar)]))),Variance <- cbind(Variance,rbind(mean(vari),sd(vari)))) # Bootstrap estimates and errors
		ifelse(length(outran!=0),Range <- cbind(Range,rbind(mean(ran[-c(outran)]),sd(ran[-c(outran)]))),Range <- cbind(Range,rbind(mean(ran),sd(ran))))       # id.
		skewness.b.sple <- c(skewness.b.sple,skewness(vari))
		
		MeanBootCentroid[tps,] <- apply(centri,2,mean)				###### ADDED BY MARTIN HUGHES 13/04/2015
		SDBootCentroid<- c(SDBootCentroid,apply(centri,2,sd))		###### ADDED BY MARTIN HUGHES 13/04/2015
		
		}
	list(Variance=Variance, Range=Range,skewness.b.sple=skewness.b.sple,MeanBootCentroid=MeanBootCentroid,SDBootCentroid=SDBootCentroid)
	}

analysis <- DtT(MS,PA,1000) 		 

#########	14/04/2015	possible code for Dmorpho as in HuangEA15 equation but applied to the centroids of the morphospace and time slices rather than comparing families and species
MeanBootCentroid<- analysis$MeanBootCentroid
SDBootCentroid <- analysis$SDBootCentroid
TotspCentroid<-apply(MS,2,mean)

DmrphTS<-c()
for (i in 2:ncol(PA)) {

	dmrphts<-c()
	for (j in 1:ncol(MeanBootCentroid)) {
		
		dmrphts<-c(dmrphts,		((MeanBootCentroid[i,j]-MeanBootCentroid[i-1,j])/diff(range(MS[,j])))^2)
	}

	DmrphTS<-c(DmrphTS,sqrt(sum(dmrphts,na.rm=TRUE)))
	
}
	
	DmrphTS<-c(0,DmrphTS)


	
Complete.DmrphTS <- matrix(99,nrow=1,length(All.PA.cond.intvls))

for (i in 1:length(DmrphTS)) {															
	v <- which(All.PA.cond.intvls==colnames(PA)[i])
	Complete.DmrphTS[,v] <- DmrphTS[i]
	Complete.DmrphTS[Complete.DmrphTS==99]<-NA
	}	
	
	
pdf(file=paste(Timespan,author.name,"DmrphTS.pdf",sep="_") ,height=11.7, width=16.5) 		# write disparity results (SoV/R, rarefaction curves, stratigrapahy, morphospace) plots to a pdf (portrait)
par(mar= c(5, 4, 4, 2) + 3.5) 

plot(x=c(1:ncol(Complete.PA.cond)),Complete.DmrphTS,ylab="DmorphoCentroid",xlab="Stratigraphic stage",xaxt="n",mgp=c(5,1.5,0),type="l")
points(x=c(1:ncol(Complete.PA.cond)),Complete.DmrphTS)
staxlab(side=1,labels=colnames(Complete.PA.cond),nlines=2,top.line=0,line.spacing=2, srt=45)
abline(v=(1:ncol(Complete.PA.cond))+.5,col="grey",lty=3)

dev.off()

D.mor.mat<-rbind(colnames(Complete.PA.cond), Complete.DmrphTS)
rownames(D.mor.mat)<-c("Stage","DmrphTS")
write.table(D.mor.mat,"a.txt",sep="\t")






