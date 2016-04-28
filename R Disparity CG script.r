
#####################################
#									#
#			Read in data			#
#									#
#####################################

# If characters are ordered, list which ones are in the brackets, otherwise the default is NO ORDERING 
OCHAR<-c()


## ...CODE TO ACCESS LIBRARIES
source("DISPARITY_functions.R")																		# most of the major functions are in this file, version produces the distance matrix via Wills 98 instead of original foote 93
library(calibrate) 																					# needed for a function for producing alternative strat data
library(zoo)																						# needed for funtion na.approx that is used to estimate values between points in rarefaction curves
library(plotrix)																					# needed for function staxlab which can be used to rotate labels and stack them)	
library(moments)																					# needed for the skewness function that is calculated for each time bin (skewness of the bootstrap replicate distribution)



###############
# READ IN FILES FROM CORRECT LOCATION INTO R FOR ANALYSIS
###############
## These find the file names in whatever working directory r is looking at
filelist<-list.files(pattern="STAGE.txt")															# find names of files in current working directory

## print(filelist) # IF USING AQUILA THIS MUST NOT BE BLOCKED
Datesfile<-filelist[1]																				# from filelist finds the presence/absence matrix
Morphologyfile<-filelist[2]																			# from filelist finds the morphology matrix
author.name<-substr(Morphologyfile,1,(nchar(Morphologyfile)-12))									# isolate the first part of the file name (ie removing the "_M_STAGE.txt" part)

m<-read.table(Morphologyfile)																		# read in the morphology matrix
M<-as.matrix(m[,-c(1)]) 																			# make as a matrix AND REMOVES FIRST COLLUMN IE THE TAXA NAMES
rownames(M)<-(m[,1]) 																				# make rownames vector and adds back into matrix
pa<-read.table(Datesfile,header=T)																	# read in the stratigraphic range data (stage level)
timetable<-read.table("TimeBinLookup.txt",header=T,sep="\t") 										# read in look-up file for geological timescale table
PA<-as.matrix(pa[,-c(1)])  																			# make as a matrix AND REMOVES FIRST COLLUMN IE THE TAXA NAMES
rownames(PA)<-(pa[,1]) 																				# make rownames vector and adds back into matrix

Orig.PA<-PA[,-1]																					# the original PA (matrix formated) before condensing or culling minus the condense collumn of pa
Orig.PA.Tax.Div <- apply(Orig.PA,2,sum)																# the taxonomic diversity of Orig.PA

#############
# below shows the code that will condense the data matrices PA and M if the condense
# list is shorter than the original after unquie characters are used
############
modal <- function(x) {
        a <- unique(x)				
        b <- rep(NA,length(a))
        for (i in 1:length(a)) b[i] <- length(which(x==a[i]))
        m <- a[which(b==max(b))]
}

con1 <- PA[,1]																						# condense1 is the 1st collumn of PA ie the condense column of original pa matrix (2nd collumn in that)
PA <- PA[,-c(1)]																					# PA is then everything except that condense column	
con2 <- unique(con1)																				# condense2 is the unique values only of the condense column


## if con2 is NOT the same as con1 ie shorter due to unique characters only being used
## the code will rebuild the matrices PA and M so that they are now the correct condensed size.

if (length(con2)!=length(con1)) {
con3 <- c()																							# condense3 is showing the places along condense2 (ie the new size of the matrix) where condensing will occur and how many taxa will be condensed
for (i in 1:length(con2)) con3 <- c(con3,length(which(con1==con2[i])))
con4 <- which(con3>1)																				# condense4 is showing where along the line in the new matrix the condensed taxa will be located 
con5<-cumsum(con3)
revised<-M[con5,]
revisedPA <- matrix(NA,max(con2),ncol(PA))
for (i in con4) {
        tbc <- M[which(con1==i),]																	# tbc= To be condensed
        new <- c()																					# new = new matrix
        for (j in 1:ncol(tbc)) {
                c <- modal(tbc[,j])
                if (length(c)>1) c <- sample(c,1,replace=FALSE)
                new <- c(new,c)
        }
        revised[i,] <- new
        revisedPA[i,] <- PA[which(con1==i)[1],]
}

## the code here gives the letter k the locations for which condense3 equals one
for (k in which(con3==1)) {
	who <- which(con1==k)
	revisedPA[k,] <- PA[who,]
}
M <- revised
PA <- revisedPA
}



## This code will rename the new matrix with the names
nm <- pa[,1]			#nm is just name but short	
naming <- c()
for (i in con2) {
	naming <- c(naming,as.character(nm[which(con1==i)[1]]))
}
## naming - vector used to give new matrix the correct taxon name from original matrix, the name taken from groups of taxa to be condensed is just the name of the 1st taxon of that group
rownames(PA)<- naming
colnames(PA)<-names(pa[,-c(1:2)])
rownames(M) <- naming


## First a distance matrix is produced 	
DistM<-MD(M,OCHAR) 																					# distance calculation for discrete characters
ResPCO<-PCO(DistM$D,"Cailliez")																		# PCO from the data (use "Cailliez" or "Lingoes") 
MS <- ResPCO$pco 																					# makes the morphospace coords from PCO a new object for use in analysis.
dissim <- c(as.dist(DistM$D))																		# dissimilarity between taxa (based on distances)
dispco <- c(dist(MS))																				# distances for taxa in PCO

																					# this is the PA matrix after condensing (missing stages still missing) column but before culling taxa that are impossible to place due to missing data
PA.cond.Tax.Div <- apply(PA,2,sum)																# the PA matrix used for taxonomic diversity prior to condensing or other editing
PA.cond<-PA


Stage <- PA																					# the PA matrix is then used to produce time bins of different SPECIFIED durations via the time bin converter source file	
Epoch<-TimeBinConverter(PA,"Epoch")
Period<-TimeBinConverter(PA,"Period")

Timeseq <- c("St","Ep","Pe")


#####################################
# This is a giant loop, it will run the analysis and 
# plot the data for each type of time partition,
#####################################
for (Timespan in Timeseq) {
	if (Timespan=="St") PA <- Stage
	if (Timespan=="Ep") PA <- Epoch
	if (Timespan=="Pe") PA <- Period



## takes the original stratigraphy  (stage level) and removes the taxa that were removed from the kill list
if (length(DistM$kill)!=0) {		
	PA<-PA[-c(DistM$kill),]
	}


## checks for intervals of time that are now empty because taxa have been killed via above code
if (length(which(colSums(PA)==0))>0) {
	removed_time_intervals <- which(colSums(PA)==0)													# removes intervals of time that = 0 because taxa are removed from them
	PA <- PA[,-(removed_time_intervals)]															# modifies the PA matrix to reflect this removed time
	}																								# if matrix does not contain intervals with no taxa, nothing happens to the PA matrix 

##########################	
# actual disparity analysis
##########################	
Taxonomic.Diversity <- apply(PA,2,sum)
analysis <- DtT(MS,PA,1000) 			 


# for stage
 Variance <- analysis$Variance
 Range <- analysis$Range
 Variance[which(Variance==0)] <- NA																	# makes 0 values NA
 Range[which(Range==0)] <- NA																		# makes 0 values NA
 
	
#function calculates the variance for each time bin with taxa resampled (same amount as bootstrap) from entire morphospace (ie expected variance/range of random sample)
TD <- apply(PA,2,sum) ; NMDISP_SoV <-c()
for (i in 1:ncol(PA)) {
        hm <- TD[i]
        nmdisp_SoV <- c()
        for (boo in 1:1000) {
                sple <- sample(1:nrow(MS),hm,replace=T)
				if (length(sple)==1) sple <- c(sple,sple)
                nmdisp_SoV <- c(nmdisp_SoV,SoV(MS[sple,]))											# this line looks at SoV
                }
        NMDISP_SoV <- cbind(NMDISP_SoV,rbind(mean(nmdisp_SoV),sd(nmdisp_SoV)))
		
		NMDISP_SoV[which(NMDISP_SoV==0)] <- NA
		}
	
# same as above but for range (also rarefies)
NMDISP_SoR <-c()
for (i in 1:ncol(PA)) {
        hm <- TD[i]
		hm <- min(TD) 																				# this uses minimum value for time bin to calculate disparity 
		if (hm==1) hm <- 2 																			# if the lowest value is one, make it 2 (cant compute range with only a single taxon!)
        nmdisp_SoR <- c()
        for (boo in 1:1000) { 																		# bootstrap using taxa from entire morphospace
                sple <- sample(1:nrow(MS),hm,replace=T)
				if (length(sple)==1) sple <- c(sple,sple)
                nmdisp_SoR <- c(nmdisp_SoR,SoR(MS[sple,]))											# this line is what makes it look at SoR
                }
        NMDISP_SoR <- cbind(NMDISP_SoR,rbind(mean(nmdisp_SoR),sd(nmdisp_SoR)))
		NMDISP_SoR[which(NMDISP_SoR==0)] <- NA
		}

###	calculates the rarefaction curve data
rec <- 0
BIG <- array(c(NA),dim=c(2,max(apply(PA,2,sum)),ncol(PA)))											# an array that will contain all of the rarefraction curves
timescan<- as.numeric(which(apply(PA,2,sum)>2)) 													# R will find all the time bins that are NOT 1 (ie only 1 taxon and compute rarefaction curves)
for (i in timescan) {	
	ms <- MS[which(PA[,i]==1),]
	if (nrow(ms)>50) { big <-(RarCrv_big(ms,500));skipseq <- big$skipseq							# if there are more than 50 taxa in any timebin, it will use RarCrv_big and sample 
		} else { big <-(RarCrv(ms,500));skipseq <- 2:nrow(ms) }										# 50 of those taxa to generate Rarefaction curve, otherwise it will sample all taxa
	
	big <- as.matrix(big$RarCrv)
	rec2 <- max(apply(big,2,sum))
	if (rec2>rec) rec <- rec2
	BIG[1:2,skipseq,i] <- big
	}	
	

#####################################
### Centre of Gravity of a clade	-	only applied to stage level

# This first set of code is to calculate the centre of gravity for the taxonomic diversity curve
# find the values of each stage in the PA condensed matrix (PA.cond) prior to any culling (PA.cull) - but includes all PA stages between oldest and youngest intervals of time



PA.cond.YS<-tail(colnames(PA.cond), n=1)															# PA.cond.YS - PA.cond youngest stage
PA.cond.OS<-head(colnames(PA.cond), n=1)															# PA.cond.OS - PA.cond oldest stage
PA.cond.redTT<- timetable[c(which(timetable[,3]==PA.cond.OS):which(timetable[,3]==PA.cond.YS)),]	# PA.cond.redTT - reduced timetable which only includes the intervals between oldest and youngest stages in main timetable table

if (Timespan=="St") {	
	All.PA.cond.intvls<-PA.cond.redTT[,3]															# all the stage names from the reduced timetable
	All.PA.cond.intvls<-as.character(All.PA.cond.intvls)											# removed the levels from the stage names vector
	Complete.PA.cond <- matrix(0,nrow(Stage),length(All.PA.cond.intvls))							# Complete.PA.cond - the revised PA.cond matrix, which first is filled with 0 but then includes all stages with observations ie ones with and without taxa but filled only with 0
	
		for (i in 1:length(colnames(Stage))) {														# Loop fills the stages of the revised PA.cond matrix that contain taxa
				p <- which(All.PA.cond.intvls==colnames(Stage)[i])
				Complete.PA.cond[,p] <- Stage[,i]
		}		
	intvl.myrs<-PA.cond.redTT[,2]														# the duration of each stage from the reduced timetable
	intvl.stage.myrs<-intvl.myrs
}	

	
if (Timespan=="Ep") {	
	All.PA.cond.intvls<-unique(PA.cond.redTT[,5])													# all the epoch names from the reduced timetable
	All.PA.cond.intvls<-as.character(All.PA.cond.intvls)
	Complete.PA.cond <- matrix(0,nrow(Epoch),length(All.PA.cond.intvls))	
		
		for (i in 1:length(colnames(Epoch))) {														# Loop fills the stages of the revised PA.cond matrix that contain taxa
				p <- which(All.PA.cond.intvls==colnames(Epoch)[i])
				Complete.PA.cond[,p] <- Epoch[,i]
		}
	intvl.epoch.myrs<-matrix(NA,nrow=1,ncol=length(All.PA.cond.intvls)) 						# isolate the number of myrs within each epoch
	colnames(intvl.epoch.myrs)<-All.PA.cond.intvls
	for (i in 1:length(All.PA.cond.intvls)) 	intvl.epoch.myrs[,i]<-sum(PA.cond.redTT[which(PA.cond.redTT[,5]==All.PA.cond.intvls[i]),2])
	intvl.myrs<- intvl.epoch.myrs[1,]													# turn this into a vector
}


if (Timespan=="Pe") {	
	All.PA.cond.intvls<-unique(PA.cond.redTT[,6])													# all the period names from the reduced timetable
	All.PA.cond.intvls<-as.character(All.PA.cond.intvls)
	Complete.PA.cond <- matrix(0,nrow(Period),length(All.PA.cond.intvls))		
		for (i in 1:length(colnames(Period))) {														# Loop fills the stages of the revised PA.cond matrix that contain taxa
				p <- which(All.PA.cond.intvls==colnames(Period)[i])
				Complete.PA.cond[,p] <- Period[,i]
		}
	intvl.period.myrs<-matrix(NA,nrow=1,ncol=length(All.PA.cond.intvls)) 						# isolate the number of myrs within each epoch
	colnames(intvl.period.myrs)<-All.PA.cond.intvls
	for (i in 1:length(All.PA.cond.intvls)) 	intvl.period.myrs[,i]<-sum(PA.cond.redTT[which(PA.cond.redTT[,6]==All.PA.cond.intvls[i]),2])
	intvl.myrs<- intvl.period.myrs[1,]													# turn this into a vector
}


Complete.PA.cond.tax.div <- colSums(Complete.PA.cond)												# the taxonomic divesity in the revised PA.cond matrix
colnames(Complete.PA.cond)<-All.PA.cond.intvls														# renaming the columns to have correct names
rownames(Complete.PA.cond)<-rownames(PA.cond)														# renaming the rows to have correct names

COG.Tax<- CG.taxo (Complete.PA.cond,	intvl.myrs)													# The taxonomic diversity Centre of gravity outputs


# This code is to calculate the centre of gravity for the morphological diversity curves (Variance ONLY)
if (length(DistM$kill)!=0) {																		# This keeps the number of PA stages the same but the taxa that were removed as a result of too much missing data are taken away
	Complete.PA.cull<-Complete.PA.cond[-c(DistM$kill),]
	} else {	Complete.PA.cull<-Complete.PA.cond	}

COG.Var<- CG.morpho(Complete.PA.cull,MS,intvl.myrs)													# The sum of variance Centre of gravity outputs




sym.early.high.test<-Max.morpho(PA,MS)

bt.strp.quantiles.mx<-quantile(sym.early.high.test$reDiffmx,probs=c(0.025,0.975)) # for the max
bt.strp.quantiles.mn<-quantile(sym.early.high.test$reDiffmn,probs=c(0.025,0.975)) # for the min



if (bt.strp.quantiles.mx[1]<0 & bt.strp.quantiles.mx[2]>0) {
		sym.e.high.or.bell<-	"e.high"
		} else {	
		sym.e.high.or.bell<-	"concave" 
}


if (bt.strp.quantiles.mn[1]<0 & bt.strp.quantiles.mn[2]>0) {
		sym.e.low.or.bell<-	"e.low"
		} else {	
		sym.e.low.or.bell<-	"convex" 
}
	

# The is the matrix containing the COG for diversity, SoV and the values that can be then used to isolate those symmetrical lines that are flat (ie early max disparity - bottom heavy!)
COG.matrix<-t(cbind(COG.Tax,COG.Var))																# A matrix combining the different centre of gravity results

COG.result<-cbind(COG.matrix,c(NA,sym.e.high.or.bell),c(NA,sym.e.low.or.bell),c(NA,sym.early.high.test$CV))
colnames(COG.result)[10:12]<-c("sym.early.high.?","sym.early.low.?","CoeffVar")


##################
### Two bin tests
### The following culmilnates in three sets of results, with two pvalues : one for the first 2 bins (early) and one for the last two bins (late)
### TEST1 - observed disparity (bootstrapped) equals the maximum
##################

two.bin.test.table<-matrix(NA,1,6)
colnames(two.bin.test.table)<- c("ely.test1", "lt.test1","ely.test2", "lt.test2","ely.test3","lt.test3")
rownames(two.bin.test.table)<-c("pvalue")
t.b.max.test<-THE.INTEGRATED.GERBER.WILLS.TEST(MS, Complete.PA.cond)

two.bin.test.table[1,1:2]<- c(t.b.max.test[[1]],t.b.max.test[[2]])
two.bin.test.table[1,3:4]<-c(t.b.max.test[[3]],t.b.max.test[[4]])
two.bin.test.table[1,5:6]<-c(t.b.max.test[[5]],t.b.max.test[[6]])

write.table(two.bin.test.table,file=paste(Timespan,author.name,"twobintest_.txt",sep="_"), sep="\t")
two.bin.test.table
	
#####################################
### VARIANCE and RANGE VS "DIVERSITY" at the Diversification phase		- applied only to the stage level	


### Diversification phase the jablonski way:
### 1.find the highest value for diversity and disparity
### 	a. for diversity this is the final value if there are multiple peaks
### 	b. for disparity this is the final value if there are multiple peaks	-	
### 2.whichever occurs last temporally is the end of the diversification phase
### 3.plot these graphs

Max_Diversity<-which(Complete.PA.cond.tax.div==max(Complete.PA.cond.tax.div))						# finds which stages contain the maxium number of OTUs and makes a vector giving the name of the stage, AND the value of how far along the taxonomic.diversity vector it is (NOT what the stage value is!!)
if (length(Max_Diversity)>1) { 																		# if there are two or more stages with a max value, it takes the last one
	Max_Diversity<-tail(Max_Diversity,n=1)
	}


### These  next 2 loops use the All.PA.cond.stages vector to construct a new matrix for disparity Variance/Range results
### just like for the centre of gravity PA matrices,
### it puts the values of range and varaiance (and their error bars) in the correct temporal location in the new matrices called
### Complete.Var.cul (all stages that have variance which was produced from condensed/culled PA matrix)
### Complete.Rng.cull (same as above but for range)
Complete.Var.cull <- matrix(0,nrow(Variance),length(All.PA.cond.intvls))
colnames(Variance)<-colnames(PA)
for (i in 1:length(colnames(Variance))) {															
	v <- which(All.PA.cond.intvls==colnames(Variance)[i])
	Complete.Var.cull[,v] <- Variance[,i]
	Complete.Var.cull[Complete.Var.cull==0]<-NA
	}
	
Complete.Rng.cull <- matrix(0,nrow(Range),length(All.PA.cond.intvls))
colnames(Range)<-colnames(PA)
for (i in 1:length(colnames(Range))) {															
	r <- which(All.PA.cond.intvls==colnames(Range)[i])
	Complete.Rng.cull[,r] <- Range[,i]
	Complete.Rng.cull[Complete.Rng.cull==0]<-NA
	}


Max_Variance<-which(Complete.Var.cull[1,]==max(Complete.Var.cull[1,],na.rm=T))						# finds which stages contain the maxium variance and makes a vector giving the numerical location of the stage (which corresponds with counting left to right the number of stages in taxnonomic diversity)
if (length(Max_Variance)>1) { 																		# if there are two or more stages with a max value, it takes the last one
	Max_Variance<-tail(Max_Variance,n=1)
	}

Max_Range<-which(Complete.Rng.cull[1,]==max(Complete.Rng.cull[1,],na.rm=T))							# finds which stages contain the maximum range and makes a vector giving the numerical location of the stage (which corresponds with counting left to right the number of stages in taxnonomic diversity)
if (length(Max_Range)>1) { 																			# if there are two or more stages with a max value, it takes the last one
	Max_Range<-tail(Max_Range,n=1)
	}

	
#####################################
#									#
#		Plotting the data			#
#									#	
#####################################


strat_name<-colnames(Complete.PA.cond) 																# gives vector of names of the stages of PA matrix that has all stages (including the stages with no fossils) and has been condensed
Time <- 1:ncol(Complete.PA.cond)																	# gives the number of time intervals in PA matrix
Time.myrs1<-PA.cond.redTT[,7]																		# this is used to put geological timescale on graph
Time.myrs2<-PA.cond.redTT[,8]	
	

### R will decide what to save output as depending on whether you have specified ordered characters or not	
if (length(OCHAR)!=0) {
pdf(file=paste(Timespan,author.name,"DispRslt_ORD.pdf",sep="_") ,height=11.7, width=16.5) 		# write disparity results (SoV/R, rarefaction curves, stratigrapahy, morphospace) plots to a pdf (portrait)
} else { pdf(file=paste(Timespan,author.name,"DispRslt_UNORD.pdf",sep="_"),height=11.7, width=16.5) }

par(mar= c(5, 4, 4, 2) + 3.5) 																		# make margins larger (just change the +3 to whatever you want)
layout(matrix(c(1:6), 2, 3, byrow = TRUE))															# use this to set up a window in which each graph will fit in (specify the number of graphs it will hold and over how many rows/columns



#####################################
###plots the axis 1 & 2 of pco coords morphospace through time
plot(MS,asp=1,xlab="PCO Axis 1", ylab="PCO Axis 2", main="Morphospace",cex=1.5,cex.axis=1.5,cex.lab=2,cex.main=2.5,mgp=c(5,1.5,0))
abline(h=0, v=0, col="red")
#x11()



#####################################
### plots the TIME RANGES of taxa in a graph - for stages it plots all, including the stages with no fossils in them

plot(0,0,xlim=c(1,ncol(Complete.PA.cond)),ylim=c(1,nrow(Complete.PA.cond)),main="Stratigraphic Range",xlab="Time Bins",ylab="Taxa",cex.axis=1.5,cex.lab=2,cex.main=2.5,mgp=c(5,1.5,0),xaxt="n")
staxlab(side=1,labels=colnames(Complete.PA.cond),nlines=2,top.line=0,line.spacing=2, srt=45)
abline(v=(1:ncol(Complete.PA.cond))+.5,col="grey",lty=3)
cl <- rep(1,nrow(Complete.PA.cond))
cl[DistM$kill] <- 2
for (i in 1:nrow(Complete.PA.cond)) segments(min(which(Complete.PA.cond[i,]==1))-.5,i, max(which(Complete.PA.cond[i,]==1))+.5,i,lwd=3,col=cl[i])



#####################################
### plots "Taxonomic DIVERSITY" before any condensing or culling, after condensing and then condensing + culling

if (Timespan=="St") {
	Time.myrs1<-PA.cond.redTT[,7]
	plot(Time.myrs1,Complete.PA.cond.tax.div,xlab="Time Mya",main="OTU Diversity - All",cex=1.5,cex.axis=1.5,cex.lab=2,cex.main=2.5,mgp=c(5,1.5,0),pch=1,xlim=c(max(Time.myrs1),min(Time.myrs1)))
	lines(Time.myrs1,Complete.PA.cond.tax.div)
}

if (Timespan=="Ep") {
	first.stage<-head(PA.cond.redTT[match(All.PA.cond.intvls,PA.cond.redTT[,5]),7],n=1)
	last.stage<-tail(PA.cond.redTT[,7],n=1)-tail(PA.cond.redTT[,2],n=1)
	Time.myrs1<-PA.cond.redTT[match(All.PA.cond.intvls,PA.cond.redTT[,5]),7] - (intvl.myrs/2)
	plot(Time.myrs1,Complete.PA.cond.tax.div,xlab="Time Mya",main="OTU Diversity - All",cex=1.5,cex.axis=1.5,cex.lab=2,cex.main=2.5,mgp=c(5,1.5,0),pch=1,xlim=c(first.stage,last.stage))
	lines(Time.myrs1,Complete.PA.cond.tax.div)
}

if (Timespan=="Pe") {
	first.stage<-head(PA.cond.redTT[match(All.PA.cond.intvls,PA.cond.redTT[,6]),7],n=1)
	last.stage<-tail(PA.cond.redTT[,7],n=1)-tail(PA.cond.redTT[,2],n=1)
	Time.myrs1<-PA.cond.redTT[match(All.PA.cond.intvls,PA.cond.redTT[,6]),7] - (intvl.myrs/2)
	plot(Time.myrs1,Complete.PA.cond.tax.div,xlab="Time Mya",main="OTU Diversity - All",cex=1.5,cex.axis=1.5,cex.lab=2,cex.main=2.5,mgp=c(5,1.5,0),pch=1,xlim=c(first.stage,last.stage))
	lines(Time.myrs1,Complete.PA.cond.tax.div)
}


#####################################
### plot VARIANCE VS TIME
### there are two approaches to this. IF its plotting at stage level, it will use the top 2 plot schemes, so it will show the
### stages with no fossils
### and for all other time binning (sub/epoch/period) it will do it like scripts V15.

# This rebuild the NMDISP_SoV so the values are in the correct temporal locations 
Complete.NMDISP.SoV <- matrix(0,nrow(NMDISP_SoV),length(All.PA.cond.intvls))						
colnames(NMDISP_SoV)<-colnames(PA)
for (i in 1:length(colnames(NMDISP_SoV))) {															
	v2 <- which(All.PA.cond.intvls==colnames(NMDISP_SoV)[i])
	Complete.NMDISP.SoV[,v2] <- NMDISP_SoV[,i]
	Complete.NMDISP.SoV[Complete.NMDISP.SoV==0]<-NA
	}

# This rebuild the NMDISP_SoR so the values are in the correct temporal locations	
Complete.NMDISP.SoR <- matrix(0,nrow(NMDISP_SoR),length(All.PA.cond.intvls))
colnames(NMDISP_SoR)<-colnames(PA)
for (i in 1:length(colnames(NMDISP_SoR))) {															
	r2 <- which(All.PA.cond.intvls==colnames(NMDISP_SoR)[i])
	Complete.NMDISP.SoR[,r2] <- NMDISP_SoR[,i]
	Complete.NMDISP.SoR[Complete.NMDISP.SoR==0]<-NA
	}

# Gives the axis in million yrs (but only do-able for STAGE level) 
if (Timespan=="St") {
	Time.myrs1<-PA.cond.redTT[,7]																		# this is used to put geological timescale on graph
	Time.myrs2<-PA.cond.redTT[,8]	
	###plot VARIANCE VS TIME
	plot(Time.myrs2,Complete.Var.cull[1,],ylim=c(0,(1.5*max(apply(Complete.Var.cull,2,sum,na.rm=T)))),pch=15,ylab="Total variance", main="Sum of Variances",xlab="Time MYRS",cex.axis=1.5,cex.lab=2,cex.main=2.5,mgp=c(5,1.5,0),xlim=c(max(Time.myrs1),min(Time.myrs2)-1))
	lines(Time.myrs2,Complete.Var.cull[1,],lwd=2)
	segments(Time.myrs2,Complete.Var.cull[1,]-Complete.Var.cull[2,],Time.myrs2,Complete.Var.cull[1,]+Complete.Var.cull[2,])
	points(Time.myrs2,Complete.NMDISP.SoV[1,],pch=16, col=2,lwd=2)													#points plotted over top of SoV graphs
	points(Time.myrs2,Complete.NMDISP.SoV[1,]+Complete.NMDISP.SoV[2,],col=2, pch=16)
	points(Time.myrs2,Complete.NMDISP.SoV[1,]-Complete.NMDISP.SoV[2,],col=2,pch=16)
	lines(Time.myrs2,Complete.NMDISP.SoV[1,],col=2,lwd=2)															#lines plotted over top of SoV graphs
	lines(Time.myrs2,Complete.NMDISP.SoV[1,]+Complete.NMDISP.SoV[2,],col=2)
	lines(Time.myrs2,Complete.NMDISP.SoV[1,]-Complete.NMDISP.SoV[2,],col=2)
	lines(Time.myrs2[which(Complete.Var.cull[1,]!="NA")],disparityCRVfun(MS,PA),col="green")	# The actual observed disparity (for variances
	###plot RANGE VS TIME
	plot(Time.myrs2,Complete.Rng.cull[1,],ylim=c(0,(1.5*max(apply(Complete.Rng.cull,2,sum,na.rm=T)))),pch=15,ylab="Total range",main="Sum of Ranges",xlab="Time MYRS",cex.axis=1.5,cex.lab=2,cex.main=2.5,mgp=c(5,1.5,0),xlim=c(max(Time.myrs1),min(Time.myrs2)))
	lines (Time.myrs2,Complete.Rng.cull[1,],lwd=2) 
	segments(Time.myrs2,Complete.Rng.cull[1,]-Complete.Rng.cull[2,],Time.myrs2,Complete.Rng.cull[1,]+Complete.Rng.cull[2,])
	points(Time.myrs2,Complete.NMDISP.SoR[1,],pch=16,col=2,lwd=2)													#points plotted over top of SoR graphs
	points(Time.myrs2,Complete.NMDISP.SoR[1,]+Complete.NMDISP.SoR[2,],col=2,pch=16)
	points(Time.myrs2,Complete.NMDISP.SoR[1,]-Complete.NMDISP.SoR[2,],col=2,pch=16)
	lines(Time.myrs2,Complete.NMDISP.SoR[1,],col=2,lwd=2)															#lines plotted over top of SoR graphs
	lines(Time.myrs2,Complete.NMDISP.SoR[1,]+Complete.NMDISP.SoR[2,],col=2)
	lines(Time.myrs2,Complete.NMDISP.SoR[1,]-Complete.NMDISP.SoR[2,],col=2)
}

if (Timespan=="Ep") {
	first.stage<-head(PA.cond.redTT[match(All.PA.cond.intvls,PA.cond.redTT[,5]),7],n=1)
	last.stage<-tail(PA.cond.redTT[,7],n=1)-tail(PA.cond.redTT[,2],n=1)
	Time.myrs1<-PA.cond.redTT[match(All.PA.cond.intvls,PA.cond.redTT[,5]),7] - (intvl.myrs/2)
	###plot VARIANCE VS TIME
	plot(Time.myrs1,Complete.Var.cull[1,],xlim=c(first.stage,last.stage),ylim=c(0,(1.5*max(apply(Complete.Var.cull,2,sum,na.rm=T)))),pch=15,ylab="Total variance", main="Sum of Variances",xlab="Time Mya",cex.axis=1.5,cex.lab=2,cex.main=2.5,mgp=c(5,1.5,0)) 
	lines(Time.myrs1,Complete.Var.cull[1,],lwd=2)
	segments(Time.myrs1,Complete.Var.cull[1,]-Complete.Var.cull[2,],Time.myrs1,Complete.Var.cull[1,]+Complete.Var.cull[2,])
	points(Time.myrs1,Complete.NMDISP.SoV[1,],pch=16, col=2,lwd=2)													#points plotted over top of SoV graphs
	points(Time.myrs1,Complete.NMDISP.SoV[1,]+Complete.NMDISP.SoV[2,],col=2, pch=16)
	points(Time.myrs1,Complete.NMDISP.SoV[1,]-Complete.NMDISP.SoV[2,],col=2,pch=16)
	lines(Time.myrs1,Complete.NMDISP.SoV[1,],col=2,lwd=2)															#lines plotted over top of SoV graphs
	lines(Time.myrs1,Complete.NMDISP.SoV[1,]+Complete.NMDISP.SoV[2,],col=2)
	lines(Time.myrs1,Complete.NMDISP.SoV[1,]-Complete.NMDISP.SoV[2,],col=2)
	lines(Time.myrs1[which(Complete.Var.cull[1,]!="NA")],disparityCRVfun(MS,PA),col="green")	# The actual observed disparity (for variances
	###plot RANGE VS TIME
	plot(Time.myrs1,Complete.Rng.cull[1,],xlim=c(first.stage,last.stage),ylim=c(0,(1.5*max(apply(Complete.Rng.cull,2,sum,na.rm=T)))),pch=15,ylab="Total range",main="Sum of Ranges",xlab="Time Mya",cex.axis=1.5,cex.lab=2,cex.main=2.5,mgp=c(5,1.5,0))
	lines (Time.myrs1,Complete.Rng.cull[1,],lwd=2) 
	segments(Time.myrs1,Complete.Rng.cull[1,]-Complete.Rng.cull[2,],Time.myrs1,Complete.Rng.cull[1,]+Complete.Rng.cull[2,])
	points(Time.myrs1,Complete.NMDISP.SoR[1,],pch=16,col=2,lwd=2)													#points plotted over top of SoR graphs
	points(Time.myrs1,Complete.NMDISP.SoR[1,]+Complete.NMDISP.SoR[2,],col=2,pch=16)
	points(Time.myrs1,Complete.NMDISP.SoR[1,]-Complete.NMDISP.SoR[2,],col=2,pch=16)
	lines(Time.myrs1,Complete.NMDISP.SoR[1,],col=2,lwd=2)															#lines plotted over top of SoR graphs
	lines(Time.myrs1,Complete.NMDISP.SoR[1,]+Complete.NMDISP.SoR[2,],col=2)
	lines(Time.myrs1,Complete.NMDISP.SoR[1,]-Complete.NMDISP.SoR[2,],col=2)
}

if (Timespan=="Pe") {
	first.stage<-head(PA.cond.redTT[match(All.PA.cond.intvls,PA.cond.redTT[,6]),7],n=1)
	last.stage<-tail(PA.cond.redTT[,7],n=1)-tail(PA.cond.redTT[,2],n=1)
	Time.myrs1<-PA.cond.redTT[match(All.PA.cond.intvls,PA.cond.redTT[,6]),7] - (intvl.myrs/2)
	###plot VARIANCE VS TIME
	plot(Time.myrs1,Complete.Var.cull[1,],xlim=c(first.stage,last.stage),ylim=c(0,(1.5*max(apply(Complete.Var.cull,2,sum,na.rm=T)))),pch=15,ylab="Total variance", main="Sum of Variances",xlab="Time Mya",cex.axis=1.5,cex.lab=2,cex.main=2.5,mgp=c(5,1.5,0)) 
	lines(Time.myrs1,Complete.Var.cull[1,],lwd=2)
	segments(Time.myrs1,Complete.Var.cull[1,]-Complete.Var.cull[2,],Time.myrs1,Complete.Var.cull[1,]+Complete.Var.cull[2,])
	points(Time.myrs1,Complete.NMDISP.SoV[1,],pch=16, col=2,lwd=2)													#points plotted over top of SoV graphs
	points(Time.myrs1,Complete.NMDISP.SoV[1,]+Complete.NMDISP.SoV[2,],col=2, pch=16)
	points(Time.myrs1,Complete.NMDISP.SoV[1,]-Complete.NMDISP.SoV[2,],col=2,pch=16)
	lines(Time.myrs1,Complete.NMDISP.SoV[1,],col=2,lwd=2)															#lines plotted over top of SoV graphs
	lines(Time.myrs1,Complete.NMDISP.SoV[1,]+Complete.NMDISP.SoV[2,],col=2)
	lines(Time.myrs1,Complete.NMDISP.SoV[1,]-Complete.NMDISP.SoV[2,],col=2)
	lines(Time.myrs1[which(Complete.Var.cull[1,]!="NA")],disparityCRVfun(MS,PA),col="green")	# The actual observed disparity (for variances
	###plot RANGE VS TIME
	plot(Time.myrs1,Complete.Rng.cull[1,],xlim=c(first.stage,last.stage),ylim=c(0,(1.5*max(apply(Complete.Rng.cull,2,sum,na.rm=T)))),pch=15,ylab="Total range",main="Sum of Ranges",xlab="Time Mya",cex.axis=1.5,cex.lab=2,cex.main=2.5,mgp=c(5,1.5,0))
	lines (Time.myrs1,Complete.Rng.cull[1,],lwd=2) 
	segments(Time.myrs1,Complete.Rng.cull[1,]-Complete.Rng.cull[2,],Time.myrs1,Complete.Rng.cull[1,]+Complete.Rng.cull[2,])
	points(Time.myrs1,Complete.NMDISP.SoR[1,],pch=16,col=2,lwd=2)													#points plotted over top of SoR graphs
	points(Time.myrs1,Complete.NMDISP.SoR[1,]+Complete.NMDISP.SoR[2,],col=2,pch=16)
	points(Time.myrs1,Complete.NMDISP.SoR[1,]-Complete.NMDISP.SoR[2,],col=2,pch=16)
	lines(Time.myrs1,Complete.NMDISP.SoR[1,],col=2,lwd=2)															#lines plotted over top of SoR graphs
	lines(Time.myrs1,Complete.NMDISP.SoR[1,]+Complete.NMDISP.SoR[2,],col=2)
	lines(Time.myrs1,Complete.NMDISP.SoR[1,]-Complete.NMDISP.SoR[2,],col=2)
}

#####################################
### plot the rare faction curves for each time bin using the taxa within that bin
plot(0,0, typ='n',xlim=c(0,max(apply(PA,2,sum))),ylim=c(0,rec),xlab="Number of Taxa sampled", ylab="Range", main="Rarefaction Curve: Stage Level",cex.axis=1.5,cex.lab=2,cex.main=2.5,mgp=c(5,1.5,0))
for (i in timescan) {
	lines(na.approx(BIG[1,,i]),col=i)
	lines(na.approx(BIG[1,,i]+BIG[2,,i]),col=i,lty=3)
	lines(na.approx(BIG[1,,i]-BIG[2,,i]),col=i,lty=3)
	mx_value <- max(BIG[1,,i],na.rm=T)							 									# last 3 lines get lines that are plotted, find the final point of line and plot the time bin name
	x_max <- which(BIG[1,,i]==mx_value) 
	text(x_max,mx_value,i)
}
dev.off()



#####################################
### morphospace timeslices

if (length(OCHAR)!=0) {
pdf(file=paste(Timespan,author.name,"Morphsp_ORD.pdf",sep="_"),height=11.7, width=16.5) # write morphosapce plots to a pdf (portrait)
} else { pdf(file=paste(Timespan,author.name,"Morphsp_UNORD.pdf",sep="_"),height=11.7, width=16.5) } # name of file _morphospace time slices_version of scritp_ordered or unordered


if (ncol(PA)<=5) {	
	par(mar= c(5, 4, 4, 2) + 3.5) 	# make margins larger (just change the +3 to whatever you want)
	layout(matrix(c(1:ncol(PA)), 1, 1, byrow = TRUE))
} else {
	par(mar= c(5, 4, 4, 2) + 3.5) 	# make margins larger (just change the +3 to whatever you want)
	layout(matrix(c(1:ncol(PA)), 2, 3, byrow = TRUE))
}

for (i in 1:ncol(PA)) {
	time.names<-colnames(PA)
	if (length(c(DistM$kill))==0) { taxa.labels<-rownames(M) 
	} else { taxa.labels<-rownames(M[-DistM$kill,]) }
	rownames(MS)<-taxa.labels
	plot(MS[,1:2],pch=20,col="grey",cex=1.5,cex.axis=1.5,cex.lab=2,cex.main=2.5,mgp=c(3.5,1.5,0),main=(time.names[i]))
	points(MS[which(PA[,i]==1),1],MS[which(PA[,i]==1),2],ps=100,pch=19,col="red")
	textxy(MS[which(PA[,i]==1),1],MS[which(PA[,i]==1),2],labs=labels(MS[which(PA[,i]==1),1]),cx=0.5, dcol = "black", m = c(0, 0))
}
dev.off()



#####################################
### output the COG in a table
write.table(COG.result,file=paste(Timespan,author.name,"COG_.txt",sep="_"), sep="\t")

	
#####################################
### save the workspace image ie all stuff imputted into R is there (functions, results) so you can recall any part of the above code
if (length(OCHAR)!=0) { save.image(file=paste(Timespan,author.name,"DispWksp_ORD.RData",sep="_"))
} else { save.image(file=paste(Timespan,author.name,"DispWksp_UNORD.RData",sep="_")) } 			# name of file_DisparityWorkspace_version of scritp_ordered or unordered		

}																									# End of Loop for stage, sub/epoch, period	


#####################################
###The End :D


#### Output terms that are useful to be saved

# DistM - Distance matrix produced from discrete morphological character matrix						# write.table(DistM$D,file=paste(Timespan,author.name,"PCO_coords.txt",sep="_"), sep="\t")
# MS - Morphospace produced from DistM using the PCO function										# write.table(MS,file=paste(Timespan,author.name,"PCO_coords.txt",sep="_"), sep="\t")
# analysis - analysis results producing the sums of ranges and sums of variances from the MS, using user specified bootstraps to make error bars
# BIG - the values to get rarefaction curves for each time bin
# COG.Tax - the centre of gravity results for taxonomic diversity	
# COG.Var - the centre of gravity results for sum of variances	
# CoG.matrix - matrix combining the 2 Centre of gravity values										# write.table(COG.matrix,file=paste(Timespan,author.name,"COG.txt",sep="_"), sep="\t")





	






