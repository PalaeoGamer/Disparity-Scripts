#####################################
#	This section are the functions	#
#####################################


### Function to make distance matrix from character data 
 MD <- function(X,ochar,equal.w=TRUE) {				# X:morphology matrix, ochar: list of ordered characters, equal.w: equal weight is TRUE
	irr <- c()
	for (i in 1:nrow(X)) irr <- c(irr,length(na.omit(X[i,])))	# This loop reads down rows, searching for taxa that are just NA (maybe due to alot of missing data and poly morphs or...... the stupidity of the researcher who made the matrix!)
	remtax <- length(which(irr==0))		# remtax is removed taxa where they are totally NA for morphology matrix (eg only ? and polymorphic data etc)
	if (remtax!=0) {
		print("Uninformative taxa removed",quote=FALSE)
		print(as.numeric(which(irr==0)))
		X <- X[-c(which(irr==0)),]	#Matrix has been modified then to not include those taxa that are just NA (taxa that are all 0's, 1's etc are NOT affected)
		}
		w <- apply(X,2,max,na.rm=T) - apply(X,2,min,na.rm=T)	# This line is looking for characters that dont vary eg column of 1's etc, (discludes NA from the search), w is the range between max and min value in character state for each character 
	bad <- which(w==0)								# "bad" value search 1: Those columns that dont vary ie have a value of 0 (due to process in line above) are labled bad
	for (i in 1:ncol(X)) {							# "bad" value search 2: looks in columns for unique characters eg column of NA's, 1's etc and removed these also
		if (length(unique(X[,i]))==1) {bad <- c(bad,i)}
		}
	bad <- sort(unique(bad)) # to stop it outputting columns picked up twice (once for each "bad" value search) it uses the unique function again
		
		if (length(bad!=0)) {
		print("Uninformative characters not included:", quote=FALSE) # Uniformative characters include columns of one character state (with or without NA's)
		print(as.numeric(bad))								
	}												# AN ERROR MESSAGE WILL APPREAR BUT THIS DOES NOT MATTER AS IT IS TAKEN CARE OFF IN THE FUNCTION WHEN COLUMNS OF UNINFORMATIVE CHARACTERS ARE REMOVED
	preochar <- rep(0,ncol(X))	# the ordering is placed on matrix, a matrix is made where all characters have 0 
	preochar[ochar] <- 1	#characters not ordered are given weight 1 within this matrix of 0
	if (length(bad)!=0) X <- X[,-c(bad)]	#the bad columns are now removed
	D <- diag(0,nrow(X))							# produces matrix of same dimensions as morphology matrix filled with zeros
	if (length(bad)!=0) preochar <- preochar[-c(bad)]	# the pre-ordered characters have bad columns removed
	ochar <- which(preochar==1)							#list of charaters that are ordered (ie have value 1)
	if (length(bad)!=0) w <- w[-c(bad)]
	no <- c(1:dim(X)[2])[-c(ochar)]							#make vector from 1 to number of dimensions in morphology matrix (no.of collumns) - those characters that are unordered (in no=NO ORDER)
	w[no] <- 1			#Sets all w range values of characters in each collumn to 1 if they are unordered.
	pw <- combn(nrow(X),2)			#pairwise combinations of all values of morphology matrix
	kill <- c()						# kill are taxa that are removed because they are too incomplete to calculate distance between them and any other taxon in the distance matrix

### The actual distance calculation is here (below)
	for (i in 1:ncol(pw)) {
		d <- (abs(X[pw[1,i],] - X[pw[2,i],]))	#for each pairwise combination of taxa: taxa A character states - taxa B character states (only need absolute value ie dont care about +ve or -ve signs)
		d[no][which(d[no]>1)] <- 1	
		d <- d*(1/w)	# the absolute values of taxa A states - taxa B states are all then multiplied by: 1 divided by max range for that state
		d2 <- sum(na.omit(d))/(length(na.omit(d)))	# omiting characters that cannot have their distance calculated	give the distance value d2	
		d <- c(na.omit(d),rep(d2,length(d)-length(na.omit(d))))
		# add the distance value of d2 to the end of d, ommiting NA's
		d <- sqrt(sum(d^2))
		if (is.nan(d)==TRUE) kill <- c(kill,c(pw[,i]))
		D[pw[1,i],pw[2,i]] <- d
		}
		DD <- D+t(D)
### The actual distance calculation is here (above)

	if (length(kill)!=0) {
		killmat <- cbind(unique(kill),rep(0,length(unique(kill))))		#make an table of the taxa that cannot have their distance calculated
		for (i in 1:dim(killmat)[1]) killmat[i,2] <- length(which(kill==killmat[i,1]))
		order <- rank(killmat[,2],ties.method="f")		#rank these taxa in the order of times they cannot have thier distances calculated with other taxa
		killmat[order,2] <- killmat[,2]
		killmat[order,1] <- killmat[,1]
		i <- dim(killmat)[1]
		while (length(which(is.nan(D)==T))!=0) {
			D[c(killmat[i,1]),] <- 999
			D[,c(killmat[i,1])] <- 999
			i <- i-1
		}
		DD <- D+t(D)
		kill <- which(D[,1]==999)
		if (length(kill)!=0) {
			print("Minimum number of taxa removed to produce distance matrix:",quote=FALSE)
			print(as.numeric(kill))
			D <- D[-c(kill),-c(kill)]
		}
	}
	D <- D+t(D)
	list(DD=DD,D=D,kill=kill,remtax=which(irr==0))
	}
	
###	Ordination of distance scores via Principal coordinate analysis
PCO <- function(D,correction=c("Lingoes","Cailliez")) {
	n <- nrow(D)
	correc <- match.arg(correction)
	gcm <- function(D) {	# Gower-centered matrix
		A <- -.5*D^2
		ai <- matrix(apply(A,1,mean),n,n)
		am <- matrix(mean(A),n,n)
		d1 <- A-ai-t(ai)+am	
		}
	eps <- sqrt(.Machine$double.eps)
	Ev <- eigen(gcm(D))$values
	Ev[c(which(abs(Ev)< eps))] <- 0
	nEv <- min(Ev)	# smallest eigenval.
	
	if (nEv<0) {
		if (correc=="Lingoes") {
			D <- sqrt((D^2)+2*abs(nEv)) - sqrt(2*diag(abs(nEv),n))
			print("Lingoes applied")
			}
	
		if (correc=="Cailliez") {
			zero <- diag(0,n) ; id <- diag(1,n)
			del1 <- gcm(D)
				AA <- -.5*D
				ai <- matrix(apply(AA,1,mean),n,n)
				am <- matrix(mean(AA),n,n)
			del2 <- AA-ai-t(ai)+am
			mat <- cbind(rbind(zero,-id),rbind(2*del1,-4*del2))
			c2 <- max(Re(eigen(mat)$values))
			D <- D+c2-diag(c2,n)
			print("Cailliez applied")
			}
		}
		
	d1 <- gcm(D)
	S <- eigen(d1)
	S$values[abs(S$values)<eps] <- 0
	list(pco=S$vectors%*%diag(sqrt(S$values)),Eval=S$values,PCAEval=S$values/(nrow(D)-1))
	}

### Disparity functions for Sum of vairances and Sum of Ranges, use all PCO axes
SoV <- function(X) sum(diag(cov(X)))
SoR <- function(X) sum(apply(X,2,range)[2,]-apply(X,2,range)[1,])

### Calculate Taxonomic diversity at each time bin
Taxonomic.Diversity <- apply(PA,2,sum)
Error.TD <- sqrt(Taxonomic.Diversity)

### Disparity analysis using SoV and SoR
DtT <- function(MS,PA,boo) {
	skewness.b.sple <- c()			# the skewness of the distribution of the bootstrap resample at everybin
	rar <- min(apply(PA,2,sum))			# the diversity value the Sum of Ranges is rarefied to
	if (rar==1) rar <- 2
	Variance <- c() ; Range <- c()
	
	MeanBootCentroid<- matrix(NA,nrow=ncol(PA),ncol=ncol(MS))		###### ADDED BY MARTIN HUGHES 13/04/2015
	SDBootCentroid<- c()		###### ADDED BY MARTIN HUGHES 13/04/2015
	
	for (tps in 1:ncol(PA)) {			# a loop for computing disparity for each time interval
		occ <- which(PA[,tps]==1)
		if (length(occ)==1) {occ <- c(occ,occ) ; rar <- 2 }
		vari <- c() ; ran <- c()
		centri<- matrix(NA,ncol=ncol(MS),nrow=boo)		###### ADDED BY MARTIN HUGHES 13/04/2015
		for (i in 1:boo) {
			b.sple <- sample(occ,length(occ),replace=T)
			r.sple <- sample(occ,rar,replace=F)
			vari <- c(vari,SoV(MS[b.sple,]))		# records disparity values of boostrap replicates
			ran <- c(ran,SoR(MS[r.sple,]))		# same but for range-based estimates
			centri[i,]<- apply(MS[b.sple,],2,mean)	###### ADDED BY MARTIN HUGHES 13/04/2015
			}
		outvar <- which(vari==0) ; outran <- which(ran==0)
		ifelse(length(outvar!=0),Variance <- cbind(Variance,rbind(mean(vari[-c(outvar)]),sd(vari[-c(outvar)]))),Variance <- cbind(Variance,rbind(mean(vari),sd(vari)))) # Bootstrap estimates and errors
		ifelse(length(outran!=0),Range <- cbind(Range,rbind(mean(ran[-c(outran)]),sd(ran[-c(outran)]))),Range <- cbind(Range,rbind(mean(ran),sd(ran))))       # id.
		skewness.b.sple <- c(skewness.b.sple,skewness(vari))
		
		MeanBootCentroid[tps,] <- apply(centri,2,mean)	###### ADDED BY MARTIN HUGHES 13/04/2015
		SDBootCentroid<- c(SDBootCentroid,apply(centri,2,sd))	###### ADDED BY MARTIN HUGHES 13/04/2015
		
		}
	list(Variance=Variance, Range=Range,skewness.b.sple=skewness.b.sple,MeanBootCentroid=MeanBootCentroid,SDBootCentroid=SDBootCentroid)
	}


###
Partial.Disparity <- function(X,gp) {
	Xc <- scale(X,scale=F) ; PD <- c()
	for (i in levels(gp)) {
		PD <- c(PD,sum(diag(Xc%*%t(Xc)/(nrow(Xc)-1))[which(gp==i)]))
		}	
	PD
	}


### Rarefaction function for SoR , boo is  number of bootstrap replicates to get standard error
### 2 versions, the RarCrv is for regular use, but if more than 50 taxa in a time bin it will use
### the RarCrv_big code insted

# rarefaction curve for time bins less than 50 taxa, it will use all of the taxa in that bin
RarCrv <- function(X,boo) {
	rCrv <- c()
	for (size in 2:nrow(X)) { 
		rRge <- c()
		for (i in 1:boo) {
			sple <- sample(1:nrow(X),size,replace=F)
			rRge <- c(rRge,SoR(X[sple,]))
			}
		rCrv <- cbind(rCrv,rbind(mean(rRge),sd(rRge)))
		}
	list(RarCrv=rCrv)
	}
	
# rarefaction for more than 50 taxa, it will use 50 taxa from the number in that timebin (equally spread not random)	
RarCrv_big <- function(X,boo) {	
	rCrv <- c()
	for (size in seq(2,nrow(X),ceiling((nrow(X)-2)/50))) { # increments can be adapted to sample size (rounded up)
		rRge <- c()
		for (i in 1:boo) {											# run the following for however many times you wish to replicate
			sple <- sample(1:nrow(X),size,replace=F)				# sample from 1 to the number of rows in X (number of taxa in morphospace that are in that time bin) without replacement
			rRge <- c(rRge,SoR(X[sple,]))							# for that sample of taxa, calculate the Sum of ranges
			}
		rCrv <- cbind(rCrv,rbind(mean(rRge),sd(rRge)))				# produce the matrix for that SoR (mean) with the standard deviation
		}
	list(RarCrv=rCrv,skipseq=seq(2,nrow(X),ceiling((nrow(X)-2)/50)))	# add the number of the time interval (based on number of in time intervals) 
	}
	
### function to produce the timebins you want
#TimeBinConverter<-function (PA,strat.level=c("Stage","SubEpoch","Epoch","Period")) {
TimeBinConverter<-function (PA,strat.level=c("Stage","Epoch","Period")) {

 
 #pa<-pa[,-c(1)]		# this removes the collumn for condensing so that it does not interfere with the time condensing
 
 empty<-matrix(0,nrow=nrow(PA), ncol=1) # make an extra column with only 0's in it (for those periods that include only an extra stage)
 
 colnames(empty)<-c("empty") # name it
 
 PA.mod<-cbind(PA,empty) # modify pa, now has empty column with 0's : important for final part of script
 
 row.numbers<-match(colnames(PA),timetable[,3]) # find the row numbers of the timetable file that match the names the stratigraphy file
 
 tt<-timetable[row.numbers,] # output a reduced timetable only including the stages from stratigraphy file
  
 tt[,6]<-as.character(tt[,6]) # turn the vaules in tt at columns 6 (periods) into class "characters" insted of factors
 tt[,5]<-as.character(tt[,5]) # turn the vaules in tt at columns 5 (epochs) into class "characters" insted of factors
#tt[,4]<-as.character(tt[,4]) # turn the vaules in tt at columns 4 (sub epochs) into class "characters" insted of factors
 tt[,3]<-as.character(tt[,3]) # turn the vaules in tt at columns 3 (stages) into class "characters" insted of factors
 
 stlvl<-match.arg(strat.level) # use this function to be able to specify things in function and therefore change the outcome
 if (stlvl=="Stage") colSumz<-split(tt[,3], tt[,3]) # split the stages into groups based on the Stages
# if (stlvl=="SubEpoch") colSumz<-split(tt[,3], tt[,4]) # split the stages into groups based on the SubEpochs
 if (stlvl=="Epoch") colSumz<-split(tt[,3], tt[,5]) # split the stages into groups based on the Epochs
 if (stlvl=="Period") colSumz<-split(tt[,3], tt[,6]) # split the stages into groups based on the Periods
 
 result<-matrix(0,nrow=nrow(PA), ncol=length(colSumz)) # make empty matrix that is the number of rows in PA, and the number of columns determined by Colsumz
 colnames(result)<-names(colSumz) # name the columns with the periods from colSumz
 rownames(result)<-rownames(PA) # name the rows with the taxa from PA
 
 for (i in 1:length(colSumz)) # loop to see if any of the periods identified have a length 1 ie 1 stage
 {
 colSumz[[i]]<-append(colSumz[[i]],"empty") # add the extra stage called empty to each
  }
 
 for (i in names(colSumz)) #eg the periods in this case jurassic, cretaceous, palaegene (but in alphabetical order)
 {	
	result[,i]<-rowSums(PA.mod[,colSumz[[i]]]) 	# sum the rows by looking at values in PA.mod if empty column needed, knowning that stages in PA.mod correpond to the periods in colSumz
 } 
  
 result <- replace(result, result > 1, 1) # changes all values that are greater than 1 to 1
 
 if (stlvl=="Stage") strat.name<-unique(tt[,3]) # order the data correctly (based on time scale,  not alphabetically)
# if (stlvl=="SubEpoch") strat.name<-unique(tt[,4])
 if (stlvl=="Epoch") strat.name<-unique(tt[,5])
 if (stlvl=="Period") strat.name<-unique(tt[,6])
 
 result<-result[,strat.name] 
}
 
 
 
### The following are functions for calculating the Centre of gravities	
# CG.taxo computes the center of gravity of taxonomic curves. Arguments are PA (the
# condensed presence-absence matrix prior to culling) and timescale (a vector of length l where
# l is the number of stages. each entry is the temporal length of the stage ie in myrs).
# for instance, timescale= [0.9,1,1.1,3, ... , 1,1.4]

CG.taxo <- function(PA,timescale) { # PA is the original, unculled PA
	ntb <- ncol(PA) ; n <- nrow(PA)
	diversity <- apply(PA,2,sum)

	accumtime <- c()														# accumulated time! simply sums the interval length from the first to the last and thus goes from 1 to ntb (ntb=number of time bins)
	for (i in 1:ntb) accumtime <- c(accumtime,sum(timescale[1:i]))
	midtime <- c()														# the point in the middle of each interval
	for (i in 2:ntb) midtime <- c(midtime,accumtime[i-1]+(accumtime[i]-accumtime[i-1])/2)
	midtime <- c(accumtime[1]/2,midtime)

	CGi <- sum(midtime)/ntb
	CGt <- sum(midtime*diversity)/sum(diversity)

	reCGt <- c()
	for ( boo in 1:1000) {
		sple <- sample(1:n,n,replace=T)
		rePA <- PA[sple,]
		rediversity <- apply(rePA,2,sum)
		rediversity[which(rediversity==0)] <- NA
		#reCGt <- c(reCGt, sum(midtime*rediversity)/sum(rediversity))
		reCGt <- c(reCGt, sum(midtime*rediversity,na.rm=T)/sum(rediversity,na.rm=T))
	}
	nCGt <- mean(reCGt)
	erCGt <- sd(reCGt)
		Foote.tst <- length(which(reCGt>CGi))
	if (Foote.tst > 500) Foote.tst <- 1000-Foote.tst
	pval <- (Foote.tst*2)/1000
	sc <- accumtime[ntb]
list(CGt=CGt, nCGt=nCGt,erCGt=erCGt,CGi=CGi, scCGt=CGt/sc, scnCGt=nCGt/sc,scerCGt=erCGt/sc,scCGi=CGi/sc,pval=pval)
}


###
# disparityCRVfun computes disparity curve (for variance only). the function is used by
# the next one - it gets disparity values again but insted of taking each interval in turn
# (and therefore independant from one another) it will sample from the whole of PA, make new PA and morphology matrix and 
# calculate disparity that way

BIG <- c()
Max.morpho <- function(PA,X) {
        # PA is the culled PA
        # X is the PCO ordination
BIG=c()
        ntb <- ncol(PA) ; n <- nrow(X)
        reDiffmx <- c()
        reDiffmn <- c()
		reDiffmx.late <- c()
		reDiffmn.late <- c()
        obs.mx <- max(disparityCRVfun(X,PA))
        obs.mn <- min(disparityCRVfun(X,PA))
        CV <- c() ; cv <- c()
        for ( boo in 1:1000) {
                sple <- sample(1:n,n,replace=T)
                rePA <- PA[sple,]
                redisparity <- disparityCRVfun(X[sple,],rePA)
                if (length(which(redisparity==0))==0) cv <- c(cv,(1+1/(4*length(redisparity)))*sd(redisparity)/mean(redisparity))
                earlyvalue <- max(redisparity[1:2])
				latevalue <- max(redisparity[(ntb-1):ntb])
                reDiffmx <- c(reDiffmx, obs.mx-earlyvalue)
                reDiffmn <- c(reDiffmn, obs.mn+earlyvalue)
                reDiffmx.late <- c(reDiffmx.late, obs.mx-latevalue)
                reDiffmn.late <- c(reDiffmn.late, obs.mn+latevalue)
                				
				BIG <- rbind(BIG,redisparity)
                }
		CV <- mean(cv)
list(reDiffmx=reDiffmx,reDiffmn=reDiffmn,CV=CV,BIG=BIG,reDiffmx.late=reDiffmx.late,reDiffmn.late=reDiffmn.late)
}


THE.INTEGRATED.GERBER.WILLS.TEST <- function(X,PA) {
	ntb <- ncol(PA)
	earlyguys <- unique(c(which(PA[,1]==1),which(PA[,2]==1)))
	lateguys <- unique(c(which(PA[,(ntb-1)]==1),which(PA[,ntb]==1)))
	if(length(earlyguys)==1) earlyguys<-c(earlyguys,earlyguys)
	if(length(lateguys)==1) lateguys<-c(lateguys,lateguys)
	earlynull <- c()
	latenull <- c()
	earlyobs <- c()
	lateobs <- c()
	for (i in 1:1000) {
		sple1 <- sample(1:nrow(X),length(earlyguys),replace=T)
		sple2 <- sample(1:nrow(X),length(lateguys),replace=T)
		earlynull <- c(earlynull,sum(diag(cov(X[sple1,]))))
		latenull <- c(latenull,sum(diag(cov(X[sple2,]))))
		earlyobs <- c(earlyobs,sum(diag(cov(X[sample(earlyguys,length(earlyguys),replace=T),]))))
		lateobs <- c(lateobs,sum(diag(cov(X[sample(lateguys,length(lateguys),replace=T),]))))
	}
earlylimit <- max(earlynull)
latelimit <- max(latenull)
#max values
early.test.max<-length(which(earlyobs>earlynull))
	if (early.test.max > 500) early.test.max <- 1000-early.test.max
	early.max.pval <- (early.test.max*2)/1000
late.test.max<-length(which(lateobs>latenull))
	if (late.test.max > 500) late.test.max <- 1000-late.test.max
	late.max.pval <- (late.test.max*2)/1000
ifelse(earlylimit>quantile(earlyobs,.95),comment.early.max<-"Early disparity does not reach maximum theorerical value (given sample size)",comment.early.max<- "-")
ifelse (latelimit>quantile(lateobs,.95), comment.late.max<-"Late disparity does not reach maximum theorerical value (given sample size)",comment.late.max<- "-")
list(earlyobs=earlyobs, lateobs=lateobs , earlynull=earlynull , latenull=latenull , earlylimit=earlylimit , latelimit=latelimit ,early.max.pval=early.max.pval , late.max.pval=late.max.pval , comment.early.max=comment.early.max , comment.late.max=comment.late.max)
}








disparityCRVfun <- function(x,pa) {
	Dcrv <- c()
	for (t in 1:ncol(pa)) {
	if ( length(which(pa[,t]==1))==1) Dcrv <- c(Dcrv,0)
	if ( length(which(pa[,t]==1))!=1) Dcrv <- c(Dcrv,sum(apply(x[which(pa[,t]==1),],2,var)))
	}
	Dcrv[which(is.na(Dcrv))] <- 0
	Dcrv
}

###
# CG.morpho computes the center of gravity of morphological curves. Arguments are PA (the
# condensed and culled PA), X (the PCO coordinates in your case) and timescale (same as above).

CG.morpho <- function(PA,X,timescale) {
	# PA is the culled PA
	# X is the PCO ordination

	ntb <- ncol(PA) ; n <- nrow(X)

	accumtime <- c()																					# accumulated time
	for (i in 1:ntb) accumtime <- c(accumtime,sum(timescale[1:i]))
	midtime <- c()
	for (i in 2:ntb) midtime <- c(midtime,accumtime[i-1]+(accumtime[i]-accumtime[i-1])/2)
	midtime <- c(accumtime[1]/2,midtime)

	CGi <- sum(midtime)/ntb
	disparity <- disparityCRVfun(X,PA)
	CGm <- sum(midtime*disparity)/sum(disparity)

	reCGm <- c()											#CoG resampled via resampled PA for morphology
	for ( boo in 1:1000) {
		sple <- sample(1:n,n,replace=T)
		rePA <- PA[sple,]									#rePA - resampled PA (whole matrix (FOOTE) not each stage seperatly (OUR way)
		redisparity <- disparityCRVfun(X[sple,],rePA)		#redisparity - the disparity calculated from this resampled PA matrix
		redisparity[which(redisparity==0)] <- NA
		reCGm <- c(reCGm, sum(midtime*redisparity,na.rm=T)/sum(redisparity,na.rm=T))
	}
	nCGm <- mean(reCGm)
	erCGm <- sd(reCGm)
	Foote.tst <- length(which(reCGm>CGi))
	if (Foote.tst > 500) Foote.tst <- 1000-Foote.tst
	pval <- (Foote.tst*2)/1000
	sc <- accumtime[ntb]
list(CGm=CGm, nCGm=nCGm,erCGm=erCGm,CGi=CGi,scCGm=CGm/sc, scnCGm=nCGm/sc,scerCGm=erCGm/sc,scCGi=CGi/sc,pval=pval)
}


# Outputs are (the list might need to be extended):
# - CGm		or		CGt	: observed CG		- the Centre of gravity of the data
# - nCGm		or		nCGt : average CG based on bootstrap (should be close to CGm)
# - erCGm		or		erCGt : error bars for nCGm (based on Foote 1993 methods)
# - CGi : CG inherent in timescale - ie what you compare observed to to understand if centre of gravity IS top/bottom/symmetrically heavy. Note that for a given dataset the CGi obtained with CG.taxo and CG.morpho should be exactly the same.
# - scCGm		or		scCGt : scaled CGm to span the [0,1] range
# - scnCGm		or		scnCGt : scaled nCGm		or		nCGt
# - scerCGm		or		scerCGt : scaled erCGm		or		erCGt
# - scCGi: scaled CGi (same logic for the output of CG.taxo)


	