#created for Perry Hooker for possible collaboration;
#2-D simulation of K function with envelopes
# R code

rm(list = ls())		#clears vectors
begin <- Sys.time()	#begins time clock

### Parameters ###
nsim = 100			#number of simulations for simulated Ks
n <- 20			#number of release site
d <- c(1,1)			#dimensions of study region
center <- c(.5,.5)	#center of cell (cell is a single point)
#band <- c(.1,.2)		#clustered band radius (used for a different contrived data set)
lambda <- 3			#exponential distribution parameter
h <- 1:40			#distances away from cell boundary

# Observed Data Construction (Contrived)
#=====================================================================
x.CO <- center[1]					#x.CO and y.CO are locations of cell surface
y.CO <- center[2]
#r <-     runif(n,band[1],band[2])		#different rs for different contrived data set
#r <-     rexp(n,lambda)
#theta <- runif(n,0,2*pi)
#x.RP.all <- r*cos(theta)+center[1]
#y.RP.all <- r*sin(theta)+center[2]
x.RP.all <- runif(n,0,1)			#this randomly places Release Points (RP) around cell
y.RP.all <- runif(n,0,1)			#thus pattern about cell is random or not clustered about cell
take <- x.RP.all >= 0  &  x.RP.all <= d[1]  &  y.RP.all >= 0 & y.RP.all <= d[2]
x.RP <- x.RP.all[take]
y.RP <- y.RP.all[take]				#needed when contrived data set allows RPs out of study region
n <- length(x.RP)					#number of release points
m <- length(x.CO)					#number of cell surface points
R <- d[1]*d[2]					#area of study region
hd <- (4*max(h))					#see vectors maxh and h.2

par(mfrow=c(1,1))					#graphic of study region
plot(x.RP,y.RP,xlim=c(0,d[1]),ylim=c(0,d[2]))
points(x.CO,y.CO,pch=16,col="red")


# Construction of Observed K and w (w is an edge correction)
#=======================================================================================
Dist <- 
((matrix(rep(x.RP),n,m,byrow=FALSE)-matrix(rep(x.CO),n,m,byrow =TRUE))^2 +
 (matrix(rep(y.RP),n,m,byrow=FALSE)-matrix(rep(y.CO),n,m,byrow =TRUE))^2)^.5		#calculates all distances from each pair of RPs and cell surface points			
Dist.min <- apply(Dist,1,min)		#takes the min distance for each RP to cell surface (applies when m>1)														
rm(Dist)					#removes Dist matrix for memory reasons											
maxh <- max(Dist.min)			 
h.2 <- 1:ceiling(maxh*hd)           #hd, maxh, and h.2 objects are use so that a series of distances with equal increments
						#needed since loop below needs integers

n.w <- 10000							#beginning of loop that estimates area of bands/shells about cell surface
w=NA									#this is needed to create w the edge correction of the K function
for(i in h.2){
x.w <- runif(n.w,center[1]-(i/hd),center[1]+(i/hd))	#drops random points about cell 
y.w <- runif(n.w,center[2]-(i/hd),center[2]+(i/hd))	#note that i is being divided to create appropriate distances
Dist <- 
((matrix(rep(x.w),n.w,m,byrow=FALSE)-matrix(rep(x.CO),n.w,m,byrow =TRUE))^2 +			#Distance matrix
 (matrix(rep(y.w),n.w,m,byrow=FALSE)-matrix(rep(y.CO),n.w,m,byrow =TRUE))^2)^.5
Dist.min.2 <- apply(Dist,1,min)  						#takes the min distance for each RP to cell surface (applies when m>1)
rm(Dist)
take <- x.w>0  &  x.w<d[1]  &  y.w>0 & y.w<d[2]				#take=1 if point lies within study region; =0 if not
num.star <- sum((take == 1) & (Dist.min.2 < i/hd) & (Dist.min.2 > (i-1)/hd))  #number of points in shell and in study region
num.tri  <- sum(              (Dist.min.2 < i/hd) & (Dist.min.2 > (i-1)/hd))	#number of points in shell in and out of study region
w[i] <- num.star/num.tri							#find proportion of band/shell that is within study region
}
w		#edge correction; proportion of ith band/shell within study region

C=NA;for(i in h.2){ C[i] = sum(Dist.min<(i/hd)) }	#i loop; back to observed data (contrived); counts num of RPs within i/hd distance from cell
IN     <- (c(C,0)-c(0,C))[1:length(C)]			#counts num of RPs within a shell/band
w.plus <- cumsum(((1-w)/w)*IN)				#if band lies outside study region (w<1) extra RPs must be added; w.plus=0 if w=1
n.e <- (C+w.plus)[max(h.2)]					#adjusted n for K function (function output scalar)
K.obs <- (R/n.e)*(C+w.plus)					#K function dependent on distance


# Construction of Simulated K (envelopes)
#=======================================================================================
K.sim = matrix(rep(NA),nsim,length(h.2))
for(j in 1:nsim){
x.sim <- runif(n,0,1)
y.sim <- runif(n,0,1)
Dist <- 
((matrix(rep(x.sim),n,m,byrow=FALSE)-matrix(rep(x.CO),n,m,byrow =TRUE))^2 +
 (matrix(rep(y.sim),n,m,byrow=FALSE)-matrix(rep(y.CO),n,m,byrow =TRUE))^2)^.5
Dist.min <- apply(Dist,1,min)  
rm(Dist)
C=NA;for(i in h.2){ C[i] = sum(Dist.min<(i/hd)) } #i loop
IN        <- (c(C,0)-c(0,C))[1:length(C)]
w.plus    <- cumsum(((1-w)/w)*IN)
n.e <- (C+w.plus)[max(h.2)]
K.sim[j,] <- (R/n.e)*(C+w.plus)  		
}
#loop drops points down many times and for each time calculates K function using same method as above including edge effect
#thus areas and proportions from above are needed

Q = matrix(NA,length(h.2),3)
for (i in h.2){ Q[i,] <- as.numeric( quantile(K.sim[,i],c(.025,.5,.975)) ) }
envelope.Low   <-  Q[,1]
envelope.Mid   <-  Q[,2]
envelope.High  <-  Q[,3]
#creates envelopes based on the quantiles of the simulations

#Graphic
#=======================================================================================	
L.obs 			 <- K.obs 		#- envelope.Mid
envelope.High.L 	<- envelope.High 	#- envelope.Mid
envelope.Low.L 	<- envelope.Low 	#- envelope.Mid
ylim = c(min(L.obs,envelope.Low.L,envelope.High.L),
		 max(L.obs,envelope.Low.L,envelope.High.L))

par(mfrow=c(1,1))
plot(h.2,L.obs,ylim=ylim,,xlim=c(0,150),type="n", xaxt="n",
		main='Modified L with Envelopes',
		xlab="Distance (h)",ylab="Modified L",
		cex=1,cex.lab=1,cex.axis=1)
points(h.2,L.obs, col=1, type="l",lwd=2)
points(h.2,envelope.High.L, col="black", type="l",lty=2)
points(h.2,envelope.Low.L, col="black", type="l",lty=2)


axis(1,at=seq(0,max(h),8),labels=seq(0,.5,.1),cex.axis=1.2)
#legend('topleft',c("Modified K","Envelopes"),lty=c(1,2),
#	,lwd=c(3,3),col=c(1,2),cex=1.2)
#abline(v=max(h))



end <- Sys.time()
end-begin
#!


