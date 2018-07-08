

setwd("/Users/jessieabbate/desktop/ebola manuscript/")

thefile <- read.table(file="ebola_cases_SL.txt",header=T)



infectedExp <- function(x){ exp(0.063*(x-169+41))}
#infectedR0 <- function(x){ exp(0.088*(x-169+41))}


####################
## Exp model (deterministic)
####################

TotInfectedExp <- function(x){exp(0.063*(x-169+41))/0.063}
TotInfectedExpmin <- function(x){exp(0.011*(x-169+41))/0.011}
TotInfectedExpmax <- function(x){exp(0.127*(x-169+41))/0.127}

temp <- max(thefile[,1])-128

thefileExp <- vector(length=temp)
thefileExpmin <- vector(length=temp)
thefileExpmax <- vector(length=temp)

	thefileExp[1] <- 1
	thefileExpmin[1] <- 1
	thefileExpmax[1] <- 1




for(i in 129:max(thefile[,1]))
{
	thefileExp[i-127] <- TotInfectedExp(i) - TotInfectedExp(i-1)
	thefileExpmin[i-127] <- TotInfectedExpmin(i) - TotInfectedExpmin(i-1)	
	thefileExpmax[i-127] <- TotInfectedExpmax(i) - TotInfectedExpmax(i-1)		
}

# #Another way to do it
# for(i in 129:max(thefile[,1]))
# {
	# thefileExp[i-127] <- thefileExp[i-127-1] + infectedExp(i)-infectedExp(i-1)
# }


plot(1,xlim=c(110,260),ylim=c(0,10),xlab="Day of the year",ylab="ln(Cumulative number of cases)",pch=16,cex=0.5)
polygon(c(145,145,169,169),c(-1,13,13,-1),col = "grey90", border = NA)
points(log(thefile[,2])~thefile[,1],pch=16,cex=0.5)


lines(log(thefileExp)~c(128:max(thefile[,1])),col=1)
lines(log(thefileExpmin)~c(128:max(thefile[,1])),col=1,lty=2)
lines(log(thefileExpmax)~c(128:max(thefile[,1])),col=1,lty=2)





####################
## BD model (stochastic)
####################

runs <- 100

#create a list for the results
results <- vector("list",length=runs)

#loop on the number of runs
for(i in c(1:runs))
{


#initialise values
beta <- 1.259*0.338/33056
gamma <-  0.338
tot <- 1
S <- 33055
I <- 1
t <- 0

#initialise results vector
pop <- c(t,S,I,tot)

#loop until max time reached or extinction
while((t<=164)&(I>0))
{
	#when is the next event happenning (drawn from an exponential distribution)
	deltat <- rexp(1,(beta*S*I+gamma*I))

	#determine which of the two events happens
	if(runif(1)<=((beta*S*I)/(beta*S*I+gamma*I)))
	{
		S=S-1;
		I=I+1;
		tot=tot+1;
	}
	else{I=I-1}

	#update time and store results
	t <- t+deltat
	pop <- rbind(pop,c(t,S,I,tot))

}

#store the results of this run in the list
results[[i]] <- pop

}

plot(1,xlim=c(110,300),ylim=c(0,10),xlab="Day of the year",ylab="ln(Total number of new cases)",pch=16)
polygon(c(145,145,169,169),c(-1,13,13,-1),col = "grey90", border = NA)
points(log(thefile[,2])~thefile[,1],pch=16,cex=0.5)
for(i in 1:runs)
{
	lines(results[[i]][,1]+126,log(results[[i]][,4]),col="light blue",lwd=0.05)	
}

results1<-results



########################
# Deterministic SIR
########################

library(deSolve)


dynamics = function(beta=1.66,gamma=1/2.2,tot=1,S0=762,I0=1,tmax=20,by=0.01)
{
	require(deSolve)
	parameters = c(beta=beta,gamma=gamma)
	state = c(S=S0,I=I0,tot=I0)
	model = function(t,state,parameters)
		with(as.list(c(state,parameters)),
		{
			dS = -beta*S*I
			dI = beta*S*I - gamma*I
			dtot=beta*S*I
			list(c(dS,dI,dtot))
		})
	times = seq(0,tmax,by)
	return(ode(y=state,times=times,func=model,parms=parameters))
}

res <- dynamics(beta= (1.259*0.338/33056), gamma=0.338,S0= 33056,I0=1,tmax=164,tot=1)
resmin <- dynamics(beta= (1.0358*0.338/1103), gamma=0.338,S0= 1103,I0=1,tmax=164,tot=1)
resmax <- dynamics(beta= (1.5357*0.338/337600), gamma=0.338,S0= 337600,I0=1,tmax=164,tot=1)

res[,1] <- res[,1]+126 
resmin[,1] <- resmin[,1]+126 
resmax[,1] <- resmax[,1]+126 
head(res)

#plot(1,xlim=c(110,300),ylim=c(0,10),xlab="Day of the year",ylab="ln(Number of new cases)",pch="•")
#polygon(c(145,145,169,169),c(-1,13,13,-1),col = "grey90", border = NA)
#points(log(thefile[,2])~thefile[,1],pch="•")
lines(res[,1],log(res[,4]),col=4)
lines(resmin[,1],log(resmin[,4]),col=4,lty=4)
lines(resmax[,1],log(resmax[,4]),col=4,lty=4)




# ######################################
# # The final (pretty) stochastic figure
# ######################################

#plot(1,xlim=c(110,300),ylim=c(0,10),xlab="Day of the year",ylab="ln(Total number of new cases)",pch=16)
plot(1,xlim=c(110,300),ylim=c(0,10),xlab="Day of the year",ylab="ln(Cumulative number of cases)",pch=16)
polygon(c(145,145,169,169),c(-1,13,13,-1),col = "grey90", border = NA)
for(i in 1:runs)
{
  lines(results[[i]][,1]+126,log(results[[i]][,4]),col="green",lwd=0.05)	
}
lines(res[,1],log(res[,4]),col=4)
lines(resmin[,1],log(resmin[,4]),col=4,lty=4)
lines(resmax[,1],log(resmax[,4]),col=4,lty=4)
points(log(thefile[,2])~thefile[,1],pch=16,cex=0.5)



###################################################################################
# ################################
# # Gillespie using the package
# ################################


# library(GillespieSSA)

# #parms = c(beta= (1.256*0.337), gamma=0.346) 
# parms = c(beta= (1.259* 0.338/33056), gamma= 0.338) 

# # R0=beta/gamma => beta = 1.35*0.326

# #initcond = c(I=1) 
# initcond = c(S=33056, I=1, R=0) 

# #terms <- c("beta*I", "gamma*I")
# terms <- c("beta*I*S", "gamma*I")

# #mat <- matrix(c(+1,-1), ncol=2, nrow=1,byrow=T)
# mat <- matrix(c(-1,+1,0,
                          # 0,-1,+1), ncol=3, nrow=2,byrow=T)
# matT <- t(mat)                          

# maxiterations <- 100

# plot(1,xlim=c(110,300),ylim=c(0,10),xlab="Day of the year",ylab="ln(Number of new cases)",pch="•")
# polygon(c(145,145,169,169),c(-1,13,13,-1),col = "grey90", border = NA)
# points(log(thefile[,2])~thefile[,1],pch="•")
# for(i in 1:maxiterations)
# {
# #	temp <- ssa(initcond, terms, mat, parms, tf=124, maxWallTime=10)$data
	# temp <- ssa(initcond, terms, matT, parms, tf=164, maxWallTime=10)$data
	# popsizes <- as.numeric(as.character(temp[,3]))
	# totCases <- popsizes
	# for(i in 2:length(popsizes))
	# {
		# diff <- popsizes[i]-popsizes[i-1]
		# if(diff<0){diff=0}
		# totCases[i] <- totCases[i-1]+diff
	# }
	# thetimes <- as.numeric(as.character(temp[,1]))+126
	# lines(thetimes,log(totCases),col=4,lwd=0.05)
# }


# #lines(log(thefileR0)~c(126:max(thefile[,1])),col=4)
# #lines(log(thefileR0min)~c(126:max(thefile[,1])),col=4,lty=2)
# #lines(log(thefileR0max)~c(126:max(thefile[,1])),col=4,lty=2)
