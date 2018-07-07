###################################################################
#   Summary statistics for Monte Carlo simulations                #
#                                                                 #
#  Accompanying:                                                  #
#                                                                 #
#  Potential impact of sexual transmission of Ebola virus         #
#  on the epidemic in West Africa                                 #
#  JL Abbate, C-L Murall, H Richner, C Althaus                    #
#                                                                 #
#   December 11, 2015                                             #        
###################################################################

#required libraries
library(MASS)
library(plotrix)

#####################################################
#    without STI   
#####################################################

out    <- read.table(file="Stoch_noSTI_summaries_1000.csv",header=T)
colnames(out)


#####################################################
#    with STI   alpha = 3 months
#####################################################

outsti    <- read.table(file="Stoch_STI_3mo_summaries_1000.csv",header=T)
colnames(outsti)


#####################################################
#    with STI   alpha = 6 months
#####################################################

outsti6    <- read.table(file="Stoch_STI_6mo_summaries_1000.csv",header=T)
colnames(outsti6)

# summaries
#
# > colnames(out)
# [1] "test"               "R0"                 "TotalNumCases"     
# [4] "NoMoreEIs"        "DayLastCase"         "EpidemicPeakNew" 
# [7] "DayEpidemicPeakNew"
# 
# > colnames(outsti)
# [1] "test"               "R0"                 "eta"               
# [4] "p"                  "q"                  "alpha"                
# [7] "stiTotalNumCases"   "stiTotNumSTIcases"  "stiNoMoreEICs"
# [10]"stiDayLastCase"     "stiEpidemicPeakNew" "stiDayEpidemicPeakNew"

#####################################################
#                                                   #
#    Tail Histograms for Stochastic Runs            #
#                                                   #
#####################################################

# the data
thefile <- read.csv(file="data_WHO_SierraLeone_updated.csv",header=T)

weekIndexCase<-17  # April 23, 2014; week ends April 27, 2014

##################################
##  For No STI
##################################

#week of last case
WeekLastCase<-as.integer(as.character(out[,5]))/7
hist(WeekLastCase,breaks=20,xlim=c(1,300))

#pdf("weeklastcase_noSTI_hist.pdf", width = 4, height = 4)

lastweek<-375
xlabels <-thefile[weekIndexCase:(lastweek+weekIndexCase-1),11] #11= the date, 10 is the week number ; epiweek 17 is simulation week 1
xlength <- c(1:lastweek)
maxfreq <-c(0,150)
WeekLastCaseNoZero<-WeekLastCase
WeekLastCaseNoZero[which(out[,3]<50)]<-NA   # remove early die-outs
hist(WeekLastCaseNoZero,breaks=20,xlim=c(1,lastweek),ylim=maxfreq,xaxt="n",xlab=NULL,las=1,main="Week of Last Acute Case",col="grey") # xaxt = xaxis text
axis(1,at=seq(1,lastweek,by=25),labels=FALSE,par(tcl=-0.5))
axis(1,xlength,labels=FALSE,lty=0)
xl<-array(xlabels,lastweek)
xlabsby<-seq(1,lastweek,50)
text(x=xlabsby,par()$usr[3]-8,labels=xl[xlabsby],srt=60,adj=1,xpd=TRUE)

#dev.off()

##################################
##  For STI 3 months Conv. Period
##################################

#week of last case
WeekLastCase<-as.integer(as.character(outsti[,10]))/7
hist(WeekLastCase,breaks=20,xlim=c(1,500))

#pdf("weeklastcase_STI_3mo_hist.pdf", width = 4, height = 4)

lastweek<-375
xlabels <-thefile[weekIndexCase:(lastweek+weekIndexCase-1),11] #11= the date, 10 is the week number ; epiweek 17 is simulation week 1
xlength <- c(1:lastweek)
maxfreq <-c(0,150)
WeekLastCaseNoZero<-WeekLastCase
WeekLastCaseNoZero[which(outsti[,7]<50)]<-NA   # remove early die-outs
hist(WeekLastCaseNoZero,breaks=20,xlim=c(1,lastweek),ylim=maxfreq,xaxt="n",xlab=NULL,las=1,main="Week of Last Acute Case",col="grey") # xaxt = xaxis text
axis(1,at=seq(1,lastweek,by=25),labels=FALSE,par(tcl=-0.5))
axis(1,xlength,labels=FALSE,lty=0)
xl<-array(xlabels,lastweek)
xlabsby<-seq(1,lastweek,50)
text(x=xlabsby,par()$usr[3]-8,labels=xl[xlabsby],srt=60,adj=1,xpd=TRUE)

#dev.off()

##################################
##  For STI 6 months Conv. Period
##################################

#week of last case
WeekLastCase<-as.integer(as.character(outsti6[,10]))/7
hist(WeekLastCase,breaks=20,xlim=c(1,500))

#pdf("weeklastcase_STI_6mo_hist.pdf", width = 4, height = 4)

lastweek<-375
xlabels <-thefile[weekIndexCase:(lastweek+weekIndexCase-1),11] #11= the date, 10 is the week number ; epiweek 17 is simulation week 1
xlength <- c(1:lastweek)
maxfreq <-c(0,150)
WeekLastCaseNoZero<-WeekLastCase
WeekLastCaseNoZero[which(outsti6[,7]<50)]<-NA   # remove early die-outs
hist(WeekLastCaseNoZero,breaks=30,xlim=c(1,lastweek),ylim=maxfreq,xaxt="n",xlab=NULL,las=1,main="Week of Last Acute Case",col="grey") # xaxt = xaxis text
axis(1,at=seq(1,lastweek,by=25),labels=FALSE,par(tcl=-0.5))
axis(1,xlength,labels=FALSE,lty=0)
xl<-array(xlabels,lastweek)
xlabsby<-seq(1,lastweek,50)
text(x=xlabsby,par()$usr[3]-8,labels=xl[xlabsby],srt=60,adj=1,xpd=TRUE)

#dev.off()

#####################################################
#                                                   #
#    Summary Statistics for Stochastic Runs         #
#                                                   #
#####################################################

#####################################################
#    Cumulative (total) Number of Cases
#####################################################

tot.out<-out["TotalNumCases"]
hist(tot.out[,1])

tot.outsti<-outsti["stiTotalNumCases"]
hist(tot.outsti[,1])

tot.outsti6<-outsti6["stiTotalNumCases"]
hist(tot.outsti6[,1])

means_out<-colMeans(tot.out)
means_out
SEM.out<-apply(tot.out,2,std.error)
SEM.out

means.outsti<-colMeans(tot.outsti,na.rm=T)
means.outsti
SEM.outsti<-apply(tot.outsti,2,std.error)
SEM.outsti

means.outsti6<-colMeans(tot.outsti6,na.rm=T)
means.outsti6
SEM.outsti6<-apply(tot.outsti6,2,std.error)
SEM.outsti6

xx<-sapply(c(tot.out),as.numeric)
yy<-sapply(c(tot.outsti),as.numeric)
zz<-sapply(c(tot.outsti6),as.numeric)
wilcox.test(xx,yy,conf.int=TRUE)
wilcox.test(xx,zz,conf.int=TRUE)

#####################################################
#    Peak Size (Incidence) Number of New Cases per day
#####################################################

ps.out<-out["EpidemicPeakNew"]
hist(ps.out[,1])

ps.outsti<-outsti["stiEpidemicPeakNew"] 
hist(ps.outsti[,1])

ps.outsti6<-outsti6["stiEpidemicPeakNew"]
hist(ps.outsti6[,1])

means_out<-colMeans(ps.out)
means_out
SEM.out<-apply(ps.out,2,std.error)
SEM.out

means.outsti<-colMeans(ps.outsti,na.rm=T)
means.outsti
SEM.outsti<-apply(ps.outsti,2,std.error)
SEM.outsti

means.outsti6<-colMeans(ps.outsti6,na.rm=T)
means.outsti6
SEM.outsti6<-apply(ps.outsti6,2,std.error)
SEM.outsti6

xxx<-sapply(c(ps.out),as.numeric)
yyy<-sapply(c(ps.outsti),as.numeric)
zzz<-sapply(c(ps.outsti6),as.numeric)
wilcox.test(xxx,yyy,conf.int=TRUE)
wilcox.test(xxx,zzz,conf.int=TRUE)

#####################################################
#    Subset data to only those runs that reached 50 or more total cumulative cases
#####################################################

# subset data to exclude epidemics with under 100 cases
out<-subset(out,out["TotalNumCases"]>49)
nrow(out)

outsti<-subset(outsti,outsti["stiTotalNumCases"]>49)
nrow(outsti)

outsti6<-subset(outsti6,outsti6["stiTotalNumCases"]>49)
nrow(outsti6)

#####################################################
#    Date of last Case
#####################################################

lastcase.out<-out["DayLastCase"]
hist(lastcase.out[,1])

lastcase.outsti<-outsti["stiDayLastCase"] 
hist(lastcase.outsti[,1])

lastcase.outsti6<-outsti6["stiDayLastCase"] 
hist(lastcase.outsti6[,1])

means_out<-colMeans(lastcase.out,na.rm=T)
means_out
SEM.out<-apply(lastcase.out,2,std.error,na.rm=T)
SEM.out

means.outsti<-colMeans(lastcase.outsti,na.rm=T)
means.outsti
SEM.outsti<-apply(lastcase.outsti,2,std.error,na.rm=T)
SEM.outsti
sd.outsti<-apply(lastcase.outsti,2,sd,na.rm=T)
sd.outsti

means.outsti6<-colMeans(lastcase.outsti6,na.rm=T)
means.outsti6
SEM.outsti6<-apply(lastcase.outsti6,2,std.error,na.rm=T)
SEM.outsti6
sd.outsti6<-apply(lastcase.outsti6,2,sd,na.rm=T)
sd.outsti6

t.test(lastcase.out,lastcase.outsti)
t.test(lastcase.out,lastcase.outsti6)

xxxx<-sapply(c(lastcase.out),as.numeric)
yyyy<-sapply(c(lastcase.outsti),as.numeric)
zzzz<-sapply(c(lastcase.outsti6),as.numeric)
wilcox.test(xxxx,yyyy,conf.int=TRUE)
wilcox.test(xxxx,zzzz,conf.int=TRUE)

sum(lastcase.outsti>730,na.rm=T) # 2 years
sum(lastcase.outsti>0,na.rm=T)

sum(lastcase.outsti6>730,na.rm=T) # 2 years
sum(lastcase.outsti6>0,na.rm=T)


#####################################################
#    Date of Epidemic peak Incidence (max # new cases per day)
#####################################################

pt.out<-out["DayEpidemicPeakNew"]
hist(pt.out[,1])


pt.outsti<-outsti["stiDayEpidemicPeakNew"]
hist(pt.outsti[,1])

pt.outsti6<-outsti6["stiDayEpidemicPeakNew"]
hist(pt.outsti6[,1])

means_out<-colMeans(pt.out,na.rm=T)
means_out
SEM.out<-apply(pt.out,2,std.error,na.rm=T)
SEM.out

means.outsti<-colMeans(pt.outsti,na.rm=T)
means.outsti
SEM.outsti<-apply(pt.outsti,2,std.error,na.rm=T)
SEM.outsti
sd.outsti<-apply(pt.outsti,2,sd,na.rm=T)
sd.outsti

means.outsti6<-colMeans(pt.outsti6,na.rm=T)
means.outsti6
SEM.outsti6<-apply(pt.outsti6,2,std.error,na.rm=T)
SEM.outsti6
sd.outsti6<-apply(pt.outsti6,2,sd,na.rm=T)
sd.outsti6

t.test(pt.out,pt.outsti)

xxxxx<-sapply(c(pt.out),as.numeric)
yyyyy<-sapply(c(pt.outsti),as.numeric)
zzzzz<-sapply(c(pt.outsti6),as.numeric)
wilcox.test(xxxxx,yyyyy)
wilcox.test(xxxxx,zzzzz)


