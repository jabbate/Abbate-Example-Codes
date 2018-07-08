#######################################################
#
#   Use Rarefaction Curve to Investigate Associations
#
#######################################################

########################################################
#        Set Workspace & Libraries               
########################################################
########################################################
rm(list=ls(all=TRUE))

setwd("~/Dropbox/Work_Labs/Roche/NacherThaiCohorts/Nacher_Severity/")

# Libraries & Functions
# pour le calcul d'enveloppe de confiance ‡ 95 %
library(boot)
library(plotrix)
library(lattice)
# Fonction utilisÈe :
source("FctTestScreenENV.txt")

########################################################
#        Decide on Associations to Investigate               
########################################################
########################################################

# Define the Question:

al <- 1      # Ascaris lumbercoides
ww <- 1      # Trichuris trichiura (whipworm)
hw <- 1      # Ancylostoma duodenale or Necator americanus  (hookworm)
ss <- 1      # Strongyloides stercoralis

ov <- 0      # ov? too few positive!

Sexes  <-"all.sexes"          # "all.sexes", "males",    or "females"

Cases  <-"compare.severe"     # "all.cases", "compare.severe" is only hyperparasitemic or cerebral",
                              # "only.hyper", "only.cerebral", or "only.simple")

cm <- 0      # cerebral malaria
s2 <- 0      # simple (0) vs. severe (1)

runasa<- 1   # if 1, run full ASA along with MCA and AFC; if 0, don't run ASA

members<-c(cm,s2,al,ww,hw,ss,ov)
names(members)<-c("cm","s2","al","ww","hw","ss","ov")
dis<-c(which(members==1))
n<-length(dis)
clip<-c(Sexes,Cases)

########################################################
#        Set up Rarefaction Parameters & Res Vector               
########################################################
########################################################

# X-values for rarefaction curve (% data subset)
j<-seq(.1,1,.1)

# number of iterations per X value
reps<-20

# define results vector with x value plus the number of y results from tt
# number of pathogens
n

# number of combinations
  ncombos<-2^n
  y<-c(1:ncombos,"X")
  ly<-length(y)
  yres<-matrix(y,reps,ly, byrow = TRUE)
  zzres<-yres

#
run_rarefaction<-1   # 1 if running the simulations
########################################################
#        Run the ASA               
########################################################
########################################################
if(run_rarefaction==1){ 
  
  set.seed(100) 
  
  for(jj in 1:length(j)){
    rarex<-j[jj]
    ############################
    # Initiation pour les tirages alÈatoires
     
    for(i in 1:reps){
      
      ############################
      # DonnÈes initiales :
      Xdat <- read.csv("Nacher_Severity_DataFull.csv", sep=",")
      dim(Xdat)
      names(Xdat)

      ############################
      # subset data 
       Xdat<-subset(Xdat,Xdat$malariacase!='NA')
      # Xdat<-subset(Xdat,Xdat$sex!='NA')
      # Xdat<-subset(Xdat,Xdat$age!='NA')
      # Xdat<-subset(Xdat, Xdat$age<21)
      # Xdat<-subset(Xdat, Xdat$age>14)
      # Xdat<-subset(Xdat,Xdat$ankylost_hw==0)
      # generate ASA variable list
      asa.vars<-c()
      
      # conditions
      
      # Sexes
      if(Sexes=="males"){Xdat<-subset(Xdat,Xdat$sex=='homme')}
      if(Sexes=="females"){Xdat<-subset(Xdat,Xdat$sex=='femme')}
      
      # Cases
      if(Cases=="compare.severe"){Xdat<-subset(Xdat,Xdat$malariacase=="hyperparasitemic" | Xdat$malariacase=="cerebral")}
      if(Cases=="only.hyper"){Xdat<-subset(Xdat,Xdat$malariacase=="hyperparasitemic")}
      if(Cases=="only.cerebral"){Xdat<-subset(Xdat,Xdat$malariacase=="cerebral")}
      if(Cases=="only.simple"){Xdat<-subset(Xdat,Xdat$malariacase=="simple")}
      if(Cases=="all.cases"){Xdat<-subset(Xdat,Xdat$malariacase!='NA')}
      
      # Re-code
      Xdat$cm<-NA
      Xdat$cm[which(Xdat$malariacase=="hyperparasitemic")]<-0
      Xdat$cm[which(Xdat$malariacase=="cerebral")]<-1
      
      Xdat$s2<-NA
      Xdat$s2[which(Xdat$malariacase=="hyperparasitemic" |Xdat$malariacase=="cerebral")]<-1
      Xdat$s2[which(Xdat$malariacase=="simple")]<-0
                  

      Xdat$PAal<-NA
      Xdat$PAal[which(Xdat$ascaris==0)]<-0
      Xdat$PAal[which(Xdat$ascaris>0)]<-1
      
      Xdat$PAtt<-NA
      Xdat$PAtt[which(Xdat$trichuris_ww==0)]<-0
      Xdat$PAtt[which(Xdat$trichuris_ww>0)]<-1
      
      Xdat$PAhw<-NA
      Xdat$PAhw[which(Xdat$ankylost_hw==0)]<-0
      Xdat$PAhw[which(Xdat$ankylost_hw>0)]<-1
      
      Xdat$PAss<-NA
      Xdat$PAss[which(Xdat$anguilul_st==0)]<-0
      Xdat$PAss[which(Xdat$anguilul_st>0)]<-1
      
      Xdat$PAov<-NA
      Xdat$PAov[which(Xdat$ov==0)]<-0
      Xdat$PAov[which(Xdat$ov>0)]<-1
      
      # Presence
      
      if(cm==1){Xdat<-subset(Xdat,Xdat$cm!='NA'); asa.vars<-c(asa.vars,"cm")}else{asa.vars<-c(asa.vars,"NA")}
      if(s2==1){Xdat<-subset(Xdat,Xdat$s2!='NA'); asa.vars<-c(asa.vars,"s2")}else{asa.vars<-c(asa.vars,"NA")}
      
      if(al==1){Xdat<-subset(Xdat,Xdat$PAal!='NA'); asa.vars<-c(asa.vars,"PAal")}else{asa.vars<-c(asa.vars,"NA")}
      if(ww==1){Xdat<-subset(Xdat,Xdat$PAtt!='NA'); asa.vars<-c(asa.vars,"PAtt")}else{asa.vars<-c(asa.vars,"NA")}
      if(hw==1){Xdat<-subset(Xdat,Xdat$PAhw!='NA'); asa.vars<-c(asa.vars,"PAhw")}else{asa.vars<-c(asa.vars,"NA")}
      if(ss==1){Xdat<-subset(Xdat,Xdat$PAss!='NA'); asa.vars<-c(asa.vars,"PAss")}else{asa.vars<-c(asa.vars,"NA")}
      if(ov==1){Xdat<-subset(Xdat,Xdat$PAov!='NA'); asa.vars<-c(asa.vars,"PAov")}else{asa.vars<-c(asa.vars,"NA")}
      
      #set ASA variable list
      asa.vars<-c(asa.vars[which(members==1)])
      
      xtab.list<-c(asa.vars[1])
      for(nn in 2:n){xtab.list<-c(xtab.list,asa.vars[nn])}
      xtab.list
      
      Xtab <- as.matrix(Xdat[,xtab.list] ) 
      dim(Xtab);
      
      #rarefaction sampling (without replacemenet):
      percent<-round(rarex*dim(Xtab)[1])
      Xtab<- Xtab[sample(1:nrow(Xtab), percent,replace=FALSE),]
      
      #are there enough individuals?
      ncombos<-2^(length(asa.vars))
      too_few<-0
      
      dim(Xdat)[1]
      if(dim(Xdat)[1]<=ncombos){too_few<-1}
      
      ############################
      # RUN THE ASA
      if(too_few==1){print("Too few Individuals"); yy<-c(rep(NA,ncombos),rarex)}else{asa.vars;
        
        #the ASA
        #
        dim(Xtab);
        #
        # Struture de sauvegarde des rÈsultats :
        ResSCRENV <- FctTestScreenENV(Xtab)
        #
        # test significatif si > 0	
        ResSCRENV[[1]];	
        #
        # Pvalue			
        ResSCRENV[[2]];	
        #
        xNC <- length(ResSCRENV[[3]]);
        x <- (1:xNC);
        #     yinf <- ResSCRENV[[8]];
        #     ysup <- ResSCRENV[[7]];
        #     yobs <- ResSCRENV[[6]];
        #
        #     plot(x,yobs,type="n",ylim=range(c(0.0,yinf,ysup,yobs)),
        #          xlab="indice de combinaison",ylab="effectif",
        #          main="jeu de donnÈes BB/BR")
        #     lines(x,yinf,lty=1,col=4,lwd=2) #lower 95% 
        #     lines(x,ysup,lty=1,col=3,lwd=2) #upper 95%
        #     points(x,yobs,col=2) #observed
        valx <- x[(ResSCRENV[[4]]!=0)];
        #     for (i in 1:length(valx) ) {
        #          abline(v=valx[i],lty=3,col=1,lwd=1)
        #          }
        #
        # combinaison significatives des germes	 
        (ResSCRENV[[3]])[valx];	
        #
        # Pvalue des combinaisons significatives 
        (ResSCRENV[[5]])[valx];	
        #
        cbind(ResSCRENV[[3]],ResSCRENV[[4]],ResSCRENV[[5]]);
        # Combinaison - compris dans l'enveloppe (1) ou non (0) - P-value associÈe ‡ chaque combinaison
        cbind(ResSCRENV[[3]][valx],ResSCRENV[[6]][valx],ResSCRENV[[8]][valx],ResSCRENV[[7]][valx]);
        # Combinaison hors de l'enveloppe - obs - IC
        tt<-cbind(ResSCRENV[[3]],ResSCRENV[[6]],ResSCRENV[[8]],ResSCRENV[[7]],ResSCRENV[[4]],ResSCRENV[[5]]);
        # Combinaisons - obs - IC - significativitÈ - P-value
        tt;
        print(paste(jj,i))
        yy<-c(tt[,6],rarex) #pvalue
        zz<-c(tt[,5],rarex) #beyond envelope sig. = 1 or 0 
        }
      yres[i,]<-yy 
      zzres[i,]<-zz
    }
    if(rarex==j[1]){rareres<-yres}else{rareres<-rbind(rareres,yres)}
    if(rarex==j[1]){rarereszz<-zzres}else{rarereszz<-rbind(rarereszz,zzres)}
  }
  ############################
  #        Save ASA results along with Parameters               
  ############################
  
  rareres<-as.data.frame(rareres)
  names(rareres)[ncombos+1]<-"X"
  for(mm in 1:ncombos){
    names(rareres)[mm]<-tt[mm,1]
    }
  for(nn in 1:n){
    pathocol<-paste("Patho",nn,sep="")
    rareres[,pathocol]<-colnames(Xtab)[nn]
    }
  rareres[,"Clip"]<-paste(clip[1],clip[2],clip[3])
  rareres[,"n"]<-n
  rareres[,"j"]<-c(j)
  rareres[,"reps"]<-reps
  
  rarereszz<-as.data.frame(rarereszz)
  names(rarereszz)[ncombos+1]<-"X"
  for(mm in 1:ncombos){
    names(rarereszz)[mm]<-tt[mm,1]
  }
  for(nn in 1:n){
    pathocol<-paste("Patho",nn,sep="")
    rarereszz[,pathocol]<-colnames(Xtab)[nn]
  }
  rarereszz[,"Clip"]<-paste(clip[1],clip[2],clip[3])
  rarereszz[,"n"]<-n
  rarereszz[,"j"]<-c(j)
  rarereszz[,"reps"]<-reps
}
########################################################
#        Save Directionality of each Association              
########################################################
if(run_rarefaction==1){
  # last run of full data ASA results
  direction<-as.character(c(rep(NA,dim(rareres)[2])))
  for (qq in 1:ncombos){
    if(tt[qq,2]>ceiling((tt[qq,4]-tt[qq,3])/2)+tt[qq,3]){direction[qq]<-"Frequent"}else
      if(tt[qq,2]<floor((tt[qq,4]-tt[qq,3])/2)+tt[qq,3]){direction[qq]<-"Rare"}else{direction[qq]<-"Random"}
  }
  rareres[,1:ncombos] <- sapply(rareres[,1:ncombos],as.character) 
  rareres<-rbind(rareres,direction)
}

########################################################
#        Save Maximum P value for each Fraction              
########################################################


# to save the maximum pvalue for the most significant fraction: ######
#
# if(run_rarefaction==1){
#   max.pvals<-matrix(NA,1,ncol(rareres))
#   # Get Average Pvalue for each fraction (jj in j)
#   for (qq in 1:ncombos){
#     max.pval<-matrix(NA,1,length(j))
#     for(jj in 1:length(j)){
#       pvals<-as.numeric(as.character(rareres[which(rareres$X==j[jj]),(qq)]))
#       max.pval[1,jj]<-max(pvals,na.rm = TRUE)
#     }
#     max.pvals[qq]<-min(max.pval,na.rm = TRUE)
#   }
#   rareres<-rbind(rareres,max.pvals[1,])
# }                                             #########

# to save only the max pvalue for the full dataset (100 % data included)
if(run_rarefaction==1){
  max.pvals<-matrix(NA,1,ncol(rareres))
  # Get Max Pvalue for each fraction (jj in j)
  for (qq in 1:ncombos){
    max.pval<-matrix(NA,1,length(j))
    pval<-as.numeric(as.character(rareres[which(rareres$X==j[length(j)]),(qq)]))
    max.pvals[qq]<-max(pval,na.rm = TRUE)
  }
  rareres<-rbind(rareres,max.pvals[1,])
}

# to save the largest p-value considered significant in any run
if(run_rarefaction==1){
  
  names(rareres)
  dim(rareres)
  rareres[,1]
  ncombos
  
  sigsonly<-rarereszz[1:(length(j)*reps),1:ncombos]
  pvalsonly<-rareres[1:(length(j)*reps),1:ncombos]
  sigpvalsonly<-matrix(NA,(length(j)*reps),ncombos)
  combo.max.sig.pval<-matrix(NA,1,ncombos)
  
  # Get Max Pvalue for each fraction (jj in j)
  for (qq in 1:ncombos){
    for(rr in 1:(length(j)*reps)){
      if(sigsonly[rr,qq]==1){
      sigpvalsonly[rr,qq]<-pvalsonly[rr,qq]}
    }
    combo.max.sig.pval[qq]<-as.numeric(as.character(max(c(sigpvalsonly[,qq]),na.rm = TRUE)))
    colnames(combo.max.sig.pval)<-colnames(pvalsonly)
    }
  max.sig.pval<-as.numeric(as.character(max(sigpvalsonly,na.rm = TRUE)))
}

# save
# write.csv(rareres,file="RarefactionDataX.csv")
# write.table(rareres,file="RarefactionDataTable.csv",append=TRUE,sep=",")

########################################################
#
#   Rarefaction Curves
#
########################################################

# generate a value that describes the rarefaction curve for each association
  
########################################################
#        Set up input and output vectors               
########################################################
  
# input dataframe
if(run_rarefaction==1){out<-rareres}
if(run_rarefaction==0){out<-read.csv(file="RarefactionDataX.csv",header=TRUE);
out<-out[,2:dim(out)[2]]
n<-out[1,"n"]
reps<-out[1,"reps"]
if(n>3){vec<-as.numeric(as.character(names(table(out[,"j"])))); j<-c(vec)}else{j<-seq(.1,1,.1)}
}  


# output dataframe
ncombos<-2^n
RareRes<-as.data.frame(array(NA,c(ncombos,1)))
names(RareRes)[1]<-"Association"
RareRes[,"Association"]<-names(out[1:ncombos])
RareRes[,"Direction"]<-as.character(out[(dim(out)[1]-1),1:ncombos])
RareRes[,"MaxPvalue"]<-as.character(out[dim(out)[1],1:ncombos])

########################################################
#        Generate Rarefaction Curves for each Association               
########################################################

for (assoc in 1:ncombos)
{
  X<-as.numeric(as.character(out$X[1:(length(j)*reps)]))
  Y<-1-as.numeric(as.character(out[1:(length(j)*reps),assoc]))
  
  #print(xyplot(Y~X,type=c("p","l")))
  
  #Very Strong Associations
  # Y>0.99   / 95%
  pval<-0.99            # minimum significance level considered
  nval<-0.95            # minimum % of individuals above pval per rep set considered
  nabove<-nval*(length(Y)/length(j)) # minimum # of individuals above pval per rep set required for = 1.
  res<-c()
  for(i in 1:(length(j)))   #for each % of data tried, create a column in RareRes & calculate how many runs were > pval.
    {
    name<-paste("VSA_",j[i],sep="")
    RareRes[assoc,name]<-length(which(Y[(i*reps-reps+1):(i*reps)]>pval))
    if(RareRes[assoc,name]<nabove){RareRes[assoc,name]<-0; res<-append(res,0)}
    if(RareRes[assoc,name]>=nabove){RareRes[assoc,name]<-1; res<-append(res,1)}
     }
  RareRes[assoc,"VSA_pos"]<-sum(res,na.rm=TRUE)
  
  if(sum(res,na.rm=TRUE)>0)   #for each % of data tried, fill in RareRes with 1 or 0 for significance
  {
   sum(Y>pval)
   length(which(Y>pval))
   newX<-X[which(Y>pval)]
   tab.newX<-table(newX)
   VSA_min<-min(as.numeric(as.character(names(tab.newX[which(tab.newX>=nabove)]))))
   VSA_max<-max(as.numeric(as.character(names(tab.newX[which(tab.newX>=nabove)]))))
   VSA_length<-VSA_max-VSA_min
   RareRes[assoc,"VSA_min"]<-VSA_min
   RareRes[assoc,"VSA_max"]<-VSA_max
   RareRes[assoc,"VSA_length"]<-VSA_length
   
   if(VSA_max==1){RareRes[assoc,"VSA_hit_1"]<-1;
   
        RareRes[assoc,"Missed_VSA"]<-0        
        
        if(max(as.numeric(as.character(names(tab.newX[which(tab.newX<nabove)]))))>VSA_min)
        {
         RareRes[assoc,"VSA_skips"]<-1;
         skip_col<-max(which(tab.newX<nabove))
         max_col<-max(which(tab.newX>=nabove))
         VSA_sure_min<-min(as.numeric(as.character(names(tab.newX[(skip_col+1):max_col]))))
         VSA_sure_length<-VSA_max-VSA_sure_min
         RareRes[assoc,"VSA_sure_min"]<-VSA_sure_min
         RareRes[assoc,"VSA_sure_length"]<-VSA_sure_length
        }else{
             RareRes[assoc,"VSA_skips"]<-0
             RareRes[assoc,"VSA_sure_min"]<-VSA_min
             RareRes[assoc,"VSA_sure_length"]<-VSA_length
              }

   }else{RareRes[assoc,"VSA_hit_1"]<-0;
   
         RareRes[assoc,"Missed_VSA"]<-1;
         RareRes[assoc,"VSA_sure_min"]<-VSA_min
         RareRes[assoc,"VSA_sure_max"]<-VSA_max
         
       if(max(as.numeric(as.character(names(tab.newX[which(tab.newX<nabove)]))))>VSA_min)
           {
           RareRes[assoc,"VSA_skips"]<-1;
           skip_col<-max(which(tab.newX<nabove))
           max_col<-max(which(tab.newX>=nabove))
           min_col<-min(which(tab.newX>=nabove))
           VSA_sure_min1<-min(as.numeric(as.character(names(tab.newX[(skip_col+1):max_col]))))
           VSA_sure_min2<-min(as.numeric(as.character(names(tab.newX[min_col:(skip_col+1)]))))
           length1<-VSA_max-VSA_sure_min1
           length2<-VSA_sure_min2-VSA_min
           VSA_sure_length<-max(length1,length2)
           RareRes[assoc,"VSA_sure_length"]<-VSA_sure_length
            }else{
                  RareRes[assoc,"VSA_skips"]<-0
                  RareRes[assoc,"VSA_sure_length"]<-VSA_length
                  }
   
    }
  }else{
      RareRes[assoc,"VSA_hit_1"]<-NA
      RareRes[assoc,"VSA_skips"]<-NA
      RareRes[assoc,"VSA_min"]<-NA
      RareRes[assoc,"VSA_max"]<-NA
      RareRes[assoc,"VSA_length"]<-NA
      RareRes[assoc,"VSA_sure_min"]<-NA
      RareRes[assoc,"VSA_sure_length"]<-NA
      RareRes[assoc,"Missed_VSA"]<-NA
         }
  
  #Strong Associations
  # Y>0.95   / 95%
  pval<-0.95
  nval<-0.95
  nabove<-nval*(length(Y)/length(j))
  res<-c()
  for(i in 1:(length(j)))   #for each % of data tried, create a column in RareRes & calculate how many runs were > pval.
  {
    name<-paste("StA_",j[i],sep="")
    RareRes[assoc,name]<-length(which(Y[(i*reps-reps+1):(i*reps)]>pval))
    if(RareRes[assoc,name]<nabove){RareRes[assoc,name]<-0; res<-append(res,0)}
    if(RareRes[assoc,name]>=nabove){RareRes[assoc,name]<-1; res<-append(res,1)}
  }
  RareRes[assoc,"StA_pos"]<-sum(res,na.rm=TRUE)
  
  if(sum(res,na.rm=TRUE)>0)   #for each % of data tried, fill in RareRes with 1 or 0 for significance
  {
    sum(Y>pval)
    length(which(Y>pval))
    newX<-X[which(Y>pval)]
    tab.newX<-table(newX)
    StA_min<-min(as.numeric(as.character(names(tab.newX[which(tab.newX>=nabove)]))))
    StA_max<-max(as.numeric(as.character(names(tab.newX[which(tab.newX>=nabove)]))))
    StA_length<-StA_max-StA_min
    RareRes[assoc,"StA_min"]<-StA_min
    RareRes[assoc,"StA_max"]<-StA_max
    RareRes[assoc,"StA_length"]<-StA_length
    
    if(StA_max==1){RareRes[assoc,"StA_hit_1"]<-1;
    
    RareRes[assoc,"Missed_StA"]<-0        
    
    if(max(as.numeric(as.character(names(tab.newX[which(tab.newX<nabove)]))))>StA_min)
    {
      RareRes[assoc,"StA_skips"]<-1;
      skip_col<-max(which(tab.newX<nabove))
      max_col<-max(which(tab.newX>=nabove))
      StA_sure_min<-min(as.numeric(as.character(names(tab.newX[(skip_col+1):max_col]))))
      StA_sure_length<-StA_max-StA_sure_min
      RareRes[assoc,"StA_sure_min"]<-StA_sure_min
      RareRes[assoc,"StA_sure_length"]<-StA_sure_length
    }else{
      RareRes[assoc,"StA_skips"]<-0
      RareRes[assoc,"StA_sure_min"]<-StA_min
      RareRes[assoc,"StA_sure_length"]<-StA_length
    }
    
    }else{RareRes[assoc,"StA_hit_1"]<-0;
    
    RareRes[assoc,"Missed_StA"]<-1;
    RareRes[assoc,"StA_sure_min"]<-StA_min
    
    if(max(as.numeric(as.character(names(tab.newX[which(tab.newX<nabove)]))))>StA_min)
    {
      RareRes[assoc,"StA_skips"]<-1;
      skip_col<-max(which(tab.newX<nabove))
      max_col<-max(which(tab.newX>=nabove))
      min_col<-min(which(tab.newX>=nabove))
      StA_sure_min1<-min(as.numeric(as.character(names(tab.newX[(skip_col+1):max_col]))))
      StA_sure_min2<-min(as.numeric(as.character(names(tab.newX[min_col:(skip_col+1)]))))
      length1<-StA_max-StA_sure_min1
      length2<-StA_sure_min2-StA_min
      StA_sure_length<-max(length1,length2)
      RareRes[assoc,"StA_sure_length"]<-StA_sure_length
    }else{
      RareRes[assoc,"StA_skips"]<-0
      RareRes[assoc,"StA_sure_length"]<-StA_length
    }
    
    }
  }else{
    RareRes[assoc,"StA_hit_1"]<-NA
    RareRes[assoc,"StA_skips"]<-NA
    RareRes[assoc,"StA_min"]<-NA
    RareRes[assoc,"StA_max"]<-NA
    RareRes[assoc,"StA_length"]<-NA
    RareRes[assoc,"StA_sure_min"]<-NA
    RareRes[assoc,"StA_sure_length"]<-NA
    RareRes[assoc,"Missed_StA"]<-NA
  }  
 
  
  #Weaker Associations
  # Y>0.90     /95%
  pval<-0.90
  nval<-0.95
  nabove<-nval*(length(Y)/length(j))
  res<-c()
  for(i in 1:(length(j)))   #for each % of data tried, create a column in RareRes & calculate how many runs were > pval.
  {
    name<-paste("WkA_",j[i],sep="")
    RareRes[assoc,name]<-length(which(Y[(i*reps-reps+1):(i*reps)]>pval))
    if(RareRes[assoc,name]<nabove){RareRes[assoc,name]<-0; res<-append(res,0)}
    if(RareRes[assoc,name]>=nabove){RareRes[assoc,name]<-1; res<-append(res,1)}
  }
  RareRes[assoc,"WkA_pos"]<-sum(res,na.rm=TRUE)
  
  if(sum(res,na.rm=TRUE)>0)   #for each % of data tried, fill in RareRes with 1 or 0 for significance if there were any (sum(res)>0)
  {
    sum(Y>pval)
    length(which(Y>pval))
    newX<-X[which(Y>pval)]
    tab.newX<-table(newX)
    WkA_min<-min(as.numeric(as.character(names(tab.newX[which(tab.newX>=nabove)]))))
    WkA_max<-max(as.numeric(as.character(names(tab.newX[which(tab.newX>=nabove)]))))
    WkA_length<-WkA_max-WkA_min
    RareRes[assoc,"WkA_min"]<-WkA_min
    RareRes[assoc,"WkA_max"]<-WkA_max
    RareRes[assoc,"WkA_length"]<-WkA_length
    
    if(WkA_max==1){RareRes[assoc,"WkA_hit_1"]<-1;
    
    RareRes[assoc,"Missed_WkA"]<-0        
    
    if(max(as.numeric(as.character(names(tab.newX[which(tab.newX<nabove)]))))>WkA_min)
    {
      RareRes[assoc,"WkA_skips"]<-1;
      skip_col<-max(which(tab.newX<nabove))
      max_col<-max(which(tab.newX>=nabove))
      WkA_sure_min<-min(as.numeric(as.character(names(tab.newX[(skip_col+1):max_col]))))
      WkA_sure_length<-WkA_max-WkA_sure_min
      RareRes[assoc,"WkA_sure_min"]<-WkA_sure_min
      RareRes[assoc,"WkA_sure_length"]<-WkA_sure_length
    }else{
      RareRes[assoc,"WkA_skips"]<-0
      RareRes[assoc,"WkA_sure_min"]<-WkA_min
      RareRes[assoc,"WkA_sure_length"]<-WkA_length
    }
    
    }else{RareRes[assoc,"WkA_hit_1"]<-0;
    
    RareRes[assoc,"Missed_WkA"]<-1;
    RareRes[assoc,"WkA_sure_min"]<-WkA_min
    
    if(max(as.numeric(as.character(names(tab.newX[which(tab.newX<nabove)]))))>WkA_min)
    {
      RareRes[assoc,"WkA_skips"]<-1;
      skip_col<-max(which(tab.newX<nabove))
      max_col<-max(which(tab.newX>=nabove))
      min_col<-min(which(tab.newX>=nabove))
      WkA_sure_min1<-min(as.numeric(as.character(names(tab.newX[(skip_col+1):max_col]))))
      WkA_sure_min2<-min(as.numeric(as.character(names(tab.newX[min_col:(skip_col+1)]))))
      length1<-WkA_max-WkA_sure_min1
      length2<-WkA_sure_min2-WkA_min
      WkA_sure_length<-max(length1,length2)
      RareRes[assoc,"WkA_sure_length"]<-WkA_sure_length
    }else{
      RareRes[assoc,"WkA_skips"]<-0
      RareRes[assoc,"WkA_sure_length"]<-WkA_length
    }
    
    }
  }else{
    RareRes[assoc,"WkA_hit_1"]<-NA
    RareRes[assoc,"WkA_skips"]<-NA
    RareRes[assoc,"WkA_min"]<-NA
    RareRes[assoc,"WkA_max"]<-NA
    RareRes[assoc,"WkA_length"]<-NA
    RareRes[assoc,"WkA_sure_min"]<-NA
    RareRes[assoc,"WkA_sure_length"]<-NA
    RareRes[assoc,"Missed_WkA"]<-NA
  }  

 #Strongest (Bonferroni-corrected) Associations (VVSA)
  # Y>0.95 /95%    
  #(bonferroni corrected: 0.5/32=0.00156 > 0.99844; 0.05/16=0.003125 > 0.996875; 0.05/8=0.00625 > 0.99375 )
  #  pval<-1-(0.05/ncombos)
  
  # using combo.max.sig.pval
  defaultpval<-0.05/ncombos
  newpvals<-combo.max.sig.pval
  newpvals[which(combo.max.sig.pval=="-Inf")]<-defaultpval
  newpvals[which(is.na(combo.max.sig.pval))]<-defaultpval
  newpvals[which(combo.max.sig.pval<defaultpval)]<-defaultpval
  
  newpval<-newpvals[assoc]
  
  pval<-1-newpval  
  nval<-0.95
  nabove<-nval*(length(Y)/length(j))
  res<-c()
  for(i in 1:(length(j)))   #for each % of data tried, create a column in RareRes & calculate how many runs were > pval.
  {
    name<-paste("VVSA_",j[i],sep="")
    RareRes[assoc,name]<-length(which(Y[(i*reps-reps+1):(i*reps)]>pval))
    if(RareRes[assoc,name]<nabove){RareRes[assoc,name]<-0; res<-append(res,0)}
    if(RareRes[assoc,name]>=nabove){RareRes[assoc,name]<-1; res<-append(res,1)}
  }
  RareRes[assoc,"VVSA_pos"]<-sum(res,na.rm=TRUE)
  
  if(sum(res,na.rm=TRUE)>0)   #for each % of data tried, fill in RareRes with 1 or 0 for significance if there were any (sum(res)>0)
  {
    sum(Y>pval)
    length(which(Y>pval))
    newX<-X[which(Y>pval)]
    tab.newX<-table(newX)
    VVSA_min<-min(as.numeric(as.character(names(tab.newX[which(tab.newX>=nabove)]))))
    VVSA_max<-max(as.numeric(as.character(names(tab.newX[which(tab.newX>=nabove)]))))
    VVSA_length<-VVSA_max-VVSA_min
    RareRes[assoc,"VVSA_min"]<-VVSA_min
    RareRes[assoc,"VVSA_max"]<-VVSA_max
    RareRes[assoc,"VVSA_length"]<-VVSA_length
    
    if(VVSA_max==1){RareRes[assoc,"VVSA_hit_1"]<-1;
    
    RareRes[assoc,"Missed_VVSA"]<-0        
    
    if(max(as.numeric(as.character(names(tab.newX[which(tab.newX<nabove)]))))>VVSA_min)
    {
      RareRes[assoc,"VVSA_skips"]<-1;
      skip_col<-max(which(tab.newX<nabove))
      max_col<-max(which(tab.newX>=nabove))
      VVSA_sure_min<-min(as.numeric(as.character(names(tab.newX[(skip_col+1):max_col]))))
      VVSA_sure_length<-VVSA_max-VVSA_sure_min
      RareRes[assoc,"VVSA_sure_min"]<-VVSA_sure_min
      RareRes[assoc,"VVSA_sure_length"]<-VVSA_sure_length
    }else{
      RareRes[assoc,"VVSA_skips"]<-0
      RareRes[assoc,"VVSA_sure_min"]<-VVSA_min
      RareRes[assoc,"VVSA_sure_length"]<-VVSA_length
    }
    
    }else{RareRes[assoc,"VVSA_hit_1"]<-0;
    
    RareRes[assoc,"Missed_VVSA"]<-1;
    RareRes[assoc,"VVSA_sure_min"]<-VVSA_min
    
    if(max(as.numeric(as.character(names(tab.newX[which(tab.newX<nabove)]))))>VVSA_min)
    {
      RareRes[assoc,"VVSA_skips"]<-1;
      skip_col<-max(which(tab.newX<nabove))
      max_col<-max(which(tab.newX>=nabove))
      min_col<-min(which(tab.newX>=nabove))
      VVSA_sure_min1<-min(as.numeric(as.character(names(tab.newX[(skip_col+1):max_col]))))
      VVSA_sure_min2<-min(as.numeric(as.character(names(tab.newX[min_col:(skip_col+1)]))))
      length1<-VVSA_max-VVSA_sure_min1
      length2<-VVSA_sure_min2-VVSA_min
      VVSA_sure_length<-max(length1,length2)
      RareRes[assoc,"VVSA_sure_length"]<-VVSA_sure_length
    }else{
      RareRes[assoc,"VVSA_skips"]<-0
      RareRes[assoc,"VVSA_sure_length"]<-VVSA_length
    }
    
    }
  }else{
    RareRes[assoc,"VVSA_hit_1"]<-NA
    RareRes[assoc,"VVSA_skips"]<-NA
    RareRes[assoc,"VVSA_min"]<-NA
    RareRes[assoc,"VVSA_max"]<-NA
    RareRes[assoc,"VVSA_length"]<-NA
    RareRes[assoc,"VVSA_sure_min"]<-NA
    RareRes[assoc,"VVSA_sure_length"]<-NA
    RareRes[assoc,"Missed_VVSA"]<-NA
  }  
  
  
  }

########################################################
#        Save Parameters               
########################################################

# save pathogens used
pathos_used<-"pathos"
for (pp in 1:n){patho_used<-out[1,(ncombos+2+pp-1)];pathos_used<-paste(pathos_used,patho_used)}
RareRes[,"Pathos"]<-pathos_used

# save parameters explored
RareRes[,"Clip"]<-out[1,"Clip"]

# save minimum pvalue for significant SCN result per combo
for(oo in 1:ncombos){
RareRes[oo,"maxSignSCNpvalue"]<-combo.max.sig.pval[oo]}

########################################################
#        Save Rarefaction Curve Results               
########################################################

RareRes
# save
# write.csv(RareRes,file="RarefactionCurvesX.csv")
# write.table(rareres,file="RarefactionCurvesTable.csv",append=TRUE,sep=",")

########################################################
#        Visualize Rarefaction Curve Results               
########################################################


#VSA
VSAres<-RareRes[ncombos:1,4:(3+length(j))] #this flips the Y-axis so '111111' is on top and '0' is on bottom
VSAres<-t(as.matrix(VSAres))  # this organizes the output from left to right (90deg rotate left will display)
colnames(VSAres)<-RareRes[ncombos:1,"Association"]
rownames(VSAres)<-names(RareRes[,4:(3+length(j))])
#image(VSAres,xlab="% data sampled",ylab="Association")
ymax<-ncombos
xmax<-length(j)
ylabs<-seq(0,(1/ncombos)*(ncombos+1),1/(ncombos-1))
xlabs<-seq(0,(1/xmax)*(xmax+1),1/(xmax-1))
image(VSAres,axes=FALSE,xlab="% data sampled",ylab="Association",col=c("red","yellow"),main="Very Strong Associations")
text(x=0, y=ylabs,labels=colnames(VSAres),cex=0.8)
text(x=xlabs, y=0-(1/xmax),labels=rownames(VSAres),srt=45,cex=0.5)
#text(x=1, y=ylabs,labels=RareRes[dim(RareRes)[1]:1,"Direction"],cex=0.8,pos=2)
sum(RareRes[,"Missed_VSA"],na.rm=TRUE)
summary<-c(RareRes[which(RareRes$VSA_min!="NA"),"Association"],
           RareRes[which(RareRes$VSA_min!="NA"),"Direction"] ,
           RareRes[which(RareRes$VSA_min!="NA"),"VSA_sure_min"],
           RareRes[which(RareRes$VSA_min!="NA"),"VSA_skips"])
pathos_used
matrix(summary,nrow=dim(RareRes[which(RareRes$VSA_min!="NA"),])[1],ncol=4,byrow=FALSE)

#StA
StAres<-RareRes[ncombos:1,(3+length(j)+10):(2+length(j)*2+10)] #this flips the Y-axis so '111111' is on top and '0' is on bottom
StAres<-t(as.matrix(StAres))  # this organizes the output from left to right (90deg rotate left will display)
colnames(StAres)<-RareRes[ncombos:1,"Association"]
rownames(StAres)<-names(RareRes[,(3+length(j)+10):(2+length(j)*2+10)])
#image(StAres,xlab="% data sampled",ylab="Association")
ymax<-ncombos
xmax<-length(j)
ylabs<-seq(0,(1/ncombos)*(ncombos+1),1/(ncombos-1))
xlabs<-seq(0,(1/xmax)*(xmax+1),1/(xmax-1))
image(StAres,axes=FALSE,xlab="% data sampled",ylab="Association",col=c("red","yellow"),main="Strong Associations")
text(x=0, y=ylabs,labels=colnames(StAres),cex=0.8)
text(x=xlabs, y=0-(1/xmax),labels=rownames(StAres),srt=45,cex=0.5)
#text(x=1, y=ylabs,labels=RareRes[dim(RareRes)[1]:1,"Direction"],cex=0.8,pos=2)
sum(RareRes[,"Missed_StA"],na.rm=TRUE)
summary<-c(RareRes[which(RareRes$StA_min!="NA"),"Association"],
           RareRes[which(RareRes$StA_min!="NA"),"Direction"] ,
           RareRes[which(RareRes$StA_min!="NA"),"StA_sure_min"],
           RareRes[which(RareRes$StA_min!="NA"),"StA_skips"])
pathos_used
matrix(summary,nrow=dim(RareRes[which(RareRes$StA_min!="NA"),])[1],ncol=4,byrow=FALSE)

#WkA
WkAres<-RareRes[ncombos:1,(2+(length(j)+10)*2):(1+length(j)*3+10*2)] #this flips the Y-axis so '111111' is on top and '0' is on bottom
WkAres<-t(as.matrix(WkAres))  # this organizes the output from left to right (90deg rotate left will display)
colnames(WkAres)<-RareRes[ncombos:1,"Association"]
rownames(WkAres)<-names(RareRes[,(2+(length(j)+10)*2):(1+length(j)*3+10*2)])
#image(WkAres,xlab="% data sampled",ylab="Association")
ymax<-ncombos
xmax<-length(j)
ylabs<-seq(0,(1/ncombos)*(ncombos+1),1/(ncombos-1))
xlabs<-seq(0,(1/xmax)*(xmax+1),1/(xmax-1))
image(WkAres,axes=FALSE,xlab="% data sampled",ylab="Association",col=c("red","yellow"),main="Weaker Associations")
text(x=0, y=ylabs,labels=colnames(WkAres),cex=0.8)
text(x=xlabs, y=0-(1/xmax),labels=rownames(WkAres),srt=45,cex=0.5)
#text(x=1, y=ylabs,labels=RareRes[dim(RareRes)[1]:1,"Direction"],cex=0.8,pos=2)
sum(RareRes[,"Missed_WkA"],na.rm=TRUE)
summary<-c(RareRes[which(RareRes$WkA_min!="NA"),"Association"],
           RareRes[which(RareRes$WkA_min!="NA"),"Direction"] ,
           RareRes[which(RareRes$WkA_min!="NA"),"WkA_sure_min"],
           RareRes[which(RareRes$WkA_min!="NA"),"WkA_skips"])
pathos_used
matrix(summary,nrow=dim(RareRes[which(RareRes$WkA_min!="NA"),])[1],ncol=4,byrow=FALSE)

#VVSA
VVSAres<-RareRes[ncombos:1,(1+length(j)*3+10*3):(1+(3+length(j)+10)*3)] #this flips the Y-axis so '111111' is on top and '0' is on bottom
VVSAres<-t(as.matrix(VVSAres))  # this organizes the output from left to right (90deg rotate left will display)
colnames(VVSAres)<-RareRes[ncombos:1,"Association"]
rownames(VVSAres)<-names(RareRes[,(2+(length(j)+10)*2):(1+length(j)*3+10*2)])
#image(VVAres,xlab="% data sampled",ylab="Association")
ymax<-ncombos
xmax<-length(j)
ylabs<-seq(0,(1/ncombos)*(ncombos+1),1/(ncombos-1))
xlabs<-seq(0,(1/xmax)*(xmax+1),1/(xmax-1))
image(VVSAres,axes=FALSE,xlab="% data sampled",ylab="Association",col=c("red","yellow"),main="VVStrong Associations")
text(x=0, y=ylabs,labels=colnames(VVSAres),cex=0.8)
text(x=xlabs, y=0-(1/xmax),labels=rownames(VVSAres),srt=45,cex=0.5)
#text(x=1, y=ylabs,labels=RareRes[dim(RareRes)[1]:1,"Direction"],cex=0.8,pos=2)
sum(RareRes[,"Missed_VVSA"],na.rm=TRUE)
summary<-c(RareRes[which(RareRes$VVSA_min!="NA"),"Association"],
           RareRes[which(RareRes$VVSA_min!="NA"),"Direction"] ,
           RareRes[which(RareRes$VVSA_min!="NA"),"VVSA_sure_min"],
           RareRes[which(RareRes$VVSA_min!="NA"),"VVSA_skips"])
pathos_used
matrix(summary,nrow=dim(RareRes[which(RareRes$VVSA_min!="NA"),])[1],ncol=4,byrow=FALSE)

#write.csv(RareRes,file="RarefactionCurvesX.csv")
combo.max.sig.pval
