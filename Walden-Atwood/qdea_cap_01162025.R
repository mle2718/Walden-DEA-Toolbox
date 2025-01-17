#This is the QDEA capacity model
#The QDEA package is available by contacting Joe Atwood
#jatwood AT montana.edu
#Author: John Walden, NOAA/NMFS/NEFSC 
#Date Modified: 06-JAN-2025
###############################################################################
rm(list=ls())
graphics.off()
###############################################################################
library(qDEAC)    
library(dplyr)
library(data.table)
library(readr)
###############################################################################
#Function to return Peer values to calculate optimal variable input levels
##############################################################################
peers_d=function(DPEERS,X){
  
  tmp2=list(nr=nrow(X),nc=nrow(X),ia=DPEERS$dmu0,ja=DPEERS$dmuz,ra=DPEERS$z,
            SMM='RCI',ZINDEX=F)
  
  # convert sparse matrix object to matrix
  peers=SM2A(tmp2)
  tmp=apply(peers,1,sum)
  pick=which(tmp==0)
  peers[cbind(pick,pick)]=1
  
  return(peers)
}
################################################################################
#Read in Data
################################################################################
df=read_csv("./data/spec_cap_distort_072024.csv", show_col_types = FALSE) #Read in data
nDMUs=nrow(df)
obs<-data.frame(1:nDMUs)
#Put inputs into a Fixed input matrix and a Variable Input Matrix
XF=as.matrix(df[,c("LEN","VHP")])
XV=as.matrix(df[,c("tda")])
# create output X and input Y and X matrices
Y<-as.matrix(df[,c("sp1_lbs", "sp2_lbs","sp3_lbs","sp4_lbs")])
X<-XF
M=ncol(Y)
N=ncol(X)
J=nrow(X)
YX<-cbind(Y,X)
###############################################################################
################################################################################
#QDEA MODEL no Bootstrap
################################################################################
###########################################################################
#Note: 3 observations are held out at each iteration i.e. qout=3/nDMUs.
tmp1=qDEAC(X,Y,X0=X,Y0=Y,orient='OUT',nqiter=4,RTS='VRS',
     qout=3/nDMUs, mcells=4, replaceA1 = F, getproject=T)
##############################################################################
eff=tmp1$effvalsq            #get efficiency scores
eff[eff<1|is.na(eff)]=1      #For those iterations that don't solve eff=1

obs<-c(1:J)                  #Put number of observations in a vector
eff<-as.data.frame(cbind(obs,eff)) #create dataframe with observations and eff
#########################################################################
#calculate optimal variable input usage
#First get peer values
PEERS<-tmp1$PEER_DATA$PEERSq
PEERS$dmuz=ifelse(!is.na(PEERS$dmuz),PEERS$dmuz,PEERS$dmu0)
PEERS$z=ifelse(!is.na(PEERS$z),PEERS$z,1)
PEERS$z=ifelse(PEERS$z <0,0,PEERS$z)
PEERS=orderBy(~dmu0+dmuz,data=PEERS)  #order by dmu0
peersq<-peers_d(PEERS,X)              #call subroutine from above

opt_xv<-peersq%*%XV                   #calculate optimal variable inputs
######################################################################
#cbind results and optimal Variable inputs and save in a dataframe
results<-as.data.frame(cbind(eff,opt_xv))  
results$cu=1/results$eff
######################################################################
