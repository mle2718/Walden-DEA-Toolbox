#This program estimates the OrderM model 
#Author: John Walden, NOAA/NMFS/NEFSC 
#Date Modified: 17-JAN-2025
rm(list=ls())
library(plyr)
library(dplyr)
library(data.table)
library(readr)
#####################################################################################
#Read in Data for DEA Model
#####################################################################################
df=read_csv("./data/spec_cap_distort_072024.csv", show_col_types = FALSE) #Read in data
nDMUs=nrow(df)
obs<-data.frame(1:nDMUs)
XF=as.matrix(df[,c("LEN","VHP")])
XV=as.matrix(df[,c("tda")])
# create output and input Y and X matrices
Y<-as.matrix(df[,c("sp1_lbs", "sp2_lbs","sp3_lbs","sp4_lbs")])
X<-XF
M=ncol(Y)
N=ncol(X)
J=nrow(X)
YX<-cbind(Y,X)
##################################################################################
####################################################################
#For output oriented model
#############################################################################
#Resample function
resample<-function(y,m){
  
  j=nrow(y)
  sample=floor(j*runif(m,min=0, max=1)+1)
  yb<-as.matrix(y[sample,])
  return(yb)
}
########################################################################
#Order-M function
eff_m<-function(x0,y0,X,Y,m,B){
  
  #xi<-as.vector(x0)
  j=nrow(X)
  xkv<-matrix(1,j,1)%*%x0
  flagx=(X<=xkv)                      #Determine observations with X<=x0
  # flagx<-as.data.frame(t(flagx))    #convert to data frame
  nx<-ncol(X)
  flag=rowSums(flagx)
  YM<-as.matrix(Y[(flag==nx),])
  
  effmk=0
  thetab=0;

  for (b in 1:B){
 
    YMb=as.matrix(resample(YM,m))   #resample is a function defined above
    mc=ncol(YMb)
    jp<-nrow(YMb)
    
    eff1=0                 #to hold results
    eff2=0                 #to hold results
    
    y0[y0==0]<-NA  #Need to handle zero valued outputs in vector
    for(i in 1:jp){
      eff1[i]<-min((YMb[i,]/y0), na.rm=TRUE)
    }
    eff2=max(eff1,1)
    res<-round(eff2,3)

    thetab[b]=res    
  }   
   
   lt=list(theta=mean(thetab),std=(sd(thetab)/sqrt(B)))
   newlist<-list(lt)
  # 
  #theta=mean(thetab)
  
  return(newlist)
}
################################################################################ 
#Calculate Order-M
################################################################################
thetab=0
cap=0
std=0
B=250  #Number of simulations
set.seed(12345)
for(j in 1:J){
  xref<-X[j,]
  yref<-Y[j,]
  thetab[j]=eff_m(xref,yref,X,Y,30,B)
  
  cap[j]=thetab[[j]]$theta
  std[j]=thetab[[j]]$std
}

ordm1<-as.data.frame(cbind(cap,std))
ordm1$cu=1/ordm1$cap
obs<-as.data.frame(obs)
ordm1<-cbind(obs,ordm1$cap,ordm1$cu)
colnames(ordm1)<-c("Obs","Cap","CU")
################################################################################