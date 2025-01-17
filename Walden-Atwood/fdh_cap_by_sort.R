#FDH estimator based on enumeration
#Author: John B. Walden, NOAA/NMFS/NEFSC
#Modified: 17-Jan-2025
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
#############################################################################
#Function for output oriented FDH
fdh<-function(x0,y0,X,Y){
  
  j=nrow(X)
  xkv<-matrix(1,j,1)%*%x0
  flagx=(X<=xkv)                  #Determine observations with X<=x0
  nx<-ncol(X)
  flag=rowSums(flagx)
  YM<-as.matrix(Y[(flag==nx),])
  jp<-nrow(YM)
  
  eff1=0                 #to hold results
  eff2=0                 #to hold results
  
  y0[y0==0]<-NA  #Need to handle zero valued outputs in vector
  for(i in 1:jp){
    eff1[i]<-min((YM[i,]/y0), na.rm=TRUE)
  }
  eff2=max(eff1,1)
  res<-round(eff2,3)
  return(res)
}
##############################################################################
#end of function
##############################################################################
res=0
#Begin loop for FDH calculation
for(j in 1:J){
  xref<-X[j,]
  yref<-Y[j,]
  
  res[j]<-fdh(xref,yref,X,Y)
}
#End Loop
#Put results in dataframe
results<-as.data.frame(cbind(obs,res))
colnames(results)<-c("Obs","Cap")
results$CU=1/results$Cap
#############################################################################