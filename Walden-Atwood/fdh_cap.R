#############################################################################
#Johansen Capacity FDH Model
#This version uses a MIP model to solve for the FDH capacity output.
#Author: John Walden, NOAA/NMFS/NEFSC
#Modified: 16-Jan-2025
############################################################################
rm(list=ls())
library(Rglpk)
library(dplyr)
library(data.table)
library(stats)
library(readr)
#####################################################################################
#Read in Data for DEA Model
#####################################################################################
df=read_csv("./data/spec_cap_distort_072024.csv", show_col_types = FALSE) #Read in data
nDMUs=nrow(df)
obs<-data.frame(1:nDMUs)
#############################################################################################
#Define Matrices
###########################################################################################
XF=as.matrix(df[,c("LEN","VHP")])
XV=as.matrix(df[,c("tda")])
# create output and input Y and X matrices
Y<-as.matrix(df[,c("sp1_lbs", "sp2_lbs","sp3_lbs","sp4_lbs")])
X<-XF
M=ncol(Y)
N=ncol(X)
J=nrow(X)
YX<-cbind(Y,X)
###################################################################################
#Build A Matrix
####################################################################################
A=matrix(0,M+N+1,J+1)
A[1:(M+N),1:J]=t(YX)
A[(M+N+1),1:J]=1.0
####################################################################################
#set objective function
obj=c(rep(0,J),1)
rest=c(rep('>=',M),rep('<=',N),'==')
max<-TRUE
objvals1=0
status1=0
peers<-matrix(0,J,J)#create matrix of zero's to hold peer weights
types1<-c(rep("B",J),"C")   #Binary values set for J constraints
res<-matrix(0,J,3)
#####################################################################################
#Solve model in Loop
#####################################################################################
for(j in 1:J){
  A[(1:M),J+1]=-A[1:M,j]
  rhs=c(rep(0,M),as.matrix(X[j,]),1) 
  sol<- Rglpk_solve_LP(obj=obj, mat=A, dir=rest, types=types1, rhs=rhs, max=max)
  peers[j,]<-sol$solution[1:J]
  if(j%%100==0|j==J)  print(paste('on dmu',j,'of',J))
  res[j,1]=j
  res[j,2]=round(sol$optimum,3)
  res[j,3]=sol$status
}
#####################################################################################
#End of loop
#####################################################################################
summary(res[,3])  #tells whether model solved at all iterations. Value of 0 means it solved

VI<-round((peers%*%XV),3)  #Determine optimal variable input levels

df2<-as.data.frame(cbind(res[,1:2],VI))
colnames(df2)[1:2]<-c("Obs","Cap")
df2$CU=1/df2$Cap
###################################################################################
##################################################################################
