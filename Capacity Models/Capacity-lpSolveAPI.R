####################################################################
#DEA Capacity Model 5-24-2022
#This program estimates the Johansen capacity model as found in the text
#Production Frontiers By Fare, Grosskopf and Lovell (1994)
#Capacity is modelled for fishing vessels in the Northwest Atlantic
#Author: John B. Walden, May 24, 2022
#Both the RGLPK and LpSolveAPI version of the program will be shown in this 
#file. 
#The optimal variable inputs needed will be calculated external to the Main LP
#model after saving the peer values from each DEA model run.
########################################################################
rm(list=ls())
####################################################################
####################################################################
library(tidyverse)
library(Rglpk)
library(lpSolveAPI)
####################################################################
#Read data into dataframe df1
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
df1<-read_csv('fishery_data.csv')   #Data file should be in Working Directory
####################################################################
# create input X data frame. Put fixed inputs into X matrix and Variable
#inputs into VX data frame.
X<-as.matrix(df1[,c("FX1","FX2")])      #FX1 and FX2 are Fixed Inputs
VX<-as.matrix(df1[,c("V1")])            #VX is the single Variable Input
#create output Y matrix
Y<-as.matrix(df1[,c("Q1", "Q2","Q3")])  #Q1, Q2, and Q3 are the Outputs. Species of Fish
####################################################################
M=ncol(Y)
N=ncol(X)
J=nrow(X)
YX=cbind(Y,X)
####################################################################
A=matrix(0,M+N+1,J+1)         #M+N+1 sets up for a VRS model
A[1:(M+N),1:J]=t(YX)
A[M+N+1,1:J]=1.0              #1's for VRS constraint
###################################################################
nr=nrow(A) 
nc=ncol(A)
obj=c(rep(0,J),1)                     #objective function coefficients
(rest=c(rep('>=',M),rep('<=',N),'='))
objtype='max'
######################################################################
LP_API=make.lp(nrow=nr,ncol=nc)
lp.control(LP_API,sense=objtype)
set.objfn(LP_API,obj)

for(i in 1:nr){
  set.row(LP_API,i,A[i,])
}

set.constr.type(LP_API,rest)
status2=0
objvals2=0
peers2<-matrix(0,J,J)                #Initiate PEERS matrix
##############################################################################
#Initiate DEA Loop
##############################################################################
for(j in 1:J){
  A[1:M,J+1]=-A[1:M,j]               #last column for outputs in A matrix is set equal to
                                     #the current observation
  rhs=c(rep(0,M),as.matrix(X[j,]),1) #Zero for outputs, row j for from input df 
                                     #and one for VRS
  set.rhs(LP_API,rhs)                
  set.column(LP_API,nc,A[,nc])       #loading revised input data into LPApi matrix
  set.objfn(LP_API,obj)
  
  (status2[j]=solve(LP_API))
  (objvals2[j]=round(get.objective(LP_API),3))
  peers2[j,]<-get.variables(LP_API)[1:J]  
  if(j%%100==0|j==J)  print(paste('on dmu',j,'of',J))
  
}
########################################################################################
#Results
#########################################################################################
e<-dea(X,Y,RTS="vrs",ORIENTATION="out")#Run dea model using Benchmarking package to test
bench<-round(e$objval,3)  
###############################################################################
summary(status2) 
summary(objvals2)
summary(bench)
summary(objvals2-bench)
################################################################################
VI2<-peers2%*%VX   #This calculates the optimal variable input levels.
summary(VI2)
