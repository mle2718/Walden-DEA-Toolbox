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
library(tidyverse)
library(Rglpk)
library(lpSolveAPI)
library(Benchmarking)
####################################################################
#Read data into dataframe df1
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
df1<-read_csv('fishery_data.csv', show_col_types = FALSE)   #Data file should be in Working Directory
####################################################################
# create input X Matrix. Put fixed inputs into X matrix and Variable
#inputs into VX Matrix.
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
peers<-matrix(0,J,J)          #create matrix of zero's to hold peer weights
#################################################################
obj=c(rep(0,J),1)                     #objective function coefficients
dir=c(rep('>=',M),rep('<=',N),'==')   #direction for signs on costraints
res=0
status=0
##############################################################################
#Initiate DEA Loop
##############################################################################
for(j in (1:J)){
  A[1:M,(J+1)]=-A[1:M,j]              #replace last column in A for m outputs 
                                      #by negative j column for output oriented 
                                      #model
  rhs=c(rep(0,M),as.matrix(X[j,]),1)  #Zero for outputs, row j for from input df 
                                      #and one for VRS
  
  sol<- Rglpk_solve_LP(obj=obj, mat=A, dir=dir, rhs=rhs, max=T) 
 
  #Solve using RGLPK
  
  peers[j,]<-sol$solution[1:J]        #store peer values for variable input calculation
  status[(j)]=sol$status              #store whether model solves at each iteration
  res[(j)]=round(sol$optimum,3)       #store results rounded to 3 decimal points
  if(j%%100==0|j==J)  print(paste('on dmu',j,'of',J))
}
VI<-peers%*%VX   #This calculates the optimal variable input levels.
###############################################################################
#Results
###############################################################################
e<-dea(X,Y,RTS="vrs",ORIENTATION="out")#Run dea model using Benchmarking package to test
bench<-round(e$objval,3)  
###############################################################################
summary(status) 
summary(res)    #Summary of efficiency Scores
summary(bench)
summary(res-bench)
#############################################################################
summary(VI)     #Summary of Variable Inputs
###############################################################################
