##########################################################################
#R Program to calculate directional distance function model
#Model is taken from "New Directions: Efficiency and Productivitys" (2004) 
#by Fare and Grosskopf.
#This program expands good outputs and contracts undesirable outputs using
#a directional distance function. Inputs are held constant in this model
#The program will build Matrices and solve with the Rglpk solver.
#Author: John B. Walden
#Date: June 14,2002
####################################################################
#First Clear any previous data stored in memory, and require lpSolveAPI and readr
rm(list=ls())
library(Rglpk)
library(tidyverse)
####################################################################
#Beginning of Data Step
####################################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#Set working directory to current directory. Data should be stored here
df1=read_csv("ddf_example.csv", show_col_types = FALSE)
###################################################################
#df1=df1[1:10,]     #For setting up small problem to see if code
                    #manipulates matrices correctly
###################################################################
#create input X matrix
X<-as.matrix(df1[,c("FX1","FX2","VX1","VX2")])
# create output Y matrix
YD<-as.matrix(df1[,c("Y1","Y2")])       #Desirable Outputs
YU<-as.matrix(df1[,"U1"])               #Undesirable Outputs
Y<-cbind(YD,YU)
YX<-cbind(Y,X)
###############################################################################
###############################################################################
J=nrow(df1)    #Number of observations
M=ncol(Y)      #Number of outputs
N=ncol(X)      #Number of Inputs                              
################################################################################
################################################################################
#Define A Matrix and constraint set
################################################################################
A<-matrix(0,(M+N+1),J+1)
A[1:(M+N),(1:J)]<-t(YX)
A[(M+N+1),(1:J)]=1
################################################################################
#End of DATA Step and setting A matrix
################################################################################
#Solve using Rglpk
##############################################################################
S=ncol(A)
res=0
status=0
obj1=c(rep(0,J),1)
dir=c(rep('>=',M),rep('<=',N),'==')
for(s in 1:(S-1)){
  
  #####################################################################  
  #set directional vectors and redefine A matix for current observation
  #################################################################### 
  dyd=as.matrix(YD[s,])
  dyu=as.matrix(YU[s,])
  A[(1:M),S]=c(-dyd,dyu)   #put direction vector in last column
                           #In this example the directional vector=observed
  ####################################################################
  #set RHS
  rhs=c(as.matrix(YD[s,]),as.matrix(YU[s,]),as.matrix(X[s,]),1)
  
  #Solve using RGLPK
  sol<- Rglpk_solve_LP(obj=obj1, mat=A, dir=dir, rhs=rhs, max=T) 
  
  status[(s)]=sol$status
  res[(s)]=round(sol$optimum,3)  
  if(s%%100==0|s==S)  print(paste('on dmu',s,'of',S-1))
}
###############################################################################
#Show results and solver status
summary(status)
summary(res)
#############################################################################
