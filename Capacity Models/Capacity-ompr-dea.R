###############################################################################
#DEA Capacity Model 5-31-2022
#This program estimates the Johansen capacity model as found in the text
#Production Frontiers By Fare, Grosskopf and Lovell (1994)
#Capacity is modelled for fishing vessels in the Northwest Atlantic
#Author: John B. Walden, May 24, 2022
#This version uses the ompr package to set up the DEA program, and then the 
#the RGLPK solver.
#The optimal variable inputs needed will be calculated external to the Main LP
#model after saving the peer values from each DEA model run.
#This program uses the OMPR package to set up the Johansen Capacity Model
##############################################################################
rm(list=ls())
#First load required packages
library(ompr)
library(dplyr)
library(tidyverse)              
library(Rglpk)
library(Benchmarking)
###############################################################################
#Read the data into df1
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
df1<-read_csv('fishery_data.csv', show_col_types = FALSE)
####################################################################
#df1<-df1[(1:10),]
# create input X Matrix. Put fixed inputs into X matrix and Variable
#inputs into VX Matrix.
X<-as.matrix(df1[,c("FX1","FX2")])      #FX1 and FX2 are Fixed Inputs
VX<-as.matrix(df1[,c("V1")])            #VX is the single Variable Input
#create output Y matrix
Y<-as.matrix(df1[,c("Q1", "Q2","Q3")])  #Q1, Q2, and Q3 are the Outputs. Species of Fish
#####################################################################
YX=cbind(Y,X)
###############################################################################
J=nrow(Y)                              #Number of observations
m<-colnames(Y)                         #Read column names of Y into m
M=length(m)
n<-colnames(X)                         #Read column names of X into m
N=length(n)
jp=1                                    #set jp=1 for initial DEA matrix
################################################################################
###############################################################################
#Set-up the DEA model using OMPR for one observation. This step sets up the A 
#matrix, rhs values and objective function to use in DEA Model
###############################################################################
model1<-MIPModel() %>%
  add_variable(lambda[j],j=1:J) %>%
  add_variable(theta) %>%
  set_objective(1*theta,"max") %>%
  add_constraint(sum_over(lambda[j]*Y[j,m], j = 1:J) >= theta*Y[jp,m], m = 1:M) %>%
  add_constraint(sum_over(lambda[j]*X[j,n], j = 1:J) <= X[jp,n], n = 1:N)%>%
  add_constraint(sum_over(lambda[j], j=1:J) ==1)

vnames1 = variable_keys(model1)                      #Extract the name of the variables for this model 
obj1=as.vector(objective_function(model1)$solution)  #Extract objective function values
names(obj1)=vnames1;                                 #set names for the objective function 
constraints1 <- (extract_constraints(model1))        #set constraints
################################################################################
#Set up A matrix and run for multiple observations
################################################################################
A = as.matrix(constraints1$matrix); 
colnames(A)=vnames1; 
(dirs1 = constraints1$sense)       #Signs for Constraints
(rhs1 = constraints1$rhs)          #RHS values for first observation
S=ncol(A)
##############################################################################
#Initialize PEER Matrix
peers<-matrix(0,(S-1),(S-1))  #S is the number of columns in A matrix, which 
                              #equals J+1. Need Peers for J observations
##############################################################################
res=0
status=0
##############################################################################
#Initiate DEA Loop
##############################################################################
for(s in 1:(S-1)){
  
  A[1:M,S]=-A[1:M,s]        #replace last column in A for m outputs 
                            #by negative j column for output oriented 
                            #model
  
  rhs=c(rep(0,M),A[(M+1):(M+N),s],1)   #Zero for outputs, row s for input 
                                       #and one for VRS
  
  sol<- Rglpk_solve_LP(obj=obj1, mat=A, dir=dirs1, rhs=rhs, max=T) 
  #Solve using RGLPK
  
  status[(s)]=sol$status
  res[(s)]=round(sol$optimum,3)
  peers[s,]<-sol$solution[1:(S-1)]
  if(s%%100==0|s==S)  print(paste('on dmu',s,'of',S-1))
}
###############################################################################
VI<-peers%*%VX   #This calculates the optimal variable input levels.
###############################################################################
e<-dea(X,Y,RTS="vrs",ORIENTATION="out")#Run dea model using Benchmarking package to test
bench<-round(e$objval,3)  
###############################################################################
summary(status)    #See if everything Solved
summary(res)       #results
summary(bench)     #Results from Benchmark package
summary(res-bench)
###############################################################################
summary(VI)
###############################################################################

