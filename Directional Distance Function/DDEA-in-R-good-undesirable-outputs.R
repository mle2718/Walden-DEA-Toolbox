##########################################################################
#R Program to calculate directional distance function model
#Model is taken from "New Directions: Efficiency and Productivitys" (2004) 
#by Fare and Grosskopf.
#This program expands good outputs and contracts undesirable outputs using
#a directional distance function. Inputs are held constant in this model
#The first program will use ompr and the Rglpk solver.
#The second program uses the same A matrix and the lpSolveAPI solver
#Author: John B. Walden
#Date: June 14,2002
####################################################################
#First Clear any previous data stored in memory, and require lpSolveAPI and readr
rm(list=ls())
library(lpSolveAPI)
library(Rglpk)
library(ompr)
library(dplyr)
library(tidyverse)
####################################################################
#Beginning of Data Step
####################################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#Set working directory to current directory. Data should be stored here
df1=read_csv("ddf_example.csv", show_col_types = FALSE)
###################################################################
#df1=df1[1:10,]
# create input X matrix
X<-as.matrix(df1[,c("FX1","FX2","VX1","VX2")])
# create output Y matrix
YD<-as.matrix(df1[,c("Y1","Y2")])       #Desirable Outputs
YU<-as.matrix(df1[,"U1"])               #Undesirable Outputs
###############################################################################
###############################################################################
J=nrow(df1)                            #Number of observations
(md<-colnames(YD))                     #Read column names of Y into m
(MD=length(md))
(mu<-colnames(YU))
(MU=length(mu))
(n<-colnames(X))                       #Read column names of X into m
(N=length(n))
 jp=1                                  #set jp=1 for initial DEA problem used
                                       #to construct the A matrix
################################################################################
###############################################################################
#Set-up the DEA model using OMPR for one observation. This step sets up the A 
#matrix, rhs values and objective function to use in DEA Model
###############################################################################
model1<-MIPModel() %>%
  add_variable(z[j],j=1:J, lb=0) %>%
  add_variable(Beta) %>%
  set_objective(1*Beta,"max") %>% 
  add_constraint(sum_over(z[j]*YD[j,md], j = 1:J) >= (1+Beta)*YD[jp,md], md = 1:MD)%>% 
  add_constraint(sum_over(z[j]*YU[j,mu], j = 1:J) >= (1-Beta)*YU[jp,mu], mu = 1:MU) %>%
  add_constraint(sum_over(z[j]*X[j,n], j = 1:J) <= X[jp,n], n = 1:N)%>%
  add_constraint(sum_over(z[j], j=1:J) ==1)

vnames1 = variable_keys(model1)                      #Extract the name of the variables for this model 
obj1=as.vector(objective_function(model1)$solution); #Extract objective function values
names(obj1)=vnames1;                                 #set names for the objective function 
constraints1 <- (extract_constraints(model1))        #set constraints
################################################################################
#Define A Matrix and constraint set
################################################################################
A = as.matrix(constraints1$matrix); 
colnames(A)=vnames1; 
dirs1 = constraints1$sense       #Signs for Constraints
rhs1 = constraints1$rhs          #RHS values for first observation
S=ncol(A) #S is needed to set up correct columns for LP
################################################################################
#End of DATA Step
################################################################################
#Solve using Rglpk
##############################################################################
res=0
status=0
for(s in 1:(S-1)){
  
  #####################################################################  
  #set directional vectors and redefine A matix for current observation
  #################################################################### 
  dyd=as.matrix(YD[s,])
  dyu=as.matrix(YU[s,])
  A[1:(MD+MU),S]=c(-dyd,dyu)
  ###################################################################
  #Set RHS
  rhs=c(as.matrix(YD[s,]),as.matrix(YU[s,]),as.matrix(X[s,]),1)
  #Solve using RGLPK
  sol<- Rglpk_solve_LP(obj=obj1, mat=A, dir=dirs1, rhs=rhs, max=T) 
  
  status[(s)]=sol$status
  res[(s)]=round(sol$optimum,3)  
  if(s%%100==0|s==S)  print(paste('on dmu',s,'of',S-1))
}
summary(status)
summary(res)
###############################################################################
###############################################################################
#Now Solve use LpAPI
#Re-use A matrix
###############################################################################
nr=nrow(A) 
nc=ncol(A)
LP_API=make.lp(nrow=nr,ncol=nc)
lp.control(LP_API,sense='max')
set.objfn(LP_API,obj1)          #obj1 comes from ompr
#Need to redefine the signs for lpSolveAPI because of different "=" vs. "=="
#############################################################################
rest=c(rep('>=',MD),rep('>=',MU), rep('<=',N),'=')
set.constr.type(LP_API,rest)

for(i in 1:nr){
  set.row(LP_API,i,A[i,])
}
#####################################################################
# run loop
objvals1=0
status1=0

for(s in 1:(S-1)){

#####################################################################  
#set directional vectors and redefine A matix for current observation
  (dyd=as.matrix(YD[s,]))
  (dyu=as.matrix(YU[s,]))
  (A[1:(MD+MU),S]=c(-dyd,dyu))
#################################################################### 
 #Next, define RHS  
  (rhs=c(as.matrix(YD[s,]),as.matrix(YU[s,]),as.matrix(X[s,]),1))
  set.rhs(LP_API,rhs) 
  set.column(LP_API,nc,A[,nc])
  set.objfn(LP_API,obj1)
    
  (status1[s]=solve(LP_API))
  (objvals1[s]=round(get.objective(LP_API),3))

  if(s%%100==0|s==S)  print(paste('on dmu',s,'of',S))
    
} # end loop   for(s in 1:S-1)
#############################################################################
#Results
#############################################################################
summary(status1)
summary(objvals1)
############################################################################
#Compare to results from RGlpk as a check
summary(res)
summary(objvals1-res)
############################################################################
