###############################################################################
#This program uses the OMPR package to set up the DEA output oriented TE model using the 
#Charnes (1981) data from the Benchmarking package
rm(list=ls())
###############################################################################
#First load required packages
library(Benchmarking)
library(ompr)
library(ompr.roi)
library(ROI.plugin.glpk)
library(dplyr)
library(Rglpk)
###############################################################################
#Read in the Charnes data into df1
data(charnes1981)
df1<-charnes1981
#df1<-df1[1:10,]
###############################################################################
J=nrow(df1)                            #Number of observations
X<-df1[,c("x1","x2","x3","x4","x5")]   #Read x variables into X data frame
Y<-df1[,c("y1","y2","y3")]             #Read y variables into Y data frame
m<-colnames(Y)                         #Read column names of Y into m
M=length(m)
n<-colnames(X)                         #Read column names of X into m
N=length(n)
e<-dea(X,Y,RTS="vrs",ORIENTATION="out")#Run dea model using Benchmarking package to test
bench<-round(e$objval,3)               #store results in data structure bench
j=1                                    #set j=1 for initial DEA matrix
################################################################################
###############################################################################
#Set-up the DEA model using OMPR for one observation. This step sets up the A 
#matrix, rhs values and objective function to use in DEA Model
###############################################################################
model1<-MIPModel() %>%
  add_variable(lambda[j],j=1:J) %>%
  add_variable(theta[j]) %>%
  set_objective(sum_over(theta[j], j = 1),"max") %>%
  add_constraint(sum_over(lambda[j]*Y[j,m], j = 1:J) >= theta[j]*Y[j,m], m = 1:M) %>%
  add_constraint(sum_over(lambda[j]*X[j,n], j = 1:J) <= X[j,n], n = 1:N)%>%
  add_constraint(sum_over(lambda[j], j=1:J) ==1)
  
(vnames1 = variable_keys(model1))                      #Extract the name of the variables for this model 
(obj1=as.vector(objective_function(model1)$solution)); #Extract objective function values
names(obj1)=vnames1;                                   #set names for the objective function 
obj1                                                   #echo names. Can be turned off later
constraints1 <- (extract_constraints(model1))          #set constraints
################################################################################
#Set up A matrix and run for multiple observations
################################################################################
A = as.matrix(constraints1$matrix); 
colnames(A)=vnames1; 
(dirs1 = constraints1$sense)       #Signs for Constraints
(rhs1 = constraints1$rhs)          #RHS values for first observation
S=ncol(A)
res=0
status=0
for(s in 1:(S-1)){

  A[1:M,S]=-A[1:M,s]         #replace last column in A for m outputs 
                             #by negative s A column for output oriented model
  
  rhs=c(rep(0,M),A[(M+1):(M+N),s],1)
  sol<- Rglpk_solve_LP(obj=obj1, mat=A, dir=dirs1, rhs=rhs, max=T) #Solve using RGLPK
  status[(s)]=sol$status
  res[(s)]=round(sol$optimum,3)
  if(s%%100==0|s==S)  print(paste('on dmu',s,'of',S-1))
  
}
summary(status)    #See if everything Solved
stopifnot(sum(status!=0)==0) #Hard stop if any entries in status are non-zero
res                #results
summary(res-bench) #comparing results against those from the Benchmark program
###############################################################################

