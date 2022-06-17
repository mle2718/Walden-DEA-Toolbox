##########################################################################
#R Program to calculate directional distance function model
#Model is taken from "New Directions: Efficiency and Productivity" (2004) 
#by Fare and Grosskopf.
#This program uses OMPR to output the initial A matrix. 
#There are three different directional vector problems shown here.
#Model #1 Expands outputs using the observed values for the directional vector
#This is equivalent to the radial output oriented DEA model
#Model #2 Uses the directional vector [1,1] to expand outputs and contract inputs.
#model #3 Uses the directional vector [1,0] to expand outputs while holding inputs
#constant
####################################################################
#First Clear any previous data stored in memory, and require lpSolveAPI and readr
rm(list=ls())
library(Rglpk)
library(ompr)
library(Benchmarking)
####################################################################
#Beginning of Data Step
####################################################################
data(charnes1981)
df1<-charnes1981
#df1<-df1[1:10,]
###################################################################
names(df1)
J=nrow(df1)
#create input X Matrix
X<-as.matrix(df1[,c("x1","x2","x3","x4","x5")])
#create output Y matrix
Y<-as.matrix(df1[,c("y1","y2","y3")])             
###################################################################
m<-colnames(Y)                         #Read column names of Y into m
M=length(m)
n<-colnames(X)                         #Read column names of X into m
N=length(n)
jp=1                                    #set j=1 for initial DEA matrix(M=ncol(Y))
###################################################################
#End of DATA Step
####################################################################
###############################################################################
#Set-up the DEA model using OMPR for one observation. This step sets up the A 
#matrix, rhs values and objective function to use in DEA Model
###############################################################################
model1<-MIPModel() %>%
  add_variable(z[j],j=1:J, lb=0) %>%
  add_variable(Beta, lb=-Inf, ub=Inf) %>%
  set_objective(1*Beta,"max") %>%
  add_constraint(sum_over(z[j]*Y[j,m], j = 1:J) >= (1+Beta)*Y[jp,m], m = 1:M) %>%
  add_constraint(sum_over(z[j]*X[j,n], j = 1:J) <= X[jp,n], n = 1:N)%>%
  add_constraint(sum_over(z[j], j=1:J) ==1)

vnames1 = variable_keys(model1)                      #Extract the name of the variables for this model 
obj1=as.vector(objective_function(model1)$solution); #Extract objective function values
names(obj1)=vnames1;                                   #set names for the objective function 
constraints1 <- (extract_constraints(model1))          #set constraints
#########################################################################################
A = as.matrix(constraints1$matrix); 
colnames(A)=vnames1; 
dirs1 = constraints1$sense       #Signs for Constraints
rhs1 = constraints1$rhs          #RHS values for first observation
S=ncol(A)                        #S is needed to set up correct columns for LP
###############################################################################
#Model #1 - Expand Outputs and Holds Inputs Constants
#Directional Vector is set equal to the Observed values
#This is equivalent to the radial output expansion
################################################################################
res=0
status=0
for(s in 1:(S-1)){
  
  A[(1:M),S]=-A[(1:M),s]       #replace last column in A for m outputs 
                                #by negative s A column for output oriented model
  ####################################################################
  #set RHS
  rhs=c(as.matrix(Y[s,]),as.matrix(X[s,]),1)
  sol<- Rglpk_solve_LP(obj=obj1, mat=A, dir=dirs1, rhs=rhs, max=T) 
  #Solve using RGLPK
  
  status[(s)]=sol$status
  res[(s)]=round(1+sol$optimum,3)  #add one to make equivalent to radial DEA
                                   #and for comparison to Benchmark
  if(s%%100==0|s==S)  print(paste('on dmu',s,'of',S-1))
}
#################################################################################
#Use Benchmarking Program to check results
#################################################################################
e<-dea(X,Y,RTS="vrs",ORIENTATION="out")
bench<-round(e$objval,3)
################################################################################
summary(status)
summary(res)
summary(bench)
summary(res-bench)
#################################################################################
#Model #2 - Expand outputs and Contract Inputs using the directional vector [Y,X]
#First, set up the ompr Model
#################################################################################
jp=1
model2=MIPModel() %>%
  add_variable(z[j], j=1:J, lb = 0)%>% 
  add_variable(Beta, lb = -Inf, ub = Inf) %>%
  set_objective(1*Beta,"max") %>%
  add_constraint(sum_over(Y[j,m] * z[j], j = 1:J) >= (1+Beta)*Y[jp,m] , m = 1:M) %>%
  add_constraint(sum_over(X[j,n] * z[j], j = 1:J) <= (1-Beta)*X[jp,n] , n = 1:N) %>%
  add_constraint(sum_over(1*z[j], j = 1:J) == 1.0)

vnames1 = variable_keys(model2)                      #Extract the name of the variables for this model 
obj1=as.vector(objective_function(model2)$solution); #Extract objective function values
names(obj1)=vnames1;                                   #set names for the objective function 
constraints1 <- (extract_constraints(model2))          #set constraints
#########################################################################################
A = as.matrix(constraints1$matrix); 
colnames(A)=vnames1; 
dirs1 = constraints1$sense       #Signs for Constraints
rhs1 = constraints1$rhs          #RHS values for first observation
S=ncol(A)
res2=0
status2=0
for(s in 1:(S-1)){
  
  A[(1:M),S]=-A[(1:M),s]            #last column in A for m outputs is negative A[,s]
  A[(M+1):(M+N),S]=A[(M+1):(M+N),s] #last column for inputs is A[(M+1):(M+N),s]
 
  rhs=c(as.matrix(Y[s,]),as.matrix(X[s,]),1)       #rhs for ddf model
  sol<- Rglpk_solve_LP(obj=obj1, mat=A, dir=dirs1, rhs=rhs, max=T) 
  #Solve using RGLPK
  
  status2[(s)]=sol$status
  res2[(s)]=round(sol$optimum,3)
  if(s%%100==0|s==S)  print(paste('on dmu',s,'of',S-1))
  
}
summary(status2)
summary(res2)
####################################################################################
#Model #3 - Expand outputs using vector [1,0] 
#################################################################################
jp=1
model3=MIPModel() %>%
  add_variable(z[j], j=1:J, lb = 0)%>% 
  add_variable(Beta, lb = -Inf, ub = Inf) %>%
  set_objective(1*Beta,"max") %>%
  add_constraint(sum_over(Y[j,m] * z[j], j = 1:J) >= Y[jp,m]+Beta , m = 1:M) %>%
  add_constraint(sum_over(X[j,n] * z[j], j = 1:J) <= X[jp,n] , n = 1:N) %>%
  add_constraint(sum_over(1*z[j], j = 1:J) == 1.0)

vnames1 = variable_keys(model3)                      #Extract the name of the variables for this model 
obj1=as.vector(objective_function(model3)$solution); #Extract objective function values
names(obj1)=vnames1;                                   #set names for the objective function 
constraints1 <- (extract_constraints(model3))          #set constraints
#########################################################################################
A = as.matrix(constraints1$matrix); 
colnames(A)=vnames1; 
dirs1 = constraints1$sense       #Signs for Constraints
#########################################################################################
S=ncol(A)
res3=0
status3=0
for(s in (1:(S-1))){

  ##################################################################################
  rhs=c(as.matrix(Y[s,]),as.matrix(X[s,]),1)       #rhs for ddf model
  #Note that the only thing which changes in this model is the RHS
  #last column in A matrix for all outputs will always be -1
  ###################################################################################
  sol<- Rglpk_solve_LP(obj=obj1, mat=A, dir=dirs1, rhs=rhs, max=T) #Solve using RGLPK
  status3[s]=sol$status
  res3[s]=round(sol$optimum,3)
  if(s%%100==0|s==S)  print(paste('on dmu',S,'of',S))
  
}
##########################################################################################
#Check this model against results from Benchmarking Package
#Use (1,1,1) for directional vector with output orientation
z<-rep(1,M)
e<-dea.direct(X,Y,DIRECT=z,RTS="vrs",ORIENTATION="out")
res4=round(e$objval,3)
###########################################################################################
summary(status3) 
summary(res3)
summary(res4)
summary(res3-res4)
##########################################################################################
