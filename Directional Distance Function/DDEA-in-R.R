##########################################################################
#R Program to calculate directional distance function model
#Model is taken from "New Directions: Efficiency and Productivitys" (2004) 
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
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
df0=read.csv("data1.csv")
###################################################################
#df0=df0[1:10,]
# create input X matrix
X=df0[,c("FX1","FX2","VX1","VX2")]
# create output Y matrix
Y=as.matrix(df0[,c("Y1")])
###################################################################
(M=ncol(Y))
(N=ncol(X))
(J=nrow(X))
YX=cbind(Y,X)
###################################################################
#End of DATA Step
####################################################################
###############################################################################
#Set-up the DEA model using OMPR for one observation. This step sets up the A 
#matrix, rhs values and objective function to use in DEA Model
###############################################################################
j=1
model1<-MIPModel() %>%
  add_variable(lambda[j],j=1:J) %>%
  add_variable(Beta[j]) %>%
  set_objective(sum_over(Beta[j], j = 1),"max") %>%
  add_constraint(sum_over(lambda[j]*Y[j,m], j = 1:J) >= (1+Beta[j])*Y[j,m], m = 1:M) %>%
  add_constraint(sum_over(lambda[j]*X[j,n], j = 1:J) <= X[j,n], n = 1:N)%>%
  add_constraint(sum_over(lambda[j], j=1:J) ==1)

(vnames1 = variable_keys(model1))                      #Extract the name of the variables for this model 
(obj1=as.vector(objective_function(model1)$solution)); #Extract objective function values
names(obj1)=vnames1;                                   #set names for the objective function 
obj1                                                   #echo names. Can be turned off later
constraints1 <- (extract_constraints(model1))          #set constraints
#########################################################################################
A = as.matrix(constraints1$matrix); 
colnames(A)=vnames1; 
(dirs1 = constraints1$sense)       #Signs for Constraints
(rhs1 = constraints1$rhs)          #RHS values for first observation
S=ncol(A)
###############################################################################
#Model #1 - Expand Outputs and Contract Inputs
#Directional Vector is set equal to the Observed values
#This is equivalent to the radial output expansion
################################################################################
res=0
status=0
for(s in 1:(S-1)){
  
 (A[(1:M),S]=-A[(1:M),s])       #replace last column in A for m outputs 
                                #by negative s A column for output oriented model
 (rhs=A[,s])                    #RHS constraints. For DDF Model, A
  sol<- Rglpk_solve_LP(obj=obj1, mat=A, dir=dirs1, rhs=rhs, max=T) #Solve using RGLPK
  status[(s)]=sol$status
  res[(s)]=1+round(sol$optimum,3)
  if(s%%100==0|s==S)  print(paste('on dmu',s,'of',S-1))
  
}
summary(status)
summary(res)
#################################################################################
#Use Benchmarking Program to check results
#################################################################################
e<-dea(X,Y,RTS="vrs",ORIENTATION="out")
bench<-round(e$objval,3)
summary(bench)
#################################################################################



















A=matrix(0,M+N+1,J+1)
A[1:(M+N),1:J]=t(YX)
A[(M+N+1),1:J]=rep(1,J)
##End of A Matrix build#############################################
obj=c(rep(0,J),1)
dir=c(rep('>=',M),rep('<=',N),'==')
res=0
status=0
for(j in (1:J)){
  (dy=as.matrix(Y[j,]))
  (dx=as.matrix(X[j,]))
  (A[1:(M+N),J+1]=c(-dy,dx))
  
  (rhs=c(as.matrix(YX[j,]),1))
  sol<- Rglpk_solve_LP(obj=obj, mat=A, dir=dir, rhs=rhs, max=T) #Solve using RGLPK
  status[(j)]=sol$status
  res[(j)]=1+round(sol$optimum,3)
  if(j%%100==0|j==J)  print(paste('on dmu',j,'of',J))
  
}
summary(status) 















obj=c(rep(0,J),1)
(rest=c(rep('>=',M),rep('<=',N),'=='))



(nr=nrow(A)); 
(nc=ncol(A))
LP_API=make.lp(nrow=nr,ncol=nc)
lp.control(LP_API,sense='max')
set.objfn(LP_API,obj)

for(i in 1:nr){
  set.row(LP_API,i,A[i,])
}

set.constr.type(LP_API,rest)
#####################################################################
# run loop
objvals1=0
status1=0

for(j in 1:J){
  (dy=as.matrix(Y[j,]))
  (dx=as.matrix(X[j,]))
  (A[1:(M+N),J+1]=c(-dy,dx))
 
  (rhs=c(as.matrix(YX[j,]),1))
  set.rhs(LP_API,rhs) 
  set.column(LP_API,nc,A[,nc])
  set.objfn(LP_API,obj)
    
  (status1[j]=solve(LP_API))
  (objvals1[j]=get.objective(LP_API))

  if(j%%100==0|j==J)  print(paste('on dmu',j,'of',J))
    
} # end loop   for(j in 1;J)

print(objvals1)

write.csv(objvals1,"~/John/joe/manual/res_out_dd.csv")
