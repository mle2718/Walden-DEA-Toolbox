##########################################################################
#R Program to calculate input oriented TE Model
#Model is taken from "Production Frontiers" (1994) by Fare, Grosskopf and Lovell
#This version calculates a CRS model
####################################################################
#First Clear any previous data stored in memory, and requires Rglpk
rm(list=ls())
library(Rglpk)
library(Benchmarking)
######################################################
#Beginning of Data Step
####################################################################
data(charnes1981)
df1<-charnes1981
###################################################################
# create input X matrix
X<-as.matrix(df1[,c("x1","x2","x3","x4","x5")])
# create output Y matrix
Y<-as.matrix(df1[,c("y1","y2","y3")])
###################################################################
(M=ncol(Y))
(N=ncol(X))
(J=nrow(X))
YX=cbind(Y,X)
###################################################################
#End of DATA Step
###################################################################
#Define A Matrix. 
A=matrix(0,M+N,J+1)
#Next, Transform YX matrix and copy to A. 
A[1:(M+N),1:J]=t(YX)
#################################################################
obj=c(rep(0,J),1)
dir=c(rep('>=',M),rep('<=',N))
res=0
status=0
for(j in (1:J)){
  (A[(M+1):(M+N),J+1]=-A[(M+1):(M+N),j]) #replace last column in A for n inputs 
                                         #by negative s A column for input 
                                         #oriented model
  
  (rhs=c(as.matrix(Y[j,]),rep(0,N)))     #rhs is being set to obs j output data, 
                                         #0 for the input data, and 1 for VRS
  sol<- Rglpk_solve_LP(obj=obj, mat=A, dir=dir, rhs=rhs, max=F) #Solve using RGLPK
  status[(j)]=sol$status
  res[(j)]=round(sol$optimum,3)
  if(j%%100==0|j==J)  print(paste('on dmu',j,'of',J))
}
################################################################################
e<-dea(X,Y,RTS="crs",ORIENTATION="in")#Run dea model using Benchmarking package to test
bench<-round(e$objval,3)               #store results in data structure bench
###############################################################################
summary(status) 
summary(res)
summary(bench)
summary(res-bench)
################################################################################
