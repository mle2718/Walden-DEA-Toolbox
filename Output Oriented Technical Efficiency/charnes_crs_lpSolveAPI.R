##########################################################################
#R Program to calculate output oriented TE Model
#Model is taken from "Production Frontiers" (1994) by Fare, Grosskopf and Lovell
#This version calculates a CRS model
####################################################################
#First Clear any previous data stored in memory, and require lpSolveAPI
rm(list=ls())
library(lpSolveAPI)
library(Benchmarking)
######################################################
#Beginning of Data Step
####################################################################
data(charnes1981)
df1<-charnes1981
names(df1)        
###################################################################
# create input X matrix
X<-df1[,c("x1","x2","x3","x4","x5")]
# create output Y matrix
Y<-df1[,c("y1","y2","y3")]
###################################################################
(M=ncol(Y))
(N=ncol(X))
(J=nrow(X))
YX=cbind(Y,X)
###################################################################
#End of DATA Step
###################################################################
#Define A Matrix. M+N rows allows for CRS. 
A=matrix(0,M+N,J+1)
#Next, Transform YX matrix and copy to A. 
A[1:(M+N),1:J]=t(YX)
#################################################################
objtype='max'
#set zero for all J colums, plus 1 for last column
obj=c(rep(0,J),1)
rest=c(rep('>=',M),rep('<=',N))
#################################################################
(nr=nrow(A)) 
(nc=ncol(A))
#############################################################################
#Set-UP LP Problem
LP_API=make.lp(nrow=nr,ncol=nc)
lp.control(LP_API,sense=objtype)

set.objfn(LP_API,obj)
set.constr.type(LP_API,rest)

for(i in 1:nr){
  set.row(LP_API,i,A[i,])   #setting up rows by reading in A matrix rows
}

objvals2=0                  #Vector to hold results
status2=0                   #Vector to hold model solve status. 0=Model solved
##############################################################################
#Start of DEA loop
##############################################################################
for(j in 1:J){
  A[1:M,J+1]=-A[1:M,j]             #replace last column in A for m outputs 
                                   #by negative j A column 
                                   #for output oriented model
 
  rhs=c(rep(0,M),as.matrix(X[j,])) #rhs is being set to 0 for outputs, observation
                                   #j input data from X
  set.rhs(LP_API,rhs) 
  set.column(LP_API,nc,A[,nc])     #loading revised input data into LPApi
  set.objfn(LP_API,obj)
  
  status2[j]=solve(LP_API)         #Solve Model
  
  objvals2[j]=round(get.objective(LP_API),3)
  
  if(j%%100==0|j==J)  print(paste('on dmu',j,'of',J))
  
} # end of loop  for(j in 1:J)
#################################################################################
#Use Benchmark to test
e<-dea(X,Y,RTS="crs",ORIENTATION="out")#Run dea model using Benchmarking package to test
bench<-round(e$objval,3)               #store results in data structure bench
################################################################################
summary(status2)
summary(objvals2)
summary(bench)
summary(objvals2-bench)

