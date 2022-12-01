##########################################################################
#R Program to calculate and bootstrap output oriented TE Model
#Second part of the program provides bootstrap results
#TE Model is taken from "Production Frontiers" (1994) by Fare, Grosskopf and Lovell
#This version calculates a VRS model
####################################################################
#First Clear any previous data stored in memory, and require lpSolveAPI and readr
rm(list=ls())
library(lpSolveAPI)
library(readr)
library(plyr)
####################################################################
#Beginning of Data Step
####################################################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
df0=read_csv("data1000.csv")
###################################################################
# create input X matrix
X<-as.matrix(df0[,c("X1","X2","X3","X4")])
# create output Y matrix
Y<-as.matrix(df0[,c("Y1")])
###################################################################
(M=ncol(Y))
(N=ncol(X))
(J=nrow(X))
YX=cbind(Y,X)
###################################################################
#End of DATA Step
###################################################################
#Define A Matrix. M+N+1 rows allows for VRS. 
A=matrix(0,M+N+1,J+1)
#Next, Transform YX matrix and copy to A. 
A[1:(M+N),1:J]=t(YX)
#Now, set the last row in the A matrix equal to 1, for all J values
#This is for VRS.
A[M+N+1,1:J]=1.0
#################################################################
#Set up and solve initial DEA Program
##################################################################
objtype='max'
#set zero for all J colums, plus 1 for last colums
obj=c(rep(0,J),1)
(rest=c(rep('>=',M),rep('<=',N),'='))
#################################################################
(nr=nrow(A)) 
(nc=ncol(A))
LP_API=make.lp(nrow=nr,ncol=nc)
lp.control(LP_API,sense=objtype)

set.objfn(LP_API,obj)
set.constr.type(LP_API,rest)

for(i in 1:nr){
  set.row(LP_API,i,A[i,])   #setting up rows by reading in A matrix rows
}

objvals2=0
status2=0
for(j in 1:J){
  A[1:M,J+1]=-A[1:M,j]
  (rhs=c(rep(0,M),as.matrix(X[j,]),1)) #rhs is being set to obs j output data, 
                                       #and J input data from X0, and 1 for VRS
  set.rhs(LP_API,rhs) 
  set.column(LP_API,nc,A[,nc]) #loading revised input data into LPApi matrix
  set.objfn(LP_API,obj)

  (status2[j]=solve(LP_API))
  (objvals2[j]=get.objective(LP_API))

  if(j%%100==0|j==J)  print(paste('on dmu',j,'of',J))
    
} # end loop   for(j in 1:J)

delete.lp(LP_API)
#################################################################
# bootstrap setup and solve using lpSolveAPI
#################################################################
set.seed(1001)
(n=J)
(m=round(sqrt(n),0))  #m is the subsample size
B=250
############################
# B=number of bootstraps
############################
statusB=matrix(NA,B,n)
boot=matrix(NA,B,n)
###########################
AB=matrix(0,M+N+1,m+1)
objB=c(rep(0,m),1)
###########################
(ncB=ncol(AB))
(nrB=nrow(AB))
LP_B=make.lp(nrow=nrB,ncol=ncB)
lp.control(LP_B,sense=objtype)
set.objfn(LP_B,objB)
set.constr.type(LP_B,rest)

for(b in 1:B){
    pickval=sample(1:J,m)
    AB[,1:m]=A[,pickval]
   
    for(i in 1:nrB){
     set.row(LP_B,i,AB[i,])
    }

   for(s in 1:J){
      (AB[,1]=A[,s])            #First column of Matrix AB is observation s From Matrix A
       AB[1:M,m+1]=-AB[1:M,1]   #for output oriented model. 
       #Last column outputs in AB = negative first column outputs in AB.
      
       rhsB=c(rep(0,M),(AB[(M+1):(N+1),1]),1)
    
       set.rhs(LP_B,rhsB) 
       set.column(LP_B,1,AB[,1])
       set.column(LP_B,ncB,AB[,ncB])
       set.objfn(LP_B,objB)
    
       statusB[b,s]=solve(LP_B)
       boot[b,s]=get.objective(LP_B)
   }# end loop for (s in 1:J)
#The next line prints to the screeen which iteration of the bootstrap
#is being calculated
  if(b%%1==0|b==B)  print(paste('on bootstrap rep',b,'of',B))  
} 
######################## end loop on b####################################################
summary(rowSums(statusB))  #test to make sure all bootstrap replications solved
                           #all values should be zero, or there was a problem
#Next part is for the bias-correction and confidence intervals
#Need the RTS used in the DEA model to do correct calculation
RTS='VRS'
# construct the "S" matrix
if(RTS=='CRS'|RTS=='crs') (beta=2/(M+N))
if(RTS=='VRS'|RTS=='Vrs') (beta=2/(M+N+1))
S=matrix(0,B,n)
  
for(s in 1:J) S[,s]=objvals2[s]-(((m/n)^beta)*(boot[,s]-objvals2[s]))

round(mean(objvals2),3) #Calculates mean TE score for initial DEA Model
round(mean(S),3)        #Calculates mean TE score across all bootstrap runs

round(quantile(S,0.025),3) #Calculates lower 2.5% tail of bootstrap scores
round(quantile(S, 0.5),3)  #Calculates median of bootstrap scores
round(quantile(S, 0.975),3) #calculates upper 2.5% tail of bootstrap scores
######################################################################################
#This next section calculates the mean bias adjusted TE score and the 5% confidence
#interval for each DMU
####################################################################################
R=matrix(0,J,4)
for (j in 1:J) {R[j,1]=objvals2[j]}
for (j in 1:J) {R[j,2]=round(mean(S[,j]),3)}
for (j in 1:J) {R[j,3]=round(quantile(S[,j],0.025),3)}
for (j in 1:J) {R[j,4]=round(quantile(S[,j],0.975),3)}

colnames(R)<-c("OTE","BCTE","OTE0.025","OTE0.975")
#####################################################################################