##########################################################################
#R Program to calculate and bootstrap output oriented TE model
#Second part of the program provides bootstrap results
#TE Model is taken from "Production Frontiers" (1994) by Fare, Grosskopf and Lovell
#This version calculates a VRS model and uses the RGLPK solver in R
#RGLPK will run slower than LpAPI, but it's syntax is easier to understand
####################################################################
#First Clear any previous data stored in memory, and require 
rm(list=ls())
library(Rglpk)
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
obj=c(rep(0,J),1)
dir=c(rep('>=',M),rep('<=',N),'==')
objvals2=0
status=0
for(j in (1:J)){
  A[1:M,(J+1)]=-A[1:M,j]         #replace last column in A for m outputs 
  #by negative j A column for output oriented model
  rhs=c(rep(0,M),A[(M+1):(M+N),j],1)       #Zero for outpus, column j for inputs 
                                           #and one for VRS
  sol<- Rglpk_solve_LP(obj=obj, mat=A, dir=dir, rhs=rhs, max=T) #Solve using RGLPK
  status[(j)]=sol$status
  objvals2[(j)]=round(sol$optimum,3)       #Hold results in objvals2
  if(j%%100==0|j==J)  print(paste('on dmu',j,'of',J))
  
}
summary(status)                           #Check to make sure everything solved
#################################################################
# Bootstrap set up and solve using RgLPK
#################################################################
set.seed(1001)
(n=J)
(m=round(sqrt(n),0)) #m is the subsample size
B=250
############################
# B=number of bootstraps
############################
statusB=matrix(NA,B,n)  #Holds results showing whether model solved
boot=matrix(NA,B,n)     #Will hold bootstrap results in Matrix of dimension (B,n)
###########################
AB=matrix(0,M+N+1,m+1)
objB=c(rep(0,m),1)
dirB=c(rep('>=',M),rep('<=',N),'==')
###########################

for(b in 1:B){
    pickval=sample(1:J,m)      #A subsample of observations numbers is chosen
    AB[,1:m]=A[,pickval]       #Subsample Data are copied from matrix A to AB

    for(s in 1:J){             
      (AB[,1]=A[,s])            #First column of Matrix AB is observation s From Matrix A
      AB[1:M,m+1]=-AB[1:M,1]    #for output oriented model. 
                                #Last column outputs = negative first column outputs.
      
     rhsB=c(rep(0,M),(AB[(M+1):(N+1),1]),1)
     sol<- Rglpk_solve_LP(obj=objB, mat=AB, dir=dirB, rhs=rhsB, max=T) #Solve using RGLPK
     
     statusB[b,s]=sol$status
     boot[b,s]=round(sol$optimum,3)
      
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
  
for(s in 1:J){
  S[,s]=objvals2[s]-(((m/n)^beta)*(boot[,s]-objvals2[s]))
}

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
 for (j in 1:J) {R[j,1]=round(objvals2[j],3)}
 for (j in 1:J) {R[j,2]=round(mean(S[,j]),3)}
 for (j in 1:J) {R[j,3]=round(quantile(S[,j],0.025),3)}
 for (j in 1:J) {R[j,4]=round(quantile(S[,j],0.975),3)}

colnames(R)<-c("OTE","BCTE","OTE0.025","OTE0.975")


