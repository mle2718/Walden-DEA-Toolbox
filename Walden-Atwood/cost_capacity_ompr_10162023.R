#####################################################################
#Program used to estimate CU using RAC and the OMPR R package.
#This is based on the manuscript "Economic measures of capacity utilization: A 
#nonparametric short-run cost function analysis" by Subhash Ray, John Walden and
#Lei Chen. European Journal of Operational Research 293 (2021) 375-387
#Modelis based on Equation #31.
# Author: John Walden, NOAA/NMFS/NEFSC
#Woods Hole, MA 02543
#July 13, 2021
#Revised August 24, 2023.
#Revised July 16, 2024.
#Revised January 17, 2025
###################################################################
####################################################################
rm(list=ls())
graphics.off()
###############################################################################
library(Rglpk)
library(doBy)
library(readr)
library(ompr)
###############################################################################
# Read in data and place in Matrices 
df=read_csv("./data/spec_cap_distort_072024.csv", show_col_types = FALSE) #Read in data

Y<-as.matrix(df[,c("sp1_lbs","sp2_lbs","sp3_lbs","sp4_lbs")])
V<-as.matrix(df[,c("tda")])
K<-as.matrix(df[,c("LEN", "VHP")])

YX=cbind(Y,V,K)
############################################################################
J=nrow(Y)                              #Number of observations
m<-colnames(Y)                         #Read column names of Y into m
M=length(m)
nv<-colnames(V)                        #Read column names of V into nv
NV=length(nv)
nf<-colnames(K)                        #Read column names of K into nf
NF=length(nf)
###########################################################################
j=1
w=0
w=as.numeric(df[j,"costda"])
r<-matrix(0,1,NF)
colnames(r)<-nf
r[1,"LEN"]=as.numeric(1149*df[1,"kpcnt2"])           #length, mean value=1149 per foot
r[1,"VHP"]=as.numeric(51*df[1,"kpcnt2"])             #hp, mean value=51 per hp
jp=1
################################################################################
#Set-up the DEA model using OMPR for one observation. This step sets up the A 
#matrix, rhs values and objective function to use in DEA Model #31 from  
#Ray, walden and Chen EJOR article
###############################################################################
model1<-MIPModel() %>%
  add_variable(mu[j], j=1:J) %>%
  add_variable(q[nv], nv=1:NV) %>%
  add_variable(z[n], n=1:NF) %>%
  add_variable(sigma)%>%
  set_objective(sum_over(w[nv]*q[nv], nv=1:NV) + sum_over(r[k]*z[k], k=1:NF),"min")%>%
  add_constraint(sum_over(mu[j]*Y[j,m], j = 1:J) >= Y[jp,m], m = 1:M)%>%
  add_constraint(sum_over(mu[j]*V[j,nv], j = 1:J) <= q[nv], nv=1:NV) %>%
  add_constraint(sum_over(mu[j]*K[j,nf], j = 1:J) <= z[nf], nf = 1:NF)%>%
  add_constraint(sum_over(mu[j], j=1:J) == sigma) %>%
  add_constraint(z[nf]== sigma*K[jp,nf], nf = 1:NF)

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
#################################################################################
max<-FALSE
objvals1=0
cu=0
status1=0
peers<-matrix(0,J,J) 
opt_fx=matrix(0,J,2)
opt_vx=0
t=0
for(j in 1:J){
  obj1["q[1]"]=as.numeric(df[j,"costda"])
  obj1["z[1]"]=as.numeric(1149*df[j,"kpcnt2"]) 
  obj1["z[2]"]=as.numeric(51*df[j,"kpcnt2"])  
  
  rhs1=c(as.matrix(Y[j,]),rep(0,NV),rep(0,NF),0,rep(0,NF)) 
  A[(M+NV+NF+2):(M+NF+NV+3),"sigma"]=as.numeric(-df[j,nf])
  sol<- Rglpk_solve_LP(obj1, A, dirs1, rhs1, max=max)
  if(j%%100==0|j==J)  print(paste('on dmu',j,'of',J))
  objvals1[j]=round(sol$optimum,4)
  cu[j]=round(sol$solution[J+4],3)
  t[j]=round(1/cu[j],2)
  peers[j,]<-sol$solution[1:(J)]
  status1[j]=sol$status
  opt_vx[j]<-sol$solution[J+NV]
  opt_fx[j,]<-sol$solution[(J+NV+1):(J+NV+2)]
}

opt_vx<-round(opt_vx,2)
opt_fx<-as.data.frame(round(opt_fx,2))
colnames(opt_fx)<-c("opt_len","opt_vhp")

final<-as.data.frame(cbind(df,objvals1,cu,t,opt_vx,opt_fx))

final$obs_cost=(final$tda*final$costda)+(final$LEN*1149*final$kpcnt2)+
  (final$VHP*51*final$kpcnt2)
