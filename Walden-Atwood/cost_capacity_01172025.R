#####################################################################
#Program used to estimate CU using Ray Average Cost.
#This is based on the manuscript "Economic measures of capacity utilization: A 
#nonparametric short-run cost function analysis" by Subhash Ray, John Walden and
#Lei Chen. European Journal of Operational Research 293 (2021) 375-387
#Modelis based on Equation #31.
# Author: John Walden NOAA/NMFS/NEFSC
#July 13, 2021
#Revised August 24, 2023.
#Revised July 16, 2024
###################################################################
####################################################################
rm(list=ls())
graphics.off()
###############################################################################
library(Rglpk)
library(doBy)
library(readr)
###############################################################################
#Read in Data
df=read_csv("./data/spec_cap_distort_072024.csv", show_col_types = FALSE) #Read in data
###############################################################################
Y<-df[,c("sp1_lbs","sp2_lbs","sp3_lbs","sp4_lbs")]
V<-as.matrix(df[,c("tda")])
K<-as.matrix(df[,c("LEN", "VHP")])
#############################################################################
knames<-colnames(K)    
m<-colnames(Y)         #Read column names of Y into m
M=length(m)
nv<-colnames(V)        #Read column names of X into m
NV=length(nv)
nf<-colnames(K)        #Read column names of K into nf
NF=length(nf)
###########################################################################
L=ncol(V)
N=ncol(K)
YX=cbind(Y,V,K)
J=nrow(YX)

vnames<-matrix("NA",1,J)
for(j in 1:J){vnames[1,j]<-paste0("mu",j)}

################################################################################
#Define A Matrix 
A=matrix(0,M+L+N+1+N,J+L+N+1)
#Next, Transform YX matrix and copy to A. This is a VRS model 
A[1:(M+L+N),1:J]=t(YX)
A[(M+L),(J+1)]=-1   #VX 1
A[(M+L+(N-1)),(J+2)]=-1   #FX 1
A[(M+L+N),J+3]=-1   #FX 2
A[(M+L+N+1),(1:J)]=1
A[(M+L+N+1),(J+L+N+1)]=-1
A[(M+L+N+2),(J+2)]=1
A[(M+L+N+2),(J+L+N+1)]=as.numeric(-df[j,"LEN"])
A[(M+L+N+3),(J+3)]=1
A[(M+L+N+3),(J+L+N+1)]=as.numeric(-df[j,"VHP"])
colnames(A)<-c(vnames,nv,knames,"sigma")
##############################################################################
######################################################################
#This next code segment calculates the cost model 
#LP Solver Rglpk
#########################################################################
(dir=c(rep('>=',M),rep('<=',L),rep('<=',N),rep('==',3)))
max<-FALSE
objvals1=0
cu=0
t=0
status1=0
peers<-matrix(0,J,J) 
opt_fx=matrix(0,J,2)
opt_vx=0
for(j in 1:J){
  vcost=df[j,"costda"]
  lcost=df[j,"kpcnt2"]*1149           #length, mean value=1149 per foot
  hpcost=df[j,"kpcnt2"]*51            #hp, mean value=51 per hp
  obj=c(rep(0,J),vcost,lcost,hpcost,0)
  rhs=c(as.matrix(Y[j,]),rep(0,L),rep(0,N),rep(0,3)) 
  A[(M+L+N+2),(J+L+N+1)]=as.numeric(-df[j,"LEN"])
  A[(M+L+N+3),(J+L+N+1)]=as.numeric(-df[j,"VHP"])
  sol<- Rglpk_solve_LP(obj, A, dir, rhs, max=max)
  if(j%%100==0|j==J)  print(paste('on dmu',j,'of',J))
  objvals1[j]=round(sol$optimum,4)
  cu[j]=round(sol$solution[J+4],3)
  t[j]=1/cu[j]
  peers[j,]<-sol$solution[1:(J)]
  status1[j]=sol$status
  opt_vx[j]<-sol$solution[J+L]
  opt_fx[j,]<-sol$solution[(J+L+1):(J+L+N)]
}
###############################################################################
cu=as.data.frame(cu)
t<-as.data.frame(t)
opt_vx<-as.data.frame(round(opt_vx,2))
colnames(opt_vx)[1]<-"opt_da"
opt_fx<-as.data.frame(round(opt_fx,2))
colnames(opt_fx)<-c("OPT_LEN","OPT_VHP")
###############################################################################
#Store values and claculate observed cost
###############################################################################
final<-as.data.frame(cbind(df,objvals1,cu,t,opt_vx,opt_fx))

final$obs_cost=(final$tda*final$costda)+(final$LEN*1149*final$kpcnt2)+
  (final$VHP*51*final$kpcnt2)
###############################################################################
