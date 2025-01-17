#DEA Model for Johansen Capacity and TE
#Author: John B. Walden, NOAA/NMFS/NEFSC
#Modified: 17-Jan-2025
############################################################################
rm(list=ls())
library(Rglpk)
library(dplyr)
library(data.table)
library(readr)
library(Benchmarking)
#####################################################################################
#Read in Data for DEA Model
#####################################################################################
df=read_csv("./data/spec_cap_distort_072024.csv", show_col_types = FALSE) #Read in data
#############################################################################################
#DEA MODEL Set Up
###########################################################################################
XF=as.matrix(df[,c("LEN","VHP")])
XV=as.matrix(df[,c("tda")])
# create output and input Y and X matrices
Y<-as.matrix(df[,c("sp1_lbs", "sp2_lbs","sp3_lbs","sp4_lbs")])
X<-XF
M=ncol(Y)
N=ncol(X)
J=nrow(X)
YX<-cbind(Y,X)
##########################################################################################
e<-dea(XF,Y,RTS="vrs", ORIENTATION = "out") #Using Benchmarking package
CAP<-as.data.frame(round(e$objval,3))       #Extract objective function value
colnames(CAP)<-"CAP"
CU<-as.data.frame(1/CAP)                    #Calculate CU
colnames(CU)<-"CU"
peers<-as.matrix(e$lambda)                  #Extract Peer Values
opt_da<-as.data.frame(round(peers%*%XV,2))  #Calculate Optimal Days Absent
colnames(opt_da)<-"OPT_DA"
final<-cbind(df$Obs,CAP,CU,opt_da)
colnames(final)[1]<-"Obs"
###########################################################################################