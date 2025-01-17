################################################################################
#DEA Model for Revenue Capacity
#Author: John Walden, NOAA/NMFS/NEFSC
#Date Modified 16-Jan-2025
################################################################################
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
df$sp1_price=(df$sp1_val/df$sp1_lbs)
df$sp2_price=(df$sp2_val/df$sp2_lbs)
df$sp3_price=(df$sp3_val/df$sp3_lbs)
df$sp4_price=(df$sp4_val/df$sp4_lbs)

prices<-as.matrix(cbind(df$sp1_price,df$sp2_price,df$sp3_price,df$sp4_price))
mean_prices<-round(colMeans(prices, na.rm = TRUE),2)
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
###################################################################################
#put mean prices in a matrix to use in the revenue model
p2<-as.matrix(t(mean_prices))
colnames(p2)<-c("sp1","sp2","sp3","sp4")
pnames<-colnames(p2)
###################################################################
#Use Benchmarking Package to get revenue max
##################################################################
r<-revenue.opt(X,Y,p2,RTS="VRS")
peers<-as.matrix(r$lambda)
opt_vx<-peers%*%XV
opt_y<-r$yopt
rev<-r$rev
df$trev=(df$sp1_lbs*p2[1])+(df$sp2_lbs*p2[2])+(df$sp3_lbs*p2[3])+(df$sp4_lbs*p2[4])
final<-as.data.frame(cbind(df$Obs,df$trev,rev,opt_vx))

colnames(final)<-c("obs","obs_rev","opt_rev","opt_da")
final$rev_cu=final$obs_rev/final$opt_rev
#####################################################################################

