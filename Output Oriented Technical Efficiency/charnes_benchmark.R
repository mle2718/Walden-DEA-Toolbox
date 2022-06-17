#DEA Program to demonstrate using CHARNES data with Benchmark Program
#This program solves an output oriented VRS model
###############################################################################
rm(list=ls()) 
library(Benchmarking)
###############################################################################
data(charnes1981)
df1<-charnes1981
###############################################################################
X<-as.matrix(df1[,c("x1","x2","x3","x4","x5")])
Y<-as.matrix(df1[,c("y1","y2","y3")])
e<-dea(X,Y,RTS="vrs",ORIENTATION="out")
bench<-round(e$objval,3)
summary(bench)
###############################################################################
