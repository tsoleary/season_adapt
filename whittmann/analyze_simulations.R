#!/usr/bin/Rscript --vanilla
args <- commandArgs(TRUE)

source("analysis_funcs.R")
library(foreach)

numreps <- 10

#parameters
inputfilenamebase <- "input"
outputfilenamebase <- "output"
allinfofilenamebase <- "allinfo"
fixedfilenamebase <- "fixedfile"
resultfilebase <- "results"

param1 <- 1
N <- 1000
nl <- 100
d <- 0.7
cap <- 2*N

ytreatments <- c(0.5,1,2,4)
gtreatments <-  c(1,seq(2,20,2))
Y <- rep(ytreatments,rep(length(gtreatments),length(ytreatments)))
G <- rep(gtreatments,length(ytreatments))

gridvec <- seq(0,1,0.04)
mutrate <- 10^(-4)
numscenarios <- length(Y)

stabcoefs <- foreach (treatment=1:numscenarios,.combine='cbind') %dopar%{
    wintergenerations <- G[treatment]
    summergenerations <- G[treatment]
    #if (wintergenerations>1){
        cyclelength <- wintergenerations + summergenerations
        numyears <- 1000
        numgen <- cyclelength*numyears
        D <- read.table(paste("output",treatment,"seasonalmutationcounts.txt",sep="_"),header=FALSE)
        res <- stabilitycoefficient(D,N,wintergenerations,summergenerations,mutrate)
        coef <- res[[1]]
        return(c(mean(coef),sd(coef),res[[2]]))
    #}
    #else {
    #    return(rep(NA,length(gridvec)+2))
    #}
}
effective.s  <- foreach (treatment=1:numscenarios,.combine='cbind') %dopar%{
    wintergenerations <- G[treatment]
    summergenerations <- G[treatment]
    #if (wintergenerations>1){
        cyclelength <- wintergenerations + summergenerations
        numyears <- 1000
        numgen <- cyclelength*numyears
        D <- read.table(paste("output",treatment,"seasonalmutationcounts.txt",sep="_"),header=FALSE)
        res <- effectives(D,N,wintergenerations,summergenerations)
        coef <- res[[1]]
        means <- res[[2]]
        return(c(mean(coef),sd(coef),res[[2]]))
    #}
    #else {
   #     return(rep(NA,length(gridvec)+2))
    #}
}

propexpected <- foreach (treatment=1:numscenarios,.combine='cbind') %dopar%{
    wintergenerations <- G[treatment]
    summergenerations <- G[treatment]
    #if (wintergenerations>1){
        cyclelength <- wintergenerations + summergenerations
        numyears <- 1000
        numgen <- cyclelength*numyears
        D <- read.table(paste("output",treatment,"seasonalmutationcounts.txt",sep="_"),header=FALSE)
        prop <- propexpecteddirection(D,N,wintergenerations,summergenerations)
        return(c(mean(prop),sd(prop)))
    #}
    #else {
    #    return(rep(NA,2))
    #}
}

table14 <- data.frame(simexperiment=rep(14,numscenarios),d=rep(d,numscenarios),c=rep(param1,numscenarios),y=Y,cap=rep(cap,numscenarios),N=rep(N,numscenarios),numloci=rep(nl,numscenarios),g=G,stabilitycoefmean=stabcoefs[1,],stabilitycoefsd=stabcoefs[2,],effectivesmean=effective.s[1,],effectivessd=effective.s[2,],propexpectedmean=propexpected[1,],propexpectedsd=propexpected[2,])
write.table(table14,"resultstable_14.txt",row.names=FALSE,col.names=TRUE,quote=FALSE)

meantable1 <- stabcoefs[3:nrow(stabcoefs),]
meantable2 <- effective.s[3:nrow(effective.s),]
write.table(meantable1,"meantable14_1.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(meantable2,"meantable14_2.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
