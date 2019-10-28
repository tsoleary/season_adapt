#!/usr/bin/Rscript --vanilla
args <- commandArgs(TRUE)

library(foreach)
library(doMC)


numreps <- 10
numcores <- 1
registerDoMC(numcores)

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


numneutralloci <- 10
selcoefsummer<-1
selcoefconst <- 0

foreach (treatment=1:length(Y)) %dopar%{
    wintergenerations <- G[treatment]
    summergenerations <- G[treatment]
    cyclelength <- wintergenerations + summergenerations
    numyears <- 1000
    numgen <- cyclelength*numyears
    inputfilename <- paste(inputfilenamebase,"_",treatment,".txt",sep="")
    sink(file=inputfilename)
    cat(paste(outputfilenamebase,"_",treatment," #filenamebase\n",sep=""))
    cat(paste(treatment,"#seed\n"))
    cat(paste(numreps,"#numreps\n"))
    cat(paste(wintergenerations,"#winter generations\n"))
    cat(paste(summergenerations,"#summer generations\n"))
    cat(paste(N,"#Nmin\n"))
    cat(paste(N,"#Nmax\n"))
    cat(paste(numyears,"#numyears\n"))
    cat("1 #number of seasonally selected traits\n")
    cat(paste(nl,"#number of seasonally selected loci\n"))
    cat("0 #number of constantly selected loci\n")
    cat(paste(numneutralloci,"#number of neutral loci\n"))
    cat("0.0001 #mutation rate at seasonally selected loci\n")
    cat("0 #mutation rate at constantly selected loci\n")
    cat("0.0001 #mutation rate at neutral loci\n")
    cat("0 #dfe\n")
    cat(paste(selcoefsummer,"#seasonal selection coefficient\n"))
    cat(paste(selcoefconst,"#constant selection coefficient\n"))
    cat("p #fitness function\n")
    cat("0 #dde\n")
    cat(paste(d,"#trait dominance\n"))
    cat(paste(param1,"#fitness parameter 1\n"))
    cat(paste(Y[treatment],"#fitness parameter 2\n"))
    cat("m #combine traits\n")
    cat(paste(cap,"#offspring number cap\n"))
    cat("2 #initialization mode\n")
    cat("1 #first sampling time\n")
    cat("1 #sampling interval\n")
    cat(paste(N,"#samplesize\n"))
    cat(paste(numyears*2,"#perturbation start\n"))
    cat("5 #perturbation run time\n")
    cat("10 #perturbation distance\n")
    cat("10 #perturbation replicates\n")
    cat("0.01 #perturbation frequency\n")
    sink()
    system(paste("./simplesim",inputfilename))
}

