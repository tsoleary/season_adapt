#This file contains code to analyze the simulation results in Wittmann et al. (2017) Segregation lift: How seasonally fluctuating selection can maintain polymorphism at many loci.

stabilitycoefficient <- function(D,N,wintergenerations,summergenerations,mutrate){#input parameters: data frame with allele-frequency trajectories, population size in the middle of winter (Nmin), number of generations per winter, number of generations per summer, also works for just 1 locus
    #determine data dimension and relevant indices
    cyclelength <- wintergenerations+summergenerations
    numreps <- length(unique(D[,1]))
    numgen <- D[nrow(D),2]#in total for a simulation, not per year
    numyears <- numgen/cyclelength #numyears should be an even number!
    coefs <- numeric(numreps)
    gridvec <- seq(0,1,0.04)
    allmeans <- matrix(NA,nrow=numreps,ncol=length(gridvec)-1)
    
    for (r in 1:numreps){
        Drep <- subset(D,D[,1]+1==r)
        #matplot(Drep[,2],Drep[,3],type="l")
        M <- matrix(1:numgen,nrow=numgen/2)

        if (wintergenerations %% 2 == 0){
            indices <- M[seq(1/2*wintergenerations,nrow(M),cyclelength),2]#middle of the cold season as reference point every year, use only data from the second half of each simulation run

                                        #compute frequency changes over one year
            before <- Drep[as.vector(indices)[1:(length(indices)-1)],3:ncol(Drep)]
            afteryear <- Drep[as.vector(indices)[2:length(indices)],3:ncol(Drep)]
        }
        else {
            indices1 <- M[floor(seq(1/2*wintergenerations,nrow(M),cyclelength)),2]
            indices2 <- M[ceiling(seq(1/2*wintergenerations,nrow(M),cyclelength)),2]

                                        #compute frequency changes over one year
            before <- (Drep[as.vector(indices1)[1:(length(indices1)-1)],3:ncol(Drep)]+Drep[as.vector(indices2)[1:(length(indices2)-1)],3:ncol(Drep)])/2
            afteryear <- (Drep[as.vector(indices1)[2:length(indices1)],3:ncol(Drep)]+Drep[as.vector(indices2)[2:length(indices2)],3:ncol(Drep)])/2
        }
        
        before[(before==0) | (before==2*N)]<-NA

        deltayear <- as.vector(afteryear-before)
        before <- as.vector(before)
        before <- before[!is.na(deltayear)]
        deltayear <- deltayear[!is.na(deltayear)]
                                        #compute average frequency changes for grid of initial frequencies
        means <- numeric(length(gridvec)-1)
        for (g in 1:(length(gridvec)-1)){
            deltas <- subset(deltayear,(before>=2*N*gridvec[g]) & (before<2*N*gridvec[g+1]))
            means[g]<-mean(deltas)/(2*N)
        }
                                        #model fitting
        midpoints <- (gridvec[1:(length(gridvec)-1)]+gridvec[2:length(gridvec)])/2
        mutationinput <- cyclelength*mutrate*(1-2*midpoints)
        meanpred <- midpoints*(1-midpoints)*(1-2*midpoints)
        meanmod <-lm(means-mutationinput ~ -1 + meanpred)
        coefs[r] <- meanmod$coefficients[1]
        allmeans[r,] <- means
    }
    return(list(coefs,colMeans(allmeans,na.rm=TRUE)))
}

stableornotsimple <- function(D,N,wintergenerations,summergenerations,mutrate){#for cases without parameter variation across loci, here we do not need to subdivide replicates, returns average and standard error or summer allele frequency across time and loci
    cyclelength <- wintergenerations+summergenerations
    numreps <- length(unique(D[,1]))
    numgen <- D[nrow(D),2]
    numyears <- numgen/cyclelength
    gridvec <- seq(0,0.5,0.025)
    midpoints <- (gridvec[1:(length(gridvec)-1)]+gridvec[2:length(gridvec)])/2
    mutationinput <- cyclelength*mutrate*(1-2*midpoints)
    nb <- length(midpoints)
    means <- matrix(NA,nrow=length(gridvec)-1,ncol=numreps)
    averagefrequencies <- numeric(numreps)
    
    for (r in 1:numreps){
        Drep <- subset(D,D[,1]+1==r)
        
        M <- matrix(1:numgen,nrow=numgen/2)
        if (wintergenerations %% 2 == 0){
            indices <- M[seq(1/2*wintergenerations,nrow(M),cyclelength),2]#middle of the year

                                        #compute frequency changes over one year
            before <- Drep[as.vector(indices)[1:(length(indices)-1)],3:ncol(Drep)]
            afteryear <- Drep[as.vector(indices)[2:length(indices)],3:ncol(Drep)]
        }
        else {
            indices1 <- M[floor(seq(1/2*wintergenerations,nrow(M),cyclelength)),2]
            indices2 <- M[ceiling(seq(1/2*wintergenerations,nrow(M),cyclelength)),2]

                                        #compute frequency changes over one year
            before <- (Drep[as.vector(indices1)[1:(length(indices1)-1)],3:ncol(Drep)]+Drep[as.vector(indices2)[1:(length(indices2)-1)],3:ncol(Drep)])/2
            afteryear <- (Drep[as.vector(indices1)[2:length(indices1)],3:ncol(Drep)]+Drep[as.vector(indices2)[2:length(indices2)],3:ncol(Drep)])/2
        }
        
        before[(before==0) | (before==2*N)]<-NA

        deltayear <- as.vector(afteryear-before)
        before <- as.vector(before)
        before <- before[!is.na(deltayear)]
        deltayear <- deltayear[!is.na(deltayear)]

        for (g in 1:nb){
            deltas <- subset(deltayear,(before>=2*N*gridvec[g]) & (before<2*N*gridvec[g+1]))
            deltas2 <- subset(deltayear,(2*N-before>=2*N*gridvec[g]) & (2*N-before<2*N*gridvec[g+1]))
            means[g,r]<-mean(c(deltas,-deltas2))/(2*N)
        }

        averagefrequencies[r] <- mean(before,na.rm=TRUE)/(2*N)
        
    }
    
    means <- means-matrix(rep(mutationinput,numreps),ncol=numreps)

    signs <- numeric(nb)
    meanofmeans <- rowMeans(means,na.rm=TRUE)
    for (b in 1:nb){
        se <- sd(means[b,],na.rm=TRUE)
        nr <- sum(!is.na(means[b,]))
        lower <- meanofmeans[b]-2*se/sqrt(nr)
        upper <- meanofmeans[b]+2*se/sqrt(nr)
        if (nr > 1){
            if (lower > 0){
                signs[b] <- 1
            }
            else if (upper < 0){
                signs[b] <- -1
            }
        }
    }
    
    b <- 0
    firstsignificant <- FALSE
    while((!firstsignificant)&(b<nb)){
        b <- b+1
        firstsignificant <- (signs[b]!=0)
    }
    stab <- (signs[b]>0)
    
    frequencyse <- sd(averagefrequencies)/(sqrt(numreps))
    return(c(stab,mean(averagefrequencies),frequencyse))
}

stableornot <- function(D,N,wintergenerations,summergenerations,mutrate){ #for cases with parameter variation across loci
    cyclelength <- wintergenerations+summergenerations
    numreps <- length(unique(D[,1]))
    numgen <- D[nrow(D),2]
    numyears <- numgen/cyclelength
    gridvec <- seq(0,0.5,0.025)
    midpoints <- (gridvec[1:(length(gridvec)-1)]+gridvec[2:length(gridvec)])/2
    mutationinput <- cyclelength*mutrate*(1-2*midpoints)
    nb <- length(midpoints)
    stabvec <- numeric(numreps)

    for (r in 1:numreps){
        Drep <- subset(D,D[,1]+1==r)
        
        M <- matrix(1:numgen,nrow=numgen/2)
        indices <- M[seq(1/4*cyclelength,nrow(M),cyclelength),2]#middle of the year

        #compute frequency changes over one year
        before <- Drep[as.vector(indices)[1:(length(indices)-1)],3:ncol(Drep)]
        afteryear <- Drep[as.vector(indices)[2:length(indices)],3:ncol(Drep)]
        before[(before==0) | (before==2*N)]<-NA

        deltayear <- as.vector(afteryear-before)
        before <- as.vector(before)
        before <- before[!is.na(deltayear)]
        deltayear <- deltayear[!is.na(deltayear)]

        #divide data into 10 subreplicates
        npr <- 10
        means <- matrix(NA,nrow=length(gridvec)-1,ncol=npr)

        pointsperpr <- floor(length(deltayear)/npr) #number of data points per subreplicate

        for (pr in 1:npr){
        #compute average frequency changes for grid of initial frequencies
            deltaspr <- deltayear[((pr-1)*pointsperpr+1):(pr*pointsperpr)]
            beforepr <- before[((pr-1)*pointsperpr+1):(pr*pointsperpr)]
            
            for (g in 1:nb){
                deltas <- subset(deltaspr,(beforepr>=2*N*gridvec[g]) & (beforepr<2*N*gridvec[g+1]))
                deltas2 <- subset(deltaspr,(2*N-beforepr>=2*N*gridvec[g]) & (2*N-beforepr<2*N*gridvec[g+1]))
                means[g,pr]<-mean(c(deltas,-deltas2))/(2*N)
            }
        }

        means <- means-matrix(rep(mutationinput,npr),ncol=npr)

        signs <- numeric(nb)
        meanofmeans <- rowMeans(means,na.rm=TRUE)
        for (b in 1:nb){
            se <- sd(means[b,],na.rm=TRUE)
            nr <- sum(!is.na(means[b,]))
            lower <- meanofmeans[b]-2*se/sqrt(nr)
            upper <- meanofmeans[b]+2*se/sqrt(nr)
            if (nr > 1){
                if (lower > 0){
                    signs[b] <- 1
                }
                else if (upper < 0){
                    signs[b] <- -1
                }
            }
        }
        
        b <- 0
        firstsignificant <- FALSE
        while((!firstsignificant)&(b<nb)){
            b <- b+1
            firstsignificant <- (signs[b]!=0)
        }
        if (signs[b]>0){
            stabvec[r] <- TRUE
        }
        else {
            stabvec[r] <- FALSE
        }
        
    }
    return(stabvec)
}

stabilitycoefficient2 <- function(D,N,wintergenerations,summergenerations,mutrate){#input parameters: data frame with allele-frequency directories, population size, number of generations per winter, number of generations per summer
    #determine data dimension and relevant indices
    cyclelength <- wintergenerations+summergenerations
    numreps <- D[nrow(D),1]+1 #because the replicate count starts from 0
    numgen <- D[nrow(D),2]
    numyears <- numgen/cyclelength
    coefs <- numeric(numreps)
    
    for (r in 1:numreps){
        Drep <- subset(D,D[,1]+1==r)
        M <- matrix(1:numgen,nrow=numgen/2)
        indices <- M[seq(1/4*cyclelength,nrow(M),cyclelength),2]#middle of the cold season as reference point every year, use only data from the second half of each simulation run

                                        #compute frequency changes over one year
        before <- Drep[as.vector(indices)[1:(length(indices)-1)],3:ncol(Drep)]
        afteryear <- Drep[as.vector(indices)[2:length(indices)],3:ncol(Drep)]
        before[(before==0) | (before==2*N)]<-NA

        deltayear <- as.vector(afteryear-before)
        before <- as.vector(before)
        before <- before[!is.na(deltayear)]
        deltayear <- deltayear[!is.na(deltayear)]
        mutationinput <- cyclelength*mutrate*(1-2*before)
        meanpred <- before*(1-before)*(1-2*before)
        meanmod <-lm(deltayear-mutationinput ~ -1 + meanpred)
        coefs[r] <- meanmod$coefficients[1]
    }
    return(coefs)
}

#works only if wintergenerations = summergenerations
effectives <- function(D,N,wintergenerations,summergenerations,mutrate){#input parameters: data frame with allele-frequency directories, population size, number of generations per winter, number of generations per summer
    #determine data dimension and relevant indices
    cyclelength <- wintergenerations+summergenerations
    numreps <- D[nrow(D),1]+1 #because the replicate count starts from 0
    numgen <- D[nrow(D),2]
    numyears <- numgen/cyclelength
    coefs <- numeric(numreps)
    gridvec <- seq(0,1,0.04)
    allmeans <- matrix(NA,nrow=numreps,ncol=length(gridvec)-1)
    
    for (r in 1:numreps){
        Drep <- subset(D,D[,1]+1==r)
        M <- matrix(1:numgen,nrow=numgen/2)
        indices <- M[seq(summergenerations,nrow(M),summergenerations),2]

    #compute frequency changes over one season and correct for direction
        before <- as.matrix(Drep[as.vector(indices)[1:(length(indices)-2)],3:ncol(Drep)])
        afterseason <- as.matrix(Drep[as.vector(indices)[2:(length(indices)-1)],3:ncol(Drep)])
        before[(before==0) | (before==2*N)]<-NA
        nb <- nrow(before)
        
        before <- c(before[seq(1,nb,2),],1-before[seq(2,nb,2),])
        afterseason <- c(afterseason[seq(1,nb,2),],1-afterseason[seq(2,nb,2),])
        deltaseason <- as.vector(afterseason-before)

        {
        if (summergenerations %% 2 == 0){
            indices <- M[seq(summergenerations+summergenerations/2,nrow(M),summergenerations),2]
            middle <- as.matrix(Drep[as.vector(indices)[1:(length(indices)-1)],3:ncol(Drep)])
        }
        else {
            indices1 <- M[seq(summergenerations+floor(summergenerations/2),nrow(M),summergenerations),2]
            middle1 <- as.matrix(Drep[as.vector(indices1)[1:(length(indices1)-1)],3:ncol(Drep)])
            indices2 <- M[seq(summergenerations+ceiling(summergenerations/2),nrow(M),summergenerations),2]
            middle2 <- as.matrix(Drep[as.vector(indices2)[1:(length(indices2)-1)],3:ncol(Drep)])
            middle <- (middle1+middle2)/2
        }
    }
        middle <- c(middle[seq(1,nrow(middle),2),],1-middle[seq(2,nrow(middle),2),])
        
        middle <- middle[!is.na(deltaseason)]
        deltaseason <- deltaseason[!is.na(deltaseason)]

    #compute average frequency changes for grid of initial frequencies
        means <- numeric(length(gridvec)-1)
        for (g in 1:(length(gridvec)-1)){
            deltas <- subset(deltaseason,(middle>=2*N*gridvec[g]) & (middle<2*N*gridvec[g+1]))
            means[g]<-mean(deltas)/(2*N)
        }

    #model fitting
        midpoints <- (gridvec[1:(length(gridvec)-1)]+gridvec[2:length(gridvec)])/2
        mutationinput <- summergenerations*mutrate*(1-2*midpoints)
        meanpred <- midpoints*(1-midpoints)
        meanmod <-lm(means-mutationinput ~ -1 + meanpred)
        coefs[r] <- meanmod$coefficients[1]
        allmeans[r,] <- means
    }
    return(list(coefs,colMeans(allmeans,na.rm=TRUE)))
}

propexpecteddirection <- function(D,N,wintergenerations,summergenerations){#input parameters: data frame with allele-frequency directories, population size, number of generations per winter, number of generations per summer
    #determine data dimension and relevant indices
    cyclelength <- wintergenerations+summergenerations
    numreps <- D[nrow(D),1]+1 #because the replicate count starts from 0
    numgen <- D[nrow(D),2]
    numyears <- numgen/cyclelength
    propexpected <- numeric(numreps)
    
    for (r in 1:numreps){
        Drep <- subset(D,D[,1]+1==r)
        M <- matrix(1:numgen,nrow=numgen/2)
        indices <- M[seq(summergenerations,nrow(M),summergenerations),seq(2,ncol(M),2)]

                                        #compute frequency changes over one season and correct for direction
        before <- as.matrix(Drep[as.vector(indices)[1:(length(indices)-2)],3:ncol(Drep)])
        afterseason <- as.matrix(Drep[as.vector(indices)[2:(length(indices)-1)],3:ncol(Drep)])
        before[(before==0) | (before==2*N)]<-NA
        deltaseason <- as.vector(afterseason-before)*rep((-1)^(1+(1:nrow(afterseason))),ncol(afterseason))
        deltaseason <- deltaseason[!is.na(deltaseason)]
        propexpected[r] <- sum(deltaseason>0)/sum(deltaseason!=0)
    }
    return(propexpected)
}

fitnesscoefvar <- function(D,wintergenerations,summergenerations){
    cyclelength <- wintergenerations+summergenerations
    numreps <- D[nrow(D),1]+1 #because the replicate count starts from 0
    numgen <- D[nrow(D),2]
    coefvar <- numeric(numreps)
    for (r in 1:numreps){
        Drep <- subset(D,D[,1]+1==r)
        Drep <- subset(Drep,Drep[,2]>numgen/2)
        coefvar[r] <- mean(Drep[,4]/Drep[,3])
    }
    return(coefvar)
}

tracking <- function(D,wintergenerations,summergenerations){
    cyclelength <- wintergenerations+summergenerations
    numreps <- D[nrow(D),1]+1 #because the replicate count starts from 0
    numgen <- D[nrow(D),2]
    fitnessmeans <- matrix(NA,nrow=cyclelength,ncol=numreps)
    fitnesscoefvars <- matrix(NA,nrow=cyclelength,ncol=numreps)
    for (r in 1:numreps){
        Drep <- subset(D,D[,1]+1==r)
        Drep <- subset(Drep,Drep[,2]>numgen/2)
        meanmatrix <- matrix(Drep[,3],nrow=cyclelength)
        coefvarmatrix <- matrix(Drep[,4]/Drep[,3],nrow=cyclelength)
        fitnessmeans[,r] <- rowMeans(meanmatrix)
        fitnesscoefvars <- rowMeans(coefvarmatrix)
    }
    return(list(fitnessmeans,fitnesscoefvars))
}


#states whether a locus is detectable according to the following criterion: change in the expected direction by at least 5% in at least half of the seasons (detecvec) + Proportion of 3-year data sets for which a locus changes by at least 5 % in the expected direction in every season (detecvecshort)
detectability <- function(D,N,wintergenerations,summergenerations){#input parameters: data frame with allele-frequency trajectories for a single locus in different replicates (input data frame has just three columns: replicate, time, allele frequency of the focal locus), population size, number of generations per winter, number of generations per summer
    #determine data dimension and relevant indices
    cyclelength <- wintergenerations+summergenerations
    numreps <- D[nrow(D),1]+1 #because the replicate count starts from 0
    numgen <- D[nrow(D),2]
    numyears <- numgen/cyclelength
    detecvec <- numeric(numreps)
    detecvecshort <- numeric(numreps)
    for (r in 1:numreps){
        Drep <- subset(D,D[,1]+1==r)
        M <- matrix(1:numgen,nrow=numgen/2)
        indices <- M[seq(summergenerations,nrow(M),summergenerations),2]

    #compute frequency changes over one season and correct for direction
        before <- as.vector(Drep[as.vector(indices)[1:(length(indices)-2)],3])
        afterseason <- as.vector(Drep[as.vector(indices)[2:(length(indices)-1)],3])
        before[(before==0) | (before==2*N)]<-NA
        nb <- length(before)
        deltaseason <- (afterseason-before)/(2*N)
        deltaseason <- deltaseason*rep(c(1,-1),length.out=nb)
        detecvec[r] <- ((sum(deltaseason>0.05)/length(deltaseason)) >= 0.5)
        windowdetect <- foreach (t=1:(nb-5),.combine='c')%do%{
            return(min(deltaseason[t:(t+5)])>=0.05)
        }
        detecvecshort[r] <- sum(windowdetect,na.rm=TRUE)/length(windowdetect)

    }
    return(rbind(detecvec,detecvecshort))
}


criticaldsingle <- function(stabcoefs,dtreatments){
    numneg <- sum(stabcoefs<0)
    if (numneg < length(dtreatments)){
        return((dtreatments[numneg]+dtreatments[numneg+1])/2)
    }
    else {
        return(NA)
    }
}

criticald <- function(stabcoefmatrix,dtreatments){#matrix with replicates in columns and d values in rows; for each replicate, the simulations for the different d values are independent, but that doesn't matter
    critd <- as.numeric(apply(stabcoefmatrix,2,function(x) criticaldsingle(x,dtreatments)))
    return(c(mean(critd,na.rm=TRUE),sd(critd,na.rm=TRUE)))
}

