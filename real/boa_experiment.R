#-----------------------------------------------------------------------------#
#                                  BOA                                        #
#                              REAL PROBLEMS                                  #
#                                                                             #
# Ciniro Ap. Leite Nametala                                                   #
# ciniro@gmail.com                                                            #
# V1.00: April-22-2021                                                        #
#                                                                             #
# R Source-code of paper "On the performance of the Bayesian Optimization     #
# Algorithm with combined scenarios of Search Algorithms and Scoring Metrics" #
# submitted to Genetic Programming and Evolvable Machines - Springer          #
#-----------------------------------------------------------------------------#

#prepare environment
rm(list=ls())
cat('\014')

#libraries
library("bnlearn")
library("plot3D")
library("stringr")
library("rgl")
library("magick")
library("data.table")
library("pheatmap")

#external functions
source("boa_functions_real.R")

#==============================================================================
#learning algorithms
algos <- c("hc","tabu")
#network scores
scores <- c("k2","aic","bic","bds")
#parameters of experiments
nexp <- 20
#==============================================================================

#parameters of selection
n.tournament <- 1 #percent
k.tournament <- 5
n.new <- 0.2 #percent of new individuals
n.explore <- 0.1 #percent of randon individuals
  
#parameters of convergence
maxgen <- 300
npop <- 500
tol <- 0.001
checkVarianceIter <- 6

#representation
ndim <- 2
linf <- -2.5
lsup <- 2.5
precisao <- 0.001
bitnec <- ceiling(log2((lsup-linf)/precisao))
casas <- contaCasasDec(precisao)
lsupbin <- RealToBin(lsup)
linfbin <- RealToBin(linf)
size <- ndim*bitnec

#parameters of plots
plotlast <- FALSE
exportConfigFile <- TRUE

#graphs
plotGraphs <- FALSE
resPlotFunc <- 0.1
plotFuncRotate <- FALSE

#select function
#"sphere","perm","rastrigin","rosenbrock","michal","schaffer"
fitfunction <- f.michal
fitname <- "michalewicz"
print(fitname)

#plot function fitness
if (plotGraphs == TRUE)
{
  x <- seq(linf,lsup,resPlotFunc)
  y <- x
  z <- matrix(NA, length(x), length(y))
  for (i in 1:length(x))
  {
    for (j in 1:length(y))
    {
      z[i,j] <- fitfunction(c(x[i],y[j]))
    }
  }
  
  if (plotFuncRotate == TRUE) {
    persp3D(x,y,z,theta=45,phi=30)
    
    rgl::persp3d(x,y,z, col = "lightblue",
                 ticktype="detailed", xlab="", ylab="", zlab="",axes=FALSE)
  }
}

#prepare the files
if (exportConfigFile == TRUE) {
  nameExperiment <- toupper(fitname)
  
  title1 <- "BAYESIAN OPTIMIZATION ALGORITHM"
  numberexp <- paste("NAME OF EXPERIMENT: ", nameExperiment, sep="")
  tExps <- paste("TOTAL OF EXPERIMENTS: ", nexp, sep="")
  sep1 <- "================================================="
  title2 <- "PARAMETERS"
  tNTournament <- paste("N TOURNAMENT: ", n.tournament, sep="")
  tKTournament <- paste("K TOURNAMENT: ", k.tournament, sep="")
  tMaxGen <- paste("MAX GENERATIONS: ", maxgen, sep="")
  tNPop <- paste("SIZE OF POPULATION: ", npop, sep="")
  tFitFunction <- paste("FIT FUNCTION: ", fitname, sep="")
  tTolerance <- paste("TOLERANCE: ", tol, sep="")
  sep2 <- "================================================="
  title3 <- "RESULTS"
  
  pathfile <- paste("results/output_", nameExperiment, sep="")
  mkdirs(pathfile)
  
  namefile_config <- "data_experiment.txt"
  allpath_config <- paste(pathfile,"/",namefile_config,sep="")
  
  namefile_geral <- "geral_experiment.csv"
  allpath_geral <- paste(pathfile,"/",namefile_geral,sep="")
  
  writeLines(c(title1,
               numberexp,
               tExps,
               sep1,
               title2,
               tNTournament,
               tKTournament,
               tMaxGen,
               tNPop,
               tFitFunction,
               tTolerance,
               sep2,
               title3), 
               allpath_config)
}

dfgeral <- NULL

#==============================================================================
for (n in 1:nexp) {
  namefile_trees <- paste("trees", n, ".txt", sep="")
  namefile_inds <- paste("ind", n, ".csv", sep="")
  namefile_stats <- paste("stats", n, ".csv", sep="")
  
  #inicializes randomly a BN
  best_bn_zero <- generate_initial_model(size)

  #generate initial population
  pop_zero <- rbn(best_bn_zero$model, n = npop)
  pop_zero <- bounds_pop(pop_zero)
  
  #get initial best cromosoma
  bestcromo_zero <- get_best(pop_zero)

  for (alg in algos) {
    for (sco in scores) {
      
      print(paste("Exec: ",alg,"-",sco))
      
      path_trees <- paste(pathfile,"/",alg,"_",sco,sep="")
      path_inds <- paste(pathfile,"/",alg,"_",sco,sep="")
      path_stats <- paste(pathfile,"/",alg,"_",sco,sep="")
      path_adjmat <- paste(pathfile,"/",alg,"_",sco,"/adjmatrices/",sep="")
      
      mkdirs(path_trees)
      mkdirs(path_inds)
      mkdirs(path_stats)
      mkdirs(path_adjmat)
      
      allpath_trees <- paste(path_trees,"/",namefile_trees,sep="")
      allpath_inds <- paste(path_inds,"/",namefile_inds,sep="")
      allpath_stats <- paste(path_stats,"/",namefile_stats,sep="")

      #inicializes randomly a BN
      best_bn <- best_bn_zero
      
      #generate initial population
      pop <- pop_zero
      
      #get initial best cromosoma
      bestcromo <- bestcromo_zero
      
      #storage stats populations by iteration
      dfstats <- matrix(data=NA,nrow=maxgen,ncol=6)
      colnames(dfstats) <- c("mean","median","sd","best","worst","score")
      
      dfinds <- NULL
      
      meanfitness <- rep(NA, maxgen)
      bestfitness <- rep(NA, maxgen)
      matpop <- bestcromo$best
      
      for (iter in 1:maxgen) {
        print(paste("Iter: ",iter," - Best: ",bestcromo$fit_best,sep=""))
        
        #let them compete
        winners <- tournament(pop, npop, k.tournament, n.tournament)
        
        #variance
        if  (iter > checkVarianceIter)
          if (length(unique(c(meanfitness[iter - 1],
                              meanfitness[iter - 2],
                              meanfitness[iter - 3],
                              meanfitness[iter - 4]))) == 1)
            break()
      
        
        #building a new BN based on winners
        best_bn <- build_bn(winners[[1]], size, best_bn$bn, alg, sco)
        write.table(amat(best_bn$bn), file=paste(path_adjmat,"iter_",iter,".csv",sep=""), row.names=FALSE, col.names=TRUE, sep = ";")
        
        #generate new individuals based on new BN
        newpop <- rbn(best_bn$model, n = ceiling(npop*n.new), winners[[1]])
        newpop <- bounds_pop(newpop)
        
        #include new individuals in pop
        pop <- replaceworse(pop, newpop)
        
        #include new randon individuals in pop
        pop_randon <- rbn(best_bn_zero$model, n = ceiling(npop*n.explore))
        pop_randon <- bounds_pop(pop_randon)
        pop <- replaceworse(pop, pop_randon)

        #evaluate best against winner
        newbestcromo <- get_best(pop)
        
        if (newbestcromo$fit_best < bestcromo$fit_best) {
          bestcromo <- newbestcromo
        }
        
        #store the stats of iteration
        lsstats <- get_stats_pop(pop)
        
        dfstats[iter,"mean"] <- round(lsstats$meanpop,2)
        dfstats[iter,"median"] <- round(lsstats$medianpop,2)
        dfstats[iter,"sd"] <- round(lsstats$sdpop,2)
        dfstats[iter,"best"] <- round(lsstats$bestpop,2)
        dfstats[iter,"worst"] <- round(lsstats$worstpop,2)
        
        meanfitness[iter] <- lsstats$meanpop
        bestfitness[iter] <- bestcromo$fit_best
        matpop <- rbind(matpop, bestcromo$best)
        
        dfinds <- rbind(dfinds, bestcromo$best)
      }
      
      write.table(dfstats[1:iter,], file=allpath_stats, row.names=FALSE, col.names=TRUE, sep = ";")
      write.table(dfinds, file=allpath_inds, row.names=FALSE, col.names=TRUE, sep = ";")
      
      write(paste(alg,"-",sco),
            file=allpath_config,
            append=TRUE)
      write(paste("BEST FIT:", bestcromo$fit_best),
            file=allpath_config,
            append=TRUE)
      write(paste("STOP:",iter),
            file=allpath_config,
            append=TRUE)
      write("--------------------------------------------",
            file=allpath_config,
            append=TRUE)
      
      results <- c(n, 
                   alg, 
                   sco, 
                   bestcromo$fit_best,
                   iter,
                   round(lsstats$meanpop,2),
                   round(lsstats$medianpop,2),
                   round(lsstats$sdpop,2))
      
      dfgeral <- rbind(dfgeral, results)
                   
      #==============================================================================
      if (plotlast == TRUE) {
        #MEAN FITNESS CONVERGENCE
        plot(meanfitness, type="l", xlab="Iterations", 
             ylab="Mean fitness", 
             main=paste(alg, "-", sco, ": Convergence - Mean fitness", sep=""), panel.first=grid(), lwd=2)
        
        abline(v=iter-1, col="blue", lty=2)
        
        #BEST FITNESS CONVERGENCE
        plot(bestfitness, type="l", xlab="Iterations",
              ylab="Best fitness",
              main="Convergence - Best fitness", panel.first=grid(), lwd=2)
        
        abline(v=iter-1, col="blue", lty=2)
        
        a <- convertCromoBinToReal(bestcromo$best)
        a <- t(apply(as.matrix(a), 2, as.numeric))
        generateContours(a,x,y,z)
        
        matpop <- apply(as.matrix(matpop), 2, as.numeric)
        rownames(matpop) <- as.character(seq(1:dim(matpop)[1]))
        pheatmap(matpop, 
                 cluster_rows = FALSE, 
                 cluster_cols = FALSE, 
                 scale = "none", 
                 color = colorRampPalette(c("white", "navy", "black"))(100), 
                 border_color = "grey60", 
                 breaks = seq(0,1,0.01), 
                 legend_breaks = seq(0,1,0.1), 
                 show_colnames = FALSE, 
                 show_rownames = FALSE, 
                 main = "Heatmap convergence")
        
        graph::plot(best_bn$bn, main=paste(alg, "-", sco, sep=""))
        print("STOP ITERATION:")
        print(iter)
        print("BEST SOLUTION:")
        print(bestcromo$best)
        print(convertCromoBinToReal(bestcromo$best))
        print("FIT BEST SOLUTION:")
        print(bestcromo$fit_best)
        print("STRING MODEL:")
        print(modelstring(best_bn$model))
      }
    }
  }
}

colnames(dfgeral) <- c("exp","alg","sco","best","iter","mean","median","sd")
write.table(dfgeral, file=allpath_geral, row.names=FALSE, col.names=TRUE, sep = ";")
  
beepr::beep(1)