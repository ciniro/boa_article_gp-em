#-----------------------------------------------------------------------------#
#                            BOA - FUNCTIONS                                  #
#                           BINARY PROBLEMS                                   #
#                                                                             #
# Ciniro Ap. Leite Nametala                                                   #
# ciniro@gmail.com                                                            #
# V1.00: April-22-2021                                                        #
#                                                                             #
# R Source-code of paper "On the performance of the Bayesian Optimization     #
# Algorithm with combined scenarios of Search Algorithms and Scoring Metrics" #
# submitted to Genetic Programming and Evolvable Machines - Springer          #
#-----------------------------------------------------------------------------#

#fitness functions

#onemax
onemax <- function(cromo) {
  n <- length(cromo)
  sizebb <- n/qtdebb
  
  fitness <- 0
  
  for(i in seq(1, n, by = sizebb)) {
    bb <- cromo[i:(i + sizebb - 1)]
    
    fitness <- fitness + sum(bb)
  }
  
  return(fitness/10)
}

#leadingones
leadingones <- function(cromo) {
  n <- length(cromo)
  sizebb <- n/qtdebb
  
  fitness <- 0
  
  for(i in seq(1, n, by = sizebb)) {
    bb <- cromo[i:(i + sizebb - 1)]
    
    fitness <- fitness + prod(bb)
  }
  
  return(fitness)
}

#trap-10
trap <- function(cromo) {
  n <- length(cromo)
  sizebb <- n/qtdebb

  ynorm <- c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0, 1)
  
  fitness <- 0
  
  for(i in seq(1, n, by = sizebb)) {
    bb <- cromo[i:(i + sizebb - 1)]
    
    ones <- sum(bb)
    
    fitbb <- 0

    fitbb <- ynorm[ones + 1]
    
    fitness <- fitness + fitbb
  }
  
  return(fitness)
}

#deceptive-10
deceptive <- function(cromo) {
  n <- length(cromo)
  sizebb <- n/qtdebb
  
  ynorm <- c(0.9999, 0.9997, 0.9992, 0.9976, 0.9933, 0.9818, 0.9503, 0.8647, 0.6321, 0, 1)
  
  fitness <- 0
  
  for(i in seq(1, n, by = sizebb)) {
    bb <- cromo[i:(i + sizebb - 1)]
    
    ones <- sum(bb)
    
    fitbb <- 0

    fitbb <- ynorm[ones + 1]
    
    fitness <- fitness + fitbb
  }
  
  return(fitness)
}

#bipolar-10
bipolar <- function(cromo) {
  n <- length(cromo)
  sizebb <- n/qtdebb

  ynorm <- c(1,0,0.25,0.5,0.75,0.8,0.75,0.5,0.25,0,1)
  
  fitness <- 0
  
  for(i in seq(1, n, by = sizebb)) {
    bb <- cromo[i:(i + sizebb - 1)]
    
    ones <- sum(bb)
    
    fitbb <- ynorm[ones + 1]
    
    fitness <- fitness + fitbb
  }
  
  return(fitness)
}

#multi
multi <- function(cromo) {
  n <- length(cromo)
  sizebb <- n/qtdebb
  
  #deceptive
  ynormf1 <- c(0.99, 0.9997, 0.9992, 0.9976, 0.9933, 0.9818, 0.9503, 0.8647, 0.6321, 0, 1)
  #bipolar
  ynormf2 <- c(1,0,0.25,0.5,0.75,0.8,0.75,0.5,0.25,0,1)

  fitness <- 0
  
  inverter <- 0
  for(i in seq(1, n, by = sizebb)) {
    bb <- cromo[i:(i + sizebb - 1)]
    
    ones <- sum(bb)
    
    fitbb <- 0
    
    if (inverter == 0) {
      fitbb <- ynormf1[ones + 1]
    } else {
      fitbb <- ynormf2[ones + 1]
    }
    
    if (inverter == 0)
        inverter <- 1
    else
      inverter <- 0
    
    fitness <- fitness + fitbb
  }
  
  return(fitness)
}

#BOA FUNCTIONS
generate_initial_model <- function(size) {
  xi <- paste("x", 1:size, sep="")
  bn <- empty.graph(nodes = xi)
  
  lv <- c("0","1")
  
  cpt <- list()
  
  for (i in 1:size) {
    cpt[[i]] <- array(c(0.5, 0.5), 
                      dim = 2, 
                      dimnames = list(x = lv))
  }
  
  names(cpt) <- xi
  
  model = custom.fit(bn, cpt)
  
  return(list(model = model, bn = bn))
}

get_best <- function(poplocalbest) {
  mpop <- apply(as.matrix(poplocalbest), 2, as.numeric)
  fits <- apply(mpop, 1, fitfunction)
  indbest <- which.max(fits)
  
  return(list(best = mpop[indbest,], fit_best = fits[indbest]))
}

get_stats_pop <- function(poplocalbest) {
  poplocalbest <- apply(as.matrix(poplocalbest), 2, as.numeric)
  fits <- apply(poplocalbest, 1, fitfunction)
  
  lsstats <- list(meanpop = mean(fits),
                  medianpop = median(fits),
                  sdpop = sd(fits),
                  bestpop = max(fits),
                  worstpop = min(fits))
  
  return(lsstats)
}

tournament <- function(poplocal, npop, k, n) {
  winners <- vector()
  fits <- vector()
  nwinners <- ceiling(npop*n)
  
  for (i in 1:nwinners) {
    ind_competitors <- sample(1:npop, k)
    winner <- get_best(poplocal[ind_competitors,])
    winnercromo <- apply(as.matrix(winner$best), 2, as.numeric)
    winners <- rbind(winners, as.numeric(winnercromo))
    fits <- c(fits, as.numeric(winner$fit_best))
  }

  return (list(winners, fits))
}

build_bn <- function(dfpop, size, previous_bn, alg, sco) {
  dfpop <- apply(dfpop, 2, as.character)
  #include levels (at least one 0 e one 1)
  dfpop <- rbind(dfpop, rep("0", size), rep("1", size))
  dfpop <- as.data.frame(dfpop)
  xi <- paste("x", 1:size, sep="")
  colnames(dfpop) <- xi
  
  bn <- NULL

  if (alg == "hc")
    bn <- hc(dfpop, score = sco, restart = 5, perturb = 3)
  else if (alg == "tabu")
    bn <- tabu(dfpop, score = sco, tabu = 5)
  
  model = bn.fit(bn, dfpop, replace.unidentifiable = TRUE)
  
  return(list(model = model, bn = bn))
}

rankpop <- function(unrankedpop) {
  unrankedpop <- apply(as.matrix(unrankedpop), 2, as.numeric)
  fits <- apply(unrankedpop, 1, fitfunction)
  return(sort(fits, index.return=TRUE)$ix)
}

replaceworse <- function(originalpop, newinds) {
  rankinds <- rankpop(originalpop)
  
  for (i in 1:dim(newinds)[1]) {
    a <- apply(as.matrix(newinds[i,]), 2, as.numeric)
    originalpop[rankinds[i],] <- a
  }
  
  return(originalpop)
}

readFile = function(filepath) {
  con = file(filepath, "r")
  value <- NULL
  while ( TRUE ) {
    line = readLines(con, n = 1, warn = FALSE)
    if ( length(line) == 0 ) {
      break
    }
    value <- line
  }
  
  close(con)
  return(value)
}

mkdirs <- function(fp) {
  if(!file.exists(fp)) {
    mkdirs(dirname(fp))
    dir.create(fp)
  }
}