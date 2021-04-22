#-----------------------------------------------------------------------------#
#                            BOA - FUNCTIONS                                  #
#                             REAL PROBLEMS                                   #
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
#sphere
#optfitness 0..
f.sphere <- function(x) 
{
  return(sum(x^2))
}

#rastrigin
#optfitness: 0
f.rastrigin <- function(x)
{
  d <- length(x)
  sum <- sum(x^2 - 10*cos(2*pi*x))
  y <- 10*d + sum
  return(y)
}

#rosenbrock
#optfitness: 0
f.rosenbrock <- function(x)
{
  d <- length(x)
  xi <- x[1:(d-1)]
  xnext <- x[2:d]
  
  sum <- sum(100*(xnext-xi^2)^2 + (xi-1)^2)
  
  y <- sum
  return(y)
}

#zakharov
#optfitness: 0..
f.zakharov <- function(x)
{
  ii <- c(1:length(x))
  sum1 <- sum(x^2)
  sum2 <- sum(0.5*ii*x)
  
  y <- sum1 + sum2^2 + sum2^4
  return(y)
}

#perm
#optim: 0..
f.perm <- function(x, b=0.5)
{
  d <- length(x)
  ii <- c(1:d)
  jj <- matrix(rep(ii,times=d), d, d, byrow=TRUE)
  
  xxmat <- matrix(rep(x,times=d), d, d, byrow=TRUE)
  inner <- rowSums((jj^ii+b)*((xxmat/jj)^ii-1))	
  outer <- sum(inner^2)
  
  y <- outer
  return(y)
}

#michalewicz
#optfitness: 0
f.michal <- function(xx, m=10)
{
  ii <- c(1:length(xx))
  sum <- sum(sin(xx) * (sin(ii*xx^2/pi))^(2*m))
  
  y <- -sum + 1.8013
  return(y)
}

#schaffer2 function
#optfitness: 0
f.schaffer <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
  
  fact1 <- (sin(x1^2-x2^2))^2 - 0.5
  fact2 <- (1 + 0.001*(x1^2+x2^2))^2
  
  y <- 0.5 + fact1/fact2
  return(y)
}

#graph functions
generateContours <- function(a, x = seq(0, 1, length.out = nrow(z)),
                             y = seq(0, 1, length.out = ncol(z)),
                             z, nlevels=30, ...) {
  ma <- max(abs(z))
  lvls <- seq(-ma, ma, length.out = nlevels)
  cols <- colorRampPalette(c("blue","white","red")) (nlevels - 1)
  filled.contour(x, y, z, 
                 main = "Optimization point", 
                 plot.axes = {axis(1); axis(2); points(x = a[1], y = a[2], pch=20, col="blue", cex=3)},
                 col=cols, levels=lvls, ...)
}

#representation functions
BinToDec <- function(x) {
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
}

DecToBin <- function(x){
  tmp <- rev(as.numeric(intToBits(x)))
  id <- seq_len(match(1,tmp,length(tmp))-1)
  tmp[-id]
}

BinToReal <- function (bit) {
  dec <- BinToDec(bit)
  dec <- linf + (precisao * dec)
  return (dec)
}

RealToBin <- function (dec) {
  bit <- (dec-linf)/precisao
  vet <- paste(round(DecToBin(bit)), collapse = "")
  vet <- str_pad(vet, bitnec, pad = "0")
  return(vet)
}

contaCasasDec <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

splitChunk <- function(x, n) {
  sst <- strsplit(x, '')[[1]]
  m <- matrix('', nrow=n, ncol=(length(sst)+n-1)%/%n)
  m[seq_along(sst)] <- sst
  apply(m, 2, paste, collapse='')
}

convertCromoBinToReal <- function(cromo){
  n <- length(cromo)
  real <- rep(NA, ndim)
  
  strcromo <- ""
  a <- apply(as.matrix(cromo), 2, as.numeric)

  strcromo <- paste(as.character(a), collapse = '')
  strcromo <- splitChunk(strcromo, bitnec)
  
  for (i in 1:ndim)
    real[i] <- BinToReal(strcromo[i])
  
  return(real)
}

convertCromoRealToBin <- function(cromo){
  n <- length(cromo)
  bin <- rep(NA, n)
  for (i in 1:n)
    bin[i] <- RealToBin(cromo[i])
  
  return(bin)
}

# ALGORITHM FUNCTIONS
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
  
  realpop <- matrix(NA, nrow = dim(poplocalbest)[1], ncol = ndim)
  
  for (i in 1:dim(poplocalbest)[1]) {
    realpop[i,] <- convertCromoBinToReal(poplocalbest[i,])
  }
  
  fits <- apply(realpop, 1, fitfunction)
  indbest <- which.min(fits)
  
  return(list(best = poplocalbest[indbest,], fit_best = fits[indbest]))
}

get_stats_pop <- function(poplocalbest) {
  realpop <- matrix(NA, nrow = dim(poplocalbest)[1], ncol = ndim)
  
  for (i in 1:dim(poplocalbest)[1]) {
    realpop[i,] <- convertCromoBinToReal(poplocalbest[i,])
  }
  
  fits <- apply(realpop, 1, fitfunction)
  
  lsstats <- list(meanpop = mean(fits),
                  medianpop = median(fits),
                  sdpop = sd(fits),
                  bestpop = min(fits),
                  worstpop = max(fits))
  
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

bounds_pop <- function(popbounds) {
  realpop <- matrix(NA, nrow = dim(popbounds)[1], ncol = ndim)
  
  for (i in 1:dim(popbounds)[1]) {
    realpop[i,] <- convertCromoBinToReal(popbounds[i,])
  }
  
  correctpop <- matrix(NA, nrow = 0, ncol = size)
  
  for (i in 1:dim(popbounds)[1]) {
    correct <- bounds_cromo(popbounds[i,])
    correctpop <- rbind(correctpop, correct)
  }
  
  return(correctpop)
}

bounds_cromo <- function(correctcromo) {
  cromoReal <- convertCromoBinToReal(correctcromo)

  for (i in 1:length(cromoReal)) {
    if (cromoReal[i] > lsup) {
      auxsup <- strsplit(lsupbin, "")
      auxsup <- as.numeric(auxsup[[1]])
      correctcromo[][(bitnec * i - bitnec + 1):(bitnec * i)] <- auxsup
    } 
    
    if (cromoReal[i] < linf) {
      auxinf <- strsplit(linfbin, "")
      auxinf <- as.numeric(auxinf[[1]])
      correctcromo[][(bitnec * i - bitnec + 1):(bitnec * i)] <- auxinf
    }
  }
  
  return(correctcromo)
}

rankpop <- function(unrankedpop) {
  unrankedrealpop <- matrix(NA, nrow = dim(unrankedpop)[1], ncol = ndim)
  
  for (i in 1:dim(unrankedpop)[1]) {
    unrankedrealpop[i,] <- convertCromoBinToReal(unrankedpop[i,])
  }
  
  fits <- apply(unrankedrealpop, 1, fitfunction)
  
  return(rev(sort(fits, index.return=TRUE)$ix))
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

printMatReal <- function(poplocalbest) {
  poplocalbest <- apply(as.matrix(poplocalbest), 2, as.numeric)
  
  realpop <- matrix(NA, nrow = dim(poplocalbest)[1], ncol = ndim)
  
  for (i in 1:dim(poplocalbest)[1]) {
    realpop[i,] <- convertCromoBinToReal(poplocalbest[i,])
  }
  
  print(realpop)
}