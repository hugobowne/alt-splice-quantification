#################################################
###set wd, define libraries
#################################################
#set your working dir. here
setwd("~/repos/alt-splice-quantification/model/scripts/")
require(plyr)
library(foreach)
library(doMC)
#################################################
###set sources, set.seed
#################################################
source("../modules/simLL_tot.R")
source("../modules/simul_tot.R")
source("../modules/plotting.R")

args <- commandArgs(trailingOnly = TRUE)


#################################################
###we define the function to run the experiment:
###the input consists of (n,k) pairs;
###the output is the posterior, the data, the runtime.
#################################################
sample_sim <- function( file , outputname , sim_num = 5, sim_frac = 10 ){
  
  #################################################
  ###load data, shape data, select data
  #################################################
  for ( i in 1:sim_num){
  load( file )
  data <- knPairs
  data <- data[,c(1,3)]
  data_pmin <- data$k2/data$n
  
  ##subsample here below:
  m <- sim_frac
  mydata <- data
  data <- mydata[sample(1:nrow(mydata), nrow(mydata)/m,
                        replace=TRUE),]
  data_pmin <- data$k2/data$n
  
  
  #################################################
  ###run simulations to get LL (max. likelihood)
  #################################################
  #set key parameters
  br <- 20 #number of bins
  N <- seq( 1.01 , 1.99, length.out = 41 ) #grid of N's
  f <- seq( 0.00 , 0.5 , length.out = 41) #grid of f's
  T <- 8 #number of simulations for each run
  
  #N <- c(1.29 , 1.47)
  #f <- c(0 , 0.06)
  
  ddf <- expand.grid( N , f ) #grid of possible values for N , f
  ddff <- do.call("rbind", replicate(T, ddf, simplify = FALSE))
  
  
  
  
  registerDoMC(cores=T) #setup
  print("cores registered.")
  print("ready, set, go!")
  # Start the clock!
  ptm <- proc.time()
  LLs <- foreach(n = 1:T) %dopar% apply( ddf ,  1 , function(x) simLL_tot(x[1] , x[2] , data , br , s = 5))
  # Stop the clock ~11645.70s
  t <- proc.time() - ptm #for drosophila, a single point for s=10, t=10 takes ~124sec.
  filename = paste(outputname , sprintf("%1.0f",i) , ".RData" , sep = "")
  save( data , ddf , LLs , t , file = filename)
  
  
  }
  
  
}

#################################################
###Now we run the experiment
#################################################
sample_sim( args[1] , args[2] )

#sample_sim("../data/zebrafish.RData" , "lmm")

