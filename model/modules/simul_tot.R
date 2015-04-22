########################################################
###this function simulates the good shit
### n = total no of reads
### m = number of simulations
### N = fluctuation probabilistic power-law exponent
### f = probability of bonafide splicing
########################################################
#data <- read.csv("pombenk.csv")
# data$pmin <- apply( data , 1 , function(x) x[2]/x[1])
# df <- data
# df.nonzero <- subset( df , df$pmin != 0)
# data <- df.nonzero

########################################################
## we define the quantile function of the power law 
##distribution and draw m RVs from the power law dist. 
##using the method of inverse transform (Smirnov) sampling
########################################################  

#N = 1.58
xmax = 1 #upper bound for power law
xmin = 0.0001 #this is also epsilon (lower bound for power law)
#b = 1 - N
qpower <- function(q , b) (q*(xmax**b-xmin**b) + xmin**b)**(1/b)

simulate_fluc1 <- function( N , data ){
  ###############################################################
  #make nice simulations i like!! and also time it.
  #
  ###############################################################
  b <- 1 - N
  ##here we create matrix s; s[1] is data$n ; s[2] is runif to determine if bonafide splicing occurs;
  ##s[3] is binomial prob. if bf splicing occurs; qpower(x[4] , b) is binom. prob. of splicing from
  ##fluctuations.
  s1 <- rbind( data$n , runif( length(data$n) )  )
  ##apply to matrix s to get vector of pmins
  P <- apply( s1 , 2 , function(x) rbinom( 1 , x[1] , qpower(x[2] , b) )/x[1] )
  #if pmin > 0.5, here we set pmin = 1 - pmin; if pmin > 1, we set pmin = 1;
  #P <- (1-P)*(P>0.5)*(P<1) + P*(P<=0.5) + (P>1)
  return(P)
}

simulate_fluc_tot <- function( N , f , data , s = 1){
  ###############################################################
  #make nice simulations i like!! and also time it.
  #
  ###############################################################
  b <- 1 - N
  ##here we create matrix s; s[1] is data$n ; s[2] is runif to determine if bonafide splicing occurs;
  ##s[3] is binomial prob. if bf splicing occurs; qpower(x[4] , b) is binom. prob. of splicing from
  ##fluctuations.
  s1 <- rbind( rep( data$n , s )  , runif( length(data$n) ) , runif( length(data$n) ) , runif( length(data$n) ) )
  ##apply to matrix s to get vector of pmins
  P <- apply( s1 , 2 , function(x) (rbinom( 1 , x[1] ,(x[2] < f)*x[3] ) + rbinom( 1 , x[1] , qpower(x[4] , b) ) )/x[1])
  #if pmin > 0.5, here we set pmin = 1 - pmin; if pmin > 1, we set pmin = 1;
  P <- (1-P)*(P>0.5)*(P<1) + P*(P<=0.5) + 0*(P>1)
}

simulate_binomial <- function( N , f , data , s = 1){
  ###############################################################
  #make nice simulations i like!! and also time it.
  #
  ###############################################################
  b <- 1 - N
  ##here we create matrix s; s[1] is data$n ; s[2] is runif to determine if bonafide splicing occurs;
  ##s[3] is binomial prob. if bf splicing occurs; qpower(x[4] , b) is binom. prob. of splicing from
  ##fluctuations.
  #s1 <- rbind( rep( data$n , s )  , runif( length(data$n) ) , runif( length(data$n) ) , runif( length(data$n) ) )
  ##apply to matrix s to get vector of pmins
  #P <- apply( s1 , 2 , function(x) (rbinom( 1 , x[1] ,(x[2] < f)*x[3] ) + rbinom( 1 , x[1] , qpower(x[4] , b) ) )/x[1])
  #if pmin > 0.5, here we set pmin = 1 - pmin; if pmin > 1, we set pmin = 1;
  P <- apply( data$n , 2 , function(x) rbinom( 1 , x , N))
  P <- (1-P)*(P>0.5)*(P<1) + P*(P<=0.5) + 0*(P>1)
}

# # Start the clock!
# ptm <- proc.time()
# mod <- simu1( 1.58 , 0.04 , data)
# # Stop the clock
# proc.time() - ptm