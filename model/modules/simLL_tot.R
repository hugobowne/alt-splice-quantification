simLL_binomial <- function( N , f , data , br , s = 1){
  values <- as.matrix(  simulate_binomial( N , f , data , s ) )
  histinfo <- hist( values , breaks = br , plot = FALSE )
  histinfo$breaks <- histinfo$breaks[ 1:length( histinfo$breaks) - 1]
  data_pmin <- data$k2/data$n
  yo <- apply( as.matrix(data_pmin) , 1 , function(x)
    log( histinfo$density[ histinfo$breaks == max( histinfo$breaks[ histinfo$breaks <= x ] ) ]/100 ))
  return( sum(yo) )
}

simLL_tot <- function( N , f , data , br , s = 1){
  values <- as.matrix(  simulate_fluc_tot( N , f , data , s ) )
  histinfo <- hist( values , breaks = br , plot = FALSE )
  histinfo$breaks <- histinfo$breaks[ 1:length( histinfo$breaks) - 1]
  data_pmin <- data$k2/data$n
  yo <- apply( as.matrix(data_pmin) , 1 , function(x)
    log( histinfo$density[ histinfo$breaks == max( histinfo$breaks[ histinfo$breaks <= x ] ) ]/100 ))
  return( sum(yo) )
}

medians <- function( vec , LLs , Ns ){
  pm1 <- rep( NA , length( vec ) )
  for (i in seq( 1 , length( vec ))){
    set <- LLs[ Ns == Ns[ i ] ]
    pm1[i] <-  median( set )
  }
  return( pm1 )
}

fit_gauss <- function( vec , pm1 , plot = FALSE){
  # Standardize the cumulative probabilities
  prob <- cumsum( pm1 ); prob <- prob / prob[ length(pm1 )]
  # Compute sum of squared residuals to a fit
  f <- function(q){
    res <- pnorm( vec , q[1], q[2]) - prob
    sum(res * res)
  }
  # Find the least squares fit
  coeff <-(fit <- nlm(f, c(1.2, 1)))$estimate
  if (plot == TRUE){
  # Plot the fit
  plot( vec , prob)
  curve(pnorm( x , coeff[1], coeff[2]), add=TRUE)
  }
  return( coeff )
}