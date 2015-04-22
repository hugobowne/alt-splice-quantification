##plotting functions

cdfplots <- function( a , b){
  clas <- c(rep("model",length( a )) , rep( "data" , length( b ))) #column vector to build dataframe
  df <- data.frame( pmin = c( a , b ), class = clas) ##dataframe with pmin & class (model/experimental)
  
  #Ecdf within species
  
  df.species <- ddply(df, .(class), summarize,
                      pmin1 = unique(pmin),
                      ecdf = ecdf(pmin)(unique(pmin))) ##build a dataframe with ecdf of model/data
  
  ggplot(df.species, aes(pmin1, ecdf, color = class)) + geom_step() +
    theme(legend.justification=c(1,0), legend.position=c(1,0))#this is the plot!
}

cdfplots_shade <- function( a , b ){
  clas <- c(rep("model",length( a )) , rep( "data" , length( b ))) #column vector to build dataframe
  df <- data.frame( pmin = c( a , b ), class = clas) ##dataframe with pmin & class (model/experimental)
  
  #Ecdf within species
  
  df.species <- ddply(df, .(class), summarize,
                      pmin1 = unique(pmin),
                      ecdf = ecdf(pmin)(unique(pmin))) ##build a dataframe with ecdf of model/data
  
  ggplot(df.species, aes(pmin1, ecdf, color = class)) + 
    geom_ribbon(aes(ymin=ecdf + 0.001, ymax=ecdf-0.001),
                alpha=0.2) +
    geom_step() +
    theme(legend.justification=c(1,0), legend.position=c(1,0))#this is the plot!
}

