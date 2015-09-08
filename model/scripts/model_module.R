M1 <- function( data , N , xlim = c(0,0.5) , logy =TRUE){
  #data <- knPairs
  #set.seed(71)
  
  df <- data
  dat <- as.matrix( df$k2/df$n )
  
  N <- sum( df$k2)/sum( df$n)
  #f <- par$Var2
  model <- simulate_binomial( N , 0 ,data)
  
  model <- as.data.frame(model)
  model$type <- 'model'
  colnames(model)[1] <- 'pmin'
  colnames(dat)[1] <- 'pmin'
  
  dat <- as.data.frame(dat)
  dat$type <- 'data'
  data <- dat
  
  #and combine into your new data frame vegLengths
  pmin1 <- rbind(data , model )
  
  #now make your lovely plot
  if (logy == TRUE){
    ggplot(pmin1, aes(pmin, colour = type)) + geom_density(alpha = 0.2 ) + coord_cartesian(xlim=xlim) +
      scale_y_log10() + annotation_logticks() + theme(legend.position= "none",
                                                      #axis.text.y=element_blank() , axis.text.x=element_blank(),
                                                      axis.line=element_blank() , axis.title.x=element_blank() , 
                                                      axis.title.y=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1))} else
                                                      {ggplot(pmin1, aes(pmin, fill = type)) + geom_density(alpha = 0.2 ) + coord_cartesian(xlim=xlim)}
  return( pmin1)
  
  
}

M2 <- function( data ,  par , xlim = c(0,0.5) , logy =TRUE){
  #data <- knPairs
  #set.seed(71)
  
  df <- data
  dat <- as.matrix( df$k2/df$n )
  
  N <- sum( df$k2)/sum( df$n)
  #f <- par$Var2
  model <- simulate_M2( data)
  
  model <- as.data.frame(model)
  model$type <- 'model'
  colnames(model)[1] <- 'pmin'
  colnames(dat)[1] <- 'pmin'
  
  dat <- as.data.frame(dat)
  dat$type <- 'data'
  data <- dat
  
  #and combine into your new data frame vegLengths
  pmin1 <- rbind(data , model )
  
  #now make your lovely plot
  if (logy == TRUE){
    ggplot(pmin1, aes(pmin, colour = type)) + geom_density(alpha = 0.2 ) + coord_cartesian(xlim=xlim) +
      scale_y_log10() + annotation_logticks() + theme(legend.position= "none",
                                                      axis.text.y=element_blank() , axis.text.x=element_blank(),
                                                      axis.line=element_blank() , axis.title.x=element_blank() , 
                                                      axis.title.y=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1))} else
                                                      {ggplot(pmin1, aes(pmin, fill = type)) + geom_density(alpha = 0.2 ) + coord_cartesian(xlim=xlim)}
  return( pmin1)
  
  
}

M3 <- function( data , par , xlim = c(0,0.5) , logy =TRUE){
  #data <- knPairs
  #set.seed(71)
  
  df <- data
  dat <- as.matrix( df$k2/df$n )
  N <- par$Var1
  f <- par$Var2
  model <- simulate_fluc_tot( N , f , data )
  
  model <- as.data.frame(model)
  model$type <- 'model'
  colnames(model)[1] <- 'pmin'
  colnames(dat)[1] <- 'pmin'
  
  dat <- as.data.frame(dat)
  dat$type <- 'data'
  data <- dat
  
  #and combine into your new data frame vegLengths
  pmin1 <- rbind(data , model )
  
  #now make your lovely plot
  if (logy == TRUE){
    ggplot(pmin1, aes(pmin, colour = type)) + geom_density(alpha = 0.2 ) + coord_cartesian(xlim=xlim) +
      scale_y_log10() + annotation_logticks() + theme(legend.position= c( 0.85 , 0.85),
                                                      axis.text.y=element_blank() , axis.text.x=element_blank(),
                                                      axis.line=element_blank() , axis.title.x=element_blank() , 
                                                      axis.title.y=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1))} else
                                                      {ggplot(pmin1, aes(pmin, fill = type)) + geom_density(alpha = 0.2 ) + coord_cartesian(xlim=xlim)}
  
  return( pmin1)
}

add_metrics <- function( mdf , i , model, model_name){
  #i <- 5
  fn <- paste("../data/" , filenames[i] , ".RData" , sep = "")
  load( fn )
  
  name <- c(filenames[i] , exper[i])
  alpha <- subset(alphas , species ==exper[i])$alpha
  col <- c( 'Var1' , 'Var2' )
  par <- as.data.frame( cbind(alpha, 0 ))
  colnames( par ) <- col
  
  dm <- model( knPairs , par , xlim = c(0 , 0.5), logy = TRUE)
  
  d <- subset(dm, type == 'data')
  m <- subset(dm, type == 'model')
  
  #stat_density(as.matrix(m))
  
  #dens <- density(d$pmin ) #density of data and model
  #plot( dens )
  
  dhist <- hist(d$pmin , breaks = seq(0,0.5,0.01) ) #data hist
  mhist <- hist(m$pmin , breaks = seq(0,0.5,0.01)) #model hist
  
  #mo <- model
  
  dif <- sum((dhist$counts - mhist$counts)^2/(dhist$counts + 1)^2)/length(dhist$counts)
  mdf <- rbind( mdf , c( name[1] , "dif" , model_name , dif ))
  
  diflog <- sum(log(dhist$density/mhist$density)^2)/length(dhist$density)
  mdf <- rbind( mdf , c( name[1] , "diflog" , model_name , diflog ))
  
  kl <- sum(dhist$density*log(dhist$density/mhist$density))
  mdf <- rbind( mdf , c( name[1] , "kl" , model_name , kl ))
  
  return( mdf )
  
}


add_rsq <- function( mdf , i , model, model_name){
  #i <- 5
  fn <- paste("../data/" , filenames[i] , ".RData" , sep = "")
  load( fn )
  
  name <- c(filenames[i] , exper[i])
  alpha <- subset(alphas , species ==exper[i])$alpha
  col <- c( 'Var1' , 'Var2' )
  par <- as.data.frame( cbind(alpha, 0 ))
  colnames( par ) <- col
  
  dm <- model( knPairs , par , xlim = c(0 , 0.5), logy = TRUE)
  
  d <- subset(dm, type == 'data')
  m <- subset(dm, type == 'model')
  
  #stat_density(as.matrix(m))
  
  #dens <- density(d$pmin ) #density of data and model
  #plot( dens )
  
  dhist <- hist(d$pmin , breaks = seq(0,0.5,0.01) , plot = FALSE) #data hist
  mhist <- hist(m$pmin , breaks = seq(0,0.5,0.01) , plot = FALSE) #model hist
  
  #mo <- model
  
  dc <- dhist$counts
  mc <- mhist$counts
  dme <- mean( dc )
  mme <- mean( mc )
  
  
  denom <- (dc - dme)^2
  num <- ( ( dc - mc ) - ( dme - mme ))^2
  denom <- (dc)^2
  num <- ( ( dc - mc ))^2
  rsq <- sum(num)/sum(denom)
  mdf <- rbind( mdf , c( name[1] , "rsq" , model_name , 1 - rsq ))
  
  
#   dif <- sum((dhist$counts - mhist$counts)^2/(dhist$counts + 1)^2)/length(dhist$counts)
#   mdf <- rbind( mdf , c( name[1] , "dif" , model_name , dif ))
#   
#   diflog <- sum(log(dhist$density/mhist$density)^2)/length(dhist$density)
#   mdf <- rbind( mdf , c( name[1] , "diflog" , model_name , diflog ))
#   
#   kl <- sum(dhist$density*log(dhist$density/mhist$density))
#   mdf <- rbind( mdf , c( name[1] , "kl" , model_name , kl ))
  
  return( mdf )
  
}