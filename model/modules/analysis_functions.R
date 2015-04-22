heat_analysis <- function( file , plot = TRUE){
  #remove stuff, set directory, load data
  #rm(list=ls(all=TRUE)) 
  setwd("~/Documents//alt_splice/automated/heat/")
  load( file )
  
  
  #define parameters necessary for heatmap
  my_palette <- colorRampPalette(c("blue" , "yellow" , "red" ))(n = 49)
  LLsf <- unlist( LLs )
  S <- 8
  ddff <- do.call("rbind", replicate(S, ddf, simplify = FALSE))
  ddff$LL <- LLsf
  
  third_q <- function( x ){
    return ( quantile( x , 0.5 ) )
  }
  
  agg1 <- aggregate( x = ddff , by=list( ddff$Var1 , ddff$Var2) , FUN=median)
  #agg1 <- subset( agg1 , Var2 == 0 )
  x <- unique( agg1$Var1 ) 
  y <- unique( agg1$Var2 ) 
  z <-  agg1$LL
  
  #dim( z ) <- c( 41 , 41 )
  
  if (plot == TRUE){
    heatmap( z , Rowv = NA , Colv = NA, labRow = NA , labCol = NA , scale = "none", col=my_palette)
    #heatmap.2( z , Rowv = NA , Colv = NA,  dendrogram='none' , symm=F,symkey=F,symbreaks=F ,
               #col = my_palette)
  
  
  #agg1 <-  subset( agg1, agg1$LL != max( LL ) )
  
  }
  return ( subset( agg1 , LL == max( LL )) )
}

compare_model_plot <- function( knPairs , par ){
  data <- knPairs
  #set.seed(71)
  
  df <- data
  dat <- as.matrix( df$k2/df$n )
  
  N <- par$Var1
  f <- par$Var2
  model <- simulate_fluc_tot( N , f , data )
  cdfplots( model , dat) 
}

compare_model_plot_data <- function( data , par ){
  #data <- knPairs
  #set.seed(71)
  
  df <- data
  dat <- as.matrix( df$k2/df$n )
  
  N <- par$Var1
  f <- par$Var2
  model <- simulate_fluc_tot( N , f , data )
  cdfplots( model , dat) 
}

compare_model_density_data <- function( data , par , xlim = c(0,0.5) , logy =TRUE){
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
  ggplot(pmin1, aes(pmin, fill = type)) + geom_density(alpha = 0.2 ) + coord_cartesian(xlim=xlim) +
    scale_y_log10()} else
    {ggplot(pmin1, aes(pmin, fill = type)) + geom_density(alpha = 0.2 ) + coord_cartesian(xlim=xlim)}
  
  
}
