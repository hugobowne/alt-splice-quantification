rm(list=ls(all=TRUE)) 
setwd("~/repos//alt-splice-quantification/model//scripts")
require( ggplot2 )


load( "~/Documents//alt_splice/automated/data/pombe.RData")
load( "~/Documents//alt_splice/automated/data/lymphoblast_san.RData")
load( "~/Documents//alt_splice/automated/data/k562_chromatin.RData")
df <- knPairs
dat <- as.matrix( df$k2/df$n )

data <- dat[dat != 0]
# h <- hist( data )
# 
# min( data )
# 
# plot( h , log = "xy")
# 
# d <- as.data.frame( data )
# 
# ggplot( d, aes(x = data)) + geom_histogram( binwidth = 0.05) + scale_x_log10() + scale_y_log10()

# x <- data
# h <- hist( x )
# plot(h$count, log = "xy")


h <- hist( data , breaks = 30)
h1 <- as.data.frame(cbind(h$breaks[1:length(h$breaks)-1] +0.01, h$counts))
ggplot( h1, aes(x = V1, y = V2)) + geom_point() + scale_x_log10() + scale_y_log10()

h1$V3 <- log10(h1$V1)
h1$V4 <- log10(h1$V2)
ggplot( h1, aes(x = V3, y = V4)) + geom_point() 


load( "~/Documents//alt_splice/automated/data/arabidopsis.RData")
df <- knPairs
dat <- as.matrix( df$k2/df$n )
data <- dat[dat != 0]

h <- hist( data , breaks = 30)
h2 <- as.data.frame(cbind(h$breaks[1:length(h$breaks)-1] +0.01, h$counts))
h2$V3 <- log10(h2$V1)
h2$V4 <- log10(h2$V2)
ggplot( h2, aes(x = V3, y = V4)) + geom_point() + geom_smooth(method='lm',formula=y~x)

load( "~/Documents//alt_splice/automated/data/lymphoblast_all.RData")
df <- knPairs
dat <- as.matrix( df$k2/df$n )
data <- dat[dat != 0]

h <- hist( data , breaks = 30)
h3 <- as.data.frame(cbind(h$breaks[1:length(h$breaks)-1] +0.01 , h$counts))
h3$V3 <- log10(h3$V1)
h3$V4 <- log10(h3$V2)
ggplot( h3, aes(x = V3, y = V4)) + geom_point() + geom_smooth(method='lm',formula=y~x)

h1$type <- 'k562'
h2$type <- 'arabidopsis'
h3$type <- 'lymphoblast'
d <- rbind( h1 , h2 , h3)

ggplot( d, aes(x = V3, y = V4 , colour = type)) + geom_point( size = 5 ) +
  geom_smooth(method='lm',formula=y~x) +
  theme( legend.position= c( 0.85 , 0.75) , 
         axis.text.y=element_blank() , axis.text.x=element_blank(),
         axis.line=element_blank() , axis.title.x=element_blank() , 
         axis.title.y=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0))

filename <- "../plots/power_hist.eps"
ggsave( file = filename )


plot( h , log="y", type='h')


#OR



cdf_log_log <- function( file ){
  load( file )
  df <- knPairs
  dat <- as.matrix( df$k2/df$n )
  
  data <- dat[dat != 0]
  P <- ecdf( data )
  
  #plot( P )
  
  x = seq( 0 , 0.5 , by = 0.0001)
  x = unique( data )
  Pf <- as.data.frame( cbind( x , P( x ) ))
  
  ggplot( Pf , aes( x = x , y = V2)) + geom_point()+ scale_x_log10() + scale_y_log10() +
    annotation_logticks(base = 10) +
    theme( axis.text.y=element_blank() , axis.text.x=element_blank(),
           axis.line=element_blank() , axis.title.x=element_blank(),
           axis.title.y=element_blank())
  
  
  
}

cdf_log_log("~/Documents//alt_splice/automated/data/chlamy_2.RData")



filenames <- c('arabidopsis' , 'celegans' , 'chlamy_2' , 'drosophila' , 'k562_chromatin',
               'k562pap' , 'lymphoblast_cambodian' , 'lymphoblast_maya' ,
               'lymphoblast_mbuti' , 'lymphoblast_mozabite' , 'lymphoblast_pathan' ,
               'lymphoblast_san' , 'lymphoblast_yakut' , 'mouse_8_weeks' ,
               'mouse_brain' , 'mouse_CNS' , 'rat_brain' , 'zebrafish' , 'gallus' ,
               'human_brain' , 'lizard' , 'xenopus' , 'pombe')


for (i in 1:length(filenames)){
  fn <- paste("../data/" , filenames[i] , ".RData" , sep = "")
  load( fn )
  
  
  save_file <- paste( "../power_plots/" , filenames[i] , ".eps" , sep = "")
  
  yoidel <- cdf_log_log( fn )
  ggsave(file= save_file)
}

#########
##doing this stuff on same plot
##############

cdf_info <- function( file , name){
  load( file )
  df1 <- knPairs
  dat <- as.matrix( df1$k2/df1$n )
  data <- dat[dat != 0]
  P <- ecdf( data )
  
  #plot( P )
  
  x = seq( 0 , 0.5 , by = 0.0001)
  x = unique( data )
  Pf1 <- as.data.frame( cbind( x , P( x ) ))
  Pf1$species <- name
  return( Pf1 )
  
}

Pf1 <- cdf_info( "../data/human_brain.RData" , 'brain')
Pf2 <- cdf_info( "../data/k562pap.RData" , 'k562')
Pf3 <- cdf_info( "../data/chlamy_2.RData" , 'chlamy')

df <- rbind( Pf1 , Pf2 , Pf3 )

ggplot( df , aes( x = x , y = V2)) + geom_point( aes( colour=species) )+ 
  scale_x_log10() + scale_y_log10() +
  theme( legend.position= c( 0.85 , 0.75) , 
         axis.text.y=element_blank() , axis.text.x=element_blank(),
         axis.line=element_blank() , axis.title.x=element_blank() , 
         axis.title.y=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1))

filename <- "../plots_05_08/mult_experiments_pow.eps"
ggsave( file = filename )



subset(df, df$x < 1)






Pf1 <- cdf_info( "../data/k562pap.RData" , 'k562')
Pf2 <- cdf_info_sub( "../data/k562pap.RData" , 'k562_sub_100' , 100)
Pf3 <- cdf_info_sub( "../data/k562pap.RData" , 'k562_sub_500' , 500)
Pf4 <- cdf_info_sub( "../data/k562pap.RData" , 'k562_sub_10000' , 1000)

df <- rbind( Pf1 , Pf2 , Pf3 , Pf4 )

ggplot( df , aes( x = x , y = V2)) + geom_point( aes( colour=species) )+ 
  scale_x_log10() + scale_y_log10() +
  theme( legend.position= c( 0.85 , 0.75) , 
         axis.text.y=element_blank() , axis.text.x=element_blank(),
         axis.line=element_blank() , axis.title.x=element_blank() , 
         axis.title.y=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1))

filename <- "../plots_05_08/k562_pow.eps"
ggsave( file = filename )

load( "../data/k562pap.RData")

cdf_info_sub <- function( file , name , m){
  load( file )
  df1 <- knPairs
  df1 <- subset(df1 , n > m)
  dat <- as.matrix( df1$k2/df1$n )
  data <- dat[dat != 0]
  P <- ecdf( data )
  
  #plot( P )
  
  x = seq( 0 , 0.5 , by = 0.0001)
  x = unique( data )
  Pf1 <- as.data.frame( cbind( x , P( x ) ))
  Pf1$species <- name
  return( Pf1 )
  
}
