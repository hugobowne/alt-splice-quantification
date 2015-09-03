###nice code

#################################################
###set wd, define libraries
#################################################
rm(list=ls(all=TRUE)) 
setwd("~/repos//alt-splice-quantification/model//scripts")
#setwd("~/splice/complete_cluster/")
require(plyr)
require(ggplot2)
library(foreach)
library(doMC)
library( dplyr )
#################################################
###set sources, set.seed
#################################################
source("../modules/simLL_tot.R")
source("../modules/simul_tot.R")
source("../modules/plotting.R")
source("../modules/analysis_functions.R")


# file = ("../data_grids/zebrafish_grid1.RData")
# par <- heat_analysis(file, TRUE)
#par

alphas <- read.csv( "alpha_mean_maxLL.csv")

exper <- c( 'ar' , 'ce' , 'ch' , 'dr' , 'kc' , 'kp', 'lc' , 'lm' , 'lmb' , 'lmo' , 'lp' ,
            'ls' , 'ly' , 'm8' , 'mb' , 'mc' , 'rb' , 'zf', 'gallus' , 'human' ,
            'lizard' , 'xenopus' , 'pombe' )

filenames <- c('arabidopsis' , 'celegans' , 'chlamy_2' , 'drosophila' , 'k562_chromatin',
               'k562pap' , 'lymphoblast_cambodian' , 'lymphoblast_maya' ,
               'lymphoblast_mbuti' , 'lymphoblast_mozabite' , 'lymphoblast_pathan' ,
               'lymphoblast_san' , 'lymphoblast_yakut' , 'mouse_8_weeks' ,
               'mouse_brain' , 'mouse_CNS' , 'rat_brain' , 'zebrafish' , 'gallus' ,
               'human_brain' , 'lizard' , 'xenopus' , 'pombe')



###initialize data frame of interest
###
mdf <- data.frame(transcriptome = NA , metric = NA , model = NA , value = NA )
#colnames( mdf ) <- c("transcriptome", "value")
i <- 5
add_metrics( mdf , i , M3)


#####
# Start the clock!
ptm <- proc.time()

for (i in 1:23){
  for (model_name in c( 'M1' , 'M2' , 'M3' )){
    for (j in 1:20){
      mdf <- add_metrics( mdf , i , get(model_name) , model_name)
      print( i )
    }
    
  }
}
i
# for (i in 1:23){
#   mdf <- add_metrics( mdf , i , M1)
#   #print(i)
# }
# 
# mdf1 <- mdf
# i
# for (i in 1:23){
#   for (j in 1:20){
#     mdf <- add_metrics( mdf , i)
#     print( c(i,j) )
#     }
# }

# Stop the clock
proc.time() - ptm

mdf <- add_metrics( mdf , 5)


save(mdf1 , file = "measures1.RData")

### stats by transcriptome and measure

#mean( log(dhist$density/mhist$density)^2 )

testdf <- group_by(mdf , transcriptome , metric)

aa <- summarise(testdf , mean(v) , sd(v))
aa$se <- aa$sd/sqrt(20)

m <- ggplot( aa , aes( x = "transcriptome" , y = 'mean(v)' ) )

m + geom_point()

mdf$v <- as.numeric(mdf$value)

metr = "diflog"
df1 <- subset( aa , metric == metr )
df1$m <- df1$`mean(v)`
m <- ggplot( df1 , aes( x = transcriptome , y = m ) )
limits <- aes(ymax = m + se, ymin = m-se )
m + geom_point() + geom_errorbar( limits ) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

