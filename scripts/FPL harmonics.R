#### HOUSE KEEPING ####

# Avoid scientific notation
options(scipen=999)

# Get current directory and set it as your working directory
home.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(home.dir)

# Source the "DLM functions.R" script
source('DLM functions FPL.R')

# Source the "Custom get.Gt and get.Ft functions.R" script
source('Custum get.Gt and get.Ft functions.R')


# Read in the data (of all farms)
data_selected <- readRDS("data_selected.rds")

data_selected$Month <- as.numeric(format(as.Date(data_selected$SlaughterDate), 
                                             "%m"))
  # not to be confused with month without capital, 
  # which is numbered over all years

# Define the name of the variables we want to consider
name <- 'FootpadLesions'
# FPL are not transformed


#### SYSTEMATICALLY TEST ALL NUMBERS OF HARMONIC WAVES ####

# Define the angular frequency we will use for defining our harmonic waves
w_monthly <- (2*pi)/12

# Create an empty data frame to store the results in
out_all <- data.frame()

# Prepare the plotting window
par(mfrow=c(2,3))
par(mar=c(2,2,2,2))

# NOTICE: Here we are testing the utility of using between 1 and 6 waves
#         This is because we have 12 months per year
for(n_waves in 1:6){
  
  # Aggregate and plot the data for the relevant data variable against the relevant time variable
  agg <- aggregate(data_selected[,name], 
                   by=list(data_selected$Month), 
                   FUN=mean, 
                   na.rm=TRUE)
  colnames(agg) <- c('Month', name)
  plot(agg, type='b', xlab='Month', ylab=name, main=paste('n_waves:', n_waves) )
  
  
  # First, make a linear model describing the relevant variable in in average water usage and the diurnal pattern given the hour of the day
  # NOTICE: this function comes from the "Custom get.Gt and get.Ft functions.R" script
  a <- make.LM.wHarmonics(Data = agg, 
                          w = w_monthly, 
                          n_waves, 
                          trend = TRUE, 
                          relevant.name = name, 
                          time.var.name = 'Month', 
                          plot.it=FALSE, 
                          remove.zeros=FALSE, 
                          round.by=3)
  
  # Extract the relevant information from the linear model
  LM <- a$LM
  Month <- as.data.frame(x=agg$Month)
  colnames(Month) <- 'Month'
  S <- summary(LM)
  adj_R_squared <- S$adj_R_squared
  R_squared <- S$R_squared
  
  # Save the extracted information to the out_all data frame for later processing
  out <- cbind('n_waves'=n_waves, 
               'R_squared'=R_squared, 
               'adj_R_squared'=adj_R_squared)
  out_all <- rbind(out_all, out)
  
  # Add this model's predictions per htime unit to the plot of the aggregated data
  pred <- predict(object = LM, newdata = Month)
  lines(pred~agg$Month, col='red')
  
}

# Test: how well does the lm fit the data?
# six example farms

a <- make.LM.wHarmonics(Data = agg, w = w_monthly, n_waves = 1, trend = FALSE, 
                          relevant.name = name, time.var.name = 'Month', 
                          plot.it=FALSE, remove.zeros=FALSE, round.by=3)
LM <- a$LM

# choose random farmhouses
for(n in c(1, 40, 80, 120, 160, 200)){
  house_set <- subset(data_selected, 
                      data_selected$Farmhouse == data_selected$Farmhouse[n] )

  pred <- predict(object = LM, newdata = house_set)

  plot(house_set$FootpadLesions, type="b")
  lines(pred, col="red")
}



#### Find the best number of waves ####
par(mfrow=c(1,1))
par(mar=c(4,4,4,4))
plot(out_all$adj_R_squared ~ out_all$n_wavesaves, 
     type='b', 
     col='blue', 
     ylim=c(0,1),  
     pch=16, 
     main='',
     ylab='Adjusted R^2', 
     xlab='Number of waves')


n_waves <- 1

par(mfrow=c(1,1))
par(mar=c(6,6,6,6))
plot(agg, type='b', xlab='Month', ylab='FPL score', main= NA )



#### Make plot with the final number of waves ####

# Make and plot the monthly mean values
agg <- aggregate(data_selected[,name], 
                 by=list(data_selected$Month),
                 FUN=mean, 
                 na.rm=TRUE)
colnames(agg) <- c('Month', name)

par(mfrow=c(1,1))
plot(agg, type='b', 
     xlab='Month', 
     ylab=name, 
     main=paste(name, n_waves, sep=' | ')
     )

# Make a linear model describing the relevant variable in average FPL and the 
# yearly pattern given the month
a <- make.LM.wHarmonics(Data = agg, w = w_monthly, n_waves, trend = TRUE,
                        relevant.name = name, time.var.name = 'Month',
                        plot.it=FALSE, remove.zeros=FALSE, round.by=3)

LM <- a$LM

lines(predict(LM)~agg$Month, col='red')




