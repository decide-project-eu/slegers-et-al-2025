# Harmonics per weather variable

source('DLM functions FPL.R')
source('Custum get.Gt and get.Ft functions.R')
library(ggplot2)
library(dplyr)

data_selected <- readRDS("data_selected.rds")
data_selected <- data_selected %>%
  dplyr::select(86:114, 313, 315) %>%
  na.omit() %>%
  rename(TotalPrecipitation = WeeklyTotalRain)


# most interesting variables
weather_names <- c("Temp05",
                   "WindSpeed09",
                   "Humidity01",
                   "TotalPrecipitation")



# look at the variation between months

par(mfrow = c(2,2))
for(name in weather_names){
  name_i <- grep(pattern = name, x = colnames(data_selected))
  actualname <- colnames(data_selected)[name_i]
    
  agg <- aggregate(x = data_selected[,name_i], 
                   by=list(data_selected$Month), 
                   FUN=mean, 
                   na.action="omit")
  plot(agg,type="b", main=actualname, xlab = "Month")
  
}



# plot deciles of one weather parameter at a time

par(mfrow=c(3,3))

for(name in weather_names){
  name_i_all <- grep(pattern = name, x = colnames(data_selected))
  for(name_i in name_i_all){
    actualname <- colnames(data_selected)[name_i]
    
    agg <- aggregate(x = data_selected[,name_i], 
                     by=list(data_selected$Month), 
                     FUN=mean, 
                     na.action="omit")
    plot(agg,type="b", main=actualname)
  }
  
}


# Rainfall per year
avg_data <- data_selected %>% group_by(HatchYear, Month) %>%
        summarise(Rain = mean(WeeklyTotalRain),
                  Wind = mean(WindSpeed05))

ggplot(data = avg_data,
       aes(x = Month,
           y = Rain,
           group = factor(HatchYear),
           color = factor(HatchYear))
       ) +
  geom_line(linewidth = 1.5) +
  theme_bw()


# Wind speed per year
ggplot(data = avg_data,
       aes(x = Month,
           y = Wind,
           group = factor(HatchYear),
           color = factor(HatchYear))
       ) +
  geom_line(linewidth = 1.5) +
  theme_bw()



# For each variable (decile) of each weather parameter, fit 1-6 waves
par(mfrow=c(2,3))
par(mar=c(2,4,2,2))

w_monthly <- (2*pi)/12

out_all <- data.frame()

# loop through weather parameters
for(name in weather_names){
  name_i_all <- grep(pattern = name, x = colnames(data_selected))
  # loop through weather variables (deciles)
  for(name_i in name_i_all){
    actualname <- colnames(data_selected)[name_i]
    
    # loop through number of waves
    
    for(n_waves in 1:6){
  
      # Aggregate and plot the data for the relevant data variable against month
      agg <- aggregate(data_selected[,name_i], 
                       by=list(data_selected$Month), 
                       FUN=mean, 
                       na.rm=TRUE)
      colnames(agg) <- c('Month', actualname)
      plot(agg, 
           type='b', 
           xlab='Month', 
           ylab=actualname, 
           main=paste('n_waves:', n_waves)
           )
      
      
      # First, make a linear model describing the relevant variable
      # NOTICE: this function comes from the "Custom get.Gt and get.Ft functions.R" script
      a <- make.LM.wHarmonics(Data = agg, 
                              w = w_monthly, 
                              n_waves, 
                              trend = FALSE, 
                              relevant.name = actualname, 
                              time.var.name = 'Month', 
                              plot.it=FALSE, 
                              remove.zeros=FALSE, 
                              round.by=3)
      
      # Extract the relevant information from the linear model
      LM <- a$LM
      Month <- as.data.frame(x=agg$Month)
      colnames(Month) <- 'Month'
      S <- summary(LM)
      adj.r.squared <- S$adj.r.squared
      r.squared <- S$r.squared
      
      # Save the extracted information to the out_all data frame for later processing
      out <- cbind('n_waves'=n_waves, 'r.squared'=r.squared, 'adj.r.squared'=adj.r.squared)
      out_all <- rbind(out_all, out)
      
      # Add this model's predictions per htime unit to the plot of the aggregated data
      pred <- predict(object = LM, newdata = Month)
      lines(pred~agg$Month, col='red')
  
    }
  }
  
}


# pick best for each variable
n_waves = c(rep(1, 9),  #temp
        rep(1, 9),  # windspeed
        rep(1,9), #hum
        0 #rain
        )


# plot deciles with best number of waves for one parameter at a time

i <- 0
par(mfrow=c(3,3))
par(mar=c(2,4,2,2))

for(name in weather_names){
  name_i_all <- grep(pattern = name, x = colnames(data_selected))
  # loop through weather variables (deciles)
  for(name_i in name_i_all){
      i <- i+1
      actualname <- colnames(data_selected)[name_i]
    
    # Aggregate and plot the data for the relevant data variable against month
      agg <- aggregate(data_selected[,name_i], 
                       by=list(data_selected$Month), 
                       FUN=mean, 
                       na.rm=TRUE)
      colnames(agg) <- c('Month', actualname)
      plot(agg, type='b', xlab='Month', ylab=actualname, main=actualname)
      
      
      # First, make a linear model describing the relevant variable
      # NOTICE: this function comes from the "Custom get.Gt and get.Ft functions.R" script
      a <- make.LM.wHarmonics(Data = agg, 
                              w = w_monthly, 
                              n_waves[i], 
                              trend = FALSE, 
                              relevant.name = actualname, 
                              time.var.name = 'Month', 
                              plot.it=FALSE, 
                              remove.zeros=FALSE, 
                              round.by=3)
      
      # Extract the relevant information from the linear model
      LM <- a$LM
      Month <- as.data.frame(x=agg$Month)
      colnames(Month) <- 'Month'
      
      # Add this model's predictions per htime unit to the plot of the aggregated data
      pred <- predict(object = LM, newdata = Month)
      lines(pred~agg$Month, col='red')
  }
}
