library(dplyr)

source('DLM functions.R')

data_selected <- readRDS("data_selected_notlagged.rds")
data_selected <- data_selected %>%
  rename_with(~ sub("^TotalPrecipitation", "WeeklyTotalRain", .x),
              starts_with("TotalPrecipitation_"))

weather_names <- c("Temp", 
                   "WindSpeed", 
                   "Humidity",
                   "WeeklyTotalRain",
                   "WeeklyMeanWindDirection")


## for each house, make a df
##    for each variable, make a MA and collect residuals (all separately!)
# After the loop:
##    make one QQ plot per variable
##    plot them 3 x 3 so that all deciles of one parameter are in one view

house_list <- readRDS('house_list.rda')

# make empty list to store residuals
names <- c("res_temp",
           "res_windspeed",
           "res_humidity",
           "res_rain",
           "res_winddir")

# rename variables
all_res <- list()
for (n in names) {
  if ((n == "res_winddir") | (n =="res_rain")) {
    all_res[[n]] <- c(0)
  } else {
    for (number in 1:9) {
      vector_name <- paste0(n, "_", number)
      all_res[[vector_name]] <- c(0)
    }
  }
}


# apply MA to house timeline. Do this for every weather variable

par(mfrow=c(2,1))
for(house in house_list){
  # select data
  df <- data_selected[data_selected$Farmhouse == house, ]
  
  # get data ready
  df <- df %>%
    dplyr::select(86:114) %>% #only weather of week 1
    mutate(
    across(c(where(is.character)), as.numeric)
  ) %>%
    na.omit()
    
  
  # loop through weather variables week 1 (column 17-53)
  for (i in 1:ncol(df)){
    Y <- df[,i]
    
    # From the "DLM functions.R" script, use the "moving.function" function
    MA <- moving.function(x = Y, n = 5, FUN = mean)
    residuals <- na.omit(Y - MA)
    
    all_res[[i]] <- c(all_res[[i]], residuals)
   
  }
  
  rm(df)
}


# QQ plot and R squared (all weather variables, 3x3 plots)

par(mfrow=c(3,3))

for (res in names(all_res)){
  Mean <- mean(all_res[[res]])
  SD <- sd(all_res[[res]])
  a <- sort(all_res[[res]])
  b <- sort(rnorm(n = length(all_res[[res]]), mean = Mean, sd = SD))
  
  LM <- lm(b~a)
  r.squared <- summary(LM)$r.squared
  print(paste(res, 'r.squared =', r.squared))
  
  plot(b ~ a,
       xlab = "residuals",
       ylab = "normal",
       main = res)
  abline(coef = c(0,1))
  text(
    x = min(a), y = max(b) - (1/10 * max(b)),
    labels = bquote(R^2 == .(format(r.squared, digits = 4))),
    pos = 4,
    cex = 1.2
  )
  
}


##################################

# Transform rain

par(mfrow=c(3,1))
hist(as.numeric(data_selected$WeeklyTotalRain))
hist(log(as.numeric(data_selected$WeeklyTotalRain) + 1))
hist(sqrt(as.numeric(data_selected$WeeklyTotalRain)))


res_rain <- c()
res_rain_sqrt <- c()

for(house in house_list){
  # open df
  df <- data_selected[data_selected$Farmhouse == house & !is.na(data_selected$WeeklyTotalRain),] %>%
    mutate(across(c(where(is.character)), as.numeric))
  
  Y <- df$WeeklyTotalRain
  sqrtY <- sqrt(df$WeeklyTotalRain)
  
  # From the "DLM functions.R" script, use the "moving.function" function
  MA <- moving.function(x = Y, n = 5, FUN = mean)
  residuals <- na.omit(Y - MA)
  
  MA2 <- moving.function(x = sqrtY, n = 5, FUN = mean)
  residuals2 <- na.omit(sqrtY - MA2)
  
  res_rain <- c(res_rain, residuals)
  res_rain_sqrt <- c(res_rain_sqrt, residuals2)
  
}

par(mfrow=c(1,2))

# QQ plot and R squared
Mean <- mean(res_rain)
SD <- sd(res_rain)
a <- sort(res_rain)
b <- sort(rnorm(n = length(res_rain), mean = Mean, sd = SD))
par(mfrow=c(2,1))
plot(b~a,
     xlab = "residuals",
     ylab = "normal",
     main = "TotalPrecipitation")
abline(coef = c(0,1))
LM <- lm(b~a)
r.squared <- summary(LM)$r.squared
text(
    x = min(a), y = max(b) - (1/6 * max(b)),
    labels = bquote(R^2 == .(format(r.squared, digits = 4))),
    pos = 4,
    cex = 1.2
)
print(paste('r.squared =', r.squared))

# QQ plot and R squared
Mean <- mean(res_rain_sqrt)
SD <- sd(res_rain_sqrt)
a <- sort(res_rain_sqrt)
b <- sort(rnorm(n = length(res_rain_sqrt), mean = Mean, sd = SD))
LM <- lm(b~a)
r.squared <- summary(LM)$r.squared
print(paste('r.squared =', r.squared))
plot(b~a,
     xlab = "residuals",
     ylab = "normal",
     main = "sqrt(TotalPrecipitation)")
abline(coef = c(0,1))
text(
    x = min(a), y = max(b) - (1/6 * max(b)),
    labels = bquote(R^2 == .(format(r.squared, digits = 4))),
    pos = 4,
    cex = 1.2
)