library(dplyr)

source('DLM functions FPL.R')

# Read in the data
data_selected <- readRDS("data_selected_notlagged.rds")


############## FOOTPAD LESIONS ##################
# Goal: see if the residuals of a moving average are Normally distributed
# Loop through each farm house, collect residuals

house_list <- readRDS('house_list.rda')
all_res_c <- c()
for(house in house_list){
  # open df
  df <- data_selected[data_selected$Farmhouse == house, ]
  
  Y <- df$FootpadLesions

  # From the "DLM functions.R" script, use the "moving.function" function
  MA <- moving.function(x = Y, n = 5, FUN = mean)
  residuals <- na.omit(Y - MA)
  
  all_res_c <- c(all_res_c, residuals)
  
  rm(df)
}

# QQ plot and R squared
Mean <- mean(all_res_c)
SD <- sd(all_res_c)
a <- sort(all_res_c)
b <- sort(rnorm(n = length(all_res_c), mean = Mean, sd = SD))
plot(b~a,
     xlab = "residuals",
     ylab = "normal",
     main = "residuals FPL no transformation")
abline(coef = c(0,1))
LM <- lm(b~a)
r.squared <- summary(LM)$r.squared
print(paste('r.squared =', r.squared))


###########
# LOG transform

# house_list is already loaded
all_res <- c()
for(house in house_list[1:10]){
  # open df
  df <- data_selected[data_selected$Farmhouse == house, ]
  
  Y <- log(df$FootpadLesions + 1)
  plot(Y, type='b', xlab = 'cycle', ylab='Footpad lesions')

  # From the "DLM functions.R" script, use the "moving.function" function
  MA <- moving.function(x = Y, n = 5, FUN = mean)
  residuals <- na.omit(Y - MA)
  
  plot
  
  all_res <- c(all_res, residuals)
  
  rm(df)
}

# QQ plot and R squared
Mean <- mean(all_res)
SD <- sd(all_res)
a <- sort(all_res)
b <- sort(rnorm(n = length(all_res), mean = Mean, sd = SD))
plot(b~a,
     xlab = "residuals",
     ylab = "normal",
     main = "Residuals FPL log-transformed")
abline(coef = c(0,1))
LM <- lm(b~a)
r.squared <- summary(LM)$r.squared
print(paste('r.squared =', r.squared))


###########
# SQRT transform

# house_list is already loaded
all_res <- c()
for(house in house_list){
  # open df
  df <- data_selected[data_selected$Farmhouse == house, ]
  
  Y <- sqrt(df$FootpadLesions)

  # From the "DLM functions.R" script, use the "moving.function" function
  MA <- moving.function(x = Y, n = 5, FUN = mean)
  residuals <- na.omit(Y - MA)
  
  all_res <- c(all_res, residuals)
  
  rm(df)
}

hist(all_res)


# QQ plot and R squared
Mean <- mean(all_res)
SD <- sd(all_res)
a <- sort(all_res)
b <- sort(rnorm(n = length(all_res), mean = Mean, sd = SD))
plot(b~a,
     xlab = "residuals",
     ylab = "normal",
     main = "Residuals FPL sqrt-transformed")
abline(coef = c(0,1))
LM <- lm(b~a)
r.squared <- summary(LM)$r.squared
print(paste('r.squared =', r.squared))
