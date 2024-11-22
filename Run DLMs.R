
## This script is for training + running the DLM with either EM algorithm or discount factor


#### HOUSE KEEPING ####

library(dplyr)
library(ggplot2)

# Avoid scientific notation
options(scipen=999)


### PREPARE DATA ###

# Get relevant functions
source('DLM functions FPL.R')
source('Custum get_Gt and get_Ft functions.R')
source('FUNCTIONS.R')

# Read in the data 
data_month <- readRDS("data_month.rds")

# sqrt-transform rainfall
data_month$WeeklyTotalRain <- sqrt(data_month$WeeklyTotalRain)

# - define the meta data
metadata_columns = c("FarmIdentification", "Farmhouse", 
                      "Flock", "MonthTotal", "Month",
                      "Switched", "Type2", "Breed",
                      "AntibioticsWeek1",
                      "MeanAgeAtSlaughter", "EndFlockSize", "FlockSizeDiff",
                      "DaysBetweenFlocks", "FractionThinned",
                      "Province", "Score", "BuildingAngle", "NumberOfHouses",
                      "MeanNDVTiterLag1", "AntibioticsLag1", "MortalityLag1",                 
                      "WindDir", "AngleDiffWk1", "FeedPrice")
# Month is the hatch date

# Define the variable we especially care about
main_variable <- "FootpadLesions"


# A function for calculating RMSE from a vector of forecast errors
RMSE <- function(x){
  return(sqrt(mean((x)^2, na.rm = TRUE)))
}



##### 20-fold CV ####

# 1 get list of farms, assigning number 1-20
set_seed(42)
N <- length(unique(data_month$FarmIdentification))
farms <- unique(data_month$FarmIdentification)
data_month$Fold <- NA
for(i in 1:20){
  farmhouses <- sample(x = farms, 
                       size = min(0.05*N, length(farms)),
                       replace = FALSE)
  farmhouses.i <- which(data_month$FarmIdentification %in% farmhouses)
  data_month$Fold[farmhouses.i] <- i
  farms <- farms[-which(farms %in% farmhouses)]
  print(length(farms))
}
table(data_month$Fold)


df <- data_month
df <- df[order(df$FarmIdentification, df$Farmhouse, df$MonthTotal), ]



####################################################
#### RUN SPECIFIC DLM
####################################################


# Select which type of model to run
# model <- "uni_nowave"
# model <- "uni_wave"
model <- "multi_wave"
# model <- "multi_nowave"

method <- "EM"
# method <- "discount"


#### Set up the hyper-parameters for the DLM #### 

# Define relevant parameters for the DLM 

# The variables to be included in the multi-variate DLM
if (model == "multi_nowave" | model == "multi_wave"){
  relevant_names <- c("FootpadLesions", "Temp05", "Humidity01",
                      "WeeklyTotalRain", "WindSpeed09")
} else if (model == "uni_nowave" | model == "uni_wave"){
  relevant_names <- c("FootpadLesions")
}


# The number of harmonic waves used to describe each variable
if (model == "uni_nowave"){
  n_waves = c(0)
} else if (model == "uni_wave"){
  n_waves = c(1)
} else if (model == "multi_nowave"){
  n_waves = c(0, 1, 1, 0, 1) # FPL and rain have no wave
} else if (model == "multi_wave"){
  n_waves =c (1, 1, 1, 0, 1) # rain has no wave
}


# Booleans defining whether a linear trend should be added to the model 
# for each variable
trend = rep(FALSE, length(relevant_names))

# The angular frequencies (omega) for the harmonic waves describing each variable
w = rep((x = (2*pi)/12), length(relevant_names))

# The name of the time-variable
time_var = rep('Month', length(relevant_names))


# Make an empty data frame to store the results in
out_all <- data.frame()
et_complete <- data.frame()
agg_all <- data.frame()

start_time <- Sys.time()


res_out_all <- data.frame() 

# Now iterate over the farms in this farm 
# - we are doing per-farm cross-validation, so each farm will be the test set 
# once, while the remaining farms will be the learning set

for(fold in 1:20){
  # Get the test set and learning set
  test_set <- subset(df, df$Fold == fold)
  learning_set <- subset(df, df$Fold != fold)
  
  #  Standardize learning set and test set based on learning set values
  # par(mfrow=c(length(relevant_names),3))
  for(name in relevant_names){
    # hist(learning_set[[name]], main=name)
    SD <- sd(na.omit(learning_set[[name]]))
    Mean <- mean(na.omit(learning_set[[name]]))
    if(name == "FootpadLesions"){
      SDFootpadLesions <- sd(na.omit(learning_set[[name]]))
      MeanFootpadLesions <- mean(na.omit(learning_set[[name]]))
    }
    
    learning_set[[name]] <- (learning_set[[name]] - Mean)/SD
    # hist(learning_set[[name]], main=name)
    
    test_set[[name]] <- (test_set[[name]] - Mean)/SD
    # hist(test_set[[name]], main=name)
  }
  
  Threshold <- (80 - MeanFootpadLesions)/SDFootpadLesions

  # Define the initial prior distribution
  a <- get_mu0_and_C0(data = learning_set,
                      relevant_names,
                      time_var,
                      trend, 
                      n_waves,
                      w)
  
  Gt <- get_Gt()
  mu0 <- a$mu0
  C0 <- a$C0
  # Optional: update mu0 to correspond with the start time of the data
  # mu0 <- update.mu0(mu0, Gt, start_time=1)
  # this is not necessary in this case because all farm house time series start 
  # at the correct month, instead of month 1
  
  
  if(method == "EM"){
      # Run the EM-algorithm to estimate the observational and system variance
      start_time <- Sys.time()
      set_seed(111)
      EM_out <- runEM_earlyStopping(data = learning_set, 
                                    stratify.by = NA,
                                    Spline.list = NA,
                                    identifyer = 'Farmhouse',
                                    V0=NA, W0=NA,
                                    C0.list = list(C0),
                                    mu0.list = list(mu0),
                                    no.better.limit = 2, 
                                    time.var = time_var, 
                                    relevant_names,
                                    round.by = 4)
      print(Sys.time() - start_time)
      
      # Extract the mu0, C0, V, and W from their respective lists in the output 
      # from the EM-agorithm
      mu0 <- EM_out$mu0.list[[1]]
      C0 <- EM_out$C0.list[[1]]
      V <- EM_out$V.list[[1]]
      W <- EM_out$W.list[[1]]
      steps <- EM_out$total.steps
      
      saveRDS(object = mu0, file = paste0('learning_values/', 
                                          model, '_EM_mu0_', 
                                          fold, '.rds',
                                          sep = ""))
      saveRDS(object = C0, file = paste0('learning_values/', 
                                         model, '_EM_C0_', fold, '.rds',
                                          sep = ""))
      saveRDS(object = V, file = paste0('learning_values/', 
                                        model, '_EM_V_', fold, '.rds',
                                          sep = ""))
      saveRDS(object = W, file = paste0('learning_values/', 
                                        model, '_EM_W_', fold, '.rds',
                                          sep = ""))
      saveRDS(object = steps, file = paste0('learning_values/', 
                                            model, '_EM_steps_', fold, '.rds',
                                          sep = ""))
  }
  
  
  if(method == "discount"){
      # Get V
      par(mfrow = c(2,3))
      V <- get_V(data = learning_set,
                 identifyer = "Farmhouse",
                 stratify.by = NA,
                 time.var = time_var,
                 relevant_names = relevant_names)$V_1
    
      # Get the discount factor (delta)
      start_time <- Sys.time()
      set_seed(111)
      delta <- optimize.delta(
        deltas = seq(from=0.5, to=1, by=0.01),
        data = learning_set,
        identifyer = "Farmhouse",
        stratify.by = NA,
        mu0.list = list(mu0),
        C0.list = list(C0),
        V.list = list(V),
        relevant_names, 
        Spline.list = NA,
        time.var = time_var
      )[[1]]
      
      print(Sys.time() - start_time)
      
      saveRDS(object = mu0, file = paste0('learning_values/', 
                                          model, '_discount_mu0_', fold, '.rds',
                                          sep = ""))
      saveRDS(object = C0, file = paste0('learning_values/', 
                                         model, '_discount_C0_', fold, '.rds',
                                          sep = ""))
      saveRDS(object = V, file = paste0('learning_values/', 
                                        model, '_discount_V_', fold, '.rds',
                                          sep = ""))
      saveRDS(object = delta, file = paste0('learning_values/', 
                                            model, '_discount_delta_', fold, '.rds',
                                          sep = ""))
      saveRDS(object = time, file = paste0('learning_values/', 
                                           model, '_discount_time_', fold, '.rds',
                                          sep = ""))
  }
  
  
  # Apply the DLM to the test set- one house at a time!

  et_test_set <- data.frame()
  e_back_test_set <- data.frame()

  house_count <- 0
  backtransformed_data <- data.frame()

  par(mfrow=c(1,1))

  houses <- unique(test_set$Farmhouse)
  for(house in houses){
    gc()

    # Calculate the progress, so we can see how far along we are
    progress <- round(which(houses == house) / length(houses)*100, 2)
    print(paste('Applying DLM to all houses in test set', '|', progress, '%'))

    # Get the house set
    house_set <- subset(test_set, test_set$Farmhouse == house)

    if(method == "discount"){
      res <- runDLM(data = house_set,
                    mu0 = mu0,
                    C0 = C0,
                    V = V,
                    W = NA,
                    adjust.W = TRUE,
                    delta = delta,
                    relevant_names,
                    Spline.list = NA,
                    time.var = time_var)
      }

    if(method == "EM"){
      res <- runDLM(data=house_set,
                    mu0 = mu0,
                    C0 = C0,
                    V = V,
                    W = W,
                    adjust.W = TRUE,
                    delta = NA,
                    relevant_names,
                    Spline.list = NA,
                    time.var = time_var)
      }


    # save the DLM output to data frame
    res_out <- extract.res(res)
    
    # Change column names: remove wave count because we only have 1 wave per variable
    colnames(res_out)[8:9] = c('MtWaveCosTemp', 'MtWaveSinTemp')
    colnames(res_out)[11:12] = c('MtWaveCosHum', 'MtWaveSinHum')
    colnames(res_out)[15:16] = c('MtWaveCosWind', 'MtWaveSinWind')

    metadata <- house_set[ , metadata_columns]

    # Get observed values (Yt)
    Yt_all <- data.frame(t(sapply(res$Yt,c)))
    if(nrow(Yt_all) == 1){
      Yt_all <- t(Yt_all)
    }
    #back-transform forecast errors
    Yt_back_all <- Yt_all[ , 1]
    Yt_back_all <- Yt_back_all * SDFootpadLesions + MeanFootpadLesions

    # Get forecasts (predicted values) (Ft)
    ft_all <- data.frame(t(sapply(res$ft,c)))
    if(nrow(ft_all) == 1){
      ft_all <- t(ft_all)
    }
    #back-transform forecasts
    ft_back_all <- ft_all * SDFootpadLesions + MeanFootpadLesions
    ft_back_all <- ft_back_all[ , 1]

    metadata$Yt_back_all <- Yt_back_all
    metadata$ft_back_all  <- ft_back_all
    metadata$Fold <- fold
    res_out <- cbind(metadata, res_out)

    weather_names <- c('Temp05',
            'Humidity01',
            'WeeklyTotalRain',
            'WindSpeed09',
            'MtTemp05',
            'MtWaveCosTemp',
            'MtWaveSinTemp',
            'MtHumidity01',
            'MtWaveCosHum',
            'MtWaveSinHum',
            'MtWeeklyTotalRain',
            'MtWindSpeed09',
            'MtWaveCosWind',
            'MtWaveSinWind',
            'FtTemp05',
            'FtHumidity01',
            'FtWeeklyTotalRain',
            'FtWindSpeed09',
            'EtTemp05',
            'EtHumidity01',
            'EtWeeklyTotalRain',
            'EtWindSpeed09',
            'UtTemp05',
            'UtHumidity01',
            'UtWeeklyTotalRain',
            'UtWindSpeed09')
    weather_cols <- which(colnames(res_out) %in% weather_names)


    # make a version without NA flocks, add back later and sort by month
    res_out_na <- res_out %>%
      filter(is.na(Flock))

    res_out <- res_out %>%
      filter(!is.na(Flock)) %>%
    # lag the weather columns, so that the weather is of the current flock
      mutate_at(vars(weather_cols), ~lag(., n = 1)) %>%
      bind_rows(res_out_na) %>%
      arrange(MonthTotal)

    res_out_all <- rbind(res_out_all , res_out)

    # Add to backtransformed data
    backtransformed_data <- rbind(
      backtransformed_data,
      cbind(
        ft_back_all,
        Yt_back_all,
        house
      )
    )

    # Get back-transformed forecast errors

    e_back_all <- Yt_back_all - ft_back_all
    # Add the meta data for this house
    e_back_all <- cbind(house_set[,metadata_columns], e_back_all)
    e_back_test_set <- rbind(e_back_test_set, e_back_all)


    # Get the forecasts errors (untransformed) for this house

    e_all <- data.frame(t(sapply(res$et,c)))
    if(nrow(e_all) == 1){
      e_all <- t(e_all)
    }
    colnames(e_all) <- paste('Et', relevant_names, sep='')
    rownames(e_all) <- NULL
    # Add the meta data for this house
    e_all <- cbind(house_set[,metadata_columns], e_all)
    et_test_set <- rbind(et_test_set, e_all)


    # Plot 15 houses with DLM per fold
    if (house_count %in% c(1:15)){
      house_set$FPL_back <- house_set$FootpadLesions * SDFootpadLesions + MeanFootpadLesions
      filename <- paste('plot_DLM_per_house/', model,
                        '_', method, '_', fold, '_', house_count, '.pdf',
                        sep = "")
      pdf(filename)
      plot(house_set$FPL_back,
           type = "b",
           main = paste('house from fold', fold),
           ylim = c(-10, 235)
           )
      lines(ft_back_all, col='red')
      lines(Yt_back_all, col='blue')
      switched_i <- which(house_set$Switched == 1)
      if(length(switched_i > 0)){
        # rug(switched_i, col='red', lwd=3)
        abline(v=switched_i, col='blue', lty=2)
      }
      dev.off()

      # also in the plot viewer
      plot(house_set$FPL_back,
           type = "b",
           main = paste('house from fold', fold),
           ylim = c(-10, 235)
           )
      lines(ft_back_all, col='red')
      lines(Yt_back_all, col='blue')
      if(length(switched_i > 0)){
        # rug(switched_i, col='red', lwd=3)
        abline(v=switched_i, col='blue', lty=2)
      }
    }
    house_count <- house_count + 1

  }

  # Calculate the overall RMSE for this test set
  RMSE_test_set <- RMSE(et_test_set[,paste0('et_', main_variable)])
  RMSE_test_set_back <- RMSE(e_back_test_set$e_back_all)

  # Look at how the RMSE depends on the number of observations we have so far
  # seen of a given house
  agg <- aggregate(x = e_back_test_set[,'e_back_all'],
                   by=list(e_back_test_set$Flock),
                   FUN=RMSE)
  colnames(agg) <- c('Flock', 'RMSE')

  agg_month <- aggregate(x = e_back_test_set[,'e_back_all'],
                   by=list(e_back_test_set$MonthTotal),
                   FUN=RMSE)
  colnames(agg_month) <- c('MonthTotal', 'RMSE')

  # Identify when the model has (on average) adapted to the specific data set
  Mean_RMSE <- mean(agg$RMSE, na.rm = TRUE)
  SD_RMSE <- sd(agg$RMSE, na.rm = TRUE)
  N_Flocks_adapted <- agg$Flock[which(agg$RMSE < (Mean_RMSE + 2 * SD_RMSE) )[1]]

  Mean_RMSE_m <- mean(agg_month$RMSE, na.rm = TRUE)
  SD_RMSE_m <- sd(agg_month$RMSE, na.rm = TRUE)
  N_Months_adapted <- agg_month$MonthTotal[which(agg_month$RMSE < (Mean_RMSE_m + 2 * SD_RMSE_m))[1]]

  # plot and save per flock
  filename <- paste('plot_RMSE/RMSE_flock_', model, '_', method, '_', fold, '.pdf',
                    sep = "")
  pdf(filename)
  par(mfrow = c(1,1))
  plot(agg,
       type='b',
       main=paste('fold:', fold)
       )
  abline(v=N_Flocks_adapted, col='red', type=2)
  dev.off()

  # plot and save per month
  filename <- paste('plot_RMSE/RMSE_month_', model, '_', method, '_', fold, '.pdf',
                    sep = "")
  pdf(filename)
  par(mfrow = c(1,1))
  plot(agg_month,
       type='b',
       main=paste('fold:', fold)
       )
  abline(v=N_Months_adapted, col='red', type=2)
  dev.off()

  # also in plot view
  plot(agg,
       type='b',
       main=paste('fold:', fold)
       )
  abline(v=N_Flocks_adapted, col='red', type=2)

  plot(agg_month,
       type='b',
       main=paste('fold:', fold)
       )
  abline(v=N_Months_adapted, col='red', type=2)


  backtransformed_data$ft_back_all_bin <- 0
  backtransformed_data$ft_back_all_bin[which(backtransformed_data$ft_back_all >= 80)] <- 1

  backtransformed_data$Yt_back_all_bin <- 0
  backtransformed_data$Yt_back_all_bin[which(backtransformed_data$Yt_back_all >= 80)] <- 1

  Sensitivity <- length(
    which(
      backtransformed_data$Yt_back_all_bin == 1 & backtransformed_data$ft_back_all_bin == 1
      )
    ) / length(
      which(
        backtransformed_data$Yt_back_all_bin == 1
        )
      )
  Specificity <- length(
    which(
      backtransformed_data$Yt_back_all_bin == 0 & backtransformed_data$ft_back_all_bin == 0)
    )/length(which(backtransformed_data$Yt_back_all_bin == 0))

  # Add the results to the output
  out_all <- rbind(out_all,
                   cbind('TestSet'= fold,
                         'RmseTestSet'= RMSE_test_set,
                         'RmseTestSetBack'= RMSE_test_set_back,
                         'NFlocksAdapted'= N_Flocks_adapted,
                         'Sensitivity'=Sensitivity,
                         'Specificity'=Specificity))

  # Also make a complete dataset of all predictions
  et_complete <- rbind(et_complete,
                       cbind(e_back_test_set,
                             'TestSet' = fold,
                             'FPL' = Yt_back_all))

  # And a dataset of all agg RMSE per flock
  agg_all <- rbind(agg_all,
                   cbind(
                     'TestSet' = fold,
                     agg)
                   )
}


# Let's look at the results
View(out_all)

mean(out_all$RMSE_test_set_back)
mean(out_all$Sensitivity)
mean(out_all$Specificity)

filename <- paste("outputs/out_", model, '_', method, '.rds',
                  sep = "")
saveRDS(out_all, file = filename)


# Save the complete dataset
filename <- paste("outputs/et_complete_", model, '_', method, '.rds',
                  sep = "")
saveRDS(et_complete, file = filename)

# Save aggregate RMSE per flock
filename <- paste("outputs/agg_all_", model, '_', method, '.rds',
                  sep = "")
saveRDS(agg_all, file = filename)


# See how long it all took
Diff.time <- difftime(time1 = Sys.time(), time2 = start_time, units='min')
print(Diff.time)

filename <- paste("outputs/DiffTime_", model, '_', method, '.rds',
                  sep = "")
saveRDS(Diff.time, file = filename)

filename <- paste("outputs/res_out_all_", model, '_', method, '.rds',
                  sep = "")
saveRDS(object = res_out_all, file = filename)


# Plot RMSE per flock for all folds simultaneously

ggplot(agg_all,
       aes(Flock, RMSE)) +
  geom_line(aes(group = Test.set), size = 0.5) +
  geom_smooth(se = FALSE, size = 2) +
  labs(title = paste(model, method)) +
  theme_classic() +
  theme(text = element_text(size = 20),
        axis.title = element_text(size=22)
        )

