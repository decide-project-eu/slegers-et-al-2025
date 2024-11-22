
start.time <- Sys.time()

library(randomForest)
library(pROC)
library(splitstackshape)
library(ggplot2)
library(dplyr)
library(ranger)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(multcompView)


#### A function for seeing the seperation of the classes achieved by the trained model
see.seperation <- function(obs.n.pred, Main=NA){
  
  # Plot a histogram showing the distribution of predictions given the true state
  # plot bar chart for age with Q1-Q3 intervals around the median
  summaries <- as.data.frame(as.data.frame(aggregate(obs.n.pred$Pred, by=list(obs.n.pred$Obs), FUN=summary)))
  labels <- summaries$Group.1
  values <- as.data.frame(summaries$x)
  Medians <- values$Median
  lower <- values$`1st Qu.`
  upper <- values$`3rd Qu.`
  
  Medians <- rbind(Medians)
  colnames(Medians) <- labels
  BP <- barplot(Medians, xaxt = "n", beside = TRUE,  ylim=c(0,1), col=c( 'grey80', 'grey30'), space=c(0.2,0), las=2, main=Main)
  # - ad rotated labels
  arrows(x0=BP, y0=lower, x1 = BP, y1 = upper, code = 3, angle = 90, length = 0.10)
  text(x = 1:length(Medians)-0.5,
       y = -0.05,
       labels = colnames(Medians),
       xpd = NA,
       ## Rotate the labels by 35 degrees.
       srt = 0,
       adj = 1)
  
  # Get the range of the medians as a proxy for how well the groups are seperated
  
  
  return(diff(range(Medians)))
  
}

# Function to create a new vector with positions of zeros within their runs
get_zero_positions <- function(vector) {
  # Identify runs of consecutive zeros
  runs <- rle(vector == 0)
  
  # Initialize a vector to store the positions
  position_vector <- rep(0, length(vector))
  
  # Assign positions within runs to the output vector
  current_position <- 1
  for (i in seq_along(runs$values)) {
    if (runs$values[i]) {
      position_vector[current_position:(current_position + runs$lengths[i] - 1)] <- 1:runs$lengths[i]
      current_position <- current_position + runs$lengths[i]
    } else {
      current_position <- current_position + runs$lengths[i]
    }
  }
  
  return(position_vector)
}

# A function for applying the ROC curve to the data and extracting the relevant results
get.ROC.results <- function(observation_test, prediction_test, plot.it = TRUE, main=NA){
  
  # Get the ROC curve of RF applied to the test set, and extract the auc
  roc.out.test <- roc(response = observation_test, predictor = prediction_test)
  auc.test <- round( roc.out.test$auc , 2)
  
  # Get the three vectors "thresholds", "sensitivities", and "specificities"
  thresholds <- round( roc.out.test$thresholds , 3)
  sensitivities <- round( roc.out.test$sensitivities , 2)
  specificities <- round( roc.out.test$specificities , 2)
  
  # Calculate the MMA for all thresholds
  MMA.all <- (sensitivities+specificities)/2
  
  # Find the maximum MMA value
  MMA.max <- max(MMA.all)
  best.i <- which(MMA.all == MMA.max)

  #Find the threshold which results in the highest MMA
  best.threshold <- thresholds[best.i]
  
  # # Identify the threshold(s) which yields a sensitivity as close to and above 0.8 as possible
  # a <- sensitivities - 0.8
  # a <- a[which(a >= 0)]
  # best.i <- which(a == min(a[which(a >= 0)]))
  # best.threshold <- thresholds[best.i]

  # Find the sensitivity and specificity associated with the maximum MMA
  Se.best <- sensitivities[best.i]
  Sp.best <- specificities[best.i]

  # Plot the ROC curve
  if(plot.it == TRUE){
    plot(roc.out.test, main=main)
    abline(h=Se.best, lty=2)
    abline(v=Sp.best, lty=2)
    text(x = 0.2, y=0.2, paste("AUC:", auc.test))
  }
  
  return(as.data.frame(cbind('AUC'=auc.test,
               'MMA.max'=MMA.max,
               'Threshold'=best.threshold,
               'Se'=Se.best,
               'Sp'=Sp.best)))
  
}



### UNI NOWAVE Discount DLM
# -----------------------------------------------------------------------------

#### Get data, if not already loaded ####
Data <- readRDS("outputs/res_out_all_uni_nowave_discount.rds")
Data <- na.omit(Data)

Data$FPLLag1 <- lag(Data$Yt_back_all)
Data$MtFPLLag1 <- lag(Data$MtFootpadLesions)
Data$UtFPLLag1 <- lag(Data$UtFt.j)

Data$Target <- 0
Data$Target[which(Data$Yt_back_all > 80)] <- 1
Data$Target <- as.factor(Data$Target)

Data <- na.omit(Data)

nrow(Data)

meta_names <- c("Type2", "Month",
                "AntibioticsWeek1",
                "MeanAgeAtSlaughter", "EndFlockSize", "FlockSizeDiff",
                "DaysBetweenFlocks", "FractionThinned",
                "Province", "Score", "BuildingAngle", "NumberOfHouses",
                "MeanNDVTiterLag1", "AntibioticsLag1", "MortalityLag1",                 
                "WindDir", "AngleDiffWk1", "FeedPrice", "Temp05", 
                "Humidity01", "WeeklyTotalRain", "WindSpeed09")

prev_obs_names <- c("FPLLag1")
raw_names <- c("FPLLag1")
dlm_names <- c('MtFPLLag1', 'UtFPLLag1') #standardized mean and forecast error


#### Repeat the procedure 20 times ####

# Set the seed for reproducability
set.seed(seed = 1225)
SEEDS <- sample(1:99999999, size = 50, replace = FALSE)
out_all <- data.frame()


for(fold in unique(Data$Fold)){

  #fold <- 1
  gc()
  
  test_set <- subset(Data, Data$Fold == fold)
  training_set <- subset(Data, Data$Fold != fold)

  observation_test <- test_set$Target

  # --------------------------------------------------------------------------------------

  #### Test each names collection on their own (in combination with the meta data names) ####

  name_collections <- c('prev_obs_names', 'raw_names', 'dlm_names')
  for(name_collection in name_collections){

    # Define relevant_names
    if(name_collection == 'prev_obs_names'){
       relevant_names <- c(get(name_collection)) # without metadata
    }else{
       relevant_names <- c(meta_names, get(name_collection)) # "dim", ut.names, mt.names, raw_names
    }
   
    # Make relevant_names into a string with plus signs between the original elements
    relevant_names <- paste(relevant_names, collapse = '+')

    # Make the formula
    f.ut <- as.formula(paste('Target ~ ' , relevant_names) )

    #Print it to console
    print(f.ut)


    start.time.RF <- Sys.time()
    # Train random forest model
    set.seed(42)
    
    RF_ranger <- ranger(formula = f.ut,
                       data = training_set,
                       probability = TRUE,
                       importance= 'impurity',
                       num.trees=500
    )
    
    print(Sys.time() - start.time.RF)

    # Predict on the test set
    prediction_test <- predict(object = RF_ranger, data= test_set)$predictions
    
    # Get the column with the positive probability for FPL score >80
    prediction_test <- prediction_test[ , '1']

    # Make the ROC curve and get relevant results from it
    out <- get.ROC.results(observation_test, prediction_test, plot.it = TRUE)

    # Add to out_all
    out <- cbind(
      'Fold' = fold, 
      'Model' = name_collection,
      out
    )
    # - we only care about the AUC right now, so the first row is as good as any
    out <- out[1,]
    out_all <- rbind(out_all, out)
    print(out)

  }
}

print(Sys.time()-start.time)

View(out_all)

# Calculate accuracy for each fold, using fold Prevalence
prev_data <- Data %>% 
  dplyr::group_by(Fold) %>% 
  dplyr::summarize(Prevalence = mean(as.numeric(Target)) -1 )
    # -1 because converting factor to numeric becomes 1 and 2's

prev_data$Fold <- prev_data$Fold + 1 #because here it starts with 1 instead of 0

out_all <- out_all %>% 
  left_join(prev_data, by = "Fold")

out_all$Accuracy <- (out_all$Se * out_all$Prevalence) + (out_all$Sp * (1-out_all$Prevalence))
  
out_all$MMA <- (out_all$Se + out_all$Sp) / 2
out_all$PPV <- (out_all$Se * out_all$Prevalence) / 
               (
                 (out_all$Se * out_all$Prevalence) + 
                 ((1 - out_all$Sp) * (1 - out_all$Prevalence))
               )
out_all$NPV <- (out_all$Sp * (1 - out_all$Prevalence)) / 
               (
                 (out_all$Sp * (1 - out_all$Prevalence)) + 
                 ((1 - out_all$Se ) * out_all$Prevalence)
               )

saveRDS(object = out_all, file = 'RF_out_all_uni_nowave_discount.rds')

# out_all <- readRDS('RF_out_all_uni_nowave_discount.rds')



#### Compare the different models  ####

LM <- lmer(AUC ~ Model + (1|Fold), data=out_all)
anova(LM)
summary(LM)

LM_mma <- lm(MMA ~ Model + (1|Fold), data=out_all)
anova(LM_mma)
summary(LM_mma)


# Calculate mean performance per model
out_all %>%
  dplyr::group_by(Model) %>%
  summarise(auc = mean(AUC),
            mma = mean(MMA),
            accuracy = mean(Accuracy),
            se = mean(Se),
            sp = mean(Sp))


# get (adjusted) weight means per group
# NB this is equal to the mean
model_means <- emmeans(object = LM,
                       specs = "Model")
model_means

comparisons <- emmeans(LM, pairwise ~ Model, infer = TRUE)
summary(comparisons)


compared <- as.data.frame(comparisons$emmeans)
compared$RF_name <- c("Previous FPL \n score only", "Raw variables", "DLM output")
compared$RF_name <- factor(compared$RF_name, levels = CLD$RF_name)

compared %>%
  ggplot(aes(x = factor(RF_name), y = emmean)) +
  geom_bar(stat = "identity", fill = "#A6C1DE") +  
  geom_errorbar(
    aes(ymin=lower.CL, ymax=upper.CL),
    size = 1.6,
    width = 0.5
  ) +
  labs(y="AUC",
       x = "") +
  theme_classic() +
  theme(
    text=element_text(size=25),
    axis.line=element_line(size=2)
  ) +
  ggsignif::geom_signif(
    y_position = c(0.76, 0.85), xmin = c(1, 1), xmax = c(2, 3), vjust = 0.5,
    annotation= c("***", "***"),
    textsize = 14,
    tip_length = 0,
    size = 2
  ) +
  scale_y_continuous(lim = c(0, 0.95))



# Train a final model based on all data

# - Define relevant_names
relevant_names <- c(meta_names, raw_names) # "dim", ut.names, mt.names, raw_names
relevant_names_basic <- c(prev_obs_names)

relevant_names <- paste(relevant_names, collapse = '+')

# - Make the formula
f.Utbasic <- as.formula(paste('Target ~ ' , relevant_names_basic) )
f.ut <- as.formula(paste('factor(Target) ~ ' , relevant_names) )

# Train random forest model
RF_ranger_basic <-  ranger(formula = f.Utbasic, 
                   data = Data, 
                   importance= 'impurity',
                   probability = TRUE,
                   num.trees=500
)

RF_ranger <-  ranger(formula = f.ut,
                   data = Data, 
                   importance= 'impurity',
                   probability = TRUE,
                   num.trees=500
)

# TRY OUT partial dependence
library(pdp)
RF_ranger_basic %>%
  partial(pred.var = "FPLLag1") %>%
  autoplot(rug = TRUE, train = Data) + theme_classic()

RF_ranger %>%
  partial(pred.var = "FPLLag1") %>%
  autoplot(rug = TRUE, train = Data) + theme_classic()


# See the importance
importance(RF_ranger)

variable.importance <- sort(importance(RF_ranger), decreasing = TRUE)

ggplot2::ggplot(
    tibble::enframe(
        variable.importance,
        name = "variable",
        value = "importance"
    ),
    ggplot2::aes(
        x = reorder(variable, importance),
        y = importance,
        fill = importance
    )
) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::coord_flip() +
    ggplot2::theme(
      text = element_text(size = 15)
    ) +
    ylab("Variable Importance") +
    xlab("") +
    # ggtitle("Information Value Summary") +
    guides(fill = "none") +
    scale_fill_gradient(low = "red", high = "blue") +
  theme_classic()


# Also for DLM
# Train a final model based on all data, except the final data set
# - Define relevant_names
relevant_names <- c(meta_names, raw_names) # "dim", ut.names, mt.names, raw_names

# - Make relevant_names into a string with plus signs between the original elements
relevant_names <- paste(relevant_names, collapse = '+')

# - Make the formula
f.ut <- as.formula(paste('Target ~ ' , relevant_names) )

# Train random forest model
RF_ranger <-  ranger(formula = f.ut, 
                   data = Data, 
                   importance= 'impurity',
                   probability = TRUE,
                   num.trees=500
)

# See the importance
importance(RF_ranger)   # using only ranger package

variable.importance <- sort(importance(RF_ranger), decreasing = TRUE)

ggplot2::ggplot(
    tibble::enframe(
        variable.importance,
        name = "variable",
        value = "importance"
    ),
    ggplot2::aes(
        x = reorder(variable, importance),
        y = importance,
        fill = importance
    )
) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::coord_flip() +
    ggplot2::theme(
      text = element_text(size = 15)
    ) +
    ylab("Variable Importance") +
    xlab("") +
    # ggtitle("Information Value Summary") +
    guides(fill = "none") +
    scale_fill_gradient(low = "red", high = "blue") +
  theme_classic()



######################################################

### Compare MMA to original DLM

out_all <- readRDS('RF_out_all_uni_nowave_discount.RDS')
out_all <- out_all %>%
  dplyr::select(Fold, name_collection, MMA) %>%
  rename(Model = name_collection)

out_dlm <- readRDS('outputs/res_out_all_uni_nowave_discount.rds')
  
out_dlm <- out_dlm %>%
  mutate(MMA = (Sensitivity + Specificity) / 2,
         Model = "DLM"
         ) %>%
  rename(Fold = Test.set) %>%
  dplyr::select(Fold, Model, MMA)

out.complete <- rbind(out_all, out_dlm)


LM <- lmer(MMA ~ Model + (1|Fold), data = out.complete)
anova(LM)

# LM <- lm(MMA ~ Model, data=out.complete)
# anova(LM)
# summary(LM)


model_means <- emmeans(object = LM,
                       specs = "Model")
model_means
comparisons <- emmeans(LM, pairwise ~ Model, infer = TRUE)
summary(comparisons)


compared <- as.data.frame(comparisons$emmeans)
compared$RF_name <- c("DLM", "RF DLM output", 
                      "RF Previous FPL \n score only", "RF Raw variables")
compared$RF_name <- factor(compared$RF_name, levels = c("DLM",
                      "RF Previous FPL \n score only",  "RF Raw variables",
                      "RF DLM output"))

compared %>%
  ggplot(aes(x = factor(RF_name), y = emmean)) +
  geom_bar(stat = "identity", fill = "#A6C1DE") +  
  geom_errorbar(
    aes(ymin=lower.CL, ymax=upper.CL),
    size = 1.6,
    width = 0.5
  ) +
  labs(y="MMA",
       x = "") +
  theme_classic() +
  theme(
    text=element_text(size=22),
    axis.line=element_line(size=2)
  ) +
  ggsignif::geom_signif(
    y_position = c(0.68, 0.73, 0.78, 0.83, 0.88),
    xmin = c(1, 1, 1, 2, 2), 
    xmax = c(2, 3, 4, 3, 4), 
    vjust = 0.5,
    annotation= c("***", "***", "***", "***", "***"),
    textsize = 10,
    tip_length = 0,
    size = 1.5
  ) +
  scale_y_continuous(lim = c(0, 0.95))


