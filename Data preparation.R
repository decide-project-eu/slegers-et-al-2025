library(tidyr)
library(dplyr)


### READ IN DATA ###

data_all_with0 <- read.csv(
  "C:/dfFPL_withSlow.csv"
  )


### DATA TWEAKING ###

# change 'null' to NA
data_all_with0[data_all_with0 == "null"] <- NA

# DON'T filter on missing FPL
# DON'T filter 0 scores. Instead, replace them with NA
data_all_with0 <- data_all_with0 %>%
  mutate(FootpadLesions = as.numeric(FootpadLesions))

# create dataset without 0 scores
data_all <- data_all_with0
data_all$FootpadLesions <- ifelse(data_all$FootpadLesions == 0, NA,
                                  data_all$FootpadLesions)

# create variable MonthTotal with number of months starting 2017-01-01
data_all$MonthTotal <- (data_all$HatchYear - 2017) * 12 + data_all$HatchMonth

saveRDS(data_all, file = 'data_all.rds')
cols <- colnames(data_all)

nrow(data_all)
nrow(data_all %>% filter(!is.na(FootpadLesions)))
length(unique(data_all$FarmIdentification))


### SELECTING HOUSES WITH >30 SCORES ###

# list of houses with >30 scores

split_data <- split(data_all, data_all$Farmhouse)
AllFarmhouses <- names(split_data)

house_list <- c()

for(house in AllFarmhouses){
     data_house <- split_data[[house]]
     if (length(data_house$FootpadLesions[!is.na(data_house$FootpadLesions)]) >= 30){
       house_list <- c(house_list, house)
       print("saved")
     } else {
       print("less than 30 rows, skipping.\n")
     }
}

saveRDS(house_list, file = 'house_list.rda')

house_list <- readRDS('house_list.rda')

# Make dataset with selected farms
## AND SHIFT WEATHER DATA ###

data_selected_notlagged <- data_all %>%
  # Give weather variables easier names
  rename_with(~ sub("^T_", "Temp", .x), starts_with("T_")) %>%
  rename_with(~ sub("^RH_", "Rain", .x), starts_with("RH_")) %>%
  rename_with(~ sub("^FH_", "WindSpeed", .x), starts_with("FH_")) %>%
  rename_with(~ sub("^U_", "Humidity", .x), starts_with("U_")) %>%
  
  # Get selected farm houses
  filter(
    Farmhouse %in% house_list
  ) %>%
  mutate(
    Month = as.numeric(format(as.Date(HatchDate), "%m")),
  )

data_selected <- data_selected_notlagged  %>%
  mutate(
    # make numeric columns numeric
    across(
      c(where(is.character), -c(HatchDate, SlaughterDate, Type, Breed, Province)),
      as.numeric
    ),
    # add type 2 and Month column
    Type2 = ifelse(Type == "conventional", "conv", "alt"),
    Switched = ifelse(
      (SwitchedToCONV == 1 | SwitchedToSlow == 1),
      1,
      0
    )
  ) %>%
  group_by(Farmhouse) %>%
  mutate(
    # Get weather data of the next flock instead of current
    Temp05 = dplyr::lead(Temp05, order_by = MonthTotal),
    Humidity01 = dplyr::lead(Humidity01, order_by = MonthTotal),
    WeeklyTotalRain = dplyr::lead(WeeklyTotalRain, order_by = MonthTotal),
    WindSpeed09 = dplyr::lead(WindSpeed09, order_by = MonthTotal),
    Switched = dplyr::lead(Switched, order_by = MonthTotal),
    
    # Start flock counting at 1 for each farmhouse
    Flock = Flock - (min(Flock) - 1)
  ) %>%
  ungroup() 


### SOME CHECKS ###

# Check for duplicate flocks (same farm house and flock number)
duplicate_flocks <- data_selected %>%
  group_by(Farmhouse, Flock) %>%
  summarise(n = n()) %>%
  filter(n > 1) %>%
  ungroup()
nrow(duplicate_flocks)

nrow(data_selected)
length(unique(data_selected$Farmhouse))
table(data_selected$Flock)
# now there are as many 'flock 1' as there are farm houses

# all flock numbers per farm:
flocknums <- as.data.frame(table(data_selected$Farmhouse, data_selected$Flock))
flocknums <- spread(flocknums, key = Var2, value = Freq)
View(flocknums)

table(data_selected$Switched)


saveRDS(data_selected, file = 'data_selected.rds')
saveRDS(data_selected_notlagged, file = 'data_selected_notlagged.rds')


### Create dataset with rows for each Month ###

# ! this gives leading and trailing NULL flocks starting Jan 2017

data_selected <- readRDS('data_selected.RDS')

house_list <- readRDS("house_list.rda")
all_months <- expand.grid(Farmhouse = house_list,
                          HatchYear = c(2017:2021),
                          Month = c(1:12))

farm_data <- data_selected %>%
  select(FarmIdentification, Farmhouse, StationId) %>%
  distinct()


# 3) Merge with the original data to fill in missing values
data_month <- merge(all_months,
                     data_selected,
                     by = c("Farmhouse", "HatchYear", "Month"),
                     all.x = TRUE)

# Fill in empty farm and weather station rows (known from farmhouse)
# Fill in empty month (recalculate)
data_month <- merge(data_month,
                       farm_data,
                       by = c("Farmhouse"),
                       all.x = TRUE) %>%
              select(
                -FarmIdentification.x,
                -StationId.x
              ) %>%
              rename(
                FarmIdentification = FarmIdentification.y,
                StationId = StationId.y
              ) %>%
              mutate(
                MonthTotal = (HatchYear - 2017) * 12 + Month
              ) %>%
              select(
                FarmIdentification,
                Farmhouse,
                HatchYear,
                Month,
                MonthTotal,
                Flock,
                FootpadLesions,
                Switched,
                everything()
              )

View(data_month)
saveRDS(data_month, file = "data_month.rds")                


# How many houses start at month 1, 2, 3...?
Flock1 <- data_month %>%
  filter(Flock == 1)

table(Flock1$MonthTotal)




#####################################

### DESCRIPTIVE STATISTICS OF SELECTED FARMS ###

# if necessary, load data
# house_list <- readRDS("house_list.rda")
# data_selected <- readRDS("data_selected.rds")

# number of houses
length(house_list)

# number of farms
length(unique(data_selected$FarmIdentification))

# number of flocks
nrow(data_selected)

#distribution of FPL scores
hist(data_selected[(data_selected$Type == "conventional"),]$FootpadLesions,
     breaks = 30,
     main = "FPL scores of selected farms",
     xlab = "FPL score")

mean(data_selected$FootpadLesions, na.rm = TRUE)
sd(data_selected$FootpadLesions, na.rm = TRUE)
quantile(data_selected$FootpadLesions,
         probs = c(0,0.25,0.5,0.75,1),
         na.rm = TRUE)

data_selected$Target <- ifelse(data_selected$FootpadLesions <= 80, 0, 1)
prop.table(table(data_selected$Target))
table(data_selected$Target)

# number of alternative flocks
alt <- data_selected %>% filter((Type == "medium")|(Type == "slow-growing"))
n_alt <- nrow(alt)
n <- nrow(data_selected)
n_alt
n_alt/n

# FPL alt flocks
hist(alt$FootpadLesions,
     breaks = 30,
     main = "FPL scores of selected farms",
     xlab = "FPL score")


# Describing broiler types per farm house:
ntypes <- data_selected %>%
  group_by(Farmhouse) %>%
  summarize(n_types = n_distinct(Type2),
            all_types = list(Type2),
            alt_after_conv = any(Type2 == "alt" & lag(Type2) == "conv", 
                                 na.rm = TRUE),
            conv_after_alt = any(Type2 == "conv" & lag(Type2) == "alt", 
                                 na.rm = TRUE),
            only_conv = !any(Type2 == "alt", na.rm = TRUE),
            only_alt = !any(Type2 == "conv", na.rm = TRUE)
            ) %>%
  ungroup()

nrow(ntypes)
onetype <- ntypes %>% filter(n_types == 1)

# number of houses with conv and alt flocks
nrow(ntypes %>% filter(n_types == 2))

# number of houses going from conv to alt
nrow(ntypes %>% filter(alt_after_conv == TRUE))

# number of houses going from alt to conv
nrow(ntypes %>% filter(conv_after_alt == TRUE))

nrow(ntypes %>% filter(only_conv == TRUE))
nrow(ntypes %>% filter(only_alt == TRUE))


# Missing data
NA_count1 <- data_selected %>% summarise(across(everything(), ~ sum(is.na(.))))
NA_count <- NA_count1[-1] %>% 
  t() %>% as.data.frame() %>% setNames(NA_count1[,1]) %>%
  rename(NAs = '0') %>%
  filter(NAs != 0)
View(NA_count)

# Per farmhouse number of houses and missing data
NA_perhouse <- data_selected %>% 
  group_by(Farmhouse) %>% 
  summarise(Total = n(),
            NA_FPL = sum(is.na(FootpadLesions)),
            NA_AngleDiff = sum(is.na(AngleDiff_w1)),
            NA_Rain = sum(is.na(WeeklyTotalRain))
  )
View(NA_perhouse)

# other variables
summary(data_selected$DaysBetweenFlocks)
summary(data_selected$EndFlockSize)
sd(data_selected$EndFlockSize)
summary(data_selected$MeanAgeAtSlaughter)
sd(data_selected$MeanAgeAtSlaughter)


###########################################
# Compare to complete data

mean(data_all$FootpadLesions)
sd(data_all$FootpadLesions)
quantile(data_all$FootpadLesions, probs = c(0,0.25,0.5,0.75,1))

# number of farms
length(unique(data_all$FarmIdentification))

# number of alternative flocks in total
alt <- data_all %>% filter((Type == "medium")|(Type == "slow-growing"))
n_alt <- nrow(alt)
n <- nrow(data_all)
n_alt
n_alt/n

data_all$Type2 <- ifelse(data_all$Type %in% c('medium', 'slow-growing'),
                              'alt',
                              'conv')

ntypes <- data_all %>%
  group_by(Farmhouse) %>%
  summarize(n_types = n_distinct(Type2),
            all_types = list(Type2),
            alt_after_conv = any(Type2 == "alt" & lag(Type2) == "conv", 
                                 na.rm = TRUE),
            conv_after_alt = any(Type2 == "conv" & lag(Type2) == "alt", 
                                 na.rm = TRUE)
            ) %>%
  ungroup()

# number of houses with only conv flocks
conv_houses <- ntypes[grepl("conv", ntypes$all_types), ] %>%
  filter(n_types == 1)
nrow(conv_houses)

# number of houses with conv and alt flocks
nrow(ntypes %>% filter(n_types == 2))

# number of houses going from conv to alt
nrow(ntypes %>% filter(alt_after_conv == TRUE))

# number of houses going from alt to conv
nrow(ntypes %>% filter(conv_after_alt == TRUE))

# length of production cycle and downtime
hist(data_all$MeanAgeAtSlaughter)
summary(data_all$MeanAgeAtSlaughter)
summary(as.numeric(data_all$DaysBetweenFlocks))

hist(data_selected$MeanAgeAtSlaughter)
summary(data_selected$MeanAgeAtSlaughter)
summary(as.numeric(data_selected$DaysBetweenFlocks))

# number of houses per farm
houses_per_farm <- data_all %>%
  group_by(FarmIdentification) %>%
  summarize(n_houses = mean(as.numeric(NumberOfHouses))) %>%
  ungroup()
# note: this is average number of houses through all years, so can be a decimal number

summary(houses_per_farm$n_houses)
hist(houses_per_farm$n_houses)
