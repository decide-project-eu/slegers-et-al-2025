library(dplyr)
library(corrplot)
library(ggplot2)

df <- readRDS("Data_selected_notlagged.rds")

#Create variables
df$Month <- as.numeric(format(as.Date(df$HatchDate), 
                                             "%m"))

# Weather variables of all weeks of production cycle (week 1-4 and slaughter)
cor_data <- df %>%
  dplyr::select(86:230) %>%
  mutate(
    across(c(where(is.character)), as.numeric)
  )
  

temp_data <- cor_data[ , grepl( "Temp05", names(cor_data) ) ] %>% na.omit()
hum_data <- cor_data[ , grepl( "Humidity01", names(cor_data) ) ]  %>% na.omit()
windspeed_data <- cor_data[ , grepl( "WindSpeed09", names(cor_data) ) ]  %>% na.omit()
rain_data <- cor_data[ , grepl( "WeeklyTotalRain", names(cor_data) ) ]  %>% na.omit()
week1_data <- cor_data[ , !grepl( "_w|_S|WindDir", names(cor_data) ) ]  %>% na.omit()

par(mfrow = c(1,1))

# Correlations between variables FIRST WEEK
matrix_w1 <- cor(week1_data)
corrplot(matrix_w1, 
         type = "lower", 
         method = "number",
         number.cex=0.6)

# Correlations between different weeks
mat <- cor(temp_data)
corrplot(mat, 
         type = "lower", 
         method = "number",
         number.cex=1)

mat <- cor(hum_data)
corrplot(mat, 
         type = "lower", 
         method = "number",
         number.cex=1)

mat <- cor(rain_data)
corrplot(mat, 
         type = "lower", 
         method = "number",
         number.cex=1)

mat <- cor(windspeed_data)
corrplot(mat, 
         type = "lower", 
         method = "number",
         number.cex=1)


# CHECK NORMALITY OF DIFFERENT DECILES (Week 1)

# Humidity_01 the best?
par(mfrow=c(3,3))
for(i in 1:9){
  var <- paste('Humidity0', i, sep = "")
  hist(
    cor_data[, var], 
    main = paste(c('Humidity0', i), sep = "")
  )
}


### Correlation coefficients in a list

# Functions for correlation matrix
library(Hmisc)

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}


# Correlations between all weather variables
res2 <- rcorr(as.matrix(cor_data))
cor_matrix <- flattenCorrMatrix(res2$r, res2$P)
View(cor_matrix)
