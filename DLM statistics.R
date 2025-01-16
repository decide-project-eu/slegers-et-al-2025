### Statistics ###
## Analysis of DLM results: compare different methods ##

library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)

## Note: this is without the uni_wave_discount

filenames_all <- c(
  'outputs/out_uni_nowave_EM.rds',
  'outputs/out_uni_wave_EM.rds',
  'outputs/out_multi_wave_EM.rds',
  'outputs/out_multi_nowave_EM.rds',
  'outputs/out_uni_nowave_discount.rds',
  'outputs/out_multi_wave_discount.rds',
  'outputs/out_multi_nowave_discount.rds'
  )

out <- data.frame()
for(name in filenames_all){
  
  # Make a combined dataframe of all results
  d <- readRDS(name)
  
  # Add a column with the complete name of the model
  d$ModelComplete <- sub('^outputs/out_', '', name)
  # Remove '.rds' at the end
  d$Model <- sub('\\.rds$', '', d$Model)
  
  d$Method <- ifelse(grepl('EM', name), 'EM',
                     'discount')
  d$Var <- ifelse(grepl('uni', name), 'univar',
                     'multivar')
  d$Wave <- ifelse(grepl('nowave', name), 'nowave',
                     'wave')
  d$Wave <- ifelse(grepl('nowave', name), 'nowave',
                     'wave')
  d$Model <- paste(d$Var, d$Wave, sep = "_")
  d$Model <- factor(d$Model,
                    levels = c("univar_nowave", "univar_wave",
                    "multivar_nowave", "multivar_wave"))
  out <- rbind(out, d)
  
}

# rename columns
names(out)[names(out) == 'RMSETestSetBack '] <- 'RMSE'
names(out)[names(out) == 'Sensitivity'] <- 'Se'
names(out)[names(out) == 'Specificity'] <- 'Sp'


# Add prevalence, accuracy and MMA

# Load a dataset with all farm houses
Data <- readRDS("outputs/res_out_all_multi_nowave_EM.rds")
Data$Target <- 0
Data$Target[which(Data$Yt_back_all > 80)] <- 1

# calculate prevalence for each fold
prevalence_data <- Data %>%
  dplyr::group_by(Fold) %>% 
  dplyr::summarize(Prevalence = mean(Target)) %>%
  dplyr::rename(TestSet = Fold)

# add prevalence data to output
out <- out %>% left_join(prevalence_data, by = "TestSet")
# calculate accuracy and MMA with prevalence, Se and Sp
out$Accuracy <- (out$Se * out$Prevalence) + (out$Sp * (1-out$Prevalence))
out$MMA <- (out$Se + out$Sp)/2

out$ModelComplete <- factor(out$ModelComplete)
out$ModelComplete <- relevel(out$ModelComplete, ref = "uni_nowave_EM.rds")


# mean MMA per model

MeanMMA <- out %>%
  group_by(ModelComplete) %>%
  summarise(MMA = mean(MMA))
MeanMMA


# Test: is there a difference between model performance?

rmse_lm <- lm(RMSE ~ ModelComplete,
                  data = out)
anova(rmse_lm)
summary(rmse_lm)

rmse_lmer <- lmer(RMSE ~ ModelComplete + (1|TestSet),
                  data = out)
anova(rmse_lmer)
summary(rmse_lmer)

library(emmeans)
library(multcomp)
library(multcompView)

ls = lsmeans(rmse_lmer, pairwise ~ ModelComplete)
S <- summary(ls)$contrasts
p <- S$p.value
p <- p[which(p < 0.05)]
range(p)
summary(ls)$contrasts


ls_cld <- multcomp::cld(object = ls,
                       adjust = "Tukey",
                       Letters = letters,
                       alpha = 0.05,
                       type = "response",
                       decreasing = TRUE)
ls_cld


se_lmer <- lmer(Se ~ ModelComplete + (1 | TestSet),
                  data = out)
anova(se_lmer)


sp_lmer <- lmer(Sp ~ ModelComplete + (1 | TestSet),
                  data = out)
anova(sp_lmer)


mma_lmer <- lmer(MMA ~ ModelComplete  + (1 | TestSet),
                  data = out)
anova(mma_lmer)


# PLOT TO COMPARE PERFORMANCE OF DLMS

# names for in the plot:
modelnames <- c(
  "univariate,\nno FPL wave", 
  "univariate,\nFPL wave",
  "multivariate,\nno FPL wave", 
  "multivariate,\nFPL wave"
  )

# RMSE (with statistical differences)
ggplot(data = out,
       aes(x = Model, y = RMSE, fill = Method)) +
  geom_boxplot(width = 0.9, linewidth = 0.75) +
  theme_classic() +
  theme(legend.position = c(0.8, 0.8),
        text=element_text(size=18),
        axis.line=element_line(linewidth=1),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_x_discrete(labels = modelnames) +
  scale_fill_manual(values = c("#2264A9", "#A6C1DE")) +
  ggsignif::geom_signif(
    y_position = c(150, 140, 140, 150, 160, 170), 
    xmin = c(0.8, 1.3, 2.1, 2.1, 2.1, 2.1), 
    xmax = c(1.9, 1.9, 2.7, 3.3, 3.7, 4.3), 
    vjust = 0.5,
    annotation= c("***", "***", "***", "***", "***", "**"),
    textsize = 6,
    tip_length = 0,
    size = 1
  ) +
  xlab("A")

# Se
ggplot(data = out,
       aes(x = Model, y = Se, fill = Method)) +
  geom_boxplot() +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size=18),
        axis.line=element_line(linewidth=1)
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_fill_manual(values = c("#2264A9", "#A6C1DE"))

# Sp
ggplot(data = out,
       aes(x = Model, y = Sp, fill = Method)) +
  geom_boxplot() +
  theme_classic() +
  theme(legend.position = "None",
        text=element_text(size=18),
        axis.line=element_line(size=1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_fill_manual(values = c("#2264A9", "#A6C1DE"))

# MMA (with statistical differences)
ggplot(data = out,
       aes(x = Model, y = MMA, fill = Method)) +
  geom_boxplot(width = 0.9, linewidth = 0.75) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.15),
        text=element_text(size=18),
        axis.line=element_line(linewidth=1),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_x_discrete(labels = modelnames) +
  scale_fill_manual(values = c("#2264A9", "#A6C1DE")) +
  xlab("B")


boxplot(out$RMSE ~ out$Model,
        outline = FALSE,
        las = 2) # without outliers!
boxplot(out$Se ~ out$Model,
        las = 2)
boxplot(out$Sp ~ out$Model,
        las = 2)


# plot sensitivity against RMSE
ggplot(out,
       aes(x = RMSE, y = Se, group = Model, color = Model)) +
  geom_point() +
  theme_classic()

# plot sensitivity against specificity
plot(out$Se ~ out$Sp)