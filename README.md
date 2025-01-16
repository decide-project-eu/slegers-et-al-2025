# fpl-prediction-models
Dynamic linear models (DLMs) and random forest models to predict footpad lesion scores of broiler flocks.

## Data

### Flock production and health data
This project uses a dataset that is also used in Slegers et al. (2024) ([link to article](https://doi.org/10.1016/j.psj.2024.104197)). The data processing steps can be found [in this repository](https://github.com/decide-project-eu/slegers-et-al-2024).

### Weather data and building direction
Weather variables were extracted from the [website of KNMI](https://www.knmi.nl/nederland-nu/klimatologie/uurgegevens) (Royal Dutch Meteorological Institute). [This Github repository](https://github.com/decide-project-eu/BuildingDetection) shows how the wind data can be extracted from the KNMI website using an API, and how the direction of the building can be determined using OpenStreetMap.

### Feed price
Feed prices were obtained from the [Agri&Food portal](https://agrimatie.nl/agrimatieprijzen/) of Wageningen University & Research (Dutch)

### Biosecurity score
The biosecurity questionnaire of AVINED which was used in this project can be found [on their website](https://www.avined.nl/themas/bedrijfsmanagement/bedrijfshygiene) (Dutch).

## R scripts
The order of the scripts in this project is:
1) [Data preparation](/scripts/Data_preparation.R) to select farm houses and make the dataset suitable for DLM
2) Data checks: [Normality of FPL data](scripts/FPL_check_gaussian.R), [optimal number of harmonics in FPL](scripts/FPL harmonics.R), [correlations between weather variables](scripts/weather_var_correlations.R), [Normality of weather variables](scripts/weather_var_check_gaussian.R), and the [optimal number of harmonics for weather variables](scripts/weather_var_harmonics.R)
3) [Run DLMs](scripts/Run_DLMs.R)
4) [DLM statistics](scripts/DLM_statistics.R) to analyse results
5) [RF and DLM outputs](scripts/RF_and_DLM_outputs.R) to run the random forest models and compare to best DLM

These scripts source the files [scripts/DLM_functions_FPL.R](scripts/DLM_functions_FPL.R) and [scripts/Custum_getGt_and_getFt_functions.R](scripts/Custum_getGt_and_getFt_functions.R)

The following R packages are used: tidyr, dplyr, ggplot2, corrplot, Hmisc (optional), ggsignif (optional), ranger, lme4, lmerTest, emmeans, multcomp, multcompView
