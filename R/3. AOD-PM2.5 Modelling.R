#####################################################################
#### This script combines the data from the previous two scripts and
#### fits a mixed effect model to make PM2.5 predictions. These 
#### predictions are made across Victoria from 11/11/2019 to
#### 25/01/2020 and then are seperated into a bushfire and a 
#### background component.
#####################################################################

set.seed(2020)

library(sf)
library(tidyverse)
library(lubridate)
library(raster)
library(dtplyr)
library(lme4)
library(lmerTest)
library(MuMIn)
library(spgwr)
library(GWmodel)
library(dtplyr)
library(gghighlight)
library(units)
## Load Data
## Create merged data set. 
#### Join AOD & copernicus meterological data by air quality station and id. Set zeroes in PM2.5 data as a very small number 
combined <- pm25_aod.matched %>% 
  mutate(aq_id = as.character(aq_id)) %>%
  left_join(cop_sf_1.sf.matched, by = c('aq_id','date')) %>%
  filter(layer != 0,PM25 != 0) %>%
  remove_missing()

# Model Fitting ------------

#### Fit a linear fixed effects model w/ day-specific slopes and intercepts

linear_me <- lmer(log(PM25) ~ (1+log(layer)|date) + t2m + d2m + layer * d2m + layer * t2m, data = combined)

# Model Prediction -------------
#### Load all AOD Data-
aod_all_areas <- aod_all_matched_imputed
metero_sat_dat_1 <- cop_sf_1.sf
#### Aggregate copernicus dataset 1 to daily averages.
metero_sat_dat_1.agg <- metero_sat_dat_1 %>%
  lazy_dt() %>%
  mutate(lat = unlist(map(metero_sat_dat_1$geometry,1)),
         long = unlist(map(metero_sat_dat_1$geometry,2))) %>%
  group_by(Date,lat,long) %>%
  summarise(u10 = mean(u10),
            d2m = mean(d2m),
            t2m = mean(t2m)) %>%
  as_tibble() %>%
  st_as_sf(coords = c('lat','long'),
           crs = crs(metero_sat_dat_1)) %>%
  st_transform(crs = crs(aod_all_areas))
#### Aggregate copernicus dataset 2 to daily averages.

#### create vector of unique dates
dates <- unique(cop_sf_1.sf$Date)
store_7A <- vector(mode = 'list',
                   length = length(dates))

for(i in seq_along(dates)){
  ### Find date i
  date = dates[[i]]
  ### Select AOD data for date i
  aod_all_areas.day <-aod_all_areas[aod_all_areas$date == date,]
  aod_all_areas.day.sp <- as(aod_all_areas.day,Class = 'Spatial')

  ### Select meterological data for date i
  metero_sat_dat_1.day <- metero_sat_dat_1.agg[metero_sat_dat_1.agg$Date == date,]
  ### Join imputed AOD data and meterological data
  joined <- aod_all_areas.day.sp %>%
    st_as_sf() %>%
    select(var1.pred) %>%
    st_transform(crs = crs(metero_sat_dat_1.day)) %>%
    st_join(metero_sat_dat_1.day,
            st_nearest_feature,
            left = TRUE)
  store_7A[[i]] <- joined
  print(i)
}
### Create merged dataset
aod_all_matched_stage2 <- do.call(rbind,store_7A)

# Predict New Data --------------
## Finalise prediction data
stage_2_data <- aod_all_matched_stage2 %>%
  rename(date = Date,layer = var1.pred) %>%  
  filter(layer != 0) %>%
  as_tibble() %>%
  filter(date %in% unique(combined$date)) %>%
  filter(date > as.Date('2019-11-10'))
## Make predictions
predictions <- predict(linear_me,newdata = stage_2_data)
## Add predictions to dataframe
stage_2_data$predictions_raw <- predictions

# Seperate into bushfire and background components --------------
## Calculate 99th percentile statewide thresholds
fire_thresholds <- aq_hour %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(date,aq_id) %>%
  summarise(PM25 = mean(pm25,na.rm=TRUE)) %>% 
  filter(date < as.Date('2019-02-01')) %>%
  as_tibble() %>%
  mutate(month = month(date)) %>%
  filter(month %in% c(11,12,1)) %>%
  group_by(month) %>%
  summarise(quantile_99 = quantile(log(PM25),probs = 0.99,na.rm=T))
## Calculate medians (net of )
fire_medians <-aq_hour %>%
  mutate(date = as.Date(datetime)) %>%
  filter(date < as.Date('2019-02-01')) %>%
  group_by(date,aq_id) %>%
  summarise(PM25 = mean(pm25,na.rm=TRUE)) %>%
  mutate(month = month(date),
         fire_driven = case_when(month == 1  & log(PM25) > fire_thresholds[fire_thresholds$month == 1,]$quantile_99 ~ 1,
                                 month == 11 & log(PM25) > fire_thresholds[fire_thresholds$month == 11,]$quantile_99 ~ 1,
                                 month == 12 & log(PM25) > fire_thresholds[fire_thresholds$month == 12,]$quantile_99 ~ 1,
                                 TRUE ~ 0),
         predictions_adj = ifelse(fire_driven == 1,NA,log(PM25))) %>%
  group_by(month) %>%
  summarise(median = median(log(PM25) ,na.rm=TRUE))

## Seperate estimated PM2.5 into a background and bushfire component
pm25_seperated.sf <- stage_2_data %>%
  mutate(month = month(date),
         fire_driven = case_when(month == 1    & predictions_raw > fire_thresholds[fire_thresholds$month == 1,]$quantile_99 ~ 1,
                                 month == 11   & predictions_raw > fire_thresholds[fire_thresholds$month == 11,]$quantile_99 ~ 1,
                                 month == 12   & predictions_raw > fire_thresholds[fire_thresholds$month == 12,]$quantile_99 ~ 1,
                                 month == 1    & predictions_raw > fire_thresholds[fire_thresholds$month == 1,]$quantile_99 ~ 1,
                                 month == 11   & predictions_raw > fire_thresholds[fire_thresholds$month == 11,]$quantile_99 ~ 1,
                                 month == 12   & predictions_raw > fire_thresholds[fire_thresholds$month == 12,]$quantile_99 ~ 1,
                                 TRUE ~ 0),

         fire_pm25 = case_when(month == 1    & fire_driven == 1 ~ predictions_raw - fire_medians[fire_medians$month == 1 ,]$median,
                               month == 11   & fire_driven == 1 ~ predictions_raw - fire_medians[fire_medians$month == 11,]$median,
                               month == 12   & fire_driven == 1 ~ predictions_raw - fire_medians[fire_medians$month == 12,]$median,
                               month == 1    & fire_driven == 1 ~ predictions_raw - fire_medians[fire_medians$month == 1 ,]$median,
                               month == 11   & fire_driven == 1 ~ predictions_raw - fire_medians[fire_medians$month == 11,]$median,
                               month == 12   & fire_driven == 1 ~ predictions_raw - fire_medians[fire_medians$month == 12,]$median,
                               TRUE ~ 0),
         background_pm25 = case_when(month == 1    & fire_driven == 1 ~ fire_medians[fire_medians$month == 1 ,]$median,
                                     month == 11   & fire_driven == 1 ~ fire_medians[fire_medians$month == 11,]$median,
                                     month == 12   & fire_driven == 1 ~ fire_medians[fire_medians$month == 12,]$median,
                                     month == 1    & fire_driven == 1 ~ fire_medians[fire_medians$month == 1 ,]$median,
                                     month == 11   & fire_driven == 1 ~ fire_medians[fire_medians$month == 11,]$median,
                                     month == 12   & fire_driven == 1 ~ fire_medians[fire_medians$month == 12,]$median,
                                     TRUE ~ predictions_raw)) %>%
  st_as_sf()

### Cross Validation
store_cv <- vector(mode = 'list', length = length(unique(combined$aq_id)))
for(i in 1:length(unique(combined$aq_id))){
  test <- combined %>%
    filter(aq_id == unique(combined$aq_id)[i])
  train <- combined %>%
    filter(aq_id != unique(combined$aq_id)[i])
  test <- test[test$date %in% train$date,]
  linear_me_cv <- lmer(log(PM25) ~ (1+log(layer)|date) + t2m + d2m + layer * d2m + layer * t2m, data = train)
  pred <- (predict(linear_me_cv,newdata = test) %>% as_tibble())
  joined <- cbind(test$aq_id,as.Date(test$date),log(test$PM25),pred) 
  store_cv[[i]] <- joined
  print(i)
}

cv <- do.call(rbind,store_cv) 
colnames(cv) <- c('aq_id','date','actual','predict')



