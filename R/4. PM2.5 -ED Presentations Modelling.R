#####################################################################
#### This script creates the hospital catchment areas, matches the
#### dependent variablse for our hospital data and fits our main model.
#### It also runs the simulations for our policy analysis.
#####################################################################

set.seed(2020)

library(absmaps)
library(sf)
library(tidyverse)
library(fable)
library(dtplyr)
library(mgcv)
library(MASS)
library(tsibble)
library(gam)
library(data.table)
library(bbplot)
library(gam)
library(gghighlight)
library(ggsci)
library(cowplot)
select <- dplyr::select
### Matching estimated PM2.5 to hospitals

##### Step #1 - Create catchment areas

stage_2_data.sf <- pm25_seperated.sf

##### Hospital Processing -------
# Load Hospital Location Data
projcrs <- 4283
hospital_sf <- read_csv("Hospital Data here") %>%
  st_as_sf(coords = c('lon','lat'),
           crs = projcrs)
# Download SA3 maps from ABS
#state     <- load_absmaps(area = "state", year = 2016,saveDirectory = "~/Documents/ETC4860/Data/SA3 Shapefiles",download = T)
vicmap <- state %>% filter(state_name_2016 == "Victoria")
# Read hospital data
hospital_presentations <- read_csv("~/Documents/ETC4860/Data/Hospital Data/vic_broad_replicate.csv")
# Select respiratory data by day,hospital and age band
hosp_1 <- hospital_sf %>%
  lazy_dt() %>%
  left_join(hospital_presentations,by=c("hospital")) %>%
  filter(diag == 'Respiratory') %>%
  select(date,hospital,age,n) %>%
  as_tibble()

# Fill in missing dates
hospital_spatiotemporal <- hosp_1 %>%
  as_tsibble(key = c(hospital,age),index = date,validate=FALSE) %>%
  tsibble::fill_gaps(n = 0L) %>% 
  as_tibble() %>%
  left_join(hospital_sf, by = c('hospital')) %>%
  filter(date %in% stage_2_data.sf$date) %>%
  select(date,hospital,geometry,age,n)

### Find all the dates
dates <- unique(hospital_spatiotemporal$date)

### List to store for list output
store_hospital <- vector(mode = 'list',length = length(dates))
for(i in 1:length(dates)){
  # Total Data set
  hospital.filtered <-  hospital_spatiotemporal %>%
    filter(date == dates[[i]]) %>%
    st_as_sf()
  # 0 to 17
  hospital.filtered.child <-  hospital_spatiotemporal %>%
    filter(date == dates[[i]],
           age == "0 to 17") %>%
    st_as_sf()
  # 18 to 64
  hospital.filtered.adult <-  hospital_spatiotemporal %>%
    filter(date == dates[[i]],
           age == "18 to 64") %>%
    st_as_sf()
  # 65 and over
  hospital.filtered.old <-  hospital_spatiotemporal %>%
    filter(date == dates[[i]],
           age == "65 and over") %>%
    st_as_sf()
  # Regressors
  pred.filtered <-  stage_2_data.sf %>%
    filter(date == dates[[i]]) %>%
    st_as_sf() %>%
    st_transform(crs(hospital.filtered$geometry))
  # Regressors + 0 to 17
  joined.child <- st_join(pred.filtered,
                          hospital.filtered.child,
                    join = st_nearest_feature,
                    left = TRUE) 
  # Regressors + 18 to 64
  joined.adult <- st_join(pred.filtered,
                          hospital.filtered.adult,
                          join = st_nearest_feature,
                          left = TRUE)
  # Regressots + 65 and over
  joined.old <- st_join(pred.filtered,
                        hospital.filtered.old,
                        join = st_nearest_feature,
                        left = TRUE) 
  # Join all age data
  joined <- rbind(joined.child,
        joined.adult,
        joined.old)
  ### Add remaining hospitals
  hospital.filtered.remainder <- hospital.filtered %>%
    filter(hospital %in% unique(joined$hospital) ==FALSE) %>%
    st_join(pred.filtered,
            join = st_nearest_feature,
            left = TRUE) %>%
    select(hospital,age,fire_pm25,background_pm25,d2m,n,t2m,date = date.x) %>%
    as_tibble()
  # Get rid of geometry
  hospital.filtered.remainder$geometry <- NULL
  joined$geometry <- NULL
  # Summarise data
  joined_sum <- joined %>%
    group_by(hospital,age) %>%
    summarise(fire_pm25 = mean(fire_pm25,na.rm=T),
              t2m = mean(t2m,na.rm=T),
              background_pm25 = mean(background_pm25,na.rm=T),
              d2m = mean(d2m,na.rm=T),
              n = mean(n,na.rm=T)) %>%
    mutate(date = dates[[i]]) %>%
    rbind(hospital.filtered.remainder)
  # Add the day's data to the storage list
  store_hospital[[i]] <-  joined_sum %>% as_tibble()
  print(i)
}

# Convert list to dataframe
stage_2_final_data <- do.call(rbind,store_hospital)

stage_2_final_data.sf <- stage_2_final_data %>%
  mutate(date_num = as.numeric(date) - min(as.numeric(date)),
         dow = weekdays(date),
         zero = ifelse(n == 0,1,0),
         id = paste0(hospital," - ",age),
         n = as.integer(n)) %>%
mutate(age = as_factor(age),
         hospital = as_factor(hospital)) %>%
  filter("Hospital filtering here") %>% 
  as.data.frame() 
names(stage_2_final_data.sf$fire_pm25) <- NULL
# Fit model

test_model_feglm <- glm.nb(n ~ hospital+age+(age*hospital)*fire_pm25+(age*hospital)*background_pm25+age*d2m+age*t2m,
                              data = stage_2_final_data.sf) 

# Counterfactual dataset where fire pm2.5 is zero
stage_2_counterfactual.sf <- stage_2_final_data.sf %>% mutate(fire_pm25 = 0)
##### Counterfactual Simulation for chapter 4 ---------

# Predicted mean under true data
mean_actual <- predict.glm(test_model_feglm,newdata = stage_2_final_data.sf,type='response')
# Predicted mean under counterfactual data
mean_counterfactual <- predict.glm(test_model_feglm,newdata = stage_2_counterfactual.sf,type='response')
# Simulate from real distribution
simulation_data <- replicate(rnegbin(mean_actual,theta = test_model_feglm$theta),n = 10000)
# Simulate from counterfactual distribution
simulation_data.alternative <- replicate(rnegbin(mean_counterfactual,theta = test_model_feglm$theta),n = 10000)

# Aggregate real simulation results  
sim_agg <- simulation_data %>% 
  as_tibble() %>%
  cbind(stage_2_final_data.sf$age,
        stage_2_final_data.sf$hospital,
        stage_2_final_data.sf$date,
        stage_2_final_data.sf$fire_pm25) %>%
  rename(age = `stage_2_final_data.sf$age`,
         hospital = `stage_2_final_data.sf$hospital`,
         date = `stage_2_final_data.sf$date`,
         fire_pm25 = `stage_2_final_data.sf$fire_pm25`) %>%
  pivot_longer(cols= starts_with("V"),names_to = 'sim_number',values_to='n')

# Aggregate counterfactual results
sim_agg.alternative <- simulation_data.alternative %>% 
  as_tibble() %>%
  cbind(stage_2_final_data.sf$age,
        stage_2_final_data.sf$hospital) %>%
  rename(age.alt = `stage_2_final_data.sf$age`,
         hospital.alt = `stage_2_final_data.sf$hospital`) %>%
  pivot_longer(cols= starts_with("V"),names_to = 'sim_number.alt',values_to='n.alt')

### Join simulation results
sim_agg.joined <- sim_agg.alternative %>%
  cbind(sim_agg) %>%
  select(age,date,hospital,sim_number,n,n.alt,fire_pm25)
