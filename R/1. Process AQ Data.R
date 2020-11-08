#####################################################################
#### This script process AOD & PM2.5 data. It creates two data sets 
##### - one statewide grid of AOD, and one AOD-PM2.5 matched set
#####################################################################
set.seed(2020)
library(sf)
library(tidyverse)
library(lubridate)
library(raster)
library(dtplyr)
library(absmaps)
library(gstat)
select <- dplyr::select
crs_global <- 4283
### Process Data  -----------

#### Load AOD data

##### Read in 47 um raster list (this was processed on a Windows VM) 

AOD_47_2020 <- readRDS("Data/AOD 2015-2019/aod_list_47_2020.rda")
AOD_47_2019 <- readRDS("Data/AOD 2015-2019/aod_list_47_2019.rda")

##### Read in 47 um raster file names (this was processed on a Windows VM) 

file_list_2020 <- readRDS("Data/AOD 2015-2019/file_list_2020.rda")
file_list_2019 <- readRDS("Data/AOD 2015-2019/file_list_2019.rda")
##### Find number of days in sample

len_2020 <- nchar(file_list_2020[[1]])
len_2019 <- nchar(file_list_2019[[1]])

##### Pull raw dates from file list

dates_raw_2020 <- substr(file_list_2020,len_2020-35,len_2020-29)
dates_processed_2020 <- strptime(dates_raw_2020, format="%Y%j", tz="UTC")

##### Format raw dates correctly

dates_raw_2019 <- substr(file_list_2019,len_2019-35,len_2019-29)
dates_processed_2019 <- strptime(dates_raw_2019, format="%Y%j", tz="UTC")

##### Combine 2019 & 2020 Data

aod_all <- c(AOD_47_2019[as.Date(dates_raw_2019, format="%Y%j", tz="UTC") > as.Date('2019-10-31')],AOD_47_2020)
dates_all <- c(dates_processed_2019[as.Date(dates_raw_2019, format="%Y%j", tz="UTC") > as.Date('2019-10-31')],dates_processed_2020)

#### Load air quality (PM2.5 data)

load("Data/Air Quality/aq_hourly.Rda")

##### Download a map of Victoria

state2016     <- load_absmaps(area = "state", year = 2016,saveDirectory = "Data/SA3 Shapefiles",download = T) %>% filter(state_name_2016 == "Victoria")

##### Create new aq station location file with global (4283) rather than domestic (4326) CRS

aq_proj <- aq_sf %>%
  st_transform(crs_global) 

##### Create sf file with average daily PM2.5 in Victoria in sample period. 
aq_final <- aq_hour %>%
  lazy_dt() %>%
  filter(state == "Vic") %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(date,aq_id) %>%
  summarise(PM25 = mean(pm25,na.rm=TRUE)) %>%
  as_tibble() %>%
  right_join(aq_proj,by = "aq_id") %>%
  as_tibble() %>%
  filter(year(date) %in% c('2019','2020')) %>%
  st_as_sf() 

##### Create total AOD grid 

#### Create storage list
store <- vector(mode = 'list',length = length(aod_all))
#### Iterate across every date
for(i in 1:length(aod_all)){
  raster <- aod_all[[i]]
  #### Aggregate from 1km to 10k grid
  raster.aggregate <- aggregate(raster,fact = 10)
  #### Convert to sp (this isn't very good but the only way I can move to raster to sf seems to be through sp)
  raster.sp <- as(raster.aggregate,'SpatialPolygonsDataFrame')
  #### Convert to sf
  raster.sf <- st_as_sf(raster.sp)
  #### Change CRS and ensure grid is in Victoria
  raster.sf.filtered <- raster.sf %>%
    st_transform(crs_global) %>%
    st_join(distinct(state2016), left = FALSE)
  #### Attach date
  raster.sf.filtered.dated <- raster.sf.filtered %>%
    mutate(date = as.Date(dates_all[i]))
  store[[i]] <- raster.sf.filtered.dated
  print(i)
}
#### Create final data set
aod_all_areas <- do.call(rbind,store)

## Fill in Missing Data ------

#### create vector of unique dates
dates <- unique(aod_all_areas$date)
store_impute <- vector(mode = 'list',
                   length = length(dates))

### Complete grid for imputation
grd.sf <- unique.data.frame(aod_all_areas$geometry)
grd.sp <-as(grd.sf,Class = 'Spatial')

for(i in seq_along(dates)){
  ### Find date i
  date = dates[[i]]
  ### Select AOD data for date i
  aod_all_areas.day <-aod_all_areas[aod_all_areas$date == date,] %>% select(var1.pred = layer)
  ### Convert to sf
  aod_all_areas.day.sp <- as(aod_all_areas.day,Class = 'Spatial')
  ### Check if there is any missing data. If there is not, add date column,
  ### add to store_impute and move on. 
  if(length(aod_all_areas.day.sp$var1.pred) == length(grd.sf)){
    store_impute[[i]] <- aod_all_areas.day %>%
      mutate(date = dates[[i]])
                 print(i)
                 
  } else {
  ### Create missing data grid (sp)
  missing <- as(grd.sf[!(grd.sf %in% aod_all_areas.day$geometry)],Class = 'Spatial')
  ### Filling in missing AOD data with idw
  aod_all_areas.day.imp <- gstat::idw(var1.pred ~ 1,aod_all_areas.day.sp,newdata = missing,idp=2.0,nmax = 100)
  ### Join imputed AOD data and meterological data
  joined <- aod_all_areas.day.imp %>%
    st_as_sf() %>%
    select('var1.pred') %>%
    st_transform(crs_global) %>%
    rbind(aod_all_areas.day) %>%
    mutate(date = dates[[i]])
  store_impute[[i]] <- joined
  print(i)
  }
}
aod_all_matched_imputed <- do.call(rbind,store_impute)

### Match AOD to PM2.5 data -------
store_matched <- vector(mode = 'list',length = length(aod_all_matched_imputed))
for(i in 1:length(aod_all)){
  ### Filter data by date
  store_matched.filtered <- aod_all_matched_imputed %>% 
    filter(date == as.Date(dates_all[[i]])) %>%
    st_transform(crs_global)
  aq.filtered <- aq_final %>% filter(date == as.Date(dates_all[[i]]))
  #### Match AOD & AQ
  raster.matched <- aq.filtered %>%
    st_join(store_matched.filtered,left=T)
  store_matched[[i]] <- as_tibble(raster.matched)
  print(i)
}

pm25_aod.matched<- do.call(rbind,store_matched) %>%
  mutate(aq_id = as.character(aq_id)) %>%
  select(-date.y) %>%
  rename(date = date.x,layer = var1.pred)

