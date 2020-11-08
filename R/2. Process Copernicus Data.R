#####################################################################
#### This script process ECMWF meteorological data. It creates one 
#### data set. 
#####################################################################

set.seed(2020)
library(tidyverse)
library(bomrang)
library(dtplyr)
library(sf)
library(raster)
library(tidync)
library(lubridate)

# New ---------------

### File names
file1 <- 'Data/data_final/meterological_data.nc'

### Read Data (File 1)
#### Load as a brick to get file names
brick_1 <-brick(file1)
#### Load as nc to get a workable file
sf_1.nc <- tidync(file1)
#### Convert nc to tibble
sf_1 <- sf_1.nc %>% hyper_tibble()
#### Pulling dates is painful. Names are in a X.YYYY.MM.DD.HH.MM format. We reformat the dates here.
names <- unique(brick_1@data@names)
year <- substr(names,2,5)
month <- substr(names,7,8)
day <- substr(names,10,11)
hour      <- substr(names,13,14)
minute    <- substr(names,16,17)
date_time <- as_datetime(paste0(year,'-',month,'-',day,' ',hour,":",minute,':00'))

#### Add names
sf_1.dated <- sf_1 %>%
  filter(time > 1050430) %>%
  rename(time_num = time)
sf_1.dated$time <- date_time[as.Date(date_time) > as.Date('2019-10-31')][sf_1.dated$time_num - min(sf_1.dated$time_num)+1]

### Convert to sf
sf_1.geom <- st_as_sf(sf_1.dated,
                 coords = c('longitude',
                            'latitude'),
                 crs = 'WGS84')
### Filter out values that are not in Victoria
sf_1.filtered <- sf_1.geom %>%
  st_transform(crs_global) %>%
  st_filter(state2016$geometry,join = st_contains)
### Filter Dates ()
sf_1.filtered <- sf_1.filtered %>%
  mutate(Date = as.Date(time)) %>%
  filter(Date < as.Date('2020-02-01'),
         Date > as.Date('2019-10-31'))

### Average to daily data
cop_sf_1 <- sf_1.filtered %>%
  as_tibble() %>%
  mutate(id = as.character(geometry)) %>%
  group_by(Date,id) %>%
  mutate(u10   = mean(u10),
         d2m   = mean(d2m),
         t2m   = mean(t2m),
         skt   = mean(skt))

### Data Set 1
cop_sf_1.sf <- cop_sf_1 %>%
  ungroup() %>%
  st_as_sf() %>%
  st_transform(st_crs(aq_sf))

### Create matching index between ground data and satellite data
match_index_1 <- st_nearest_feature(aq_sf,cop_sf_1.sf)
### Apply matching index to get the grids that align with ground stations
filter_list_1 <- cop_sf_1.sf$id[match_index_1]
### Select data
cop_sf_1.sf.matched <- cop_sf_1.sf %>%
  as_tibble() %>%
  lazy_dt() %>%
  filter(id %in% filter_list_1) %>%
  group_by(Date,id) %>%
  summarise(u10   = mean(u10),
            d2m   = mean(d2m),
            t2m   = mean(t2m),
            skt   = mean(skt)) %>%
  as_tibble() 
### Add air quality names
match <- data.frame(aq_sf$aq_id,filter_list_1)
colnames(match) <- c('aq_id','id')
cop_sf_1.sf.matched <- cop_sf_1.sf.matched %>% 
  left_join(match)%>%
  rename(date = Date) %>%
  mutate(aq_id = as.character(aq_id))
