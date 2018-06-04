# DATA CHECKER
#
#********************************************************************************************
# Author: Julia Thorpe
# Written for AWEAR Case Studies, January-June 2018, as part of an ongoing PhD
# project on Engineering Systems Design in Healthcare at DTU in collaboration
# with Rigshospitalet-Glostrup

# This script imports and plots data from a specified participant and
# timeframe. This is used to regularly check that the data is being recorded and
# looks as expected.

# SETUP -------------------------------------------------------------------------------------------
# Set the drive, load packages and functions.

setwd("M:/PhD_Folder/CaseStudies/Data_analysis/source")

# Load required packages:
# (Note: had trouble installing packages which I resolved by using .libpaths() to find where packages are installed and manually moving the packages from wherever R temporarily installed them to there)

# For mobility pack (downloaded from Git)
library(DataCombine)
library(zoo)
library(plyr)

library(data.table)
library(dtplyr)
library(dplyr)
library(magrittr)
library(lubridate)
library(ggplot2)
library(reshape2)
library(maptools)
library(plotly)
library(maps)
library(mapproj)
library(sp)
library(caTools)
library(geosphere) # added in v3

# Load custom functions:
source("JRT_utils.R")
source("git_mobility.R")
source("jrt_mobility.R")

# Define constants:
folder <- "../../Data_dumps/dump_current_analysis/" # path to folder that holds multiple .csv files, downloaded from nightingale webportal
not.in.use <- c("bluetooth","hardware_info","wearable","wifi")
to_plot <- c("activity", 
             "battery", 
             "exps", 
             "location", 
             "screen", 
             "steps")
users <- list(julia = "93a6d31c-e216-48d8-a9a2-f2f72362548d",
              dean = "b1316280-38a6-45e1-9bb9-7afb2a1a2a96",
              per = "7bf5fec3-f46f-419e-9573-001fe9b47d81",
              P01FA = "6321f7ef-a958-44ad-b8e5-5aa04bc004e1",
              P03JJ = "e9f44eb5-8962-4894-83c6-783025c6eaea")

# SETTINGS -------------------------------------------------------------------------------------------
# Select user and period of interest:

# Test users:
userid <- users$julia

# set the period of interest
d.start <- as.POSIXct("2018-02-11") # yyyy-mm-dd 
d.stop <- as.POSIXct("2018-02-18")

# IMPORT AND RESTRUCTURE DATA -------------------------------------------------------------------------------------------

# Notes:
# x = not used (currently not planned for analysis)
# tba = to be added
# 
# 1. "activity.csv" -- in use, coded activity types
# 2. "battery.csv" -- in use
# 3. "bluetooth.csv" -- x
# 4. "calllog.csv" -- tba
# 5. "experience_sampling.csv" -- in use, answers to daily self-reports
# 6. "hardware_info.csv" -- x
# 7. "location.csv" -- in use, GPS data from watch and phone
# 8. "screen.csv" -- in use, screen on/off transitions (phone)
# 9. "sms.csv" -- tba
# 10. "step_count.csv" -- in use, step counts from watch and phone
# 11. "wearable.csv" -- not in use, used to be stepcount from watch (empty in dump from 2017-11-28)
# 12. "wifi.csv" -- x

datasets.all <- get.data(folder, not.in.use) %>% 
  lapply(restructure, userid, d.start, d.stop)
datasets.all <- Filter(function(x) !is.null(x)[1],datasets.all) # remove any null dataframes

remove(d.start, d.stop, folder, not.in.use, userid, users)

# VISUALISE DATA-------------------------------------------------------------------------------------------

#show_plots(datasets.all, to_plot)

# MOBILITY -------------------------------------------------------------------------------------

# Notes from reading the thesis:
# - what is our sampling frequency for GPS data?
# - 

# References: see files downloaded from GitHub. These use time and distance
# between points against thresholds to assign points to "stay events", based on
# work in an article by Zheng.

# -Clean up whole stay event section and put into utils where it makes sense
# -Read through code for detecting stay events and write up logic
# -Create my own based on above with improvements (eg in numbering of stays and trajectories)

# -For each day: 
#  .. Number of trips/stays
#  .. Time spent out of / at home
#  .. Distance covered in trajectories
#  .. Plot the day: over time, plot a bar with "home", "transit", "location A", "location B" etc
#  .. Work out how to annotate with the logbooks

location<-datasets.all$location # Get location dataset
home <- findhome(location,"lat","lon")

# location.projected <- mapproject(location$lon, location$lat,projection = "orthographic") # project GPS data onto 2D plane (currently not in use, need to look into it)

gps.log <- select(location, lat, lon, timestamp, intervals, dates, times) # Get GPS log file
#gps.traj <- filter(gps.log,timestamp>=as.POSIXct("2018-02-11"), timestamp<=as.POSIXct("2018-02-12")) # reduce to single day
gps.traj <- gps.log #for entire period in dataset (instead of single day)

# TODO: create own code for this based on paper
mobility_stay_test<- stayevent(gps.traj, coor = c("lon","lat"), time = "timestamp", dist.threshold = 60,
                               time.threshold = 10, time.units = "mins")

# create group column in main dataset with results from stay event detection above 
gps.traj %<>% mutate(stayeventgroup=mobility_stay_test$stayeventgroup)  %<>% 
  mutate(stay=(!is.na(stayeventgroup))*1)  %<>% # create binary "stay" group 
  mutate(staygo_group=cumsum(c(0,abs(diff(stay))))+1) # allocate a group number for all stay and go events (counts up)

# Stay events: get centroids and times
detach(package:plyr) #need to find out whether I need plyr at all
stay.info <- gps.traj %>%
  filter(stay==1) %>%
  group_by(staygo_group) %>%
  summarize(c.lat = mean(lat),c.lon = mean(lon),arrive = min(timestamp),depart=max(timestamp))

# Classify stays as "home" or "other" based on distance
#home_threshold <- 16 # in meters; think this was based on how far away 2 places could be within same 4 decimals level of GPS info used for home calculation
home_threshold <- 25 # 
p1<-rev(home) # reverses elements to get lat and lon in correct order
p2<-stay.info[,c("c.lon","c.lat")]

stay.info %<>% mutate(dist2home=distGeo(p1, p2, a=6378137, f=1/298.257223563)) %>% # create "home" group column (TRUE/FALSE) in stay.info summary
  mutate(home=dist2home<home_threshold)

remove(mobility_stay_test,p1,p2)


# TODO: get google map in the background using ggmap and ggplotly, see 2.2.4 here: http://plotly-book.cpsievert.me/maps.html
# TODO: make traces have proper names

# Plot stay centroids and home:
plot_ly(stay.info,
        x=~c.lon,
        y=~c.lat,
        color = ~home,
        type = "scatter",
        mode = "markers") %>%
add_trace(x=home$lon4, y=home$lat4, name='detected home', color = NA, marker = list(color = 'red'))

# Calculate time spent at home and/or out of home

# To "fill in" nights at home, set arrive/depart times to midnight for the first and last home stays respectively
N<-nrow(stay.info)
if(stay.info$home[1]) stay.info$arrive[1] <- round(stay.info$arrive[1],"days")
if(stay.info$home[N]) stay.info$depart[N] <- round(stay.info$depart[N],"days")

# Get all "home" stay event durations
stay.info %<>% mutate(durations=as.numeric(difftime(depart,arrive,units="hours")))

# Sum all home durations in hours
time.home<-sum(stay.info$durations[stay.info$home])

# PUWYLO:
#TODO: create function to do this by day so it can be applied to dataset of any period
# maybe: group_by day, then by event type (home, stay out, journey), then caculate durations


home.groups <- (stay.info %>% filter(home) %>% select(staygo_group))[[1]]
out.groups <- (stay.info %>% filter(!home) %>% select(staygo_group))[[1]]

gps.traj %<>% mutate(event_type=NA)

gps.traj$event_type[gps.traj$staygo_group %in% home.groups]="stay_home"
gps.traj$event_type[gps.traj$staygo_group %in% out.groups]="stay_out"
gps.traj$event_type[is.na(gps.traj$stayeventgroup)]="go"

# gps.traj %<>% cbind(select(location,dates,times))

plot_ly(gps.traj,
      x=~dates %>% format.Date(format="%d/%m"),
      y=~times %>% format.Date(format="%H:%M"),
      type="scatter",
      mode="markers",
      color=~event_type,
      #colors="Set1"
      colors=c("orange","grey","black")) %>%
  layout(yaxis = list(title = 'Time of day'),
         xaxis = list(title = 'Date'))


# todo: look at the "mobility boundaries" from the assessment. Should I be calculating this?
# todo: create time/day plot with colours showing whether it is "home","out" or "go" events


  

# Plot results overlaid over GPS trace for all datapoints:
plot_ly(gps.traj,
        x=~lon,
        y=~lat,
        type = "scatter",
        mode = "lines+markers",
        marker = list(size = 5),
        line = list(width = 1),
        color = I('grey40')) %>%
  add_trace(x=~lon,
            y=~lat,
            type = "scatter",
            mode = "markers",
            marker = list(size = 10),
            color = ~as.factor(stayeventgroup),
            colors = "Set1")

# my attempt at the signal processing method
dist <- c(0, distGeo(gps.traj[c("lon","lat")],a=6378137, f=1/298.257223563)) #append 0 to keep same length (at the first point, no distance is covered)
cum_dist<-cumsum(dist)


#...need to implement rest of method, so far just the distance signals...

gps.traj %<>% cbind(dist,cum_dist)


# Plot the stay events on the cumulative distanc plot for a visual comparison
plot_ly(gps.traj,
        x=~timestamp,
        y=~cum_dist,
        #name="cumulative distance (m)",
        type = "scatter",
        mode = "lines",
        color = I('black')
) %>%
  add_trace(x=~timestamp,
            y=~cum_dist,
            type = "scatter",
            mode = "markers",
            marker = list(size = 10),
            color = ~as.factor(stayeventgroup),
            colors = "Set3" 
  )


# LIFESPACE (DIST/AREA METRICS) ----

# Minimum convex polygon
# Reference: http://mgritts.github.io/2016/04/02/homerange-mcp/
source("metrics/mcp.R") #TODO: move to top

xy <- select(location, lat, lon)
mcp.results <- get_mcp(xy, quantile=.99)

mcp <- mcp.results$mcp
mcp.poly <- mcp.results$mcp.poly
centroid <- mcp.results$centroid

# Plot results
plot_ly(xy, 
        x=~lon, 
        y=~lat, 
        type="scatter", 
        mode = "markers",
        name = "location data",
        marker = list(size = 2)
)%>%
  add_trace(x=~centroid[2], y=~centroid[1], mode = "markers", marker = list(size = 10), name = 'centroid') %>%
  add_trace(x=~mcp$lon, y=~mcp$lat, mode = "markers+lines", name = 'mcp')


# Standard Deviation Elipse
# A function exists in python, need to find or write for R, eg using this explanation:
# http://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-statistics-toolbox/h-how-directional-distribution-standard-deviationa.htm