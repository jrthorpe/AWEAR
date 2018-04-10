# DATA CHECKER
#
#********************************************************************************************
# Author: Julia Thorpe
# Written for AWEAR Case Studies, January-May 2018, as part of an ongoing PhD
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

# Define constants:
folder <- "../../Data_dumps/dump_current_analysis/" # path to folder that holds multiple .csv files, downloaded from nightingale webportal
not.in.use <- c("bluetooth","hardware_info","wearable","wifi")
activity.selection <- c("Still", "Foot", "Vehicle", "Tilting", "Bicycle") # main activity types
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
d.start <- as.POSIXct("2018-02-10") # yyyy-mm-dd 
d.stop <- as.POSIXct("2018-02-1")

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

# QUALITY CHECK -------------------------------------------------------------------------------------------

# Is data coming from watch, phone or both?
#TODO: move into separate function
for (i in 1:length(datasets.all)){
  #browser()
  j <- datasets.all[[i]]
  if(!is.null(j)){
    #print(distinct(j,funf_version))
    cat("Data in ",names(datasets.all)[i],"comes from ", distinct(j,dsource)[,1],"\n");
    cat("Last reading in",names(datasets.all)[i],"is on",capture.output(max(j$timestamp)),"\n");
  }
}

#lapply(restructure, userid, d.start, d.stop)

# VISUALISE DATA-------------------------------------------------------------------------------------------

show_plots(datasets.all, to_plot)

# BEHAVIOURAL METRICS -------------------------------------------------------------------------------------------

# (THIS WILL MOVE INTO SEPARATE SCRIPT)

# Notes from reading the thesis:
# - what is our sampling frequency for GPS data?
# - 

# I: Mobility 
location<-datasets.all$location
location.projected <- mapproject(location$lon, location$lat,
                                 projection = "orthographic")


# Detect trajectories
# Detect stay points / hotspots
# Detect home

# TRAJECTORIES ----

# Version T1: Based on time and distance threshold
gps.log <- select(location, lat, lon, timestamp, intervals) %>% 
  arrange(timestamp)

# split into trajectories
threshold.t <- 10/60 # 10 minutes in hours
gps.traj <- mutate(gps.log,breakpoints=intervals>threshold.t) %>%
  mutate(traj_id=c(1,cumsum(gps.traj$breakpoints[-1])+1))

plot_ly(gps.traj,
        x=~lon,
        y=~lat,
        type = "scatter",
        mode = "lines+markers",
        color = ~as.factor(traj_id))


# Calculated directly on location data: 
# -> minimum convex polygon (area)
# -> ...


# -> Action range
# -> Standard deviation elipse 
# -> Distance covered (including for trips only - later)
# -> ...


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
  
#testing plot_geo
# (deleted, need to start over)

# Standard Deviation Elipse
# A function exists in python, need to find or write for R, eg using this explanation:
# http://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-statistics-toolbox/h-how-directional-distribution-standard-deviationa.htm




#- Calculating points of interest (or stops or stays) in the data: -#

# References: see files downloaded from GitHub. These use time and distance
# between points against thresholds to assign points to "stay events", based on
# work in an article by Zheng.




# Number of trips
# Number of hotspots
# Time spent out of / at home:
  
  # Get home location as statistical mode of all lat lon data for a user

loc.counts <- count(location,lat,lon) %>% arrange(desc(n))
home <- select(loc.counts[1,],c(lat,lon))
#NOTE: the above method does absolutely not work to detect home! (See Julia's
#data, where Glostrup was detected as home). This could be because as long as
#google maps is on, the system gets lots of data, whereas at home there won't be
#many recordings. Need to detect clusters and select most dense or something
#like that.

# plot of unique locations in dataset with marker size based on number of occurences
plot_ly(loc.counts, 
        x = ~lon, 
        y = ~lat, 
        type = "scatter", mode = 'markers',
        marker = list(size = ~n))

plot_ly(loc.counts, 
        x = ~lon, 
        y = ~lat, 
        type = "scatter", mode = 'markers',
        marker = list(size = ~n))

  # Classify all points as "at home" if within 100m of home location, "out of home" if not
  # Separate into bouts with periods etc...

# II: Activity:

# Activity bouts (active and sedentary)

# Stepcounts

## FOR TESTING

# for looking at drifting timestamp and watch data
plot_ly(datasets.all$step_count, x=~index, y = ~timestamp, type="scatter", mode="markers", color=~funf_version)


# trying out code found on Github:
source("git_mobility.R")
library("DataCombine")
library("zoo")
library("plyr")
mobility_stay<- stayevent(mobility, coor = c("lon","lat"), time = "datetime", dist.threshold = 100,
                 time.threshold = 30, time.units = "mins", groupvar = "id")

mobility_stay_test<- stayevent(gps.log, coor = c("lon","lat"), time = "timestamp", dist.threshold = 100,
                          time.threshold = 30, time.units = "mins")
plot_ly(mobility_stay_test,
        x=~lon,
        y=~lat,
        type = "scatter",
        mode = "markers",
        color = ~as.factor(stayeventgroup))
# next task:
# plot all points in light grey then "add trace" to plot stayevents in colour

