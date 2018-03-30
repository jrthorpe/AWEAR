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
d.start <- as.POSIXct("2018-02-01") # yyyy-mm-dd 
d.stop <- as.POSIXct("2018-03-16")

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
for (i in 1:length(datasets.all)){
  #browser()
  j <- datasets.all[[i]]
  if(!is.null(j)){
    #print(distinct(j,funf_version))
    cat("Data in ",names(datasets.all)[i],"comes from ", distinct(j,dsource)[,1],"\n");
    cat("Last reading in",names(datasets.all)[i],"is on",capture.output(max(j$timestamp)),"\n");
  }
}

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
# -- Life Space --

# Minimum convex polygon
# Reference: http://mgritts.github.io/2016/04/02/homerange-mcp/

# calculate distances from all points to centroid of the location data
xy <- select(location, lat, lon)
centroid <- apply(xy,2,mean)
distances <- sqrt(((xy[, 1] - centroid[1])^2) + ((xy[, 2] - centroid[2])^2))

# get subset of points within specified quantile of distances
indx <- 1:length(distances)
percentages <- indx[distances <= quantile(distances, .999)]
xy.subset <- xy[percentages, ]

# get minimum convex polygon
mcp.points <- chull(xy.subset[, 1], xy.subset[, 2]) # index of points that lie on mcp
mcp <- xy.subset[mcp.points,] # coords of mcp
mcp <- rbind(mcp[nrow(mcp), ], mcp) # repeat last point to close shape

mcp.poly<-Polygon(mcp) # creates a polygon object with area attribute (access uing @area)

plot_ly(xy, 
        x=~lon, 
        y=~lat, 
        type="scatter", 
        mode = "markers",
        name = "location data",
        marker = list(size = 2)
        #color = ~dsource
        )%>%
  #add_trace(x=~centroid[1], y=~centroid[2], mode = "markers", marker = list(size = 10), name = 'centroid') %>%
  add_trace(x=~mcp$lon, y=~mcp$lat, mode = "lines", name = 'mcp') %>%
  add_trace(x=~rxy$lon, y=~rxy$lat, mode = "markers", marker = list(size = 2))

rxy<-xy[sample(nrow(xy), 1000), ]


#testing plot_geo
# (deleted, need to start over)

# Standard Deviation Elipse

# Action Range

# Distance covered (including for trips only)

# -- Other --
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


