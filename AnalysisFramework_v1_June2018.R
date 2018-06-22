# AWEAR: BEHAVIOURAL ANALYSIS FRAMEWORK
#
#********************************************************************************************
# Author: Julia Thorpe
# Written for AWEAR Case Studies, January-June 2018, as part of an ongoing PhD
# project on Engineering Systems Design in Healthcare at DTU in collaboration
# with Rigshospitalet-Glostrup

# This script ...

# SETUP -------------------------------------------------------------------------------------------

# Load packages and functions.
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
library(geosphere) # added in v3 of data checker
library(mapview) # added in v1 of analysis framework

# For DBSCAN clustering
library(fpc)
library(dbscan)
library(factoextra)

# Load custom functions:
source("M:/PhD_Folder/CaseStudies/Data_analysis/source/JRT_utils.R")
source("M:/PhD_Folder/CaseStudies/Data_analysis/source/jrt_mobility.R")
source("M:/PhD_Folder/CaseStudies/Data_analysis/source/metrics/mcp.R")

# Define constants:
folder <- "M:/PhD_Folder/CaseStudies/Data_dumps/dump_current_analysis/" # path to folder that holds multiple .csv files, downloaded from nightingale webportal
not.in.use <- c("activity", "battery", "experience_sampling", "screen", "step_count", #for when only using location data
                "bluetooth","hardware_info","wearable","wifi","calllog","sms")
to_plot <- c("activity", # from data checker, for debugging purposes only (visualise all datasets)
             "battery",
             "exps",
             "location",
             "screen",
             "steps")
users <- list(julia = "93a6d31c-e216-48d8-a9a2-f2f72362548d",
              dean = "b1316280-38a6-45e1-9bb9-7afb2a1a2a96",
              luna = "4afe61c5-c5cc-4de4-8f50-9499057668ad",
              per = "7bf5fec3-f46f-419e-9573-001fe9b47d81",
              verena = "a0c74a88-b293-4b31-92c0-967502b28132",
              agzam = "52152dda-0f54-46d1-affa-6dea85f840dc",
              anja = "36f9d061-c5e9-4f30-91b2-83351ff40288",
              nina = "b12888e9-ee4f-42bc-9494-673e5a005e0e",
              #P01FA = "6321f7ef-a958-44ad-b8e5-5aa04bc004e1",
              P03JJ = "e9f44eb5-8962-4894-83c6-783025c6eaea",
              P06SS = "f9f24838-c844-42d4-8343-b20ebdd220f3",
              P07MG = "4fbfddd0-a346-41c8-be8d-f8804b5068d3",
              P08UH = "a85f299e-7a09-4ba7-bc19-8a200c2686c2",
              P10JL = "d05fa984-8d3b-4405-b417-211d1a3f50d6",
              P13NB = "75cc5240-7805-4550-aeae-873df9984710")

# SETTINGS ==========================================
# Select user and period of interest:

# select user:
userid <- users$dean

# set the period of interest:
d.start <- as.POSIXct("2018-02-26") # yyyy-mm-dd 
d.stop <- as.POSIXct("2018-03-28")

# IMPORT AND RESTRUCTURE DATA ------------------------

datasets.all <- get.data(folder, not.in.use) %>% 
  lapply(restructure, userid, d.start, d.stop)
datasets.all <- Filter(function(x) !is.null(x)[1],datasets.all) # remove any null dataframes

remove(d.start, d.stop, folder, not.in.use, userid, users)

# VISUALISE DATA --------------------------------------

#show_plots(datasets.all, to_plot)

# MOBILITY  ==========================================

#** Moblity Setup: create datasets and variables required for mobility calculations -------

# define variables:
loc.accuracy <- 25 # threshold for accuracy of location points in meters
dT <- 5  # delta T, time window in minutes
dD <- 100 # delta D, diagonal distance boundary in meters
time.threshold.stay <- 10 # minimum duration of a stay, in minutes
time.threshold.go <- 5 # cut-off for filtering out "go" events to/from same location, in minutes
dist.threshold <- 30 # distance in meters within which two centriods belong to same stay location
doi <- as.POSIXct("2018-03-09") # day of interest

# preprocessing steps
gps.log.pre <- datasets.all$location %>% filter(accuracy<=loc.accuracy) # get rid of data points with low accuracy
keep.ratio <- nrow(gps.log.pre)/nrow(datasets.all$location)

# GPS log file
gps.log <- gps.log.pre %>% 
  select(lat, lon, timestamp, intervals.alt, dates, times) %>% 
  filter(timestamp>=doi, timestamp<=(doi %m+% days(8))) # get timeframe of just the day of interest
remove(gps.log.pre)

# calculate home coordinates based on location data
home <- find_home(gps.log,"lat","lon")

#** Extract trajectories: get series of stay/go events ("mobility traces") ====

# # For a single day:
# # append trajectory information to the gps log file, including if the points is
# # a stay/go and location ID (note: all "go" points are assigned location ID of 0)
# singleday <- gps.log %>% filter(timestamp>=doi, timestamp<=(doi %m+% days(1)))
# gps.traj.day <- get_trajectories(df=singleday,
#                             dT = dT,
#                             dD = dD,
#                             T.stay = time.threshold.stay,
#                             T.go = time.threshold.go,
#                             dist.threshold = dist.threshold)
# 
# # summarise all trajectory events
# traj.summary.day <- gps.traj.day %>% group_by(traj.event) %>%
#   summarize(T.start = min(timestamp), is.stay = median(is.stay), loc.id = mean(loc.id)) %>%
#   mutate(T.end = c(T.start[-1],max(gps.traj.day$timestamp))) %>%
#   mutate(durations = difftime(T.end,T.start, units = "mins"))


# For mulitday sets:
gps.traj <- gps.log %>% group_by(dates) %>% do(get_trajectories(.,
                                                             dT = dT,
                                                             dD = dD,
                                                             T.stay = time.threshold.stay,
                                                             T.go = time.threshold.go,
                                                             dist.threshold = dist.threshold))

traj.summary <- gps.traj %>% group_by(dates,traj.event) %>%
  summarize(T.start = timestamp[1], T.end = timestamp[length(timestamp)], is.stay = median(is.stay), loc.id = mean(loc.id)) %>%
  mutate(T.end = c(T.start[-1],T.end[length(T.end)])) %>%
  mutate(durations = difftime(T.end,T.start, units = "mins"))

# saveRDS(traj.summary,"M:/PhD_Folder/CaseStudies/Data_analysis/output/dean_summary.Rds")
# write.csv(traj.summary,"M:/PhD_Folder/CaseStudies/Data_analysis/output/dean_summary.csv")

#** Plot results for visual confirmation ----

# basic plot
plot_ly(gps.traj,
        x=~lon,
        y=~lat,
        type = "scatter",
        mode = "markers",
        color = ~as.factor(loc.id)) #or traj.event

# locations/events overlaid over GPS trace for all datapoints:
plot_ly(gps.traj,
        x=~lon, y=~lat,
        type = "scatter",
        mode = "lines+markers",
        color = I('grey40')) %>%
  add_trace(x=~lon, y=~lat,
            type = "scatter",
            mode = "lines+markers",
            color = ~as.factor(traj.event), #or loc.id
            colors = "Set1")

# map in background:
mappoints <- gps.traj
coordinates(mappoints) <- ~ lon + lat
proj4string(mappoints) <- "+init=epsg:4326"
mapview(mappoints, zcol = "traj.event", burst = TRUE, map.types = "OpenStreetMap") #or loc.id

# METRICS CALCULATION ================================================

# -For each day/week: 
#  .. Number of trips/stays
#  .. Time spent out of / at home
#  .. Distance covered in trajectories
#  .. "Action Range"
#  .. Plot the day: over time, plot a bar with "home", "transit", "location A", "location B" etc
#  .. Work out how to annotate with the logbooks

#** Spatial (Lifespace distance/area metrics) -----------------------

# Minimum convex polygon
# Reference: http://mgritts.github.io/2016/04/02/homerange-mcp/

xy <- select(gps.log, lat, lon)
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



# Action Range: (mean and max per period)
# Straight line distance between home and most distal point of a journey

# For each trajectory event (including both stays and moves):
# Measure distance from each point to home
# In trajectory summary, store max
# Could just do this for moves, since this is already done with centroids for the stays


#** Temporal -----------------------



# get centroids of all stay locations:
stay.centroids <- gps.traj %>% filter(loc.id > 0) %>% group_by(loc.id) %>% 
  summarize(c.lon = mean(lon),c.lat = mean(lat))

# identify location id for home
stay.centroids %<>% mutate(dist.home = distm(x=rev(home), 
                                            y=select(stay.centroids,c.lon,c.lat), 
                                            fun=distGeo))

cat("closest stay location to detected home is: ",
    round(min(stay.centroids$dist.home),2),"meters" )

N.stay <- stay.centroids %>% filter(dist.home > min(dist.home)) %>% nrow()
N.go <- traj.summary %>% filter(loc.id == 0) %>% nrow()

# Time outside of home
home.id <- stay.centroids %>% filter(dist.home == min(dist.home)) %>%
  select(loc.id) %>% as.numeric()

T.not.home <- 
  sum((traj.summary %>% filter(loc.id != home.id))$durations) %>% as.numeric() # duration of all events out of home
T.go <- 
  sum((traj.summary %>% filter(loc.id == 0))$durations) %>% as.numeric() # duration of all "go" events
T.stay.out <- 
  sum((traj.summary %>% filter((loc.id != home.id) & (loc.id > 0)))$durations) %>% as.numeric() # duration of all "stays" outside home

# DEBUGGING AREA ----

# Issue (2): Anomolies breaking up stays so that not detected based on time
# interval not being long enough Solution: one cause is low accuracy points that
# can be filtered out, e.g. by removing all points with accurace less than 25m.
# To justify this number, can show density distributions plots.

# Get relevant data to see problem
testInheritedMethods()<-datasets.all$location %>% select(lat, lon, timestamp, intervals.alt, dsource,accuracy) %>%
  filter(timestamp>=as.POSIXct("2018-02-14"), timestamp<=as.POSIXct("2018-02-15"))

# Show distributions indicating that selecting 30m as an upper limit will keep most of the data
hist(test$accuracy,breaks=100) # histogram of gps data accuracies
d <- density(test$accuracy) # returns the density data 
plot(d) # plots the results

# Attempts at using group_by and do for getting trajectories by date

test1 <- gps.log %>% group_by(dates) %>% do(data.frame(get_trajectories(.,
                                                                        dT = dT,
                                                                        dD = dD,
                                                                        T.stay = time.threshold.stay,
                                                                        T.go = time.threshold.go,
                                                                        dist.threshold = dist.threshold)))

test2 <- gps.log %>% group_by(dates) %>% do(get_stays(., dT, dD))

# End of DEBUGGING AREA ---

# OTHER ATTEMPTS ----

# working on multiday sets:

test2 <-  lapply(gps.log.bydates, get_trajectories,
                 dT = dT,
                 dD = dD,
                 T.stay = time.threshold.stay,
                 T.go = time.threshold.go,
                 dist.threshold = dist.threshold)

gps.log.bydates <- split(gps.log, gps.log$dates) #%>% 
for(d in gps.log.bydates){
  
  traj.bydate <- get_trajectories(df=d,
                                  dT = dT,
                                  dD = dD,
                                  T.stay = time.threshold.stay,
                                  T.go = time.threshold.go,
                                  dist.threshold = dist.threshold)
  # summarise all trajectory events
  summary.bydate <- traj.bydate %>% group_by(traj.event) %>%
    summarize(T.start = min(timestamp), is.stay = median(is.stay), loc.id = mean(loc.id)) %>%
    mutate(T.end = c(T.start[-1],max(gps.traj$timestamp))) %>%
    mutate(durations = difftime(T.end,T.start, units = "mins"))
}