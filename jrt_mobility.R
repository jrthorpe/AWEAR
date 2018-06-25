# HOME LOCATION DETECTION ----
find_home <- function(df,lat,lon){
  # Get "mode" of all points after dropping one decimal place (4 places)
  # https://en.wikipedia.org/wiki/Decimal_degrees
  
  # df containing latitude and longitude with more than 4 decimal places
  # lat: name of latitude column
  # lon: name of longitude column
  #browser()
  locations4 <- data.frame(lat4=round(df[,lat],4),lon4=round(df[,lon],4))
  loc.counts <- locations4 %>% count(lat4,lon4, sort=TRUE)
  home <- as.numeric(select(loc.counts[1,],c(lon4,lat4)))
  
  return(home)
}

update_home <- function(df,home,dist.threshold){
  # df containing latitude, longitude and "is.stay" column
  # home: current home estimation (based on find_home function)
  # dist.threshold: distance in meters from current home coords within which a point is considered "home"

  # update home location based on all stays that are home (stay points close to current home estimation)
  df.tmp <- ungroup(df) %>% mutate(homedist = distGeo(home,
                                                            ungroup(df) %>% select(lon,lat),
                                                            a=6378137, f=1/298.257223563))
  home.updated <- df.tmp %>% filter(is.stay==1,homedist<dist.threshold) %>%
    select(lon, lat) %>% summarise_all(mean) %>% as.numeric()
  
  return(home.updated)
}

# TRAJECTORY EXTRACTION ("MOBILITY TRACES") // FULL STAY DETECTION ALGORITHM  ----
get_trajectories <- function(df,
                            dT,
                            dD,
                            T.stay,
                            T.go,
                            dist.threshold){
  ##* Initial identification of stays ---------
  #browser()
  traj <- get_stays(df, dT, dD) 
  
  ##* Spatial clustering to identify stay locations ---------
  
  # cluster all "stay" data points using density-based clustering (DBSCAN method)
  clusters <- cluster_spatial(traj,dist.threshold)
  
  # append info to gps.traj as location ID variable
  traj %<>% mutate(loc.id=clusters) %>%
    mutate(traj.event=get_events(loc.id))
  traj %<>% mutate_cond(loc.id == 0, is.stay = 0) # uupdate is.stay column where outliers are now "go" points
  
  # summarise the trajectories according to "stay" and "go" events, indicating locations
  traj.summary <- traj %>% group_by(traj.event) %>%
    summarize(T.start = min(timestamp), is.stay = min(is.stay), loc.id = mean(loc.id)) %>%
    mutate(T.end = c(T.start[-1],max(traj$timestamp))) %>%
    mutate(durations = difftime(T.end,T.start, units = "mins"))
  
  ##* Temporal clustering ---------
  #browser()
  # merge together stay or go events that are interrupted briefly 
  merge.temporal <- merge_temporal(traj.summary,T.go=T.go, T.stay=T.stay)
  
  #browser()
  if(nrow(merge.temporal)>0){
  # assign all filtered events the corresponding stay location id (NOTE: loc.id is zero for "go" data)
  for(m in 1:nrow(merge.temporal)){
    tmp <- traj$traj.event==merge.temporal[m,"traj.event"]
    traj[tmp,"loc.id"] <- merge.temporal[m,"loc.id"]
  }
    
  # update is.stay column where stays become goes or vice versa based on new loc.id
  traj %<>% mutate(is.stay = ifelse(loc.id == 0, 0, 1)) 
  
  # update traj.events based on revised loc.id column
  traj %<>% mutate(traj.event=get_events(loc.id))
  }
  
  # --- end of mobility traces algorithm --- 

  
  return(traj)
}


# STAY DETECTION: DISTANCE & TIME ----
get_stays <- function(df, deltaT, deltaD){
  # Detect stay points by measuring distance across the area covered in time windows of deltaT minutes
  # df:
  # deltaT:
  # deltaD:
  
  #TODO: create logic for if the gps.log dataframe has other names for timestamp, lat and lon timestamp=NULL,lat=NULL,lon=NULL,
  
  #df<-gps.log #for debugging
  #distances <- vector(length = nrow(df)) #debugging
  #t.start <- df$timestamp #debugging
  
  is.stay <- vector(mode = "logical", length = nrow(df))
  i <- 1
  while (i <= nrow(df)){
    
    # get all GPS data points in a timeframe of deltaT
    points <- df %>% 
      filter(timestamp>=timestamp[i], timestamp<=timestamp[i]+deltaT*60) %>%
      select(lat, lon)
    #t.start[i:(i+nrow(points)-1)] <- df$timestamp[i] #debugging
    
    # diagonal distance across bounding box enclosing the selected points
    bob <- data.frame(bbox(SpatialPoints(points)))
    ddistance <- distGeo(bob$min, bob$max, a=6378137, f=1/298.257223563) #diagonal distance across bounding box
    #distances[i:(i+nrow(points)-1)] <- ddistance #debugging
    
    # if distance below threshold all points are "stay", otherwise current point is "go"
    if(ddistance <= deltaD){ 
      is.stay[i:(i+nrow(points)-1)] <- TRUE
      i <- (i+nrow(points)) # go to next point after current stay set
    } else{
      is.stay[i] <- 0
      i <- (i+1) # go to next point
    }
    
  }
  
  # Create trajectories data frame indicating if stay/go and with stay/go event numbers 
  stays <- cbind(df,is.stay) 
  
  return(stays)
  
}

get_events <- function(dat){
  # explain function...

  # Inputs explained:
  # dat: column from trajectories data frame that describes movements, e.g. is.stay or loc.id
  
  #dat <- gps.traj[,split.by]
  change.detect <- c(0,diff(dat)!=0)
  events <- cumsum(change.detect)+1
  
  return(events)
}




# MERGE STAYS: SPATIAL PROXIMITY ----

cluster_spatial <- function(df, dist.threshold){
  # explain function...
  # Note: used distance matrix as an input to ensure the use of distances that are meaningful for GPS data

  # Inputs explained:
  # gps.traj: required inputs are lon,lat,is.stay
  # dist.threshold:

  points <- df %>% filter(is.stay==1) %>% select(lon,lat) # get GPS coords
  distances <- distm(points, fun=distGeo) # get a distance matrix of all points to all other points
  min.points <- 2 # specify MinPts parameter for clusterering algorithm "dbscan"

  # For selecting appropriate eps value:
  # dbscan::kNNdistplot(dist.df, k = min.points)
  # abline(h = dist.threshold, lty = 2)

  # Apply dbscan clustering using the distance matrix as an input
  db <- fpc::dbscan(distances, method="dist", eps = dist.threshold, MinPts = min.points)

  # # For visualising the clusters:
  # fviz_cluster(db, data = distances, stand = FALSE,
  #              ellipse = FALSE, show.clust.cent = FALSE,
  #              geom = "point",palette = "jco", ggtheme = theme_classic())

  # Get clusters vector for gps.traj
  
  clusters <- vector(mode="numeric",length = nrow(df)) # vector of zeros
  clusters[df$is.stay==1] <- db$cluster # cluster ID's assigned to stay rows, outliers are 0 therefore become same as "go" points.
  #clusters[gps.traj$is.stay==0] <- -1

  return(clusters)
}

# MERGE STAYS: SPATIAL PROXIMITY -- NO LONGER IN USE! ----

merge_spatial <- function(points, dist.threshold){
  # explain function...
  # points: required inputs are lat,lon
  # dist.threshold: 
  
  merges<-list()
  while (nrow(points)>1) {
    
    # select first point in the set and calculate distances to all others
    tmp <- select(points[1,],lon,lat)
    distances <- distGeo(tmp, select(points,lon,lat), a=6378137, f=1/298.257223563)
    
    # test distances against threshold
    is.close <- distances<dist.threshold
    
    # store group numbers of any close stay events
    if (sum(is.close)>1){
      merge.set <- filter(points,is.close) %>%
        select(event.id)
      merges <- c(merges, merge.set)
    } 
    
    # remove current centroid and any that are merged
    points %<>% filter(!is.close)  #(note: the current centroid is filtered out as it is close to itself)
    
  }
  
  N <- length(merges)
  merges <- setNames(merges,paste0("merge", 1:N))
  
  return(merges)
}

# MERGE STAYS: TEMPORAL PROXIMITY ----

merge_temporal <- function(traj.summary, T.go, T.stay) {
  # Merge stays that are broken up by very short "go" events; and merge "go" events that are split by very short "stays"
  # This can also be interpreted as filtering out too short go and stay events.
  
  # traj.summary: requires columns is.stay, durations, traj.segment, loc.id
  # time threshold: 
  
  
  
  # Find "go" events that: 
  # (1) are less than 5 minutes, AND
  # (2) start and end at same location
  
  #browser() 
  # Get all "go" segments below time threshold
  cond1 <- filter(traj.summary, !is.stay, durations <= T.go) %>%
    select(traj.event)
  
  # Make sure the short duration event is not first, last or only event:
  cond1 %<>% filter(traj.event > 1 & traj.event < max(traj.summary$traj.event))
  
  if(nrow(cond1)>0){
  # Get the loc.id's of stays before and after each of the short "go" events
  stay.before <-
    filter(traj.summary, traj.event %in% (cond1$traj.event - 1)) %>% select(loc.id)
  stay.after <-
    filter(traj.summary, traj.event %in% (cond1$traj.event + 1)) %>% select(loc.id)
  
  # Test whether the "go" starts and ends at the same location
  cond2 <- cond1[stay.before == stay.after, ]
  
  # Append the loc ID that the identified "go" points should be assigned to
  go.remove <- cbind(cond2, stay.before[stay.before == stay.after, ])
  } else {
    go.remove <- cond1
    }
  
  # Find "stay" events of duration below a time threshold and convert to "go":
  cond3 <- filter(traj.summary, is.stay, durations <= T.stay) %>%
    select(traj.event)
  if(nrow(cond3)>0){
  stay.remove <- cbind(cond3,loc.id=0)
  } else {
    stay.remove <- cond3
    }
  
  # Append the filtered stays to the filtered gos
  merge.temporal <- rbind(go.remove,stay.remove)
  
  return(merge.temporal)
  
}

# MINIMUM CONVEX POLYGON ----
get_mcp <- function(locations, quantile) {
  # calculates minimum convex polygon from GPS data
  #
  
  #browser() #for debugging
  
  # make sure the input looks as expected (set of lat and lon data)
  
  # calculate distances from all points to centroid of the location data
  centroid <- apply(locations,2,mean)
  distances <- sqrt(((locations[, 1] - centroid[1])^2) + ((locations[, 2] - centroid[2])^2))
  
  # get subset of points within specified quantile of distances
  indx <- 1:length(distances)
  percentages <- indx[distances <= quantile(distances, quantile)]
  locations.subset <- locations[percentages, ]
  
  # get minimum convex polygon
  mcp.points <- chull(locations.subset[, 1], locations.subset[, 2]) # index of points that lie on mcp
  mcp <- locations.subset[mcp.points,] # coords of mcp
  mcp <- rbind(mcp[nrow(mcp), ], mcp) # repeat last point to close shape
  
  mcp.poly<-Polygon(mcp) # creates a polygon object with area attribute (access uing @area)
  
  return(list(mcp=mcp, mcp.poly=mcp.poly, centroid=centroid))
}

# SUMMARISE TRAJECTORIES ----
summarise_trajectories <- function(gps.traj, dist.threshold) {
  
  # Required columns in gps.traj: dates, traj.event, timestamp, lat, lon, is.stay, homedist, durations
  
  # summarise by trajectory segment, with new columns for action range and "is.home" flag
  traj.summary <- gps.traj %>% group_by(dates,traj.event) %>%
    summarize(T.start = timestamp[1], T.end = timestamp[length(timestamp)], 
              is.stay = median(is.stay), loc.id = mean(loc.id),
              clat = ifelse(is.stay, mean(lat), NA), clon = ifelse(is.stay, mean(lon), NA),    # stays centroids
              action.range = ifelse(is.stay,mean(homedist,na.rm=TRUE),max(homedist,na.rm=TRUE))# distance from home (furthest on move, average in stay)
    ) %>%
    mutate(T.end = c(T.start[-1],T.end[length(T.end)])) %>%
    mutate(durations = difftime(T.end,T.start, units = "mins")) %>%
    mutate(is.home = ifelse(is.stay==1 & action.range<dist.threshold, TRUE, FALSE))
  
  # get displacements between stays
  tmp.displacements <- distGeo(ungroup(traj.summary) %>% filter(is.stay==1) %>% select(clon,clat),
                               a=6378137, f=1/298.257223563)
  traj.summary[,"displacements"] <- NA
  traj.summary[traj.summary$is.stay==1,"displacements"] <- c(0,tmp.displacements); remove(tmp.displacements)
  
  # set first stay displacement of each day to zero
  traj.summary %<>% group_by(dates,is.stay) %>% 
    mutate(displacements = ifelse(is.stay==1 & row_number()==1, 0, displacements)) %>% ungroup()
  
  return(traj.summary)
}

