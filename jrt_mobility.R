



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
  home <- as.numeric(select(loc.counts[1,],c(lat4,lon4)))
  
  return(home)
}


# STAY DETECTION: DISTANCE & TIME ----
get_stays <- function(df, deltaT, deltaD){
  # Detect stay points by measuring distance across the area covered in time windows of deltaT minutes
  # df:
  # deltaT:
  # deltaD:
  
  #TODO: create logic for if the gps.log dataframe has other names for timestamp, lat and lon timestamp=NULL,lat=NULL,lon=NULL,
  
  is.stay <- vector(mode = "logical", length = nrow(gps.log)) 
  i <- 1
  while (i <= nrow(df)){
    
    # get all GPS data points in a timeframe of deltaT
    points <- df %>% 
      filter(timestamp>=timestamp[i], timestamp<=timestamp[i]+deltaT*60) %>%
      select(lat, lon)
    
    # diagonal distance across bounding box enclosing the selected points
    bob <- data.frame(bbox(SpatialPoints(points)))
    ddistance <- distGeo(bob$min, bob$max, a=6378137, f=1/298.257223563) #diagonal distance across bounding box
    
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
  gps.traj <- cbind(gps.log,is.stay) %>%
    mutate(event.id = cumsum(c(0,abs(diff(is.stay))))+1) # assign a number to each stay or go event
  
  
  return(gps.traj)
  
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

