findhome <- function(df,lat,lon){
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

