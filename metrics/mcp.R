get_mcp <- function(lon,lat, quantile) {
  # calculates minimum convex polygon from GPS data
  #
  
   #browser() #for debugging
  locations <- data.frame(lon=lon,lat=lat)
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
  
  mcp.poly.coords <- data.frame(mcp.poly@coords[,c("lon","lat")], row.names = NULL) #polygon coordinates in lat/lon
  mcp.poly.area <- areaPolygon(mcp.poly.coords, a=6378137, f=1/298.257223563)/1000000 #polygon area in km2
  
  return(mcp.poly.area)
}
