# Mobility functions from Github

displacement<- function(df, coor = NULL, groupvar = NULL){
  if (is.atomic(df)) {
    df <- data.frame(x = df)
  }
  if (is.null(coor)) {
    stop("Geographic coordinates must be supplied.")
  }
  
  if (is.null(groupvar)) {
    df2<- df[,coor]
    suppressMessages(
      df2 <- slide(df2, Var = coor[1], slideBy = 1)
    )
    suppressMessages(
      df2 <- slide(df2, Var = coor[2], slideBy = 1)
    )
    message("groupvar is not specified, treating all the data as one group.")
  }
  
  if (!is.null(groupvar)) {
    df2<- cbind(df[,coor], df[groupvar])
    suppressMessages(
      df2 <- slide(df2, Var = coor[1], GroupVar = groupvar, slideBy = 1)
    )
    suppressMessages(
      df2 <- slide(df2, Var = coor[2], GroupVar = groupvar, slideBy = 1)
    )
  }
  
  coor2<- paste0(coor, "1")
  latlon<- df2[,coor]
  latlonlag2<- df2[,coor2]
  dist<- distVincentyEllipsoid(latlon, latlonlag2, a=6378137, b=6356752.3142, f=1/298.257223563)
  
  return(dist)
}


gridcentroid <- function(df, coor = NULL, cell.size = NULL,
                         cell.units = "meters", point.id = NULL){
  if (is.atomic(df)) {
    df <- data.frame(x = df)
  }
  if (is.null(coor)) {
    stop("Geographic coordinates must be supplied.")
  }
  if (is.null(cell.size)) {
    stop("A cell size (in meters) must be provided.")
  }
  if (!cell.units == "meters") {
    stop("Cell size must be in meters")
  }
  if (is.na(point.id)){
    df$point.id <- seq(1, nrow(df))
    point.id <- "point.id"
  }
  
  grid.points <- df[c(coor, point.id)]
  
  if (nrow(grid.points) == 1){
    staypoint.id <- 1
    staypointlon <- df[1, coor[1]]
    staypointlat <- df[1, coor[2]]
    df <- cbind(df, staypoint.id, staypointlon, staypointlat)
  } else {
    points <- SpatialPoints(grid.points)
  }
}

groupdist<- function(df, coor = NULL, threshold = NULL, groupvar = NULL) {
  if (is.atomic(df)) {
    df <- data.frame(x = df)
  }
  if (is.null(coor)) {
    stop("Geographic coordinates must be supplied.")
  }
  if (is.null(threshold)) {
    stop("A threshold must be specified.")
  }
  #browser()
  if (is.null(groupvar)) {
    group<- rep(1, nrow(df))
    df2<- cbind(df[,coor], group)
    suppressMessages(
      df2 <- slide(df2, Var = coor[1], slideBy = -1)
    )
    suppressMessages(
      df2 <- slide(df2, Var = coor[2], slideBy = -1)
    )
    suppressMessages(
      df2 <- slide(df2, Var = coor[1], slideBy = 1)
    )
    suppressMessages(
      df2 <- slide(df2, Var = coor[2], slideBy = 1)
    )
    message("groupvar is not specified, treating all the data as one group.")
  }
  
  if (!is.null(groupvar)) {
    group<- rep(1, nrow(df))
    df2<- cbind(df[,coor], df[groupvar], group)
    suppressMessages(
      df2 <- slide(df2, Var = coor[1], GroupVar = groupvar, slideBy = 1)
    )
    suppressMessages(
      df2 <- slide(df2, Var = coor[2], GroupVar = groupvar, slideBy = 1)
    )
  }
  
  #browser()
  coor2<- paste0(coor, "1")
  latlon<- df2[,coor]
  latlonlag2<- df2[,coor2]
  dist2<- distVincentyEllipsoid(latlon, latlonlag2, a=6378137, b=6356752.3142, f=1/298.257223563)
  df2<- cbind(df2, dist2)
  df2$disttestgroup<- ifelse(df2$dist2<= threshold, 1, 0)
  
  df2$distgroup<- seqgroup(df2, var = "disttestgroup")
  df2$distgroup<- ifelse(df2$disttestgroup==0, NA, df2$distgroup)
  
  suppressMessages(
    df2 <- slide(df2, Var = "distgroup", slideBy = -1)
  )
  df2$distgroup<- ifelse(is.na(df2$distgroup), df2$`distgroup-1`, df2$distgroup)
  
  return(df2$distgroup)
}


grouptime<- function(df, time = NULL, units = c("auto", "secs", "mins", "hours", "days", "weeks"),
                     threshold = NULL, groupvar = NULL) {
  if (is.atomic(df)) {
    df <- data.frame(x = df)
  }
  #browser
  if (!is.POSIXct(df[time][[1]])){
    stop("Time variable must be POSIXct format.")
  }
  
  if (is.null(threshold)){
    stop("A threshold must be supplied to generate time groups.")
  }
  
  if (is.null(groupvar)){
    timediff<- as.numeric(difftime(df[nrow(df),time], df[1,time], units = units))
    df$timegroup<- ifelse(timediff>=threshold, 1, 0)
    return(df$timediff)
  }
  #browser()
  if (!is.null(groupvar)){
    df3<- ddply(df, .(get(groupvar)), function(z){
      data.frame(timediff = as.numeric(difftime(z[nrow(z),time], z[1,time]), units = units))
    })
    df3$timegroup<- ifelse(df3$timediff>=threshold, 1, 0)
    names(df3)<- c(groupvar, "timediff","timegroup")
    df3<- subset(df3, !is.na(get(groupvar)))
    df<- mergewithorder(df, df3, by=groupvar)
    
    return(df$timegroup)
  }
}


mergewithorder<- function(df1, df2, by = NULL) {
  if (is.atomic(df1)) {
    df1 <- data.frame(x = df1)
  }
  if (is.atomic(df2)) {
    df2 <- data.frame(x = df2)
  }
  if (is.null(by)) {
    stop("Unique identifier must be specified to merge.")
  }
  
  df1$sequenceid<- seq(1:nrow(df1))
  df1<- merge(df1, df2, by = by, all = T)
  df1<- df1[order(df1$sequenceid),]
  df1$sequenceid<- NULL
  
  return(df1)
}


radiusofgyration<- function(df, coor = NULL, time = NULL, time.units = c("hour", "date", "month"), groupvar = NULL) {
  if (is.atomic(df)) {
    df <- data.frame(x = df)
  }
  if (is.null(coor)) {
    stop("Geographic coordinates must be supplied.")
  }
  if (is.null(time)) {
    stop("A time object must be supplied.")
  }
  if (!is.POSIXct(df[time][[1]])){
    stop("Time variable must be POSIXct format.")
  }
  if (length(time.units)!=1) {
    stop("time.units must be a single value.")
  }
  if (!time.units %in% c("hour", "date", "month")) {
    stop("time.units must be hour, date, or month. ")
  }
  
  if (is.null(groupvar)) {
    df2<- df[,c(coor, time)]
    f<- get(time.units)
    df2[,time.units]<- f(df2[,time])
    df2<- df2[order(df2[,time]),]
    df2$group<- seqgroup(df2, var = time.units)
    message("groupvar is not provided, using time.units as grouping varialbe.")
  } else if (!is.null(groupvar)) {
    df2<- df[,c(groupvar, coor, time)]
    f<- get(time.units)
    df2[,time.units]<- f(df2[,time])
    df2<- transform(df2, newid =
                      as.numeric(interaction(get(groupvar), get(time.units), drop=TRUE)))
    df2<- df2[order(df2[,groupvar],df2[,time]),]
    df2$group<- seqgroup(df2, var = "newid")
  }
  
  df3<- ddply(df2, .(group), function(z){
    data.frame(rg = sdspatialpoints(z, coor = c("lon","lat")))
  })
  df2<- merge(df2, df3, by="group")
  return(df2$rg)
}

sdspatialpoints<- function(df, coor = NULL, ...){
  if (is.atomic(df)) {
    df <- data.frame(x = df)
  }
  if (is.null(coor)) {
    stop("Geographic coordinates must be supplied.")
  }
  
  if (!is.null(coor)) {
    meancoor<- c(mean(df[,"lon"]),mean(df[,"lat"]))
    latlon<- cbind(df[,coor])
    latlon$dist<- distVincentyEllipsoid(latlon, meancoor, a=6378137, b=6356752.3142, f=1/298.257223563)
    latlon$dist2<- latlon$dist^2
    N<- nrow(latlon)
    distsum<- sum(latlon$dist2)
    sd<- sqrt(distsum/N)
    return(sd)
  }
}

seqgroup<- function(df, var = NULL) {
  if (is.atomic(df)) {
    df <- data.frame(x = df)
  }
  #browser()
  if (is.null(var)) {
    group<- seq(1:nrow(df))
  } else {
    group<- rep(1, nrow(df))
    
    df2<- cbind(df[var], group)
    suppressMessages(
      df2 <- slide(df2, Var = var, slideBy = -1)
    )
    var2<- paste0(var,"-1")
    
    df2[1,var2]= df2[1,var]
    df2$group<- ifelse(df2[var]==df2[var2],0, cumsum(df2$group)+1)
    is.na(df2$group)<- df2$group==0
    df2[1,"group"]<- 1
    df2$group<- na.locf(df2$group,na.rm=FALSE)
    df2$group<- as.numeric(df2$group)
    
  }
  return(df2$group)
}

stayevent<- function (df, coor = NULL, time = NULL, dist.threshold = NULL, time.threshold = NULL, 
                      time.units = c("auto", "secs", "mins", "hours", "days", "weeks"), 
                      groupvar = NULL, ...) {
  df$distgroup<- groupdist(df, coor = coor, threshold = dist.threshold, groupvar = groupvar)
  df$timegroup<- grouptime(df, time = time, units = time.units, threshold = time.threshold, groupvar = "distgroup")
  df$stayeventgroup<- ifelse(!is.na(df$distgroup) & df$timegroup == 1, df$distgroup, NA)
  #browser()
  stayevents<- aggregate(cbind(get(coor[1]), get(coor[2]))~stayeventgroup, df, mean)
  names(stayevents)<- c("stayeventgroup", "stayeventlon", "stayeventlat")
  
  df<- mergewithorder(df, stayevents, by="stayeventgroup")
  df$distgroup<- NULL
  df$timegroup<- NULL
  
  return(df)
}

