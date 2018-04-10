
# Not implemented in the restructure function:
# - threshold and flags column
# - select only relevant columns

get.data <- function(folder,NIU) {
  # imports project data stored as a set of .csv files and creates a list of 
  # dataframes for use in further analysis; note that this does not take user or
  # date into account but imports ALL data
  # inputs:
  #  - folder: path to folder containing data files (character string)
  #  - NIU: vector of filenames that are not to be imported and used (each given as a character string)
  # output:
  # - list of dataframes each with data from a single csv file in "folder"  
  
  #browser() #for debugging
  
  out <- list() #initialise list to store all dataframes

  # Get a list of all .csv files (and their root names) in folder and exclude those that are not of interest
  file_list <- list.files(path=folder, pattern="*.csv")
  file_names <- tstrsplit(file_list,".", fixed = TRUE, keep=1)[[1]]
  if(!missing(NIU)) {
    file_names <- file_names [! file_names %in% NIU]
  }

  for (file in file_names){
    # import each datafile and add to the list for output
    out[[file]]<-read.csv(paste(folder, file,".csv", sep=''),header=TRUE, sep=",", stringsAsFactors=FALSE)
    cat(file,"file imported\n");
  }
  return(out)
}

restructure <- function(dat, userid, d.start, d.stop) {
  # Subsets a dataframe containing all data from a particular source (eg device
  # sensor or log) to a specified user and timeframe then restructures the 
  # columns for further processing/anaysis.
  # inputs:
  #  - dat: dataframe containing all data from a particular source (eg location or activity), including all users over all time
  #  - userid: participant code as a character string
  #  - d.start: start date of period of interest (yyyy-mm-dd)
  #  - d.stop: end date of period of interest (yyyy-mm-dd)
  # output:
  # - dat: restructured dataframe
  
  # TODO: if d.start and d.stop are same date, then there is no data in dat once filtered and the code that follows generates errors
  
  #browser() # for debugging
  cat("\n *** \n Processing file containing: ",colnames(dat)[-(1:7)], "\n\n");
  
  # Quality Check: make sure there are timestamp and user columns in the dataset:
  stopifnot("timestamp" %in% colnames(dat))
  stopifnot("user" %in% colnames(dat))
  
  # Get data for specified participant only
  if(userid %in% dat$user){
    dat %<>% filter(user==userid)
    dat$user <- NULL # removes "user" column since this is not of interest when same for entire dataset
    cat("Complete: Data extracted for specified user \n");
    } 
  else{
    cat("User ID is not in dataset, returning NULL \n");
    return(NULL)
    }

  # Extract only the period of interest:
  if(d.start <= max(dat$timestamp) & d.stop >= min(dat$timestamp)){ #checks whether the period of interest overlaps with the timeframe in the dataset
    dat %<>% 
      mutate(timestamp=as.POSIXct(timestamp,format="%Y-%m-%d %H:%M:%S")) %>% #NOTE: this just puts it into CEST format which is not correct, need to fix
      filter(timestamp>=d.start, timestamp<=d.stop) %>%
      arrange(timestamp)
    cat("Complete: Data extracted for period of interest\n");
  }
  else{
    cat("No data for the specified user within the period of interest, returning NULL\n",
        "The specified user has data recorded between ", min(dat$timestamp), "and", max(dat$timestamp), "for this dataset.\n");
    return(NULL)
  }
  

  # Add new columns that are useful for the plotting and algorithms to come:
  # to do: include inputs to the function to indicate which new columns are desired
  dat %<>% 
    mutate(intervals = c(NA,(diff(timestamp)/(60*60)))) %>% # intervals between readings in hours
    mutate(dates = as.character.Date(as.Date(dat$timestamp))) %>%  # dates
    mutate(times = as.POSIXct(strftime(dat$timestamp, format="%H:%M:%S"), format="%H:%M:%S")) %>%  # add times for plotting time on y axis (note: to get the times in the right format, they are all just assigned the same date)
    mutate(index = c(1:nrow(dat))) # index for testing purposes

  # If dat has "funf_version" column, use it create a new column indicating the data source
  if("funf_version" %in% colnames(dat)){
    tmp<-as.data.frame(tstrsplit(dat$funf_version,"-")) # "funf_version" ends in "fitness" for the phone and "wearable" for the watch
    dat$dsource<-as.character(tmp[,2]) # data source column
    dat <- dat %>% mutate(dsource = replace(dsource, dsource=="fitness","phone"))
    dat <- dat %>% mutate(dsource = replace(dsource, dsource=="wearable","watch"))
    remove(tmp)
    #dat$funf_version <- NULL
    cat("Complete: funf_version converted to dsource column indicating whether reading is from phone or watch\n");
    }
  else {
    cat("No funf_Version column in the dataset, therefore no dsource column created to indicate whether a reading is from the phone or watch\n");
  }
  #browser()
  if("activity" %in% colnames(dat)){
    #Description of activities:
    act.list<-list("0"="Vehicle",  # The device is in a vehicle, such as a car.
                   "1"="Bicycle",  # The device is on a bicycle.
                   "2"="Foot",     # The device is on a user who is walking or running.
                   "3"="Still",    # The device is still (not moving).
                   "4"="Unknown",  # Unable to detect the current activity.
                   "5"="Tilting",  # The device angle relative to gravity changed significantly.
                   "7"="Walking",  # The device is on a user who is walking.
                   "8"="Running")  # The device is on a user who is running.
    # create a new column to describe the activities
    dat$label<-unlist(act.list[as.character(dat$activity)])
    cat("Complete: column 'label' created describing activity based on activity code \n");
    
    # Keep only the entries with max confidence for each timestamp (can still result in multiple activities per timestamp if several at same time have same (max) confidence)
    before <- nrow(dat)
    dat %<>%
      group_by(timestamp) %>%
      filter(confidence == max(confidence)) %>%
      ungroup() %>%
      as.data.frame()
    after <- nrow(dat)
    cat("Complete: activity data reduced from ", before, " to ", after, " rows by keeping only the max-confidence activities for each unique timestamp \n");
  }
  
  if("plugged" %in% colnames(dat)){
    #Description of charging modes:
    bat.list<-list("0"="unplugged",
                   "1"="plugged",
                   "2"="USB plugged")
    # create a new column to describe the modes
    dat$mode<-unlist(bat.list[as.character(dat$plugged)])
    cat("Complete: column 'mode' created describing charging status based on 'plugged' variable \n");
    }

  # "flags" column to show where there are duplicates or long gaps in the data:
  # threshold <-12 # in hours, to set how long warrants a "gap" flag
  # dat <- dat %>% mutate(flag = as.factor(c("normal","duplicate","gap")[1+(intervals==0)*1+(intervals>threshold)*2]))
  
  return(dat)
}


show_plots <- function(datasets, to_plot) {
  # datasets is a list of dataframes; to_plot is a vector of the dataframes to be plotted
  custom.colors <- c("red","orange","gold","green","dodgerblue","purple","pink") #https://www.w3.org/TR/css-color-3/#svg-color
  
   # Activity:
  if("activity" %in% to_plot & "activity" %in% names(datasets)){
    activity <- (datasets.all$activity)
    
    # Aggregates by day in a stacked bar indicating activity type with colour
    activity.plot1 <- plot_ly(as.data.frame(activity %>% count(label,dates)), 
            x = ~dates %>% format.Date(format="%d/%m"), 
            y = ~n, 
            type = 'bar', color = ~label) %>%
      layout(yaxis = list(title = 'Count'),
             xaxis = list(title = 'Date'),
             margin = list(b = 100, t=100),
             barmode = 'stack',
             title = "Total activity readings per day by type")
    
    # Activities over the day
    activity.plot2 <- plot_ly(activity,
            x=~dates %>% format.Date(format="%d/%m"),
            y=~times %>% format.Date(format="%H:%M"),
            type="scatter",
            mode="markers",
            color=~label) %>%
      layout(yaxis = list(title = 'Time of day'),
             xaxis = list(title = 'Date'),
             margin = list(b = 100, t=100),
             title = "Activities over time for each day")
  }
  else{
    activity.plot1<-NULL
    activity.plot2<-NULL
  }
  
  # Battery:
  if("battery" %in% to_plot & "battery" %in% names(datasets)){
    battery<-datasets.all$battery
    battery.plot1 <- plot_ly(battery, 
            x=~timestamp, 
            y=~level, 
            type="scatter", 
            mode="markers", 
            color=~mode) %>%
      layout(margin = list(b = 100, t=100),
             xaxis = list(title = 'Time'),
             title = "Battery level and charging patterns")
  }
  else{
    battery.plot1<-NULL
  }
  
  # Experience Sampling:
  if("exps" %in% to_plot & "experience_sampling" %in% names(datasets)){
    exps<-datasets.all$experience_sampling
    exps.plot1 <- plot_ly(exps, 
             x=~dates %>% format.Date(format="%d/%m"),
             y=~times %>% format.Date(format="%H:%M"), 
             type="scatter", 
             #color=~question_id,
             mode="markers",
             symbol = ~question_id,
             symbols=1:3,
             marker = list(size = 10)) %>%
      layout(yaxis = list(title = 'Time of day'),
             xaxis = list(title = 'Date'),
             margin = list(b = 100, t=100),
             title = "Experience sampling: logged responses")
    
    # Currently unused (need a better way of visualising the answers)
    # exps.plot2 <- plot_ly(exps, 
    #                       x=~timestamp, 
    #                       y=~answer, 
    #                       type="scatter", 
    #                       color=~question_id,
    #                       mode="lines+markers",
    #                       marker = list(size = 12)) %>%
    #   layout(margin = list(b = 100, t=100),
    #          xaxis = list(title = 'Time'),
    #          title = "Answers to experience sampling questions")
    exps.plot2<-NULL
    
  }
  else{
    exps.plot1<-NULL
    exps.plot2<-NULL
  }
  
  # # Trying to get the legend to be ordered by date instead of weekday
  # plot_ly(location%>%group_by(dates),x=~lon,y=~lat,
  #         type = "scatter", mode = 'lines+markers',
  #         symbol = ~dates,
  #         symbols = shapes,
  #         #sort = FALSE,
  #         legendgroup=~dates,
  #         color = ~weekdays(as.Date(dates)))%>%
  # layout(legend=list(traceorder="grouped"))
  
  # # example to test:
  # # this works for showing individual traces by date, and indicating colour by weekday, but needs work selecting more appropriate colours
  # dummy<-data.frame(x=c(1:11),
  #                  y=runif(11),
  #                  dats = distinct(location, dates))
  # dummy$wdays <- weekdays(as.Date(dummy$dates))
  # dummy$colrs <- rep(cust.colors,length.out=11)
  # #dummy$colrs <- rep(1:7,length.out=11)
  # 
  # plot_ly(dummy,x=~x,y=~y,
  #         type = "scatter", mode = 'lines+markers',
  #         color=~dates,
  #         colors = ~colrs)
  
  # Location:
  if("location" %in% to_plot & "location" %in% names(datasets)){
    # Need ways to get rid of error points (eg some point very far away within same timeframe as others)
    location<-datasets.all$location
    location.watch<-filter(datasets.all$location, dsource == "watch")
    shapes<-rep(1:6,length.out=nrow(distinct(location, dates))) # creates a vector to assign 6 different symbols to the dates (repeating symbols)
    
    # Daily traces on blank lat/lon axes indicating day in colour:
    # plotly example: https://plot.ly/~ben_derv/0/gps/#code
    location.plot1 <- plot_ly(location,x=~lon,y=~lat,
            type = "scatter", mode = 'lines+markers',
            symbol = ~dates,
            symbols = shapes,
            color = ~weekdays(as.Date(dates))
            #colors = custom.colors
            ) %>%
      layout(margin = list(b = 100, t=100),
             title = "Location traces by day from watch and phone combined")
    location.plot2 <-  plot_ly(location.watch,x=~lon,y=~lat,
            type = "scatter", mode = 'lines+markers',
            color = ~dates,
            colors = ~custom.colors[1:length(unique(dates))])  %>%
      layout(margin = list(b = 100, t=100),
             title = "Location traces by day from watch only")
    # Look at watch vs phone data:
    location.plot3 <- plot_ly(location,x=~lon,y=~lat,
                              type = "scatter", mode = 'lines+markers',
                              color = ~dsource,
                              colors = c("cadetblue","purple")) %>%
      layout(margin = list(b = 100, t=100),
             title = "Location data from watch and phone")
  }
  else{
    location.plot1<-NULL
    location.plot2<-NULL
    location.plot3<-NULL
  }
  
  # Screen
  if("screen" %in% to_plot & "screen" %in% names(datasets)){
    screen<-datasets.all$screen
    screen.plot1 <- plot_ly(screen,
            x=~dates %>% format.Date(format="%d/%m"),
            y=~times %>% format.Date(format="%H:%M"), 
            type="scatter", 
            mode="markers",
            symbol = ~screen_on,
            symbols= c(1,4)) %>%
    layout(yaxis = list(title = 'Time of day'),
           xaxis = list(title = 'Date'),
           margin = list(b = 100, t=100),
           title = "Phone screen on/off transitions")
  }
  else{
    screen.plot1 <- NULL
  }
  
  # Steps
  if("steps" %in% to_plot & "step_count" %in% names(datasets)){
    steps<-datasets.all$step_count
    steps.watch<-filter(datasets.all$step_count, dsource == "watch")
    # to add: steps over the day
    
    
    # Total steps per day
    steps.plot1<- plot_ly(steps, 
            x=~timestamp, 
            y=~step_count, 
            type="scatter", 
            mode = "lines+markers",
            color = ~dsource,
            colors = c("cadetblue","purple")) %>%
      layout(yaxis = list(title = 'Step count'),
            xaxis = list(title = 'Date'),
            margin = list(b = 100, t=100),
            title = "Recorded stepcount for watch and phone")
  }
  else{
    steps.plot1 <- NULL
  }
  
 results<-list(activity.plot1,activity.plot2,
               battery.plot1,
               exps.plot1,exps.plot2,
               location.plot1,location.plot2,location.plot3,
               screen.plot1,
               steps.plot1)
  
  return(Filter(Negate(is.null), results))
}



# Activity data cleanup for plotting and analysis

# make it so that there is only one activity per timestamp (ie no duplicate timestamps)
# TO DO:
# check if "foot" is always with either walking or running then remove
# prioritise then remove duplicate timestamps (ie if all have same confidence then don't use "unknown" or "tilting")

grouptime<- function(df, time = NULL, units = c("auto", "secs", "mins", "hours", "days", "weeks"),
                     threshold = NULL, groupvar = NULL) {
  if (is.atomic(df)) {
    df <- data.frame(x = df)
  }
  
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


qualitychecks<- function(df){
  # df containing dsource column indicating device from which data was recorded, and timestamps
  #
  
  # check that the required columns exist:
  
  
  # device source:
  cat("Data in ", deparse(substitute(df)), "comes from ", levels(distinct(df,dsource)[,1]),"\n");
  
  # last reading:
  cat("Last reading in", deparse(substitute(df)), "is on", capture.output(max(j$timestamp)), "\n");
  
  
}
