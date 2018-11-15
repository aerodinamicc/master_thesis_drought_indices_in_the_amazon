# CWD with the fixed value of 100 mm month PET ----
cwd <- function(x, clim_year_start = 1, et = 100) {
  x <- x %>% 
    arrange(year, month) %>% 
    # compute precipitation of previous month
    mutate(prev_prec = lag(monthly_prec, default = 0))
  
  # restrict to complete year for easier blocking
  if (x$month[1] != 1) {
    n <- nrow(x)
    x <- x[which(x$month == 1)[1]:n, ]
  }
  if (tail(x$month, 1) != 12) {
    n <- nrow(x)
    x <- x[1:tail(which(x$month == 12), 1), ]
  }
  
  # create potentially altered year cycle for climatological year
  clim_year_months <- c(1:12)
  if (clim_year_start != 1) {
    clim_year_months <- clim_year_months[c(clim_year_start:12, 1:(clim_year_start - 1))]  
  }
  
  # loop in month-wise blocks
  for (i in 1:12) {
    .month <- clim_year_months[i]
    .prev_month <- clim_year_months[i - 1]
    if (.month == 1) {
      x$wd[x$month == .month] <- 0
    } else {
      wd0 <- x$wd[x$month == .prev_month] - et + x$monthly_prec[x$month == .month]
      x$wd[x$month == .month] <- ifelse(wd0 >= 0, 0, wd0)
    }
  }
  x
}


# Tidy TRMM data ----
ncdf_to_tidy <- function(nc, .var, .start, .time = "time",
                         .res = "1 day", .station = "station",
                         .lon = "lon", .lat = "lat") {
  time <- ncvar_get(nc, .time)
  time <- seq(as.Date(.start), by = .res, length.out = length(time))
  year <- year(time) %>% as.integer
  month <- month(time) %>% as.integer
  day <- day(time)
  
  var_tidy <- ncvar_get(nc, .var) %>%
    as_tibble %>%
    mutate(year = year,
           month = month,
           day = day) %>%
    gather(station, .var, starts_with("V"))
  
  landid <- tibble(station = ncvar_get(nc, .station),
                   lon = ncvar_get(nc, .lon),
                   lat = ncvar_get(nc, .lat)) %>%
    mutate(station = paste0("V", station + 1))
  
  left_join(var_tidy, landid, by = "station") %>%
    dplyr::select(lon, lat, year, month, day, .var)
}

# Refine resolution for datasets that have a more coarse one ----

refineResolution <- function(data, .var = ""){
  #a move of 0.5 degree is apllied to any of the following directions so that the each cell's resolution is refined
  east <- data
  north <- data
  northEast <- data
  
  east$lon <- east$lon + 0.5
  north$lat <- north$lat + 0.5
  northEast$lon <- east$lon
  northEast$lat <- north$lat
  
  refinedData <- data %>% 
    rbind(east) %>%
    rbind(north) %>%
    rbind(northEast)
  refinedData <- refinedData %>% arrange(lon, lat, year, month)
  return(refinedData)
}

# Crop a tibble in order to contain only Amazon relevant pixels ----
cropTibbleToAmazonExtent <- function(tableToBeCropped){
  
  #We use the TRMM extent as guidence
  lonTrmm <- ncvar_get(trmm, "lon")
  latTrmm <- ncvar_get(trmm, "lat")
  
  lonLatFromTRMM <- paste(as.character(lonTrmm), " ", as.character(latTrmm))
  
  tableToBeCropped <- tableToBeCropped %>% 
    mutate(coor = paste(as.character(lon), " ", as.character(lat))) %>%
    group_by(coor) %>%
    mutate(cellStays =  coor %in% lonLatFromTRMM) %>%
    filter(cellStays == TRUE) %>%
    ungroup() %>%
    mutate(stays = (!(lon == -76 & lat == -12.0))) %>%
    filter(stays == TRUE) %>%
    dplyr::select(-cellStays, -stays, -coor) %>%
    arrange(lon, lat, year, month)
  
  #We exclude the -76, -12 pixel because of its NA values
  
  return(tableToBeCropped)
}

# Reading in various relevant for the analysis ncdf files (PET, AET, SPEI data). A rectangle following Amazon shape is extracted initially ----
readFromNcFile <- function(ncfile, .var="", .lon = "lon", .lat = "lat", moreDims = FALSE){
  time <- ncvar_get(ncfile, "time")
  lon <- ncvar_get(ncfile, .lon)
  lat <- ncvar_get(ncfile, .lat)
  date <- 0
  #used for spei only
  correctionInDegrees <- 0
  
  #the values are the max min for both longitute and laditude - respectively
  lons <- seq(-79.5, -50.5, by= 1)
  lats <- seq(-20.5, 4.5, by = 1)
  
  if (.var == "ET_mean")
  {
    date <- as.Date(paste(as.character(time),as.character(01)), format = "%Y%m%d")
  }
  if(.var %in% c("spei", "pre", "pet")) #pre and pet refer to the CRU data
  {
    time0 <- ymd("1900-01-01")
    date <- time0 + days(time)
    correctionInDegrees <- 0.25
    lons <- seq(-79.5, -50.5, by= 0.5)
    lats <- seq(-20.5, 4.5, by = 0.5)
  }
  if (.var == "pe")
  {
    time0 <- ymd("1948-01-01")
    date <- time0 + months(time)
    lon[181:360] <- seq(from = -179.5, to = -0.5, by = 1)
  }
  
  grid <- expand.grid(lon = lons, lat = lats)
  
  #Getting the data from the NetCDF-------------------------------------------------------
  extractedData <- tibble()
  
  for (i in 1:nrow(grid)){
    amazon_lon_index <- which(lon == grid$lon[i] - correctionInDegrees)
    amazon_lat_index <- which(lat == grid$lat[i] - correctionInDegrees)
    
    #PE has one more dimension which needs to be taken care of
    startVector <- c(amazon_lon_index, amazon_lat_index, 1)
    countVector <- c(1, 1, -1)
    
    if(moreDims)
    {
      startVector <- c(amazon_lon_index, amazon_lat_index, 1, 1)
      countVector <- c(1, 1, 1, -1)
    }
    
    value <- ncvar_get(ncfile, .var,
                       start = startVector,
                       count = countVector) ####gives NA values
    singleRow <- tibble(value, date) %>% mutate(
      year = year(date),
      month = month(date),
      lon = grid$lon[i],
      lat = grid$lat[i]
    )
    extractedData <- extractedData %>% rbind(singleRow)
  }
  
  return(extractedData)
}

# compute Cwd with input from various metrices ----
computeCwd <- function(x, cwdList) {
  x <- x %>%
    group_by(lon, lat) %>%
    mutate(record_number = row_number(lon)) %>%
    ungroup() %>%
    arrange(lon, lat, year, month)
  
  for (cwdId in 1:length(cwdList)) {
    cwd <- cwdList[cwdId]
    cwdElements <- strsplit(cwd, "_") # e.g. "cwd" "preCru" "petCruAvg" "r"
    if (cwdElements[[1]][2] == "aet") {
      petVar <- cwdElements[[1]][3] # variation of the pet data - Princeton/CRU
      for (month in 1:12) {
        relevantRows <- x$month == month
        if(month == 1){
          wd <- x$aet[relevantRows] - x[relevantRows, petVar]
          wdNegative <- sapply(wd, function (value) ifelse(value >= 0, 0, value))
          x[relevantRows, cwd] <- wdNegative
        }
        else{
          prev_cwd_value <- x[x$month == month - 1, cwd]
          wd <- prev_cwd_value + x$aet[relevantRows] - x[relevantRows, petVar]
          wdNegative <- sapply(wd, function (value) ifelse(value >= 0, 0, value))
          x[relevantRows, cwd] <- wdNegative
        }
      }
    }
    else {
      prec <- cwdElements[[1]][2]
      etMetric <- cwdElements[[1]][3]
      isReset <- cwdElements[[1]][4] == "r"
      loops<- 96 #in case no resetting is applied; 8 * 12 because of the 8 year period
      
      if (isReset) {
        loops <- 12
      }
      
      for (position in 1:loops) {
        relevantRows <- 0
        
        if (isReset) {
          relevantRows <- x$month == position
        }
        else {
          relevantRows <- x$record_number == position
        }
        
        if (position == 1) {
          wd <- x[relevantRows, prec]  - x[relevantRows, etMetric]
          wdNegative <- sapply(wd, function (value) ifelse(value >= 0, 0, value))
          x[relevantRows, cwd] <- wdNegative
        }
        else {
          prev_cwd_value <- 0
          if (isReset) {
            prev_cwd_value <- x[x$month == position - 1, cwd]
          }
          else{
            prev_cwd_value <- x[x$record_number == position - 1, cwd]
          }
          
          wd <- prev_cwd_value - x[relevantRows, etMetric] + x[relevantRows, prec]
          wdNegative <- sapply(wd, function (value) ifelse(value >= 0, 0, value))
          x[relevantRows, cwd] <- wdNegative
        }
      }
    }
  }
  return(x %>% dplyr::select(-record_number, -et100))
}

# compute MCWD (based on cwd) and anomalies ----
computeMcwdAndAnomalies <- function(data, anoInScope = tibble(), anoSdInScope = tibble(), var = "")
{
  data %>% 
    group_by(lon, lat, year) %>%
    summarise(mcwd = round(min(!!as.name(var)))) -> mcwd
  
  mcwd %>% 
    group_by(lon, lat) %>%
    summarise(mcwdMean = round(mean(mcwd)), mcwdSd = round(sd(mcwd))) -> mcwd_avg
  
  mcwd%>% 
    left_join(mcwd_avg, by = c("lon", "lat")) %>% 
    mutate(mcwd_diff = mcwd - mcwdMean, mcwd_diff_sd = ifelse(mcwdSd == 0, 0, (mcwd - mcwdMean)/mcwdSd)) -> anomaly
  
  if(ncol(anoInScope) < 1)
  {
    anoInScope <- anomaly %>%
      ungroup %>%
      select(year, lon, lat, mcwd_diff)
  }
  else
  {
    anoInScope <- anoInScope %>% add_column(anomaly$mcwd_diff)
  }
  
  if(ncol(anoSdInScope) < 1)
  {
    anoSdInScope <- anomaly %>%
      ungroup %>%
      select(year, lon, lat, mcwd_diff_sd)
  }
  else
  {
    anoSdInScope <- anoSdInScope %>% add_column(anomaly$mcwd_diff_sd)
  }
  
  colnames(anoInScope)[ncol(anoInScope)] <- substring(var, 5, nchar(var))
  colnames(anoSdInScope)[ncol(anoSdInScope)] <- substring(var, 5, nchar(var))
  
  data %>% 
    group_by(lon, lat, month) %>%
    summarise(mcwdMeanByMonth = mean((!!as.name(var))), mcwdSdByMonth = sd((!!as.name(var)))) -> mcwd_monthly_avg
  #return(mcwd_monthly_avg)
  
  colName <- paste(var, "_monthlySd", sep = "")
  data <- data %>% 
    left_join(mcwd_monthly_avg, by = c("month", "lon", "lat")) %>%
    mutate(!!colName := ifelse(mcwdSdByMonth == 0, 0, ((!!as.name(var)) - mcwdMeanByMonth)/mcwdSdByMonth)) %>%
    dplyr::select(-mcwdSdByMonth, -mcwdMeanByMonth)
  
  returnedList <- list("AmazonWD" = data, "anomalies" = anoInScope, "anomaliesSd" = anoSdInScope)
  return(returnedList)
}
