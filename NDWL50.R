## Number of waterlogging days (NDWL50)
## By: H. Achicanoy
## December, 2022

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))
source('https://raw.githubusercontent.com/CIAT-DAPA/agro-clim-indices/main/_main_functions.R')

root <- 'F:/Data/'

ref <- terra::rast(paste0(root,'roi/africa.tif'))

# Soil variables
scp <- terra::rast(paste0(root,'soils/africa_scp_new.tif'))
scp <- scp %>% terra::resample(ref) %>% terra::mask(ref)
sst <- terra::rast(paste0(root,'soils/africa_ssat.tif'))
sst <- sst %>% terra::resample(ref) %>% terra::mask(ref)

# Calculate NDWL function
calc_NDWL50 <- function(yr, mn){
  outfile <- paste0(out_dir,'/NDWL50-',yr,'-',mn,'.tif')
  cat(outfile, "\n")
  if(!file.exists(outfile)){
    dir.create(dirname(outfile),F,T)
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    # Sequence of dates
    if(as.numeric(yr) > 2020 & mn == '02'){
      dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-28')), by = 'day')
    } else {
      dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    }
    # Files
    pr_fls <- paste0(pr_pth,'chirps-v2.0.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    pr_fls <- pr_fls[file.exists(pr_fls)]
    tx_fls <- paste0(tx_pth,'/Tmax.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    tx_fls <- tx_fls[file.exists(tx_fls)]
    tm_fls <- paste0(tm_pth,'/Tmin.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    tm_fls <- tm_fls[file.exists(tm_fls)]
    if(as.numeric(yr) > 2020){
      dts_chr <- as.character(dts)
      dts_chr <- gsub(pattern = yr, replacement = yrs_mpg$Baseline[yrs_mpg$Future == yr], x = dts_chr)
      dts <- as.Date(dts_chr); rm(dts_chr)
    }
    sr_fls <- paste0(sr_pth,'Solar-Radiation-Flux_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.1.nc')
    sr_fls <- sr_fls[file.exists(sr_fls)]
    # Read variables
    prc <- terra::rast(pr_fls)
    prc <- prc %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
    prc[prc == -9999] <- NA
    tmx <- terra::rast(tx_fls)
    tmx <- tmx %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
    tmx[tmx == -9999] <- NA
    tmn <- terra::rast(tm_fls)
    tmn <- tmn %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
    tmn[tmn == -9999] <- NA
    tav <- (tmx + tmn)/2
    srd <- terra::rast(sr_fls)
    srd <- srd %>% terra::crop(terra::ext(ref))
    #The line `srd <- srd/1000000` scales the solar radiation data by dividing each value by \(10^6\). This step is likely needed to convert the units of solar radiation.
    srd <- srd/1000000
    srd <- srd %>% terra::resample(x = ., y = ref) %>% terra::mask(ref)
    
    # Maximum evapotranspiration
    ETMAX <- terra::lapp(x = terra::sds(srd,tmn,tav,tmx), fun = peest)
    rm(list=c("tmn", "tmx", "tav", "srd"))
    gc(verbose=F, full=T, reset=T)
    
    # Compute water balance model
    date <- paste0(yr,'-',mn)
    if(date %in% c('1983-01','2021-01','2041-01')){
      AVAIL <<- ref
      AVAIL[!is.na(AVAIL)] <- 0
    } else {
      AVAIL <<- terra::rast(paste0(dirname(outfile),'/AVAIL.tif'))
    }
    
    eabyep_calc <- function(soilcp = scp, soilsat = ssat, avail = AVAIL,rain = prc[[1]], evap = ETMAX[[1]]){
      
      avail   <- min(avail, soilcp)
      
      # ERATIO
      percwt <- min(avail/soilcp*100, 100)
      percwt <- max(percwt, 1)
      eratio <- min(percwt/(97-3.868*sqrt(soilcp)), 1)
      
      demand  <- eratio * evap
      result  <- avail + rain - demand
      logging <- result - soilcp
      logging <- max(logging, 0)
      logging <- min(logging, soilsat)
      # runoff  <- result - logging + soilcp
      avail   <- min(soilcp, result)
      avail   <- max(avail, 0)
      # runoff  <- max(runoff, 0)
      
      out     <- list(Availability = c(AVAIL, avail),
                      # Demand       = demand,
                      # Eratio       = eratio,
                      Logging      = logging
      )
      return(out)
    }
    
    watbal <- 1:terra::nlyr(ETMAX) %>%
      purrr::map(.f = function(i){
        water_balance <- eabyep_calc(soilcp  = scp,
                                     soilsat = sst,
                                     avail   = AVAIL[[terra::nlyr(AVAIL)]],
                                     rain    = prc[[i]],
                                     evap    = ETMAX[[i]])
        AVAIL <<- water_balance$Availability
        return(water_balance)
      })
    LOGGING <- watbal %>% purrr::map('Logging') %>% terra::rast()
    # Calculate number of soil waterlogging days (if logging is above 0)
    # Note NDWL50 uses ssat * 0.5 (so this means soil is at 50% toward saturation)
    # here it suffices if the soil is above field capacity
    NDWL50  <- sum(LOGGING > (sst*0.5))
    terra::writeRaster(NDWL50, outfile)
    terra::writeRaster(AVAIL[[terra::nlyr(AVAIL)]], paste0(dirname(outfile),'/AVAIL.tif'), overwrite = T)
    
    #clean up
    rm(list=c("prc", "ETMAX", "AVAIL", "watbal", "ERATIO", "LOGGING", "NDWL50"))
    gc(verbose=F, full=T, reset=T)
  }
}

# # Historical setup
yrs <- 1983:2016
#yrs <- 2004:2016 
mns <- c(paste0('0',1:9),10:12)
#yrs <- 2012:2012
#mns <- c(paste0(10:12))


stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
names(stp) <- c('yrs','mns')
stp <- stp %>%
  dplyr::arrange(yrs, mns) %>%
  base::as.data.frame()
pr_pth <- 'F:/Data/' # Precipitation
tm_pth <- paste0(root,'Tmin') # Minimum temperature
tx_pth <- paste0(root,'Tmax') # Maximum temperature
sr_pth <- paste0(root,'solar_radiation_flux/solar_radiation_flux/') # Solar radiation
out_dir <- paste0(root,'NDWL50')
1:nrow(stp) %>%
  purrr::map(.f = function(i){
    calc_NDWL50(yr = stp$yrs[i], mn = stp$mns[i]); gc(verbose=F, full=T, reset=T)
    if (i%%5 == 0) {
      tmpfls <- list.files(tempdir(), full.names=TRUE)
      1:length(tmpfls) %>% purrr::map(.f = function(k) {system(paste0("rm -f ", tmpfls[k]))})
    }
  })


#}