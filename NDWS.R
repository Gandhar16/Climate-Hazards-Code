## Number of soil water stress days (NDWS)
## By: H. Achicanoy
## December, 2022

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))
source('https://raw.githubusercontent.com/CIAT-DAPA/agro-clim-indices/main/_main_functions.R')

root <- 'E:/Data/'

ref <- terra::rast(paste0(root, 'roi/africa.tif'))

align_rasters <- function(raster, reference) {
  # Check if extents match
  if (!isTRUE(all.equal(terra::ext(raster), terra::ext(reference)))) {
    # Adjust resolution and extent to match the reference
    raster <- terra::project(raster, crs(reference)) # Ensure CRS match
    raster <- terra::resample(raster, reference)    # Adjust to match resolution and extent
  }
  return(raster)
}


# Soil variables
scp <- terra::rast(paste0(root, 'soils/africa_scp_new.tif'))
scp <- align_rasters(scp, ref)
scp <- terra::mask(scp, ref)
sst <- terra::rast(paste0(root, 'soils/africa_ssat.tif'))
sst <- align_rasters(sst, ref)
sst <- terra::mask(sst, ref)


# Calculate NDWS function
calc_ndws <- function(yr, mn){
  outfile <- paste0(out_dir,'/NDWS-',yr,'-',mn,'.tif')
  cat(outfile, "\n")  # Debugging line to print the value of outfile
  
  if(!file.exists(outfile)){
    # Create the directory structure if it doesn't exist
    dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

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
    sr_fls <- paste0(sr_pth,'/Solar-Radiation-Flux_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.1.nc')
    sr_fls <- sr_fls[file.exists(sr_fls)]
    #if(as.numeric(yr) > 2020){
    #  dts_chr <- as.character(dts)
    #  dts_chr <- gsub(pattern = yr, replacement = yrs_mpg$Baseline[yrs_mpg$Future == yr], x = dts_chr)
    #  dts <- as.Date(dts_chr); rm(dts_chr)
    #}
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
    if (yr <= 2020) {srd <- srd/1000000}
    srd <- srd %>% terra::resample(x = ., y = ref) %>% terra::mask(ref)
    
    # Maximum evapotranspiration
    ETMAX <- terra::lapp(x = terra::sds(srd,tmn,tav,tmx), fun = peest)
    rm(list=c("tmn", "tmx", "tav", "srd"))
    gc(verbose=F, full=T, reset=T)
    
    # Compute water balance model
    date <- paste0(yr,'-',mn)
    if(date %in% c('1982-01','2021-01','2041-01','2061-01','2081-01')){
      AVAIL <<- ref
      AVAIL[!is.na(AVAIL)] <- 0
    } else {
#      avail_fl <- list.files(path=dirname(outfile), pattern="AVAIL-")
#      avail_fl <- avail_fl[grep(pattern="\\.tif", avail_fl)]
#      avail_fl <- avail_fl[length(avail_fl)]
#      AVAIL <<- terra::rast(paste0(dirname(outfile),"/", avail_fl))
#      AVAIL <<- AVAIL[[terra::nlyr(AVAIL)]]
      AVAIL <<- terra::rast('F:/Data/NDWS/AVAIL.tif')
         }
    
    eabyep_calc <- function(soilcp = scp, soilsat = ssat, avail = AVAIL,rain = prc[[1]], evap = ETMAX[[1]]){
      
      avail   <- min(avail, soilcp)
      
      # ERATIO
      percwt <- min(avail/soilcp*100, 100)
      percwt <- max(percwt, 1)
      eratio <- min(percwt/(97-3.868*sqrt(soilcp)), 1)
      
      demand  <- eratio * evap
      result  <- avail + rain - demand
      # logging <- result - soilcp
      # logging <- max(logging, 0)
      # logging <- min(logging, soilsat)
      # runoff  <- result - logging + soilcp
      avail   <- min(soilcp, result)
      avail   <- max(avail, 0)
      # runoff  <- max(runoff, 0)
      
      out     <- list(Availability = c(AVAIL, avail),
                      # Demand       = demand,
                      Eratio       = eratio
                      # Logging      = logging
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
    ERATIO <- watbal %>% purrr::map('Eratio') %>% terra::rast()
    # Calculate number of soil water stress days
    NDWS   <- terra::app(x = ERATIO, fun = function(ERATIO){ifelse(ERATIO < 0.5, 1, 0)}) %>% sum()
    terra::writeRaster(NDWS, outfile)
    terra::writeRaster(AVAIL, paste0(dirname(outfile),'/AVAIL-',yr,'-',mn,'.tif'))
    
    #clean up
    rm(list=c("prc", "ETMAX", "AVAIL", "watbal", "ERATIO", "NDWS"))
    gc(verbose=F, full=T, reset=T)
  }
}

# # Historical setup
 #yrs <- 1983:2016
yrs <- 2002:2016
mns <- c(paste0('0',1:9),10:12)
 stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
 names(stp) <- c('yrs','mns')
 stp <- stp %>%
   dplyr::arrange(yrs, mns) %>%
   base::as.data.frame()
 pr_pth <- root # Precipitation
 tm_pth <- paste0(root, 'Tmin') # Minimum temperature
 tx_pth <- paste0(root, 'Tmax') # Maximum temperature
 sr_pth <- paste0(root, 'solar_radiation_flux/solar_radiation_flux') # Solar radiation
 out_dir <- paste0(root,'NDWS')
 1:nrow(stp) %>%
   purrr::map(.f = function(i){calc_ndws(yr = stp$yrs[i], mn = stp$mns[i]); gc(verbose=F, full=T, reset=T)})
