# R options setup
g <- gc(reset = TRUE); rm(list = ls())  # Clean up the environment
options(warn = -1, scipen = 999)        # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, gtools, lubridate))
source('https://raw.githubusercontent.com/CIAT-DAPA/agro-clim-indices/main/_main_functions.R')

# Define root directory and reference raster
root <- 'F:/Data/'
ref <- terra::rast(paste0(root, 'roi/africa.tif'))

# Load and process soil variables
scp <- terra::rast(paste0(root, 'soils/africa_scp_new.tif')) %>% terra::resample(ref) %>% terra::mask(ref)
sst <- terra::rast(paste0(root, 'soils/africa_ssat.tif')) %>% terra::resample(ref) %>% terra::mask(ref)

# Path definitions for input data and output directory
pr_pth <- 'F:/Data/'
tm_pth <- paste0(root, 'Tmin/')
tx_pth <- paste0(root, 'Tmax/')
sr_pth <- paste0(root, 'solar_radiation_flux/solar_radiation_flux/')
out_dir <- paste0(root, 'NDWL0_daily')

# Function to calculate NDWL0 for a specific day
calc_ndwl0_day <- function(yr, mn, dy){
  outfile <- paste0(out_dir, '/Logging-', yr, '-', mn, '-', dy, '.tif')
  cat("Attempting to process: ", outfile, "\n")
  
  if (!dir.exists(dirname(outfile))) {
    dir.create(dirname(outfile), recursive = TRUE, showWarnings = TRUE)
  }
  
  date <- as.Date(paste0(yr, '-', mn, '-', dy))
  # Check for required files
  pr_file <- paste0(pr_pth, 'chirps-v2.0.', gsub('-', '.', x=date, fixed=TRUE), '.tif')
  tx_file <- paste0(tx_pth, 'Tmax.', gsub('-', '.', x=date, fixed=TRUE), '.tif')
  tm_file <- paste0(tm_pth, 'Tmin.', gsub('-', '.', x=date, fixed=TRUE), '.tif')
  sr_file <- paste0(sr_pth, 'Solar-Radiation-Flux_C3S-glob-agric_AgERA5_', gsub('-', '', x=date, fixed=TRUE), '_final-v1.1.nc')
  sr_file <- sr_file[file.exists(sr_file)]
  
  if (file.exists(pr_file) && file.exists(tx_file) && file.exists(tm_file) && file.exists(sr_file)) {
    cat("All necessary files exist, processing...\n")
    prc <- terra::rast(pr_file) %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
    tmx <- terra::rast(tx_file) %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
    tmn <- terra::rast(tm_file) %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
    tav <- (tmx + tmn) / 2
    srd <- terra::rast(sr_file)
    srd <- srd %>% terra::crop(terra::ext(ref))
    #The line `srd <- srd/1000000` scales the solar radiation data by dividing each value by \(10^6\). This step is likely needed to convert the units of solar radiation.
    srd <- srd/1000000
    srd <- srd %>% terra::resample(x = ., y = ref) %>% terra::mask(ref)
    ETMAX <- terra::lapp(x = terra::sds(srd, tmn, tav, tmx), fun = peest)
    
    # Setup AVAIL based on existing data or initialized state
    avail_file <- paste0(dirname(outfile), '/AVAIL.tif')
    if (file.exists(avail_file)) {
      AVAIL <<- terra::rast(avail_file)
    } else {
      AVAIL <<- ref  # Initialize based on the reference raster
      AVAIL[!is.na(AVAIL)] <- 0
    }
    
    # Compute water balance model
    water_balance <- eabyep_calc(soilcp = scp, soilsat = sst, avail = AVAIL, rain = prc, evap = ETMAX)
    AVAIL <<- water_balance$Availability
    logging <- water_balance$Logging
    
    # Save daily NDWL0 and update AVAIL
    terra::writeRaster(logging, outfile, overwrite = TRUE)
    #terra::writeRaster(AVAIL, avail_file, overwrite = TRUE)
    cat("Saved: ", outfile, "\n")
  } else {
    cat("One or more input files are missing for ", date, "\n")
  }
}

# Setup and execute processing for a range of days
yrs <- 1983:2016
mns <- c(paste0('0',1:9),10:12)
days <- 1:31

stp <- expand.grid(yrs, mns, days) %>% as.data.frame()
names(stp) <- c('yrs', 'mns', 'days')
stp <- stp %>% arrange(yrs, mns, days)

# Iterate through each date to process
purrr::map(1:nrow(stp), function(i){
  yr <- stp$yrs[i]
  mn <- stp$mns[i]
  dy <- stp$days[i]
  if (dy <= lubridate::days_in_month(as.Date(paste0(yr, '-', mn, '-01')))) {
    calc_ndwl0_day(yr = yr, mn = mn, dy = dy)
    gc(verbose = FALSE, full = TRUE, reset = TRUE)
  }
})
