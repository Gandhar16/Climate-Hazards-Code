# Set options and clean environment
g <- gc(reset = TRUE)
rm(list = ls())
options(warn = -1, scipen = 999)


library(terra)
library(tidyverse)
library(raster)


source('https://raw.githubusercontent.com/CIAT-DAPA/agro-clim-indices/main/AWCPTF.R')
suppressMessages(if(!require(pacman)){install.packages('pacman');library(pacman)} else {library(pacman)})
suppressMessages(pacman::p_load(tidyverse,terra,raster))



# Load necessary packages
if(!require(pacman)) {
  install.packages('pacman')
}
library(pacman)
p_load(tidyverse, terra, raster)

# Define paths and system-dependent roots
root <- "F:/Data/"  # Update this path to your local data directory
soils_root <- "F:/Data/SoilData/soil/"  # Update this path to your local SoilGrids data directory


# Load reference and soil depth data
ref <- terra::rast(paste0(root, "/roi/africa.tif"))
sdp <- terra::rast(paste0(root, "SoilData/gyga_af_erzd__m_1km.tif"))

# Soil data path


# Load soil property data
bld <- terra::rast(list.files(paste0(soils_root), pattern = 'bdod_.*cm_mean.tif$', full.names = TRUE) %>% sort())
cec <- terra::rast(list.files(paste0(soils_root), pattern = 'cec_.*cm_mean.tif$', full.names = TRUE) %>% sort())
cly <- terra::rast(list.files(paste0(soils_root), pattern = 'clay_.*cm_mean.tif$', full.names = TRUE) %>% sort())
phx <- terra::rast(list.files(paste0(soils_root), pattern = 'phh2o_.*cm_mean.tif$', full.names = TRUE) %>% sort())
snd <- terra::rast(list.files(paste0(soils_root), pattern = 'sand_.*cm_mean.tif$', full.names = TRUE) %>% sort())
slt <- terra::rast(list.files(paste0(soils_root), pattern = 'silt_.*cm_mean.tif$', full.names = TRUE) %>% sort())
orc <- terra::rast(list.files(paste0(soils_root), pattern = 'ocd_.*cm_mean.tif$', full.names = TRUE) %>% sort())

# Combine and process soil layers
soil <- rast(list(orc, cec, phx, snd, slt, cly, bld))
soil <- soil %>% crop(., ext(ref)) %>% resample(., ref) %>% mask(., ref)

# Extract and arrange soil data at different depth levels
crd <- ref %>% 
  terra::as.data.frame(xy = TRUE, na.rm = TRUE) %>%
  dplyr::mutate(id = row_number()) %>%
  dplyr::select(id, x, y)

soil_data <- cbind(crd, terra::extract(soil, crd[,c('x','y')]))
soil_data2 <- soil_data %>%
  tidyr::pivot_longer(names_to = 'var', values_to = 'val', -c(id, x, y)) %>%
  tidyr::separate(col = 'var', into = c('var', 'depth'), sep = "_", extra = "merge") %>%
  tidyr::pivot_wider(names_from = 'var', values_from = 'val') %>%
  dplyr::arrange(id) %>%
  dplyr::rename(CLYPPT = clay, SNDPPT = sand, SLTPPT = silt, ORCDRC = ocd, BLD = bdod, CEC = cec, PHIHOX = phh2o)

print(soil_data2)
soil_dup <- soil_data2
#debugging coode
# Check if necessary columns exist and are not entirely NA
#necessary_columns <- c("SNDPPT", "SLTPPT", "CLYPPT", "ORCDRC", "BLD", "CEC", "PHIHOX")
#missing_columns <- necessary_columns[!necessary_columns %in% names(soil_data2)]
#if (length(missing_columns) > 0) {
  #cat("Missing columns:", missing_columns, "\n")
}# else {
  # Apply transformations and create new variables only if all columns are present
  #soil_data2$ORCDRC_adj = ifelse(is.na(soil_data2$ORCDRC), NA, soil_data2$ORCDRC/10)
  #soil_data2$BLD_adj = ifelse(is.na(soil_data2$BLD), NA, soil_data2$BLD * 0.001)
  
  # Using corrected column names and adjusted variables
  #soil_data2 <- cbind(soil_data2, AWCPTF(SNDPPT = soil_data2$SNDPPT,
      #                                   SLTPPT = soil_data2$SLTPPT,
     #                                    CLYPPT = soil_data2$CLYPPT,
    #                                     ORCDRC = soil_data2$ORCDRC_adj,
   #                                      BLD = soil_data2$BLD_adj,
  #                                       CEC = soil_data2$CEC,  # Assuming CEC is correctly named
   #                                      PHIHOX = soil_data2$PHIHOX/10,
   #                                      h1 = -10, h2 = -20, h3 = -33))
}

# Print the first few rows to confirm changes
head(soil_data2)
#End debugging


#debugging
print(names(soil_data2))

# Ensure that the correct column names are used
# For example, confirm whether it's 'BLD' or 'BLDFIE', and 'CEC' or 'CECSOL'

# Assuming the correct column names are BLD and CEC and correcting the division and multiplication:
#if("BLD" %in% names(soil_data2) && "CEC" %in% names(soil_data2)) {
 # soil_data2$BLD_adj = ifelse(is.na(soil_data2$BLD), NA, soil_data2$BLD * 0.001)
  #soil_data2$CEC_adj = ifelse(is.na(soil_data2$CEC), NA, soil_data2$CEC)
  
  # Correcting the function call
  soil_data2 <- cbind(soil_data2, AWCPTF(SNDPPT = soil_data2$SNDPPT,
                                         SLTPPT = soil_data2$SLTPPT,
                                         CLYPPT = soil_data2$CLYPPT,
                                         ORCDRC = soil_data2$ORCDRC,
                                         BLD = soil_data2$BLD,
                                         CEC = soil_data2$CEC,
                                         PHIHOX = soil_data2$PHIHOX/10,
                                         h1 = -10, h2 = -20, h3 = -33))
}# else {
  #print("Check column names: BLD and CEC are not found in soil_data2")
}
#end Debugging

print(soil_data2)
#soil_data2 <- cbind(soil_data2, AWCPTF(SNDPPT = soil_data2$SNDPPT,
#                                       SLTPPT = soil_data2$SLTPPT,
#                                       CLYPPT = soil_data2$CLYPPT,
#                                       ORCDRC = soil_data2$ORCDRC,
#                                       BLD = soil_data2$BLDFIE,
#                                       CEC = soil_data2$CECSOL,
#                                       PHIHOX = soil_data2$PHIHOX/10,
#                                       h1 = -10, h2 = -20, h3 = -33))

#print(soil_data2)

# Calculate ASW in mm for each soil horizon
soil_data2$tetaFC <- soil_data2$WWP + soil_data2$AWCh3
soil_data2$AWSat <- soil_data2$tetaS - soil_data2$tetaFC

# Depth adjustments
#depths <- c("sl1", "sl2", "sl3", "sl4", "sl5", "sl6", "sl7")
#values <- c(0, 5, 15, 30, 60, 100, 200)
#names(values) <- depths
#soil_data2$depth <- as.numeric(replace(soil_data2$depth, soil_data2$depth %in% names(values), values[soil_data2$depth]))

print(soil_data2)


#debug
library(dplyr)
library(purrr)

# Function to convert depth descriptions to numeric values based on the lower bound
convert_depth <- function(depth_str) {
  if (grepl("cm_mean", depth_str)) {
    bounds <- as.numeric(unlist(strsplit(gsub("cm_mean", "", depth_str), "-")))
    return(bounds[1])  # Return only the lower bound
  } else {
    return(NA_real_)  # Return NA if the format does not match
  }
}

# Pre-process the data
soil_data2 <- soil_data2 %>%
  mutate(depth = sapply(depth, convert_depth)) %>%
  filter(!is.na(depth)) # Removing rows where depth is NA
print(soil_data2)
#end debug



# Define the soilcap_calc function
soilcap_calc <- function(x, y, rdepth=60, minval, maxval) {
  if (length(x) != length(y)) {stop("length of x and y must be the same")}
  rdepth <- max(c(rdepth, minval)) # Cross check
  rdepth <- min(c(rdepth, maxval)) # Cross-check
  wc_df <- data.frame(depth = y, wc = x)
  if (!rdepth %in% wc_df$depth) {
    wc_df1 <- wc_df[which(wc_df$depth < rdepth),]
    wc_df2 <- wc_df[which(wc_df$depth > rdepth),]
    y1 <- wc_df1$wc[nrow(wc_df1)]; y2 <- wc_df2$wc[1]
    x1 <- wc_df1$depth[nrow(wc_df1)]; x2 <- wc_df2$depth[1]
    ya <- (rdepth - x1) / (x2 - x1) * (y2 - y1) + y1
    wc_df <- rbind(wc_df1, data.frame(depth = rdepth, wc = ya), wc_df2)
  }
  wc_df <- wc_df[which(wc_df$depth <= rdepth),]
  wc_df$soilthick <- wc_df$depth - c(0, wc_df$depth[1:(nrow(wc_df) - 1)])
  wc_df$soilcap <- wc_df$soilthick * wc_df$wc
  soilcp <- sum(wc_df$soilcap) * 10 # in mm
  return(soilcp)
}

print(soilcap_calc)
# Additional soil capacity calculations
# Get root depth per pixel
rdepths <- terra::extract(x = sdp, y = base::as.data.frame(soil_data2[,c('x','y')]))
names(rdepths)[2] <- 'rdepth'
soil_data2$rdepth <- rdepths$rdepth; rm(rdepths)


#debugging
soil_data4 <- soil_data2 %>%
  dplyr::group_by(id) %>%
  dplyr::group_split() %>%  # Corrected: No argument needed here
  purrr::map(.f = function(px) {
    # Assuming rdepth needs to be constant and not extracted from px$rdepth unless you have it defined in each px
    scp  <- soilcap_calc(x = px$AWCh3, y = px$depth, rdepth = px$rdepth, minval = 45, maxval = 100)
    ssat <- soilcap_calc(x = px$AWSat, y = px$depth, rdepth = px$rdepth, minval = 45, maxval = 100)
    df <- data.frame(id = unique(px$id),
                     x  = unique(px$x),
                     y  = unique(px$y),
                     scp  = scp,
                     ssat = ssat)
    return(df)
  }) %>%
  dplyr::bind_rows()
#end debugg




soil_data4 <- soil_data2 %>%
  dplyr::group_by(id) %>%
  dplyr::group_split(id) %>%
  purrr::map(.f = function(px) {
    scp  <- soilcap_calc(x = px$AWCh3, y = px$depth, rdepth = px$rdepth, minval = 45, maxval = 100)
    ssat <- soilcap_calc(x = px$AWSat, y = px$depth, rdepth = px$rdepth, minval = 45, maxval = 100)
    df <- data.frame(id = unique(px$id),
                     x  = unique(px$x),
                     y  = unique(px$y),
                     scp  = scp,
                     ssat = ssat)
    return(df)
  }) %>%
  dplyr::bind_rows()

print(soil_data4)


sscp <- raster::rasterFromXYZ(soil_data4[,c('x','y','scp')], crs = raster::crs(ref))
ssat <- raster::rasterFromXYZ(soil_data4[,c('x','y','ssat')], crs = raster::crs(ref))

outdir <- paste0(root, '/soils')
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
writeRaster(sscp, filename = paste0(outdir, '/africa_scp_new.tif'), overwrite = TRUE)
writeRaster(ssat, filename = paste0(outdir, '/africa_ssat.tif'), overwrite = TRUE)





#debugg
library(raster)
library(sp)
#library(rgdal)

# Assuming soil_data4 is already loaded and structured as shown

# Step 1: Create Spatial Points DataFrame
coordinates(soil_data4) <- ~x+y  # Setting coordinate columns
crs(soil_data4) <- crs(ref)  # Use the CRS of the reference raster

# Step 2: Use properties from reference raster
# Create a raster template with the same extent, resolution, and CRS as the reference raster
raster_template <- raster(extent(ref), nrows=nrow(ref), ncols=ncol(ref))
crs(raster_template) <- crs(ref)

# Step 3: Rasterize both scp and ssat using the reference raster template
scp_raster <- rasterize(soil_data4, raster_template, field="scp", fun=mean)
ssat_raster <- rasterize(soil_data4, raster_template, field="ssat", fun=mean)

# Step 4: Write the Raster to files
outdir <- paste0(root, '/soils')  # Update this to your specific directory path
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
writeRaster(scp_raster, filename = paste0(outdir, '/africa_scp.tif'), format="GTiff", overwrite=TRUE)
writeRaster(ssat_raster, filename = paste0(outdir, '/africa_ssat.tif'), format="GTiff", overwrite=TRUE)

#end Debug



# Output directory and save rasters



#debug
table(names(soil_data2))  # This will give you a count of all column names
dup_names <- names(which(table(names(soil_data2)) > 1))  # Identifying duplicates
print(dup_names)  # Printing duplicate column names

