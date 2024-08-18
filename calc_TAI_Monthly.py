import os
import rasterio
import xarray as xr
import glob
import numpy as np
from affine import Affine
from rasterio.warp import reproject, Resampling
from datetime import datetime
import calendar

def get_affine_transform_from_coords(lon, lat):
    lon_res = lon[1] - lon[0]
    lat_res = lat[1] - lat[0]
    return Affine.translation(lon[0] - lon_res / 2, lat[0] - lat_res / 2) * Affine.scale(lon_res, lat_res)

def process_data_in_chunks(file_list, chunk_size, ref_raster):
    aggregated_data = np.zeros((ref_raster.height, ref_raster.width), dtype=np.float32)
    for file in file_list:
        with xr.open_dataset(file) as src:
            solar_data = src['Solar_Radiation_Flux']
            lon, lat = src['lon'].values, src['lat'].values
            src_transform = get_affine_transform_from_coords(lon, lat)
            num_chunks = solar_data.shape[1] // chunk_size + 1

            for i in range(num_chunks):
                start = i * chunk_size
                end = start + chunk_size
                solar_data_chunk = solar_data[:, start:end].to_numpy()

                if solar_data_chunk.ndim > 2:
                    solar_data_chunk = solar_data_chunk.sum(axis=0)

                out_solar_data_chunk = np.empty_like(aggregated_data)
                reproject(
                    source=solar_data_chunk,
                    destination=out_solar_data_chunk,
                    src_transform=src_transform,
                    src_crs={'init': 'EPSG:4326'},
                    dst_transform=ref_raster.transform,
                    dst_crs=ref_raster.crs,
                    resampling=Resampling.bilinear
                )

                aggregated_data += out_solar_data_chunk

    aggregated_data[np.isnan(aggregated_data)] = np.nan
    return aggregated_data

def calculate_pet(Tmax, Tmin, solar_radiation):
    Tmean = (Tmax + Tmin) / 2
    Tdiff = np.abs(Tmax - Tmin)
    Ra = solar_radiation
    
    # Check for unreasonable solar radiation values
    if np.any(Ra > 1000):  # Assuming solar radiation values should be within a realistic range
        print("Warning: High solar radiation values detected")
    
    # PET = 0.0023 * Ra * (Tmean + 17.8) * (Tdiff ** 0.5)
    PET = 0.0023  * (Tmean + 17.8) * (Tdiff ** 0.5)
    
    # Debug output for checking maximum PET values calculated
    print(f"Max PET calculated: {np.nanmax(PET)}")
    
    return PET



def calculate_tai(precip_data, pet_data):
    deficit = np.where(precip_data < pet_data, pet_data - precip_data, 0)
   
    total_pet = np.sum(np.where(precip_data < pet_data, pet_data, 0), axis=0)
    total_pet = np.where(total_pet == 0, np.finfo(float).eps, total_pet)
 
    tai = np.where(total_pet > 0, 100 * deficit / total_pet, np.nan)
    tai = np.where(tai > 5000, -1, tai)
    return tai

def load_monthly_data(file_pattern, ref_raster, aggregation='sum'):
    file_list = glob.glob(file_pattern)
    if not file_list:
        raise FileNotFoundError(f"No files found for pattern {file_pattern}")

    if aggregation == 'sum':
        data_aggregated = np.zeros((ref_raster.height, ref_raster.width), dtype=np.float32)
        for file in file_list:
            data = load_climate_data(file, ref_raster)
            data_aggregated += data
    elif aggregation == 'mean':
        data_list = [load_climate_data(file, ref_raster) for file in file_list]
        data_aggregated = np.nanmean(data_list, axis=0)

    return data_aggregated

def load_climate_data(file_path, ref_raster):
    with rasterio.open(file_path) as src:
        data = src.read(1, masked=True)
        out_data = np.empty(shape=(ref_raster.height, ref_raster.width), dtype=rasterio.float32)
        reproject(
            source=data,
            destination=out_data,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=ref_raster.transform,
            dst_crs=ref_raster.crs,
            resampling=Resampling.bilinear
        )
        out_data[out_data == src.nodata] = np.nan
        return out_data

def calc_monthly_tai(year, month, ref_raster, precip_dataset_path, tmax_dataset_path, tmin_dataset_path, solar_rad_path):
    solar_rad_pattern = os.path.join(solar_rad_path, f'*Solar-Radiation-Flux*_{year}{month:02d}*.nc')
    solar_radiation = process_data_in_chunks(glob.glob(solar_rad_pattern), 500, ref_raster)

    precip_pattern = os.path.join(precip_dataset_path, f'chirps-v2.0.{year}.{month:02d}*.tif')
    tmax_pattern = os.path.join(tmax_dataset_path, f'Tmax.{year}.{month:02d}*.tif')
    tmin_pattern = os.path.join(tmin_dataset_path, f'Tmin.{year}.{month:02d}*.tif')

    precip_data = load_monthly_data(precip_pattern, ref_raster, 'sum')
    tmax_data = load_monthly_data(tmax_pattern, ref_raster, 'mean')
    tmin_data = load_monthly_data(tmin_pattern, ref_raster, 'mean')

    pet_data = calculate_pet(tmax_data, tmin_data, solar_radiation)
    tai = calculate_tai(precip_data, pet_data)
    return tai

def save_tai_raster(tai, output_path, ref_raster):
    with rasterio.open(output_path, 'w', **ref_raster.meta) as dst:
        dst.write(tai.astype(rasterio.float32), 1)

# For leap year check
def is_valid_date(year, month, day):
    if month == 2 and day == 29:
        return calendar.isleap(year)
    return True


# Parameters
start_year = 1983
end_year = 2016
start_month = 1
end_month = 12

# File paths
ref_raster_path = "F://Data//roi//kenya_zambia.tif"
precip_dataset_path = "F://Data"
tmax_dataset_path = "F://Data//Tmax"
tmin_dataset_path = "F://Data//Tmin"
solar_rad_dataset_path = "F://Data//solar_radiation_flux//solar_radiation_flux"
output_directory = "F://Data//TAI//MonthlyData"

# Calculation loop
ref_raster = rasterio.open(ref_raster_path)
for year in range(start_year, end_year + 1):
    for month in range(start_month, end_month + 1):
        try:
            tai = calc_monthly_tai(year, month, ref_raster, precip_dataset_path, tmax_dataset_path, tmin_dataset_path, solar_rad_dataset_path)
            output_file = os.path.join(output_directory, f"TAI_{year}_{month:02d}.tif")
            save_tai_raster(tai, output_file, ref_raster)
            print(f"TAI calculation complete for {year}-{month:02d}")
        except Exception as e:
            print(f"Error on {year}-{month:02d}: {e}")
    