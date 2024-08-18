# # import os
# # import tkinter as tk
# # from tkinter import filedialog, messagebox
# # import rasterio
# # import xarray as xr
# # import glob
# # from affine import Affine
# # import numpy as np
# # from rasterio.warp import reproject, Resampling

# # def get_affine_transform_from_coords(lon, lat):
# #     lon_res = lon[1] - lon[0]
# #     lat_res = lat[1] - lat[0]
# #     return Affine.translation(lon[0] - lon_res / 2, lat[0] - lat_res / 2) * Affine.scale(lon_res, lat_res)

# # def process_data_in_chunks(file, chunk_size, ref_raster):
# #     with xr.open_dataset(file) as src:
# #         solar_data = src['Solar_Radiation_Flux']
# #         lon, lat = src['lon'].values, src['lat'].values
# #         src_transform = get_affine_transform_from_coords(lon, lat)
# #         num_chunks = solar_data.shape[1] // chunk_size + 1

# #         aggregated_data = np.zeros((ref_raster.height, ref_raster.width), dtype=np.float32)

# #         for i in range(num_chunks):
# #             start = i * chunk_size
# #             end = start + chunk_size
# #             solar_data_chunk = solar_data[:, start:end].to_numpy()

# #             if solar_data_chunk.ndim > 2:
# #                 solar_data_chunk = solar_data_chunk.sum(axis=0)

# #             out_solar_data_chunk = np.empty_like(aggregated_data)
# #             try:
# #                 reproject(
# #                     source=solar_data_chunk,
# #                     destination=out_solar_data_chunk,
# #                     src_transform=src_transform,
# #                     src_crs={'init': 'EPSG:4326'},
# #                     dst_transform=ref_raster.transform,
# #                     dst_crs=ref_raster.crs,
# #                     resampling=Resampling.bilinear
# #                 )
# #             except Exception as e:
# #                 print(f"Reprojection failed: {e}")
# #                 continue

# #             aggregated_data += out_solar_data_chunk

# #         aggregated_data[np.isnan(aggregated_data)] = np.nan

# #     return aggregated_data

# # def load_solar_radiation_data(solar_rad_path, ref_raster, year, month, day):
# #     file_pattern = os.path.join(solar_rad_path, f'*{year}{month:02d}{day:02d}*.nc')
# #     solar_rad_files = glob.glob(file_pattern)
# #     if len(solar_rad_files) == 0:
# #         raise FileNotFoundError(f"No solar radiation files found for {year}-{month:02d}-{day:02d}")

# #     chunk_size = 500
# #     daily_solar_radiation = process_data_in_chunks(solar_rad_files[0], chunk_size, ref_raster)
# #     return daily_solar_radiation

# # def load_climate_data(file_path, ref_raster):
# #     with rasterio.open(file_path) as src:
# #         data = src.read(1, masked=True)
# #         out_data = np.empty(shape=(ref_raster.height, ref_raster.width), dtype=rasterio.float32)
# #         reproject(
# #             source=data,
# #             destination=out_data,
# #             src_transform=src.transform,
# #             src_crs=src.crs,
# #             dst_transform=ref_raster.transform,
# #             dst_crs=ref_raster.crs,
# #             resampling=Resampling.bilinear
# #         )
# #         out_data[out_data == src.nodata] = np.nan
# #         return out_data

# # def calculate_pet(Tmax, Tmin, solar_radiation):
# #     Tmean = (Tmax + Tmin) / 2
# #     Tdiff = np.abs(Tmax - Tmin)
# #     Ra = solar_radiation
# #     PET = 0.0023 * Ra * (Tmean + 17.8) * (Tdiff ** 0.5)
# #     return PET

# # def calculate_tai(precip_data, pet_data):
# #     deficit = np.where(precip_data < pet_data, pet_data - precip_data, 0)
# #     total_deficit = np.sum(deficit, axis=0)
# #     total_pet = np.sum(np.where(precip_data < pet_data, pet_data, 0), axis=0)
# #     total_pet = np.where(total_pet == 0, np.finfo(float).eps, total_pet)
# #     tai = np.where(total_pet > 0, 100 * total_deficit / total_pet, np.nan)
# #     return tai

# # def calc_daily_tai(year, month, day, ref_raster, precip_dataset_path, tmax_dataset_path, tmin_dataset_path, solar_rad_path):
# #     solar_radiation = load_solar_radiation_data(solar_rad_path, ref_raster, year, month, day)

# #     precip_file = os.path.join(precip_dataset_path, f'chirps-v2.0.{year}.{month:02d}.{day:02d}.tif')
# #     tmax_file = os.path.join(tmax_dataset_path, f'Tmax.{year}.{month:02d}.{day:02d}.tif')
# #     tmin_file = os.path.join(tmin_dataset_path, f'Tmin.{year}.{month:02d}.{day:02d}.tif')

# #     precip_data = load_climate_data(precip_file, ref_raster)
# #     tmax_data = load_climate_data(tmax_file, ref_raster)
# #     tmin_data = load_climate_data(tmin_file, ref_raster)

# #     pet_data = calculate_pet(tmax_data, tmin_data, solar_radiation)

# #     tai = calculate_tai(precip_data, pet_data)

# #     return tai

# # def save_tai_raster(tai, output_path, ref_raster):
# #     with rasterio.open(output_path, 'w', **ref_raster.meta) as dst:
# #         dst.write(tai.astype(rasterio.float32), 1)

# # def run_daily_calculation():
# #     try:
# #         year = int(year_entry.get())
# #         month = int(month_entry.get())
# #         day = int(day_entry.get())
# #         ref_raster_path = ref_raster_entry.get()
# #         precip_dataset_path = precip_dataset_entry.get()
# #         tmax_dataset_path = tmax_dataset_entry.get()
# #         tmin_dataset_path = tmin_dataset_entry.get()
# #         solar_rad_dataset_path = solar_rad_dataset_entry.get()
# #         output_directory = filedialog.askdirectory()

# #         ref_raster = rasterio.open(ref_raster_path)

# #         tai = calc_daily_tai(year, month, day, ref_raster, precip_dataset_path, tmax_dataset_path, tmin_dataset_path, solar_rad_dataset_path)
# #         output_file = os.path.join(output_directory, f"TAI_{year}{month:02d}{day:02d}.tif")
# #         save_tai_raster(tai, output_file, ref_raster)

# #         messagebox.showinfo("Success", f"TAI calculation complete for {year}-{month:02d}-{day:02d}")

# #     except Exception as e:
# #         messagebox.showerror("Error", str(e))

# # # GUI Layout
# # root = tk.Tk()
# # root.title("Daily TAI Calculator")

# # tk.Label(root, text="Year:").pack()
# # year_entry = tk.Entry(root)
# # year_entry.pack()

# # tk.Label(root, text="Month:").pack()
# # month_entry = tk.Entry(root)
# # month_entry.pack()

# # tk.Label(root, text="Day:").pack()
# # day_entry = tk.Entry(root)
# # day_entry.pack()

# # tk.Label(root, text="Reference Raster:").pack()
# # ref_raster_entry = tk.Entry(root)
# # ref_raster_entry.pack()
# # tk.Button(root, text="Browse", command=lambda: select_file(ref_raster_entry)).pack()

# # tk.Label(root, text="Precipitation Data Path:").pack()
# # precip_dataset_entry = tk.Entry(root)
# # precip_dataset_entry.pack()
# # tk.Button(root, text="Browse", command=lambda: select_directory(precip_dataset_entry)).pack()

# # tk.Label(root, text="Max Temperature Data Path:").pack()
# # tmax_dataset_entry = tk.Entry(root)
# # tmax_dataset_entry.pack()
# # tk.Button(root, text="Browse", command=lambda: select_directory(tmax_dataset_entry)).pack()

# # tk.Label(root, text="Min Temperature Data Path:").pack()
# # tmin_dataset_entry = tk.Entry(root)
# # tmin_dataset_entry.pack()
# # tk.Button(root, text="Browse", command=lambda: select_directory(tmin_dataset_entry)).pack()

# # tk.Label(root, text="Solar Radiation Data Path:").pack()
# # solar_rad_dataset_entry = tk.Entry(root)
# # solar_rad_dataset_entry.pack()
# # tk.Button(root, text="Browse", command=lambda: select_directory(solar_rad_dataset_entry)).pack()

# # tk.Button(root, text="Calculate Daily TAI", command=run_daily_calculation).pack()

# # root.mainloop()






# # -------------------------------------------------------------------------- Working for some days only not eveyr day

# import os
# import rasterio
# import xarray as xr
# import glob
# from affine import Affine
# import numpy as np
# from rasterio.warp import reproject, Resampling
# from datetime import datetime, timedelta

# def get_affine_transform_from_coords(lon, lat):
#     lon_res = lon[1] - lon[0]
#     lat_res = lat[1] - lat[0]
#     return Affine.translation(lon[0] - lon_res / 2, lat[0] - lat_res / 2) * Affine.scale(lon_res, lat_res)

# def process_data_in_chunks(file, chunk_size, ref_raster):
#     with xr.open_dataset(file) as src:
#         solar_data = src['Solar_Radiation_Flux']
#         lon, lat = src['lon'].values, src['lat'].values
#         src_transform = get_affine_transform_from_coords(lon, lat)
#         num_chunks = solar_data.shape[1] // chunk_size + 1

#         aggregated_data = np.zeros((ref_raster.height, ref_raster.width), dtype=np.float32)

#         for i in range(num_chunks):
#             start = i * chunk_size
#             end = start + chunk_size
#             solar_data_chunk = solar_data[:, start:end].to_numpy()

#             if solar_data_chunk.ndim > 2:
#                 solar_data_chunk = solar_data_chunk.sum(axis=0)

#             out_solar_data_chunk = np.empty_like(aggregated_data)
#             try:
#                 reproject(
#                     source=solar_data_chunk,
#                     destination=out_solar_data_chunk,
#                     src_transform=src_transform,
#                     src_crs={'init': 'EPSG:4326'},
#                     dst_transform=ref_raster.transform,
#                     dst_crs=ref_raster.crs,
#                     resampling=Resampling.bilinear
#                 )
#             except Exception as e:
#                 print(f"Reprojection failed: {e}")
#                 continue

#             aggregated_data += out_solar_data_chunk

#         aggregated_data[np.isnan(aggregated_data)] = np.nan

#     return aggregated_data

# def load_solar_radiation_data(solar_rad_path, ref_raster, year, month, day):
#     file_pattern = os.path.join(solar_rad_path, f'*Solar-Radiation-Flux*_{year}{month:02d}{day:02d}*.nc')
#     solar_rad_files = glob.glob(file_pattern)
#     if len(solar_rad_files) == 0:
#         raise FileNotFoundError(f"No solar radiation files found for {year}-{month:02d}-{day:02d}")

#     chunk_size = 500
#     daily_solar_radiation = process_data_in_chunks(solar_rad_files[0], chunk_size, ref_raster)
#     return daily_solar_radiation

# def load_climate_data(file_path, ref_raster):
#     with rasterio.open(file_path) as src:
#         data = src.read(1, masked=True)
#         out_data = np.empty(shape=(ref_raster.height, ref_raster.width), dtype=rasterio.float32)
#         reproject(
#             source=data,
#             destination=out_data,
#             src_transform=src.transform,
#             src_crs=src.crs,
#             dst_transform=ref_raster.transform,
#             dst_crs=ref_raster.crs,
#             resampling=Resampling.bilinear
#         )
#         out_data[out_data == src.nodata] = np.nan
#         return out_data

# def calculate_pet(Tmax, Tmin, solar_radiation):
#     Tmean = (Tmax + Tmin) / 2
#     Tdiff = np.abs(Tmax - Tmin)
#     Ra = solar_radiation
#     PET = 0.0023 * Ra * (Tmean + 17.8) * (Tdiff ** 0.5)
#     return PET

# def calculate_tai(precip_data, pet_data):
#     deficit = np.where(precip_data < pet_data, pet_data - precip_data, 0)
#     # total_deficit = np.sum(deficit, axis=0)
#     total_pet = np.sum(np.where(precip_data < pet_data, pet_data, 0), axis=0)
#     total_pet = np.where(total_pet == 0, np.finfo(float).eps, total_pet)
#     # tai = np.where(total_pet > 0, 100 * total_deficit / total_pet, np.nan)
#     tai = np.where(total_pet > 0, 100 * deficit / total_pet, np.nan)

#     return tai

# def calc_daily_tai(year, month, day, ref_raster, precip_dataset_path, tmax_dataset_path, tmin_dataset_path, solar_rad_path):
#     solar_radiation = load_solar_radiation_data(solar_rad_path, ref_raster, year, month, day)

#     precip_file = os.path.join(precip_dataset_path, f'chirps-v2.0.{year}.{month:02d}.{day:02d}.tif')
#     tmax_file = os.path.join(tmax_dataset_path, f'Tmax.{year}.{month:02d}.{day:02d}.tif')
#     tmin_file = os.path.join(tmin_dataset_path, f'Tmin.{year}.{month:02d}.{day:02d}.tif')

#     precip_data = load_climate_data(precip_file, ref_raster)
#     tmax_data = load_climate_data(tmax_file, ref_raster)
#     tmin_data = load_climate_data(tmin_file, ref_raster)

#     pet_data = calculate_pet(tmax_data, tmin_data, solar_radiation)

#     tai = calculate_tai(precip_data, pet_data)

#     return tai

# def save_tai_raster(tai, output_path, ref_raster):
#     with rasterio.open(output_path, 'w', **ref_raster.meta) as dst:
#         dst.write(tai.astype(rasterio.float32), 1)

# # Parameters
# start_year = 1983   
# end_year = 1983
# start_month = 1
# end_month = 12

# # File paths
# ref_raster_path = "F://Data//roi//world.tif"
# precip_dataset_path = "F://Data"
# tmax_dataset_path = "F://Data//Tmax"
# tmin_dataset_path = "F://Data//Tmin"
# solar_rad_dataset_path = "F://Data//solar_radiation_flux//solar_radiation_flux"
# output_directory = "F://Data//TAI//NewData"

# # Calculation loop
# ref_raster = rasterio.open(ref_raster_path)
# for year in range(start_year, end_year + 1):
#     for month in range(start_month, end_month + 1):
#         start_day = 1
#         end_day = (datetime(year, month, 1) + timedelta(days=31)).day
#         for day in range(start_day, end_day + 1):
#             try:
#                 tai = calc_daily_tai(year, month, day, ref_raster, precip_dataset_path, tmax_dataset_path, tmin_dataset_path, solar_rad_dataset_path)
#                 output_file = os.path.join(output_directory, f"TAI_{year}{month:02d}{day:02d}.tif")
#                 save_tai_raster(tai, output_file, ref_raster)
#                 print(f"TAI calculation complete for {year}-{month:02d}-{day:02d}")
#             except Exception as e:
#                 print(f"Error on {year}-{month:02d}-{day:02d}: {e}")




# Edited code of above
import os
import rasterio
import xarray as xr
import glob
from affine import Affine
import numpy as np
from rasterio.warp import reproject, Resampling
from datetime import datetime, timedelta
import calendar

def get_affine_transform_from_coords(lon, lat):
    lon_res = lon[1] - lon[0]
    lat_res = lat[1] - lat[0]
    return Affine.translation(lon[0] - lon_res / 2, lat[0] - lat_res / 2) * Affine.scale(lon_res, lat_res)

def process_data_in_chunks(file, chunk_size, ref_raster):
    with xr.open_dataset(file) as src:
        solar_data = src['Solar_Radiation_Flux']
        lon, lat = src['lon'].values, src['lat'].values
        src_transform = get_affine_transform_from_coords(lon, lat)
        num_chunks = solar_data.shape[1] // chunk_size + 1

        aggregated_data = np.zeros((ref_raster.height, ref_raster.width), dtype=np.float32)

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

def load_solar_radiation_data(solar_rad_path, ref_raster, year, month, day):
    file_pattern = os.path.join(solar_rad_path, f'*Solar-Radiation-Flux*_{year}{month:02d}{day:02d}*.nc')
    solar_rad_files = glob.glob(file_pattern)
    if len(solar_rad_files) == 0:
        raise FileNotFoundError(f"No solar radiation files found for {year}-{month:02d}-{day:02d}")

    chunk_size = 500
    daily_solar_radiation = process_data_in_chunks(solar_rad_files[0], chunk_size, ref_raster)
    return daily_solar_radiation

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

def calculate_pet(Tmax, Tmin, solar_radiation):
    Tmean = (Tmax + Tmin) / 2
    Tdiff = np.abs(Tmax - Tmin)
    Ra = solar_radiation
    # PET = 0.0023 * Ra * (Tmean + 17.8) * (Tdiff ** 0.5)
    PET = 0.0023 * (Tmean + 17.8) * (Tdiff ** 0.5)
    print(f"Max PET calculated: {np.nanmax(PET)}")
    return PET

# Working PET
# def calculate_pet(Tmax, Tmin, solar_radiation):
#     Tmean = (Tmax + Tmin) / 2
#     Tdiff = np.abs(Tmax - Tmin)
#     Ra = solar_radiation
    
#     # Check for unreasonable solar radiation values
#     if np.any(Ra > 1000):  # Assuming solar radiation values should be within a realistic range
#         print("Warning: High solar radiation values detected")
    
#     # PET = 0.0023 * Ra * (Tmean + 17.8) * (Tdiff ** 0.5)
#     PET = 0.0023  * (Tmean + 17.8) * (Tdiff ** 0.5)
    
#     # Debug output for checking maximum PET values calculated
#     print(f"Max PET calculated: {np.nanmax(PET)}")
    
#     return PET



def calculate_tai(precip_data, pet_data):
    deficit = np.where(precip_data < pet_data, pet_data - precip_data, 0)
   
    total_pet = np.sum(np.where(precip_data < pet_data, pet_data, 0), axis=0)
    total_pet = np.where(total_pet == 0, np.finfo(float).eps, total_pet)
 
    tai = np.where(total_pet > 0, 100 * deficit / total_pet, np.nan)
    return tai

# Working TAI
# def calculate_tai(precip_data, pet_data):
#     deficit = np.where(precip_data < pet_data, pet_data - precip_data, 0)
#     total_pet = np.sum(np.where(precip_data < pet_data, pet_data, 0), axis=0)
#     total_pet = np.where(total_pet == 0, np.finfo(float).eps, total_pet)
#     tai = np.where(total_pet > 0, 100 * deficit / total_pet, np.nan)
#     return tai

def calc_daily_tai(year, month, day, ref_raster, precip_dataset_path, tmax_dataset_path, tmin_dataset_path, solar_rad_path):
    solar_radiation = load_solar_radiation_data(solar_rad_path, ref_raster, year, month, day)

    precip_file = os.path.join(precip_dataset_path, f'chirps-v2.0.{year}.{month:02d}.{day:02d}.tif')
    tmax_file = os.path.join(tmax_dataset_path, f'Tmax.{year}.{month:02d}.{day:02d}.tif')
    tmin_file = os.path.join(tmin_dataset_path, f'Tmin.{year}.{month:02d}.{day:02d}.tif')

    precip_data = load_climate_data(precip_file, ref_raster)
    tmax_data = load_climate_data(tmax_file, ref_raster)
    tmin_data = load_climate_data(tmin_file, ref_raster)

    pet_data = calculate_pet(tmax_data, tmin_data, solar_radiation)
    tai = calculate_tai(precip_data, pet_data)
    return tai

def save_tai_raster(tai, output_path, ref_raster):
    with rasterio.open(output_path, 'w', **ref_raster.meta) as dst:
        dst.write(tai.astype(rasterio.float32), 1)

# Parameters
start_year = 1983   
end_year = 2016
start_month = 1
end_month = 12

# File paths
ref_raster_path = "E://Data//roi//kenya_zambia.tif"
precip_dataset_path = "E://Data"
tmax_dataset_path = "E://Data//Tmax"
tmin_dataset_path = "E://Data//Tmin"
solar_rad_dataset_path = "E://Data//solar_radiation_flux//solar_radiation_flux"
output_directory = "E://Data//TAI//New"


# Calculation loop
ref_raster = rasterio.open(ref_raster_path)
for year in range(start_year, end_year + 1):
    for month in range(start_month, end_month + 1):
        # Correctly calculate the number of days in the month
        end_day = calendar.monthrange(year, month)[1]
        for day in range(1, end_day + 1):
            try:
                tai = calc_daily_tai(year, month, day, ref_raster, precip_dataset_path, tmax_dataset_path, tmin_dataset_path, solar_rad_dataset_path)
                output_file = os.path.join(output_directory, f"TAI_{year}{month:02d}{day:02d}.tif")
                save_tai_raster(tai, output_file, ref_raster)
                print(f"TAI calculation complete for {year}-{month:02d}-{day:02d}")
            except Exception as e:
                print(f"Error on {year}-{month:02d}-{day:02d}: {e}")





# # ----------------------------------------------------------------------------------------------------- Working Full



# import os
# import rasterio
# import xarray as xr
# import glob
# from affine import Affine
# import numpy as np
# from rasterio.warp import reproject, Resampling
# from datetime import datetime, timedelta
# import calendar

# def get_affine_transform_from_coords(lon, lat):
#     lon_res = lon[1] - lon[0]
#     lat_res = lat[1] - lat[0]
#     return Affine.translation(lon[0] - lon_res / 2, lat[0] - lat_res / 2) * Affine.scale(lon_res, lat_res)

# def process_data_in_chunks(file, chunk_size, ref_raster):
#     with xr.open_dataset(file) as src:
#         solar_data = src['Solar_Radiation_Flux']
#         lon, lat = src['lon'].values, src['lat'].values
#         src_transform = get_affine_transform_from_coords(lon, lat)
#         num_chunks = solar_data.shape[1] // chunk_size + 1

#         aggregated_data = np.zeros((ref_raster.height, ref_raster.width), dtype=np.float32)

#         for i in range(num_chunks):
#             start = i * chunk_size
#             end = start + chunk_size
#             solar_data_chunk = solar_data[:, start:end].to_numpy()

#             if solar_data_chunk.ndim > 2:
#                 solar_data_chunk = solar_data_chunk.sum(axis=0)

#             out_solar_data_chunk = np.empty_like(aggregated_data)
#             try:
#                 reproject(
#                     source=solar_data_chunk,
#                     destination=out_solar_data_chunk,
#                     src_transform=src_transform,
#                     src_crs={'init': 'EPSG:4326'},
#                     dst_transform=ref_raster.transform,
#                     dst_crs=ref_raster.crs,
#                     resampling=Resampling.bilinear
#                 )
#             except Exception as e:
#                 print(f"Reprojection failed: {e}")
#                 continue

#             aggregated_data += out_solar_data_chunk

#         aggregated_data[np.isnan(aggregated_data)] = np.nan

#     return aggregated_data

# def load_solar_radiation_data(solar_rad_path, ref_raster, year, month, day):
#     # file_pattern = os.path.join(solar_rad_path, f'Solar-Radiation-Flux*_{year}{month:02d}{day:02d}*.nc')
#     file_pattern = os.path.join(solar_rad_path, f'Solar-Radiation-Flux_C3S-glob-agric_AgERA5_{year}{month:02d}{day:02d}_final-v1.1.nc')
#     solar_rad_files = glob.glob(file_pattern)
#     if len(solar_rad_files) == 0:
#         raise FileNotFoundError(f"No solar radiation files found for {year}-{month:02d}-{day:02d}")

#     chunk_size = 500
#     daily_solar_radiation = process_data_in_chunks(solar_rad_files[0], chunk_size, ref_raster)
#     return daily_solar_radiation

# def load_climate_data(file_path, ref_raster):
#     with rasterio.open(file_path) as src:
#         data = src.read(1, masked=True)
#         out_data = np.empty(shape=(ref_raster.height, ref_raster.width), dtype=rasterio.float32)
#         reproject(
#             source=data,
#             destination=out_data,
#             src_transform=src.transform,
#             src_crs=src.crs,
#             dst_transform=ref_raster.transform,
#             dst_crs=ref_raster.crs,
#             resampling=Resampling.bilinear
#         )
#         out_data[out_data == src.nodata] = np.nan
#         return out_data

# def calculate_pet(Tmax, Tmin, solar_radiation):
#     Tmean = (Tmax + Tmin) / 2
#     Tdiff = np.abs(Tmax - Tmin)
#     Ra = solar_radiation
#     PET = 0.0023 * Ra * (Tmean + 17.8) * (Tdiff ** 0.5)
#     print(PET)
#     return PET


# def calculate_tai(precip_data, pet_data):
#     # Define a larger epsilon to handle small values in pet_data
#     epsilon = 1e-3  # Adjust this value as needed

#     deficit = np.where(precip_data < pet_data, pet_data - precip_data, 0)
#     tai = np.zeros_like(pet_data)

#     # Handle non-zero pet_data values
#     nonzero_indices = pet_data > epsilon
#     tai[nonzero_indices] = 100 * deficit[nonzero_indices] / pet_data[nonzero_indices]

#     # Handle small pet_data values (between 0 and epsilon)
#     small_indices = (pet_data > 0) & (pet_data <= epsilon)
#     tai[small_indices] = 100 * deficit[small_indices] / epsilon

#     return tai


# def calculate_tai(precip_data, pet_data):
#     deficit = np.where(precip_data < pet_data, pet_data - precip_data, 0)
#     total_deficit = np.sum(deficit, axis=0)
#     total_pet = np.sum(np.where(precip_data < pet_data, pet_data, 0), axis=0)
#     total_pet = np.where(total_pet == 0, np.finfo(float).eps, total_pet)
#     tai = np.where(total_pet > 0, 100 * total_deficit / total_pet, np.nan)
#     # tai = np.where(total_pet > 0, 100 * deficit / total_pet, np.nan)
# #     # New formula
# #     deficit = np.where(precip_data < pet_data, pet_data - precip_data, 0)
# #     tai = np.where(pet_data> 0, 100 * deficit / pet_data, 0)
# #     print(tai)
#     return tai


# # def calculate_tai(precip_data, pet_data):
# #     # Define a small epsilon to avoid division by zero
# #     epsilon = 1e-10
# #     deficit = np.where(precip_data < pet_data, pet_data - precip_data, 0)
# #     # Add epsilon to pet_data in the division
# #     tai = np.where(pet_data > 0, 100 * deficit / (pet_data + epsilon), 0)
# #     print(tai)
# #     return tai




# # def calculate_tai(precip_data, pet_data):
# #     # Ensure dimensions match
# #     assert precip_data.shape == pet_data.shape, "Shape mismatch between precipitation and PET data"
    
# #     # Calculate deficit for each day
# #     daily_deficit = np.where(precip_data < pet_data, pet_data - precip_data, 0)
    
# #     # Calculate total PET for each day
# #     total_pet = np.sum(np.where(precip_data < pet_data, pet_data, 0), axis=0)
    
# #     # Avoid division by zero
# #     total_pet = np.where(total_pet == 0, np.finfo(float).eps, total_pet)
    
# #     # Calculate TAI for each day
# #     daily_tai = np.where(total_pet > 0, 100 * daily_deficit / total_pet, np.nan)

# #     return daily_tai





#     # return tai
# def calc_daily_tai(year, month, day, ref_raster, precip_dataset_path, tmax_dataset_path, tmin_dataset_path, solar_rad_path):
#     try:
#         solar_radiation = load_solar_radiation_data(solar_rad_path, ref_raster, year, month, day)
#     except FileNotFoundError:
#         print(f"Solar radiation data not found for {year}-{month:02d}-{day:02d}")
#         return None

#     try:
#         precip_file = os.path.join(precip_dataset_path, f'chirps-v2.0.{year}.{month:02d}.{day:02d}.tif')
#         tmax_file = os.path.join(tmax_dataset_path, f'Tmax.{year}.{month:02d}.{day:02d}.tif')
#         tmin_file = os.path.join(tmin_dataset_path, f'Tmin.{year}.{month:02d}.{day:02d}.tif')

#         precip_data = load_climate_data(precip_file, ref_raster)
#         tmax_data = load_climate_data(tmax_file, ref_raster)
#         tmin_data = load_climate_data(tmin_file, ref_raster)
#     except FileNotFoundError as e:
#         print(f"Climate data not found for {year}-{month:02d}-{day:02d}: {e}")
#         return None

#     pet_data = calculate_pet(tmax_data, tmin_data, solar_radiation)
#     tai = calculate_tai(precip_data, pet_data)
#     return tai


# def save_tai_raster(tai, output_path, ref_raster):
#     with rasterio.open(output_path, 'w', **ref_raster.meta) as dst:
#         dst.write(tai.astype(rasterio.float32), 1)

# # Parameters
# start_year = 1983   
# end_year = 1983
# start_month = 1
# end_month = 12

# # File paths
# # ref_raster_path = "C:/Users/gkhandag/Downloads/NE1_HR_LC_SR_W/NE1_HR_LC_SR_W.tif"
# ref_raster_path = "F://Data//roi//africa.tif"
# precip_dataset_path = "F://Data"
# tmax_dataset_path = "F://Data//Tmax"
# tmin_dataset_path = "F://Data//Tmin"
# solar_rad_dataset_path = "F://Data//solar_radiation_flux//solar_radiation_flux"
# output_directory = "F://Data//TAI"


# # Debugg

# # def check_raster_extent(raster_path):
# #     with rasterio.open(raster_path) as src:
# #         print(f"Bounds: {src.bounds}")
# #         print(f"CRS: {src.crs}")
# #         print(f"Dimensions: {src.width}x{src.height}")
# #         print(f"Resolution: {src.res}")

# # # Call this function with your reference raster
# # check_raster_extent(ref_raster_path)

# # End Debugg




# # Calculation loop
# # Calculation loop


# # Calculation loop
# ref_raster = rasterio.open(ref_raster_path)
# for year in range(start_year, end_year + 1):
#     for month in range(start_month, end_month + 1):
#         start_day = 1
#         end_day = calendar.monthrange(year, month)[1]
#         for day in range(start_day, end_day + 1):
#             try:
#                 tai = calc_daily_tai(year, month, day, ref_raster, precip_dataset_path, tmax_dataset_path, tmin_dataset_path, solar_rad_dataset_path)
#                 if tai is not None:
#                     output_file = os.path.join(output_directory, f"TAI_{year}_{month:02d}_{day:02d}.tif")
#                     save_tai_raster(tai, output_file, ref_raster)
#                     print(f"TAI calculation complete for {year}-{month:02d}-{day:02d}")
#             except Exception as e:
#                 print(f"Error on {year}-{month:02d}-{day:02d}: {e}")



# # Checking the lats and lons of Solar_Radiation_Flux

# # import xarray as xr

# # # Load the NetCDF file
# # file_path = 'F://Data//solar_radiation_flux//solar_radiation_flux//Solar-Radiation-Flux_C3S-glob-agric_AgERA5_19830103_final-v1.1.nc'
# # data = xr.open_dataset(file_path)

# # # Print latitude and longitude ranges
# # lat = data['lat']
# # lon = data['lon']

# # print("Latitude range:", lat.min().values, "to", lat.max().values)
# # print("Longitude range:", lon.min().values, "to", lon.max().values)
