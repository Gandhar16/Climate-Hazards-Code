# import numpy as np
# import pandas as pd
# import rasterio
# from rasterio.mask import mask
# import os
# from shapely.geometry import box
# import geopandas as gpd
# from datetime import datetime, timedelta
# from rasterio.transform import from_bounds

# # New function for wet bulb temperature calculation
# def wet_bulb_temperature(t, rh):
#     tw = t * np.arctan(0.151977 * np.sqrt(rh + 8.313659)) + np.arctan(t + rh) - np.arctan(rh - 1.676331) + 0.00391838 * rh**1.5 * np.arctan(0.023101 * rh) - 4.686035
#     return tw

# def clip_raster_data(input_raster_path, ref_raster, top=70, left=-180, right=180, bottom=-60):
#     if not os.path.exists(input_raster_path):
#         print(f"File not found: {input_raster_path}")
#         return None

#     with rasterio.open(input_raster_path) as src:
#         geom = box(left, bottom, right, top)  # Create a bounding box with the specified extent
#         geo = gpd.GeoDataFrame({'geometry': geom}, index=[0], crs=ref_raster.crs)
#         out_image, _ = mask(src, shapes=geo.geometry, crop=True)
#         return out_image[0]

# def calc_hsh_single_month(ref_path, file_paths, output_dir, top=70, left=-180, right=180, bottom=-60):
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)

#     with rasterio.open(ref_path) as ref_raster:
#         monthly_hi = None

#         for tx_file, tm_file, rh_file in file_paths:
#             tmx = clip_raster_data(tx_file, ref_raster, top, left, right, bottom)
#             tmn = clip_raster_data(tm_file, ref_raster, top, left, right, bottom)
#             rhm = clip_raster_data(rh_file, ref_raster, top, left, right, bottom)

#             if tmx is not None and tmn is not None and rhm is not None:
#                 tav = (tmx + tmn) / 2
#                 hi = wet_bulb_temperature(tav, rhm)

#                 if monthly_hi is None:
#                     monthly_hi = hi
#                 else:
#                     monthly_hi += hi

#         if monthly_hi is not None:
#             # Compute the average for the month
#             monthly_hi /= len(file_paths)

#             # Extract month and year from the first file name
#             date_str = os.path.basename(file_paths[0][0]).split('.')[1:3]
#             year_month = '.'.join(date_str)

#             output_path = os.path.join(output_dir, f"HSH_TW_{year_month}.tif")
#             save_tif_file(monthly_hi, output_path, top, left, right, bottom, ref_raster.crs)
#             print(f"Month: {year_month} - Completed")

#         # Delete daily files to free up space
#         for tx_file, tm_file, rh_file in file_paths:
#             if os.path.exists(tx_file):
#                 os.remove(tx_file)
#             if os.path.exists(tm_file):
#                 os.remove(tm_file)
#             if os.path.exists(rh_file):
#                 os.remove(rh_file)

#         print(f"Daily files for {year_month} deleted to free up space.")

# def save_tif_file(data, file_path, top, left, right, bottom, crs):
#     height, width = data.shape
#     transform = from_bounds(left, bottom, right, top, width, height)

#     with rasterio.open(
#         file_path,
#         'w',
#         driver='GTiff',
#         height=height,
#         width=width,
#         count=1,
#         dtype=data.dtype,
#         crs=crs,
#         transform=transform,
#     ) as dst:
#         dst.write(data, 1)
#     print("TIF file saved successfully!")

# def generate_file_paths_by_month(start_date, end_date, base_dir):
#     file_paths_by_month = {}
#     current_date = start_date

#     while current_date <= end_date:
#         date_str = current_date.strftime('%Y.%m.%d')
#         tx_file = os.path.join(base_dir, 'Tmax', f'Tmax.{date_str}.tif')
#         tm_file = os.path.join(base_dir, 'Tmin', f'Tmin.{date_str}.tif')
#         rh_file = os.path.join(base_dir, 'RHum', f'RH.{date_str}.tif')

#         year_month = current_date.strftime('%Y-%m')
#         if year_month not in file_paths_by_month:
#             file_paths_by_month[year_month] = []

#         file_paths_by_month[year_month].append((tx_file, tm_file, rh_file))
#         current_date += timedelta(days=1)

#     return file_paths_by_month

# def main():
#     start_date = datetime(2013, 12, 1)
#     end_date = datetime(2016, 12, 31)
#     base_dir = "E:/Data"
#     ref_path = "D:/Download/NE1_HR_LC_SR_W/NE1_HR_LC_SR_W/NE1_HR_LC_SR_W.tif"
#     output_dir = "E:/Data/HSH_TW/Monthly"

#     file_paths_by_month = generate_file_paths_by_month(start_date, end_date, base_dir)

#     for year_month, file_paths in file_paths_by_month.items():
#         calc_hsh_single_month(ref_path, file_paths, output_dir)

# if __name__ == "__main__":
#     main()



import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask
import os
from shapely.geometry import box
import geopandas as gpd
from datetime import datetime, timedelta
from rasterio.transform import from_bounds
import gc

# New function for wet bulb temperature calculation
def wet_bulb_temperature(t, rh):
    tw = t * np.arctan(0.151977 * np.sqrt(rh + 8.313659)) + np.arctan(t + rh) - np.arctan(rh - 1.676331) + 0.00391838 * rh**1.5 * np.arctan(0.023101 * rh) - 4.686035
    return tw

def clip_raster_data(input_raster_path, ref_raster, top=70, left=-180, right=180, bottom=-60):
    if not os.path.exists(input_raster_path):
        print(f"File not found: {input_raster_path}")
        return None

    with rasterio.open(input_raster_path) as src:
        geom = box(left, bottom, right, top)  # Create a bounding box with the specified extent
        geo = gpd.GeoDataFrame({'geometry': geom}, index=[0], crs=ref_raster.crs)
        out_image, _ = mask(src, shapes=geo.geometry, crop=True)
        return out_image[0]

def calc_hsh_single_month(ref_path, file_paths, output_dir, top=70, left=-180, right=180, bottom=-60):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with rasterio.open(ref_path) as ref_raster:
        monthly_hi = None

        for tx_file, tm_file, rh_file in file_paths:
            tmx = clip_raster_data(tx_file, ref_raster, top, left, right, bottom)
            tmn = clip_raster_data(tm_file, ref_raster, top, left, right, bottom)
            rhm = clip_raster_data(rh_file, ref_raster, top, left, right, bottom)

            if tmx is not None and tmn is not None and rhm is not None:
                tav = (tmx + tmn) / 2
                hi = wet_bulb_temperature(tav, rhm)

                if monthly_hi is None:
                    monthly_hi = hi
                else:
                    monthly_hi += hi

                # Free memory used by the daily rasters
                del tmx, tmn, rhm, tav, hi
                gc.collect()

        if monthly_hi is not None:
            # Compute the average for the month
            monthly_hi /= len(file_paths)

            # Extract month and year from the first file name
            date_str = os.path.basename(file_paths[0][0]).split('.')[1:3]
            year_month = '.'.join(date_str)

            output_path = os.path.join(output_dir, f"HSH_TW_{year_month}.tif")
            save_tif_file(monthly_hi, output_path, top, left, right, bottom, ref_raster.crs)
            print(f"Month: {year_month} - Completed")

        # Free memory used by the monthly heat index raster
        del monthly_hi
        gc.collect()

def save_tif_file(data, file_path, top, left, right, bottom, crs):
    height, width = data.shape
    transform = from_bounds(left, bottom, right, top, width, height)

    with rasterio.open(
        file_path,
        'w',
        driver='GTiff',
        height=height,
        width=width,
        count=1,
        dtype=data.dtype,
        crs=crs,
        transform=transform,
    ) as dst:
        dst.write(data, 1)
    print("TIF file saved successfully!")

def generate_file_paths_by_month(start_date, end_date, base_dir):
    file_paths_by_month = {}
    current_date = start_date

    while current_date <= end_date:
        date_str = current_date.strftime('%Y.%m.%d')
        tx_file = os.path.join(base_dir, 'Tmax', f'Tmax.{date_str}.tif')
        tm_file = os.path.join(base_dir, 'Tmin', f'Tmin.{date_str}.tif')
        rh_file = os.path.join(base_dir, 'RHum', f'RH.{date_str}.tif')

        year_month = current_date.strftime('%Y-%m')
        if year_month not in file_paths_by_month:
            file_paths_by_month[year_month] = []

        file_paths_by_month[year_month].append((tx_file, tm_file, rh_file))
        current_date += timedelta(days=1)

    return file_paths_by_month

def main():
    start_date = datetime(1983, 1, 1)
    end_date = datetime(2016, 12, 31)
    base_dir = "E:/Data"
    ref_path = "D:/Download/NE1_HR_LC_SR_W/NE1_HR_LC_SR_W/NE1_HR_LC_SR_W.tif"
    output_dir = "E:/Data/HSH_TW/Monthly"

    file_paths_by_month = generate_file_paths_by_month(start_date, end_date, base_dir)

    for year_month, file_paths in file_paths_by_month.items():
        calc_hsh_single_month(ref_path, file_paths, output_dir)

if __name__ == "__main__":
    main()
