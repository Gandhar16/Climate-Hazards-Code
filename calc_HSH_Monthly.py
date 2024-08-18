# import os
# from datetime import datetime, timedelta
# import numpy as np
# import rasterio
# from rasterio.mask import mask
# from rasterio.transform import from_bounds
# import geopandas as gpd
# from shapely.geometry import box
# import pandas as pd

# # Constants for heat index calculation
# HEAT_INDEX_CONSTANTS = {
#     'c1': -8.78469475556, 'c2': 1.61139411, 'c3': 2.33854883889,
#     'c4': -0.14611605, 'c5': -0.012308094, 'c6': -0.0164248277778,
#     'c7': 2.211732e-3, 'c8': 7.2546e-4, 'c9': -3.582e-6
# }

# def heat_index(tmean, rhum):
#     constants = HEAT_INDEX_CONSTANTS
#     hi = np.where(tmean >= 25,
#                   constants['c1'] + (constants['c2'] * tmean) + (constants['c3'] * rhum) +
#                   (constants['c4'] * tmean * rhum) + (constants['c5'] * tmean**2) +
#                   (constants['c6'] * rhum**2) + (constants['c7'] * tmean**2 * rhum) +
#                   (constants['c8'] * tmean * rhum**2) + (constants['c9'] * tmean**2 * rhum**2),
#                   tmean)
#     return hi

# def clip_raster_data(input_raster_path, ref_raster, top=70, left=-180, right=180, bottom=-60):
#     if not os.path.exists(input_raster_path):
#         print(f"File not found: {input_raster_path}")
#         return None

#     with rasterio.open(input_raster_path) as src:
#         geom = box(left, bottom, right, top)  # Create a bounding box with the specified extent
#         geo = gpd.GeoDataFrame({'geometry': geom}, index=[0], crs=ref_raster.crs)
#         out_image, _ = mask(src, shapes=geo.geometry, crop=True)
#         return out_image[0]

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
#     print(f"TIF file saved successfully at {file_path}")

# def generate_file_paths(start_date, end_date, base_dir):
#     file_paths = []
#     current_date = start_date
#     while current_date <= end_date:
#         date_str = current_date.strftime('%Y.%m.%d')
#         tx_file = os.path.join(base_dir, 'Tmax', f'Tmax.{date_str}.tif')
#         tm_file = os.path.join(base_dir, 'Tmin', f'Tmin.{date_str}.tif')
#         rh_file = os.path.join(base_dir, 'RHum', f'RH.{date_str}.tif')
#         file_paths.append((tx_file, tm_file, rh_file))
#         current_date += timedelta(days=1)
#     return file_paths

# def aggregate_monthly_data(file_paths):
#     monthly_data = {}

#     for tx_file, tm_file, rh_file in file_paths:
#         date_str = os.path.basename(tx_file).split('.')[1:4]
#         date_str = '.'.join(date_str)
#         date = datetime.strptime(date_str, '%Y.%m.%d')
#         month_str = date.strftime('%Y-%m')

#         if month_str not in monthly_data:
#             monthly_data[month_str] = {'tx': [], 'tm': [], 'rh': []}

#         with rasterio.open(tx_file) as src:
#             monthly_data[month_str]['tx'].append(src.read(1))
#         with rasterio.open(tm_file) as src:
#             monthly_data[month_str]['tm'].append(src.read(1))
#         with rasterio.open(rh_file) as src:
#             monthly_data[month_str]['rh'].append(src.read(1))

#     return monthly_data

# def calculate_monthly_averages(monthly_data):
#     monthly_averages = {}

#     for month, data in monthly_data.items():
#         tx_avg = np.mean(data['tx'], axis=0)
#         tm_avg = np.mean(data['tm'], axis=0)
#         rh_avg = np.mean(data['rh'], axis=0)

#         monthly_averages[month] = {'tx': tx_avg, 'tm': tm_avg, 'rh': rh_avg}

#     return monthly_averages

# def calc_hsh_monthly(ref_path, monthly_averages, output_dir, top=70, left=-180, right=180, bottom=-60):
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)

#     with rasterio.open(ref_path) as ref_raster:
#         for month, data in monthly_averages.items():
#             tmx = data['tx']
#             tmn = data['tm']
#             rhm = data['rh']

#             tav = (tmx + tmn) / 2
#             hi = heat_index(tav, rhm)

#             output_path = os.path.join(output_dir, f"HI_{month}.tif")
#             save_tif_file(hi, output_path, top, left, right, bottom, ref_raster.crs)
#             print(f"Month: {month} - Completed")

# def main():
#     start_date = datetime(1983, 1, 1)
#     end_date = datetime(2016, 12, 31)
#     base_dir = "E:/Data"
#     ref_path = "D:/Download/NE1_HR_LC_SR_W/NE1_HR_LC_SR_W/NE1_HR_LC_SR_W.tif"
#     output_dir = "E:/Data/HSH/Final-Monthly"

#     file_paths = generate_file_paths(start_date, end_date, base_dir)
#     monthly_data = aggregate_monthly_data(file_paths)
#     monthly_averages = calculate_monthly_averages(monthly_data)
#     calc_hsh_monthly(ref_path, monthly_averages, output_dir)

# if __name__ == "__main__":
#     main()




import os
from datetime import datetime, timedelta
import numpy as np
import rasterio
from rasterio.mask import mask
from rasterio.transform import from_bounds
import geopandas as gpd
from shapely.geometry import box
import pandas as pd

# Constants for heat index calculation
HEAT_INDEX_CONSTANTS = {
    'c1': -8.78469475556, 'c2': 1.61139411, 'c3': 2.33854883889,
    'c4': -0.14611605, 'c5': -0.012308094, 'c6': -0.0164248277778,
    'c7': 2.211732e-3, 'c8': 7.2546e-4, 'c9': -3.582e-6
}

def heat_index(tmean, rhum):
    constants = HEAT_INDEX_CONSTANTS
    hi = np.where(tmean >= 25,
                  constants['c1'] + (constants['c2'] * tmean) + (constants['c3'] * rhum) +
                  (constants['c4'] * tmean * rhum) + (constants['c5'] * tmean**2) +
                  (constants['c6'] * rhum**2) + (constants['c7'] * tmean**2 * rhum) +
                  (constants['c8'] * tmean * rhum**2) + (constants['c9'] * tmean**2 * rhum**2),
                  tmean)
    return hi

def clip_raster_data(input_raster_path, ref_raster, top=70, left=-180, right=180, bottom=-60):
    if not os.path.exists(input_raster_path):
        print(f"File not found: {input_raster_path}")
        return None

    with rasterio.open(input_raster_path) as src:
        geom = box(left, bottom, right, top)  # Create a bounding box with the specified extent
        geo = gpd.GeoDataFrame({'geometry': geom}, index=[0], crs=ref_raster.crs)
        out_image, _ = mask(src, shapes=geo.geometry, crop=True)
        return out_image[0]

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
    print(f"TIF file saved successfully at {file_path}")

def generate_file_paths(start_date, end_date, base_dir):
    file_paths = []
    current_date = start_date
    while current_date <= end_date:
        date_str = current_date.strftime('%Y.%m.%d')
        tx_file = os.path.join(base_dir, 'Tmax', f'Tmax.{date_str}.tif')
        tm_file = os.path.join(base_dir, 'Tmin', f'Tmin.{date_str}.tif')
        rh_file = os.path.join(base_dir, 'RHum', f'RH.{date_str}.tif')
        file_paths.append((tx_file, tm_file, rh_file))
        current_date += timedelta(days=1)
    return file_paths

def aggregate_monthly_data(file_paths):
    monthly_data = {}

    for tx_file, tm_file, rh_file in file_paths:
        date_str = os.path.basename(tx_file).split('.')[1:4]
        date_str = '.'.join(date_str)
        date = datetime.strptime(date_str, '%Y.%m.%d')
        month_str = date.strftime('%Y-%m')

        if month_str not in monthly_data:
            monthly_data[month_str] = {'tx': [], 'tm': [], 'rh': []}

        with rasterio.open(tx_file) as src:
            monthly_data[month_str]['tx'].append(src.read(1))
        with rasterio.open(tm_file) as src:
            monthly_data[month_str]['tm'].append(src.read(1))
        with rasterio.open(rh_file) as src:
            monthly_data[month_str]['rh'].append(src.read(1))

    return monthly_data

def calculate_monthly_averages(monthly_data):
    monthly_averages = {}

    for month, data in monthly_data.items():
        tx_avg = np.mean(data['tx'], axis=0)
        tm_avg = np.mean(data['tm'], axis=0)
        rh_avg = np.mean(data['rh'], axis=0)

        monthly_averages[month] = {'tx': tx_avg, 'tm': tm_avg, 'rh': rh_avg}

    return monthly_averages

def calc_hsh_monthly(ref_path, monthly_averages, output_dir, top=70, left=-180, right=180, bottom=-60):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with rasterio.open(ref_path) as ref_raster:
        for month, data in monthly_averages.items():
            tmx = data['tx']
            tmn = data['tm']
            rhm = data['rh']

            tav = (tmx + tmn) / 2
            hi = heat_index(tav, rhm)

            year, month_num = month.split('-')
            output_path = os.path.join(output_dir, f"HI_{year}_{month_num}.tif")
            save_tif_file(hi, output_path, top, left, right, bottom, ref_raster.crs)
            print(f"Month: {month} - Completed")

def main():
    start_date = datetime(1983, 2, 1)
    end_date = datetime(1983, 12, 31)
    # end_date = datetime(2016, 12, 31)
    base_dir = "E:/Data"
    ref_path = "D:/Download/NE1_HR_LC_SR_W/NE1_HR_LC_SR_W/NE1_HR_LC_SR_W.tif"
    output_dir = "E:/Data/HSH/Final-Monthly"

    file_paths = generate_file_paths(start_date, end_date, base_dir)
    monthly_data = aggregate_monthly_data(file_paths)
    monthly_averages = calculate_monthly_averages(monthly_data)
    calc_hsh_monthly(ref_path, monthly_averages, output_dir)

if __name__ == "__main__":
    main()
