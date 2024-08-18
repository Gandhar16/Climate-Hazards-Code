# import os
# import numpy as np
# import rasterio
# from calendar import monthrange
# import pandas as pd

# def thr_hum_idx(tmax, rhum):
#     # Calculate THI
#     thi = (1.8 * tmax + 32) - ((0.55 - 0.0055 * rhum) * (1.8 * tmax - 26.8))
#     return thi

# def calc_monthly_thi(start_year, end_year, start_month, end_month, tx_pth, rh_pth, out_dir, top=70, left=-180, right=180, bottom=-60):
#     # Loop over years and months
#     for year in range(start_year, end_year + 1):
#         for month in range(start_month, end_month + 1):
#             last_day = monthrange(year, month)[1]
#             dates = pd.date_range(start=f'{year}-{month:02d}-01', end=f'{year}-{month:02d}-{last_day}', freq='D')
            
#             tx_files = [os.path.join(tx_pth, f'Tmax.{date.strftime("%Y.%m.%d")}.tif') for date in dates]
#             rh_files = [os.path.join(rh_pth, f'RH.{date.strftime("%Y.%m.%d")}.tif') for date in dates]
#             thi_sum = None
#             count = 0

#             for date, tx_file, rh_file in zip(dates, tx_files, rh_files):
#                 with rasterio.open(tx_file) as tx_src, rasterio.open(rh_file) as rh_src:
#                     tx_data = tx_src.read(1)
#                     rh_data = rh_src.read(1)
#                     thi_data = thr_hum_idx(tx_data, rh_data)
                    
#                     if thi_sum is None:
#                         thi_sum = np.zeros_like(thi_data)
                    
#                     thi_sum += thi_data
#                     count += 1

#                 # Remove data from memory
#                 del tx_data
#                 del rh_data
#                 del thi_data

#             thi_monthly = thi_sum / count if count > 0 else None
#             if thi_monthly is not None:
#                 height, width = thi_monthly.shape
#                 transform = rasterio.transform.from_bounds(left, bottom, right, top, width, height)
#                 monthly_outfile = os.path.join(out_dir, f'THI_monthly-{year}-{month:02d}.tif')
#                 os.makedirs(os.path.dirname(monthly_outfile), exist_ok=True)

#                 with rasterio.open(
#                     monthly_outfile,
#                     'w',
#                     driver='GTiff',
#                     height=height,
#                     width=width,
#                     count=1,
#                     dtype=rasterio.float32,
#                     crs='EPSG:4326',
#                     transform=transform,
#                     compress='lzw',
#                     nodata=np.nan
#                 ) as dst:
#                     dst.write(thi_monthly, 1)

# def run_script(start_year, end_year, start_month, end_month, sce_climate):
#     tx_pth = 'E:\\Data\\Tmax'
#     rh_pth = 'E:\\Data\\RHum'
#     out_dir = 'E:\\Data\\THI\\Monthly'
#     calc_monthly_thi(start_year, end_year, start_month, end_month, tx_pth, rh_pth, out_dir)

# if __name__ == "__main__":
#     start_year = 1983
#     end_year = 2016
#     start_month = 1
#     end_month = 12
#     sce_climate = "historical"
#     run_script(start_year, end_year, start_month, end_month, sce_climate)

# import os
# import numpy as np
# import rasterio
# from calendar import monthrange
# import pandas as pd
# from tqdm import tqdm

# def thr_hum_idx(tmax, rhum):
#     # Calculate THI
#     thi = (1.8 * tmax + 32) - ((0.55 - 0.0055 * rhum) * (1.8 * tmax - 26.8))
#     return thi

# def calc_thi_monthly(start_year, end_year, start_month, end_month, tx_pth, rh_pth, out_dir, ref_path, top=70, left=-180, right=180, bottom=-60):
#     # Calculate total number of months for the progress bar
#     total_months = (end_year - start_year) * 12 + (end_month - start_month + 1)
    
#     # Loop over years and months
#     with tqdm(total=total_months, desc='Overall Progress', unit='month') as pbar:
#         for year in range(start_year, end_year + 1):
#             for month in range(start_month, end_month + 1):
#                 last_day = monthrange(year, month)[1]
#                 dates = pd.date_range(start=f'{year}-{month:02d}-01', end=f'{year}-{month:02d}-{last_day}', freq='D')
                
#                 monthly_thi_sum = None
#                 daily_count = 0

#                 tx_files = [os.path.join(tx_pth, f'Tmax.{date.strftime("%Y.%m.%d")}.tif') for date in dates]
#                 rh_files = [os.path.join(rh_pth, f'RH.{date.strftime("%Y.%m.%d")}.tif') for date in dates]

#                 for date, tx_file, rh_file in zip(dates, tx_files, rh_files):
#                     with rasterio.open(tx_file) as tx_src, rasterio.open(rh_file) as rh_src:
#                         tx_data = tx_src.read(1)
#                         rh_data = rh_src.read(1)
#                         thi_data = thr_hum_idx(tx_data, rh_data)
                        
#                         if monthly_thi_sum is None:
#                             monthly_thi_sum = np.zeros_like(thi_data)
                        
#                         monthly_thi_sum += thi_data
#                         daily_count += 1

#                 monthly_thi_avg = monthly_thi_sum / daily_count

#                 height, width = monthly_thi_avg.shape
#                 transform = rasterio.transform.from_bounds(left, bottom, right, top, width, height)

#                 outfile = os.path.join(out_dir, f'THI_monthly-{year}-{month:02d}.tif')
#                 os.makedirs(os.path.dirname(outfile), exist_ok=True)

#                 with rasterio.open(
#                     outfile,
#                     'w',
#                     driver='GTiff',
#                     height=height,
#                     width=width,
#                     count=1,
#                     dtype=rasterio.float32,
#                     crs='EPSG:4326',
#                     transform=transform,
#                     compress='lzw',
#                     nodata=np.nan
#                 ) as dst:
#                     dst.write(monthly_thi_avg, 1)

#                 # Clear the monthly THI data from memory
#                 del monthly_thi_sum, monthly_thi_avg

#                 # Update the progress bar
#                 pbar.update(1)

# def run_script(start_year, end_year, start_month, end_month, sce_climate):
#     ref_path = "D:/Download/NE1_HR_LC_SR_W/NE1_HR_LC_SR_W/NE1_HR_LC_SR_W.tif"
#     tx_pth = f'E:\\Data\\Tmax'
#     rh_pth = f'E:\\Data\\RHum'
#     out_dir = 'E:\\Data\\THI\\Monthly'
#     calc_thi_monthly(start_year, end_year, start_month, end_month, tx_pth, rh_pth, out_dir, ref_path)

# if __name__ == "__main__":
#     start_year = 1983
#     end_year = 2016
#     start_month = 1
#     end_month = 12
#     sce_climate = "historical"
#     run_script(start_year, end_year, start_month, end_month, sce_climate)



import os
import numpy as np
import rasterio
from calendar import monthrange
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor

def thr_hum_idx(tmax, rhum):
    # Calculate THI
    thi = (1.8 * tmax + 32) - ((0.55 - 0.0055 * rhum) * (1.8 * tmax - 26.8))
    return thi

def process_daily_thi(tx_file, rh_file):
    with rasterio.open(tx_file) as tx_src, rasterio.open(rh_file) as rh_src:
        tx_data = tx_src.read(1)
        rh_data = rh_src.read(1)
        thi_data = thr_hum_idx(tx_data, rh_data)
    return thi_data

def calc_thi_monthly(start_year, end_year, start_month, end_month, tx_pth, rh_pth, out_dir, ref_path, top=70, left=-180, right=180, bottom=-60):
    # Calculate total number of months for the progress bar
    total_months = (end_year - start_year) * 12 + (end_month - start_month + 1)
    
    # Loop over years and months
    with tqdm(total=total_months, desc='Overall Progress', unit='month') as pbar:
        for year in range(start_year, end_year + 1):
            for month in range(start_month, end_month + 1):
                last_day = monthrange(year, month)[1]
                dates = pd.date_range(start=f'{year}-{month:02d}-01', end=f'{year}-{month:02d}-{last_day}', freq='D')
                
                monthly_thi_sum = None
                daily_count = 0

                tx_files = [os.path.join(tx_pth, f'Tmax.{date.strftime("%Y.%m.%d")}.tif') for date in dates]
                rh_files = [os.path.join(rh_pth, f'RH.{date.strftime("%Y.%m.%d")}.tif') for date in dates]

                with ThreadPoolExecutor() as executor:
                    results = list(executor.map(process_daily_thi, tx_files, rh_files))
                    
                    for thi_data in results:
                        if monthly_thi_sum is None:
                            monthly_thi_sum = np.zeros_like(thi_data)
                        monthly_thi_sum += thi_data
                        daily_count += 1

                monthly_thi_avg = monthly_thi_sum / daily_count

                height, width = monthly_thi_avg.shape
                transform = rasterio.transform.from_bounds(left, bottom, right, top, width, height)

                outfile = os.path.join(out_dir, f'THI_monthly-{year}-{month:02d}.tif')
                os.makedirs(os.path.dirname(outfile), exist_ok=True)

                with rasterio.open(
                    outfile,
                    'w',
                    driver='GTiff',
                    height=height,
                    width=width,
                    count=1,
                    dtype=rasterio.float32,
                    crs='EPSG:4326',
                    transform=transform,
                    compress='lzw',
                    nodata=np.nan
                ) as dst:
                    dst.write(monthly_thi_avg, 1)

                # Clear the monthly THI data from memory
                del monthly_thi_sum, monthly_thi_avg

                # Update the progress bar
                pbar.update(1)

def run_script(start_year, end_year, start_month, end_month, sce_climate):
    ref_path = "D:/Download/NE1_HR_LC_SR_W/NE1_HR_LC_SR_W/NE1_HR_LC_SR_W.tif"
    tx_pth = f'E:\\Data\\Tmax'
    rh_pth = f'E:\\Data\\RHum'
    out_dir = 'E:\\Data\\THI\\Monthly'
    calc_thi_monthly(start_year, end_year, start_month, end_month, tx_pth, rh_pth, out_dir, ref_path)

if __name__ == "__main__":
    start_year = 1984
    end_year = 1984
    start_month = 4
    end_month = 12
    sce_climate = "historical"
    run_script(start_year, end_year, start_month, end_month, sce_climate)
