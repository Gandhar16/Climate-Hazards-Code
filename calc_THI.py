# # -----------------------------------------------------Complete Working Code 

# import os
# import numpy as np
# import rasterio
# from calendar import monthrange
# import pandas as pd
# import tkinter as tk
# from tkinter import simpledialog
# import matplotlib.pyplot as plt
# from rasterio.plot import show

# def thr_hum_idx(tmax, rhum):
#     # Calculate THI
#     thi = (1.8 * tmax + 32) - ((0.55 - 0.0055 * rhum) * (1.8 * tmax - 26.8))
#     return thi

# def calc_thi(years, months, tx_pth, rh_pth, out_dir, ref_path):
#     # Loop over years and months
#     for year in years:
#         for month in months:
#             outfile1 = os.path.join(out_dir, f'THI_mean-{year}-{month:02d}.tif')
#             outfile2 = os.path.join(out_dir, f'THI_max-{year}-{month:02d}.tif')
#             outfile3 = os.path.join(out_dir, 'daily', f'THI_daily-{year}-{month:02d}.tif')
            
#             os.makedirs(os.path.dirname(outfile3), exist_ok=True)
#             last_day = monthrange(year, month)[1]
#             dates = pd.date_range(start=f'{year}-{month:02d}-01', end=f'{year}-{month:02d}-{last_day}', freq='D')
            
#             tx_files = [os.path.join(tx_pth, f'Tmax.{date.strftime("%Y.%m.%d")}.tif') for date in dates]
#             rh_files = [os.path.join(rh_pth, f'RH.{date.strftime("%Y.%m.%d")}.tif') for date in dates]

#             thi_daily_list = []

#             for tx_file, rh_file in zip(tx_files, rh_files):
#                 with rasterio.open(tx_file) as tx_src, rasterio.open(rh_file) as rh_src:
#                     tx_data = tx_src.read(1)
#                     rh_data = rh_src.read(1)
#                     thi_data = thr_hum_idx(tx_data, rh_data)
#                     thi_daily_list.append(thi_data)
                    
#             thi_daily_array = np.array(thi_daily_list)
#             thi_mean_array = np.mean(thi_daily_array, axis=0)
#             thi_max_array = np.max(thi_daily_array, axis=0)

#             with rasterio.open(ref_path) as ref:
#                 ref_meta = ref.meta.copy()
#                 ref_meta.update(dtype=rasterio.float32, count=1, compress='lzw', nodata=np.nan)
                
#                 with rasterio.open(outfile1, 'w', **ref_meta) as dst:
#                     dst.write(thi_mean_array.astype(rasterio.float32), 1)
#                 with rasterio.open(outfile2, 'w', **ref_meta) as dst:
#                     dst.write(thi_max_array.astype(rasterio.float32), 1)

#                 ref_meta.update(count=len(dates))
#                 with rasterio.open(outfile3, 'w', **ref_meta) as dst:
#                     for i in range(thi_daily_array.shape[0]):
#                         dst.write(thi_daily_array[i].astype(rasterio.float32), i + 1)

# def plot_thi(file_path):
#     # Plot THI data
#     with rasterio.open(file_path) as src:
#         fig, ax = plt.subplots()
#         show(src, ax=ax, title='THI Data')
#         plt.show()

# def run_script(years, months, sce_climate):
#     ref_path = 'E:\\Data\\roi\\kenya_zambia.tif'
#     tx_pth = f'E:\\Data\\Tmax\\{years[0]}'
#     rh_pth = f'E:\\Data\\RHum\\{years[0]}'
#     out_dir = 'E:\\Data\\THI'
#     calc_thi(years, months, tx_pth, rh_pth, out_dir, ref_path)
#     for year in years:
#         for month in months:
#             plot_thi(os.path.join(out_dir, f'THI_mean-{year}-{month:02d}.tif'))

# def select_years():
#     # Get years from user
#     input_years = simpledialog.askstring("Input", "Enter years (comma separated, e.g., 2020,2021,2022):",
#                                          parent=root)
#     return [int(year) for year in input_years.split(',')]

# def select_months():
#     # Get months from user
#     input_months = simpledialog.askstring("Input", "Enter months (comma separated, e.g., 1,2,3):",
#                                           parent=root)
#     return [int(month) for month in input_months.split(',')]

# root = tk.Tk()
# root.title("THI Calculation")

# # GUI components
# scenario_var = tk.StringVar(root)
# scenario_var.set("historical")
# scenario_menu = tk.OptionMenu(root, scenario_var, "historical", "future")
# scenario_menu.pack()

# run_button = tk.Button(root, text="Run", command=lambda: run_script(select_years(), select_months(), scenario_var.get()))
# run_button.pack()

# root.mainloop()









# Daily THI, and mean and max calculation with the extent as wel--------------------------------------------------------------------------

# import os
# import numpy as np
# import rasterio
# from calendar import monthrange
# import pandas as pd

# def thr_hum_idx(tmax, rhum):
#     # Calculate THI
#     thi = (1.8 * tmax + 32) - ((0.55 - 0.0055 * rhum) * (1.8 * tmax - 26.8))
#     return thi

# def calc_thi(start_year, end_year, start_month, end_month, tx_pth, rh_pth, out_dir, ref_path, top=70, left=-180, right=180, bottom=-60):
#     # Loop over years and months
#     for year in range(start_year, end_year + 1):
#         for month in range(start_month, end_month + 1):
#             outfile1 = os.path.join(out_dir, f'THI_mean-{year}-{month:02d}.tif')
#             outfile2 = os.path.join(out_dir, f'THI_max-{year}-{month:02d}.tif')
#             outfile3 = os.path.join(out_dir, 'daily', f'THI_daily-{year}-{month:02d}.tif')
            
#             os.makedirs(os.path.dirname(outfile3), exist_ok=True)
#             last_day = monthrange(year, month)[1]
#             dates = pd.date_range(start=f'{year}-{month:02d}-01', end=f'{year}-{month:02d}-{last_day}', freq='D')
            
#             tx_files = [os.path.join(tx_pth, f'Tmax.{date.strftime("%Y.%m.%d")}.tif') for date in dates]
#             rh_files = [os.path.join(rh_pth, f'RH.{date.strftime("%Y.%m.%d")}.tif') for date in dates]

#             thi_daily_list = []

#             for tx_file, rh_file in zip(tx_files, rh_files):
#                 with rasterio.open(tx_file) as tx_src, rasterio.open(rh_file) as rh_src:
#                     tx_data = tx_src.read(1)
#                     rh_data = rh_src.read(1)
#                     thi_data = thr_hum_idx(tx_data, rh_data)
#                     thi_daily_list.append(thi_data)
                    
#             thi_daily_array = np.array(thi_daily_list)
#             thi_mean_array = np.mean(thi_daily_array, axis=0)
#             thi_max_array = np.max(thi_daily_array, axis=0)

#             height, width = thi_mean_array.shape
#             transform = rasterio.transform.from_bounds(left, bottom, right, top, width, height)

#             with rasterio.open(
#                 outfile1,
#                 'w',
#                 driver='GTiff',
#                 height=height,
#                 width=width,
#                 count=1,
#                 dtype=rasterio.float32,
#                 crs='EPSG:4326',
#                 transform=transform,
#                 compress='lzw',
#                 nodata=np.nan
#             ) as dst:
#                 dst.write(thi_mean_array, 1)

#             with rasterio.open(
#                 outfile2,
#                 'w',
#                 driver='GTiff',
#                 height=height,
#                 width=width,
#                 count=1,
#                 dtype=rasterio.float32,
#                 crs='EPSG:4326',
#                 transform=transform,
#                 compress='lzw',
#                 nodata=np.nan
#             ) as dst:
#                 dst.write(thi_max_array, 1)

#             with rasterio.open(
#                 outfile3,
#                 'w',
#                 driver='GTiff',
#                 height=height,
#                 width=width,
#                 count=len(dates),
#                 dtype=rasterio.float32,
#                 crs='EPSG:4326',
#                 transform=transform,
#                 compress='lzw',
#                 nodata=np.nan
#             ) as dst:
#                 for i in range(thi_daily_array.shape[0]):
#                     dst.write(thi_daily_array[i], i + 1)

# def run_script(start_year, end_year, start_month, end_month, sce_climate):
#     ref_path = 'C:\Users\gkhandag\Downloads\NE1_HR_LC_SR_W\NE1_HR_LC_SR_W.tif'
#     tx_pth = f'F:\\Data\\Tmax'
#     rh_pth = f'F:\\Data\\RHum'
#     out_dir = 'F:\\Data\\THI\\daily'
#     calc_thi(start_year, end_year, start_month, end_month, tx_pth, rh_pth, out_dir, ref_path)

# if __name__ == "__main__":
#     start_year = 2000
#     end_year = 2016
#     start_month = 1
#     end_month = 12
#     sce_climate = "historical"
#     run_script(start_year, end_year, start_month, end_month, sce_climate)

# Daily Code----------------------------------------------

import os
import numpy as np
import rasterio
from calendar import monthrange
import pandas as pd

def thr_hum_idx(tmax, rhum):
    # Calculate THI
    thi = (1.8 * tmax + 32) - ((0.55 - 0.0055 * rhum) * (1.8 * tmax - 26.8))
    return thi

# def calc_thi(start_year, end_year, start_month, end_month, tx_pth, rh_pth, out_dir, ref_path, top=70, left=-180, right=180, bottom=-60):
#     # Loop over years and months
#     for year in range(start_year, end_year + 1):
#         for month in range(start_month, end_month + 1):
#             last_day = monthrange(year, month)[1]
#             dates = pd.date_range(start=f'{year}-{month:02d}-01', end=f'{year}-{month:02d}-{last_day}', freq='D')
            
#             tx_files = [os.path.join(tx_pth, f'Tmax.{date.strftime("%Y.%m.%d")}.tif') for date in dates]
#             rh_files = [os.path.join(rh_pth, f'RH.{date.strftime("%Y.%m.%d")}.tif') for date in dates]

#             for date, tx_file, rh_file in zip(dates, tx_files, rh_files):
#                 outfile = os.path.join(out_dir, f'THI_daily-{date.strftime("%Y-%m-%d")}.tif')
#                 os.makedirs(os.path.dirname(outfile), exist_ok=True)

#                 with rasterio.open(tx_file) as tx_src, rasterio.open(rh_file) as rh_src:
#                     tx_data = tx_src.read(1)
#                     rh_data = rh_src.read(1)
#                     thi_data = thr_hum_idx(tx_data, rh_data)

#                     height, width = thi_data.shape
#                     transform = rasterio.transform.from_bounds(left, bottom, right, top, width, height)

#                     with rasterio.open(
#                         outfile,
#                         'w',
#                         driver='GTiff',
#                         height=height,
#                         width=width,
#                         count=1,
#                         dtype=rasterio.float32,
#                         crs='EPSG:4326',
#                         transform=transform,
#                         compress='lzw',
#                         nodata=np.nan
#                     ) as dst:
#                         dst.write(thi_data, 1)
# This will save all the days file and not just 1f file for a single month
def calc_thi(start_year, end_year, start_month, end_month, tx_pth, rh_pth, out_dir, ref_path, top=70, left=-180, right=180, bottom=-60):
    # Loop over years and months
    for year in range(start_year, end_year + 1):
        for month in range(start_month, end_month + 1):
            last_day = monthrange(year, month)[1]
            dates = pd.date_range(start=f'{year}-{month:02d}-01', end=f'{year}-{month:02d}-{last_day}', freq='D')
            
            tx_files = [os.path.join(tx_pth, f'Tmax.{date.strftime("%Y.%m.%d")}.tif') for date in dates]
            rh_files = [os.path.join(rh_pth, f'RH.{date.strftime("%Y.%m.%d")}.tif') for date in dates]

            for date, tx_file, rh_file in zip(dates, tx_files, rh_files):
                outfile = os.path.join(out_dir, f'THI_daily-{date.strftime("%Y-%m-%d")}.tif')
                os.makedirs(os.path.dirname(outfile), exist_ok=True)

                with rasterio.open(tx_file) as tx_src, rasterio.open(rh_file) as rh_src:
                    tx_data = tx_src.read(1)
                    rh_data = rh_src.read(1)
                    thi_data = thr_hum_idx(tx_data, rh_data)

                    height, width = thi_data.shape
                    transform = rasterio.transform.from_bounds(left, bottom, right, top, width, height)

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
                        dst.write(thi_data, 1)





def run_script(start_year, end_year, start_month, end_month, sce_climate):
    ref_path = 'C:/Users/gkhandag/Downloads/NE1_HR_LC_SR_W/NE1_HR_LC_SR_W.tif'
    tx_pth = f'F:\\Data\\Tmax'
    rh_pth = f'F:\\Data\\RHum'
    out_dir = 'F:\\Data\\THI\\daily'
    calc_thi(start_year, end_year, start_month, end_month, tx_pth, rh_pth, out_dir, ref_path)

if __name__ == "__main__":
    start_year = 2013
    end_year = 2016
    start_month = 1
    end_month = 12
    sce_climate = "historical"
    run_script(start_year, end_year, start_month, end_month, sce_climate)
