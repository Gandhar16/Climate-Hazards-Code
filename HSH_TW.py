

# # --------------------------------------------------
# import numpy as np
# import pandas as pd
# import rasterio
# from rasterio.mask import mask
# import os
# from shapely.geometry import box
# import geopandas as gpd
# import tkinter as tk
# from tkinter import filedialog, messagebox
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
# from matplotlib.colors import ListedColormap, BoundaryNorm

# # New function for wet bulb temperature calculation
# def wet_bulb_temperature(t, rh):
#     tw = t * np.arctan(0.151977 * np.sqrt(rh + 8.313659)) + np.arctan(t + rh) - \
#          np.arctan(rh - 1.676331) + 0.00391838 * rh**1.5 * np.arctan(0.023101 * rh) - 4.686035
#     return tw

# def clip_raster_data(input_raster_path, ref_raster):
#     if not os.path.exists(input_raster_path):
#         print(f"File not found: {input_raster_path}")
#         return None

#     with rasterio.open(input_raster_path) as src:
#         bounds = ref_raster.bounds
#         geom = box(*bounds)
#         geo = gpd.GeoDataFrame({'geometry': geom}, index=[0], crs=ref_raster.crs)
#         out_image, _ = mask(src, shapes=geo.geometry, crop=True)
#         return out_image[0]

# def calc_hsh(year_start, year_end, month_start, month_end, ref_path, t_pth, rh_pth):
#     tw_avgs = []
#     for year in range(year_start, year_end + 1):
#         for month in range(month_start, month_end + 1):
#             last_day = pd.Period(f"{year}-{month}").days_in_month
#             dts = pd.date_range(start=f"{year}-{month}-01", periods=last_day)

#             t_files = [f"{t_pth}/Tmax.{date.strftime('%Y-%m-%d').replace('-', '.')}.tif" for date in dts]
#             rh_files = [f"{rh_pth}/RH.{date.strftime('%Y-%m-%d').replace('-', '.')}.tif" for date in dts]

#             with rasterio.open(ref_path) as ref_raster:
#                 temp = np.stack([clip_raster_data(f, ref_raster) for f in t_files])
#                 rh = np.stack([clip_raster_data(f, ref_raster) for f in rh_files])

#                 tw = wet_bulb_temperature(temp, rh)
#                 tw_avgs.append(np.nanmean(tw, axis=0))

#     tw_avg_combined = np.nanmean(tw_avgs, axis=0)
#     return tw_avg_combined, ref_raster.transform, ref_raster.crs

# def visualize_heat_stress(data, transform, crs):
#     new_window = tk.Toplevel(root)
#     new_window.title("Heat Stress Map Visualization")

#     colors = ["green", "yellow", "orange", "red"]
#     cmap = ListedColormap(colors)
#     upper_boundary = max(35, np.nanmax(data) + 1)  # Adjust upper boundary to ensure monotonically increasing
#     norm = BoundaryNorm([0, 25, 30, 35, upper_boundary], cmap.N)

#     fig, ax = plt.subplots()
#     cax = ax.imshow(data, cmap=cmap, norm=norm, interpolation='none')
    
#     # Create colorbar with manual label positioning
#     cbar = fig.colorbar(cax, orientation='horizontal', pad=0.14)
#     cbar.set_ticks([12.5, 27.5, 32.5, (35 + upper_boundary) / 2])
#     cbar.ax.set_xticklabels(['Mild', 'Caution', 'Danger', 'Extreme Danger'], horizontalalignment='left')
#     cbar.ax.get_xticklabels()[0].set_horizontalalignment('left')
#     cbar.ax.get_xticklabels()[-1].set_horizontalalignment('right')
#     cbar.ax.xaxis.set_ticks_position('top')
#     cbar.ax.xaxis.set_label_position('top')

#     ax.set_title("Heat Stress Map")
#     ax.set_xlabel("Longitude")
#     ax.set_ylabel("Latitude")
#     ax.grid(True)

#     canvas = FigureCanvasTkAgg(fig, master=new_window)
#     canvas_widget = canvas.get_tk_widget()
#     canvas_widget.pack(fill=tk.BOTH, expand=True)
#     canvas.draw()


# def save_tif_file(data, file_path, transform, crs):
#     with rasterio.open(
#         file_path,
#         'w',
#         driver='GTiff',
#         height=data.shape[0],
#         width=data.shape[1],
#         count=1,
#         dtype=data.dtype,
#         crs=crs,
#         transform=transform,
#     ) as dst:
#         dst.write(data, 1)
#     messagebox.showinfo("Success", "TIF file saved successfully!")

# def run_hsh_calculation():
#     try:
#         year_start = int(year_start_entry.get())
#         year_end = int(year_end_entry.get())
#         month_start = int(month_start_entry.get())
#         month_end = int(month_end_entry.get())
#         ref_path = ref_path_entry.get()
#         t_path = tx_path_entry.get()
#         rh_path = rh_path_entry.get()

#         tw_avg, transform, crs = calc_hsh(year_start, year_end, month_start, month_end, ref_path, t_path, rh_path)
#         visualize_button.config(state='normal', command=lambda: visualize_heat_stress(tw_avg, transform, crs))
#         save_button.config(state='normal', command=lambda: save_tif_file(tw_avg, save_path_entry.get(), transform, crs))
#         messagebox.showinfo("Success", "HSH Calculation Completed Successfully!")
#     except Exception as e:
#         messagebox.showerror("Error", f"An error occurred: {e}")

# def visualize_from_tif_file():
#     file_path = tif_file_path_entry.get()
#     if not file_path:
#         messagebox.showerror("Error", "Please select a TIF file to visualize.")
#         return

#     try:
#         with rasterio.open(file_path) as src:
#             data = src.read(1)
#             transform = src.transform
#             crs = src.crs

#         visualize_heat_stress(data, transform, crs)
#     except Exception as e:
#         messagebox.showerror("Error", f"An error occurred while reading the TIF file: {e}")

# def select_file(entry):
#     file_path = filedialog.askopenfilename()
#     entry.delete(0, tk.END)
#     entry.insert(0, file_path)

# root = tk.Tk()
# root.title("HSH Calculator")

# # GUI setup for inputs
# tk.Label(root, text="Year Start:").grid(row=0, column=0)
# year_start_entry = tk.Entry(root)
# year_start_entry.grid(row=0, column=1)

# tk.Label(root, text="Year End:").grid(row=1, column=0)
# year_end_entry = tk.Entry(root)
# year_end_entry.grid(row=1, column=1)

# tk.Label(root, text="Month Start:").grid(row=2, column=0)
# month_start_entry = tk.Entry(root)
# month_start_entry.grid(row=2, column=1)

# tk.Label(root, text="Month End:").grid(row=3, column=0)
# month_end_entry = tk.Entry(root)
# month_end_entry.grid(row=3, column=1)

# tk.Label(root, text="Reference Path:").grid(row=4, column=0)
# ref_path_entry = tk.Entry(root)
# ref_path_entry.grid(row=4, column=1)
# tk.Button(root, text="Browse", command=lambda: select_file(ref_path_entry)).grid(row=4, column=2)

# tk.Label(root, text="Temperature Path:").grid(row=5, column=0)
# tx_path_entry = tk.Entry(root)
# tx_path_entry.grid(row=5, column=1)
# tk.Button(root, text="Browse", command=lambda: select_file(tx_path_entry)).grid(row=5, column=2)

# tk.Label(root, text="Relative Humidity Path:").grid(row=6, column=0)
# rh_path_entry = tk.Entry(root)
# rh_path_entry.grid(row=6, column=1)
# tk.Button(root, text="Browse", command=lambda: select_file(rh_path_entry)).grid(row=6, column=2)

# # GUI setup for running calculation and visualization
# tk.Button(root, text="Run Calculation", command=run_hsh_calculation).grid(row=8, column=0, columnspan=3)
# visualize_button = tk.Button(root, text="Visualize Heat Stress Map", state='disabled')
# visualize_button.grid(row=9, column=0, columnspan=3)

# # GUI elements for saving the TIF file
# tk.Label(root, text="Save As (TIF Path):").grid(row=10, column=0)
# save_path_entry = tk.Entry(root)
# save_path_entry.grid(row=10, column=1)
# tk.Button(root, text="Browse", command=lambda: select_file(save_path_entry)).grid(row=10, column=2)
# save_button = tk.Button(root, text="Save TIF File", state='disabled')
# save_button.grid(row=11, column=0, columnspan=3)

# # GUI elements for visualizing a TIF file
# tk.Label(root, text="Visualize TIF File:").grid(row=12, column=0)
# tif_file_path_entry = tk.Entry(root)
# tif_file_path_entry.grid(row=12, column=1)
# tk.Button(root, text="Browse", command=lambda: select_file(tif_file_path_entry)).grid(row=12, column=2)
# visualize_tif_button = tk.Button(root, text="Visualize TIF", command=visualize_from_tif_file)
# visualize_tif_button.grid(row=13, column=0, columnspan=3)

# # Start the Tkinter main loop
# root.mainloop()






# -Save Daily COde -----------------------------------------------

import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask
import os
from shapely.geometry import box
import geopandas as gpd
from datetime import datetime, timedelta
from rasterio.transform import from_bounds

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

def calc_hsh_single_day(ref_path, tx_file, tm_file, rh_file, output_dir, top=70, left=-180, right=180, bottom=-60):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with rasterio.open(ref_path) as ref_raster:
        tmx = clip_raster_data(tx_file, ref_raster, top, left, right, bottom)
        tmn = clip_raster_data(tm_file, ref_raster, top, left, right, bottom)
        rhm = clip_raster_data(rh_file, ref_raster, top, left, right, bottom)

        tav = (tmx + tmn) / 2
        hi = wet_bulb_temperature(tav, rhm)

        # Extract date from file name and format it as year-month-day
        date_str = os.path.basename(tx_file).split('.')[1:4]
        date_str = '.'.join(date_str)
        # print(date_str)
        date = datetime.strptime(date_str, '%Y.%m.%d').strftime('%Y-%m-%d')
        
        output_path = os.path.join(output_dir, f"HSH_TW_{date}.tif")
        save_tif_file(hi, output_path, top, left, right, bottom, ref_raster.crs)
        print(f"Date: {date} - Completed")

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

# def save_tif_file(data, file_path, top, left, right, bottom, crs):
#     height, width = data.shape
#     transform = from_bounds(left, bottom, right, top, width, height)

#     print(f"Saving TIF file with extent: top={top}, left={left}, right={right}, bottom={bottom}")
#     print(f"Data shape: height={height}, width={width}")
#     print(f"Transform: {transform}")

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

def main():
    start_date = datetime(2012, 9, 24)
    end_date = datetime(2016, 12, 31)
    base_dir = "F:/Data"
    ref_path = "C:/Users/gkhandag/Downloads/NE1_HR_LC_SR_W/NE1_HR_LC_SR_W.tif"
    output_dir = "F:/Data/HSH_TW/Daily"

    file_paths = generate_file_paths(start_date, end_date, base_dir)

    for tx_file, tm_file, rh_file in file_paths:
        calc_hsh_single_day(ref_path, tx_file, tm_file, rh_file, output_dir)

if __name__ == "__main__":
    main()



