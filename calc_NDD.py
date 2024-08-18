# # With Visualization and GUI
import os
import calendar
import numpy as np
import rasterio
from datetime import datetime
import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors



def calc_ndd(start_year, end_year, start_month, end_month, pr_pth, out_dir, ref_raster_path):
    outfile = f"{out_dir}/NDD-{start_year}-{end_year}-{start_month:02d}-{end_month:02d}.tif"
    if not os.path.exists(outfile):
        os.makedirs(os.path.dirname(outfile), exist_ok=True)
        date_list = []
        for year in range(start_year, end_year + 1):
            for month in range(start_month, end_month + 1):
                last_day = calendar.monthrange(year, month)[1]
                date_list.extend([datetime(year, month, day) for day in range(1, last_day + 1)])

        files = [f"{pr_pth}/chirps-v2.0.{date.strftime('%Y.%m.%d')}.tif" for date in date_list]
        existing_files = list(filter(os.path.exists, files))

        with rasterio.open(ref_raster_path) as ref:
            ref_data, _ = ref.read(1), ref.transform
            ndd_array = np.zeros_like(ref_data, dtype=np.int16)
            for file in existing_files:
                with rasterio.open(file) as src:
                    src_data = src.read(1)
                    src_data[src_data == -9999] = np.nan
                    ndd_array += np.nan_to_num(src_data < 1, nan=0).astype(np.int16)
            ref_meta = ref.meta.copy()
            ref_meta.update(dtype=rasterio.int16, count=1)
            with rasterio.open(outfile, 'w', **ref_meta) as dst:
                dst.write(ndd_array, 1)
        print(f"Processed {outfile}")
    else:
        print(f"{outfile} already exists, skipping.")
    return outfile

def visualize_ndd(outfile):
    with rasterio.open(outfile) as src:
        ndd_data = src.read(1)
        plt.imshow(ndd_data, cmap='viridis')
        plt.colorbar(label='Number of Dry Days')
        plt.show()


def visualize_from_tif_file():
    file_path = tif_file_path_entry.get()
    if not file_path or not os.path.exists(file_path):
        messagebox.showerror("Error", "Please select a valid TIF file to visualize.")
        return

    try:
        with rasterio.open(file_path) as src:
            data = src.read(1)  # Ensure only the first band is read
            if data is None or np.all(data == 0):
                raise ValueError("The TIF file contains no data or only zero values.")
            plt.imshow(data, cmap='viridis')
            plt.colorbar(label='Data Values')
            plt.show()
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred while processing the TIF file: {e}")


def select_file(entry):
    file_path = filedialog.askopenfilename()
    entry.delete(0, tk.END)
    entry.insert(0, file_path)

def run_calculation():
    start_year = int(year_start_entry.get())
    end_year = int(year_end_entry.get())
    start_month = int(month_start_entry.get())
    end_month = int(month_end_entry.get())
    ref_raster_file = ref_path_entry.get()
    data_path = data_path_entry.get()
    output_path = save_path_entry.get()

    outfile = calc_ndd(start_year, end_year, start_month, end_month, data_path, output_path, ref_raster_file)
    visualize_ndd(outfile)

def visualize_all_ndd(output_files):
    for file in output_files:
        visualize_ndd(file)


# GUI Setup
root = tk.Tk()
root.title("NDD Calculator")

# Input Fields
tk.Label(root, text="Year Start:").grid(row=0, column=0)
year_start_entry = tk.Entry(root)
year_start_entry.grid(row=0, column=1)

tk.Label(root, text="Year End:").grid(row=1, column=0)
year_end_entry = tk.Entry(root)
year_end_entry.grid(row=1, column=1)

tk.Label(root, text="Month Start:").grid(row=2, column=0)
month_start_entry = tk.Entry(root)
month_start_entry.grid(row=2, column=1)

tk.Label(root, text="Month End:").grid(row=3, column=0)
month_end_entry = tk.Entry(root)
month_end_entry.grid(row=3, column=1)

tk.Label(root, text="Reference Path:").grid(row=4, column=0)
ref_path_entry = tk.Entry(root)
ref_path_entry.grid(row=4, column=1)
tk.Button(root, text="Browse", command=lambda: select_file(ref_path_entry)).grid(row=4, column=2)

tk.Label(root, text="Data Path:").grid(row=5, column=0)
data_path_entry = tk.Entry(root)
data_path_entry.grid(row=5, column=1)
tk.Button(root, text="Browse", command=lambda: select_file(data_path_entry)).grid(row=5, column=2)

# GUI elements for saving the TIF file
tk.Label(root, text="Save As (TIF Path):").grid(row=6, column=0)
save_path_entry = tk.Entry(root)
save_path_entry.grid(row=6, column=1)
tk.Button(root, text="Browse", command=lambda: select_file(save_path_entry)).grid(row=6, column=2)

# GUI elements for visualizing a TIF file
tk.Label(root, text="Visualize TIF File:").grid(row=7, column=0)
tif_file_path_entry = tk.Entry(root)
tif_file_path_entry.grid(row=7, column=1)
tk.Button(root, text="Browse", command=lambda: select_file(tif_file_path_entry)).grid(row=7, column=2)

# Buttons for running calculation, visualization, and saving
tk.Button(root, text="Run Calculation", command=run_calculation).grid(row=8, column=0, columnspan=3)
tk.Button(root, text="Visualize TIF", command=visualize_from_tif_file).grid(row=9, column=0, columnspan=3)

root.mainloop()


# New code without GUI

import os
import calendar
import numpy as np
import rasterio
from datetime import datetime

def calc_ndd(start_year, end_year, start_month, end_month, pr_pth, out_dir, ref_raster_path):
    outfile = f"{out_dir}/NDD-{start_year}-{end_year}-{start_month:02d}-{end_month:02d}.tif"
    if not os.path.exists(outfile):
        os.makedirs(os.path.dirname(outfile), exist_ok=True)
        date_list = []
        for year in range(start_year, end_year + 1):
            for month in range(start_month, end_month + 1):
                last_day = calendar.monthrange(year, month)[1]
                date_list.extend([datetime(year, month, day) for day in range(1, last_day + 1)])

        files = [f"{pr_pth}/chirps-v2.0.{date.strftime('%Y.%m.%d')}.tif" for date in date_list]
        existing_files = list(filter(os.path.exists, files))

        with rasterio.open(ref_raster_path) as ref:
            ref_data, _ = ref.read(1), ref.transform
            ndd_array = np.zeros_like(ref_data, dtype=np.int16)
            for file in existing_files:
                with rasterio.open(file) as src:
                    src_data = src.read(1)
                    src_data[src_data == -9999] = np.nan
                    ndd_array += np.nan_to_num(src_data < 1, nan=0).astype(np.int16)
            ref_meta = ref.meta.copy()
            ref_meta.update(dtype=rasterio.int16, count=1)
            with rasterio.open(outfile, 'w', **ref_meta) as dst:
                dst.write(ndd_array, 1)
        print(f"Processed {outfile}")
    else:
        print(f"{outfile} already exists, skipping.")
    return outfile

# Manually set paths and parameters
start_year = 1983
end_year = 1983
start_month = 1
end_month = 1
pr_pth = 'F:/Data'  # Example: '/data/chirps'
out_dir = 'F:/Data/NDD'  # Example: '/output'
ref_raster_path = 'F:/Data/roi/kenya_zambia.tif'  # Example: '/data/reference.tif'

# Run the NDD calculation
outfile = calc_ndd(start_year, end_year, start_month, end_month, pr_pth, out_dir, ref_raster_path)



# New code:
# import os
# import calendar
# import numpy as np
# import rasterio
# from datetime import datetime

# def calc_ndd_monthly(start_year, end_year, start_month, end_month, pr_pth, out_dir, ref_raster_path):
#     os.makedirs(out_dir, exist_ok=True)  # Ensure output directory exists
#     for year in range(start_year, end_year + 1):
#         for month in range(start_month, end_month + 1):
#             last_day = calendar.monthrange(year, month)[1]
#             date_list = [datetime(year, month, day) for day in range(1, last_day + 1)]
#             files = [f"{pr_pth}/chirps-v2.0.{date.strftime('%Y.%m.%d')}.tif" for date in date_list]
#             existing_files = list(filter(os.path.exists, files))

#             outfile = f"{out_dir}/NDD-{year}-{month:02d}.tif"
#             if not os.path.exists(outfile):  # Only process if file doesn't already exist
#                 with rasterio.open(ref_raster_path) as ref:
#                     ref_data = ref.read(1)
#                     ndd_array = np.zeros_like(ref_data, dtype=np.int16)

#                     for file in existing_files:
#                         with rasterio.open(file) as src:
#                             src_data = src.read(1)
#                             src_data[src_data == -9999] = np.nan  # Handle NoData values
#                             ndd_array += np.nan_to_num(src_data < 1, nan=0).astype(np.int16)

#                     ref_meta = ref.meta.copy()
#                     ref_meta.update(dtype=rasterio.int16, count=1)
#                     with rasterio.open(outfile, 'w', **ref_meta) as dst:
#                         dst.write(ndd_array, 1)
                
#                 print(f"Processed {outfile}")
#             else:
#                 print(f"{outfile} already exists, skipping.")

# # Manually set paths and parameters
# start_year = 1984
# end_year = 2016
# start_month = 1
# end_month = 12
# pr_pth = 'F:/Data'  # Example: '/data/chirps'
# out_dir = 'F:/Data/NDD'  # Example: '/output'
# ref_raster_path = 'F:/Data/roi/kenya_zambia.tif'  # Example: '/data/reference.tif'

# # Run the NDD calculation monthly
# calc_ndd_monthly(start_year, end_year, start_month, end_month, pr_pth, out_dir, ref_raster_path)
