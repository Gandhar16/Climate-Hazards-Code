

import os
import numpy as np
import rasterio
from datetime import datetime
from calendar import monthrange
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from rasterio.merge import merge

# Function definitions for calculate_ptot, calculate_ntx, process_scenario...


def calculate_ptot(start_year, end_year, months, pr_path, out_dir):
    for year in range(start_year, end_year + 1):
        for month in months:
            outfile = os.path.join(out_dir, f'PTOT-{year}-{month:02d}.tif')
            if not os.path.exists(outfile):
                total_precip = None
                days_in_month = monthrange(year, month)[1]
                files = [os.path.join(pr_path, f'chirps-v2.0.{year}.{month:02d}.{day:02d}.tif') for day in range(1, days_in_month + 1)]
                files = [f for f in files if os.path.exists(f)]
                if files:
                    rasters = [rasterio.open(f) for f in files]
                    data, _ = merge(rasters)
                    total_precip = np.sum(data, axis=0)  # Sum over the third dimension, which are the daily data layers
                    
                    # Closing all raster files after processing
                    for raster in rasters:
                        raster.close()

                if total_precip is not None:
                    # Writing the monthly total precipitation data to a file
                    with rasterio.open(files[0]) as src:
                        profile = src.profile
                        profile.update({
                            'count': 1,
                            'dtype': 'float32'
                        })
                        with rasterio.open(outfile, 'w', **profile) as dest:
                            dest.write(total_precip, 1)

# Main code
# def calculate_ptot(start_year, end_year, months, pr_path, out_dir):
#     outfile = os.path.join(out_dir, f'PTOT-{start_year}-{months[0]:02d}_to_{end_year}-{months[-1]:02d}.tif')
#     if not os.path.exists(outfile):
#         total_precip = None
#         for year in range(start_year, end_year + 1):
#             for month in months:
#                 days_in_month = monthrange(year, month)[1]
#                 files = [os.path.join(pr_path, f'chirps-v2.0.{year}.{month:02d}.{day:02d}.tif') for day in range(1, days_in_month + 1)]
#                 files = [f for f in files if os.path.exists(f)]
#                 if files:
#                     rasters = [rasterio.open(f) for f in files]
#                     data, _ = merge(rasters)
#                     if total_precip is None:
#                         total_precip = np.sum(data, axis=0)
#                     else:
#                         total_precip += np.sum(data, axis=0)
#                     for raster in rasters:
#                         raster.close()

#         if total_precip is not None:
#             profile = rasters[0].profile
#             profile.update({'count': 1, 'dtype': 'float32'})
#             with rasterio.open(outfile, 'w', **profile) as dest:
#                 dest.write(total_precip, 1)

# def calculate_ntx(start_year, end_year, months, threshold, tx_path, out_dir):
#     threshold_dir = os.path.join(out_dir, f'NTx{threshold}')
#     os.makedirs(threshold_dir, exist_ok=True)
#     outfile = os.path.join(threshold_dir, f'NTx{threshold}-{start_year}-{months[0]:02d}_to_{end_year}-{months[-1]:02d}.tif')
#     if not os.path.exists(outfile):
#         ntx = None
#         for year in range(start_year, end_year + 1):
#             for month in months:
#                 days_in_month = monthrange(year, month)[1]
#                 files = [os.path.join(tx_path, f'Tmax.{year}.{month:02d}.{day:02d}.tif') for day in range(1, days_in_month + 1)]
#                 files = [f for f in files if os.path.exists(f)]
#                 if files:
#                     rasters = [rasterio.open(f) for f in files]
#                     data, _ = merge(rasters)
#                     monthly_ntx = np.sum(data > threshold, axis=0)
#                     if ntx is None:
#                         ntx = monthly_ntx
#                     else:
#                         ntx += monthly_ntx
#                     for raster in rasters:
#                         raster.close()

#         if ntx is not None:
#             profile = rasters[0].profile
#             profile.update({'count': 1, 'dtype': 'uint32'})
#             with rasterio.open(outfile, 'w', **profile) as dest:
#                 dest.write(ntx, 1)
                
def calculate_ntx(start_year, end_year, months, threshold, tx_path, out_dir):
    threshold_dir = os.path.join(out_dir, f'NTx{threshold}')
    os.makedirs(threshold_dir, exist_ok=True)
    outfile = os.path.join(threshold_dir, f'NTx{threshold}-{start_year}-{months[0]:02d}_to_{end_year}-{months[-1]:02d}.tif')
    ntx_accumulated = None

    for year in range(start_year, end_year + 1):
        for month in months:
            days_in_month = monthrange(year, month)[1]
            files = [os.path.join(tx_path, f'Tmax.{year}.{month:02d}.{day:02d}.tif') for day in range(1, days_in_month + 1)]
            files = [f for f in files if os.path.exists(f)]
            if files:
                rasters = [rasterio.open(f) for f in files]
                data, _ = merge(rasters)
                monthly_ntx = np.sum(data > threshold, axis=0)
                if ntx_accumulated is None:
                    ntx_accumulated = monthly_ntx
                else:
                    ntx_accumulated += monthly_ntx
                for raster in rasters:
                    raster.close()

    if ntx_accumulated is not None:
        profile = rasters[0].profile
        profile.update({'count': 1, 'dtype': 'uint32'})
        with rasterio.open(outfile, 'w', **profile) as dest:
            dest.write(ntx_accumulated, 1)


# def process_scenario(scenario, years, months, pr_path, tx_path, out_dir, thresholds):
#     for year in years:
#         for month in months:
#             # Call the function to calculate total monthly precipitation
#             calculate_ptot(year, month, pr_path, out_dir)

#             # Call the function to calculate the number of days above temperature thresholds
#             for threshold in thresholds:
#                 calculate_ntx(year, month, threshold, tx_path, out_dir)
def process_scenario(scenario, start_year, end_year, months, pr_path, tx_path, out_dir, thresholds):
    # Call the function to calculate total precipitation for the entire duration
    calculate_ptot(start_year, end_year, months, pr_path, out_dir)

    # Call the function to calculate the number of days above temperature thresholds for the entire duration
    for threshold in thresholds:
        calculate_ntx(start_year, end_year, months, threshold, tx_path, out_dir)


def get_directory(entry):
    path = filedialog.askdirectory()
    entry.delete(0, tk.END)
    entry.insert(0, path)


# def run_processing():
#     start_year = int(start_year_entry.get())
#     end_year = int(end_year_entry.get())
#     selected_months = [int(month) for month in months_entry.get().split(',')]
#     threshold = int(threshold_entry.get())
#     pr_path = pr_path_entry.get()
#     tx_path = tx_path_entry.get()
#     out_dir = out_dir_entry.get()

#     years = list(range(start_year, end_year + 1))

#     if scenario_combobox.get() == 'Historical':
#         process_scenario('historical', years, selected_months, pr_path, tx_path, out_dir, [threshold])

#     messagebox.showinfo("Success", "Processing completed successfully.")
#     # Assuming visualization for the first month of the start year
#     visualization_filename = f'PTOT-{start_year}-{selected_months[0]:02d}_to_{end_year}-{selected_months[-1]:02d}.tif'
#     visualize_button.config(command=lambda: visualize_data(out_dir, visualization_filename))

# def run_processing():
#     start_year = int(start_year_entry.get())
#     end_year = int(end_year_entry.get())
#     selected_months = [int(month) for month in months_entry.get().split(',')]
#     threshold = int(threshold_entry.get())  # Single threshold
#     pr_path = pr_path_entry.get()
#     tx_path = tx_path_entry.get()
#     out_dir = out_dir_entry.get()

#     years = list(range(start_year, end_year + 1))

#     if scenario_combobox.get() == 'Historical':
#         # Ensure the 'thresholds' argument is a list, even if it's a single value
#         process_scenario('historical', years, selected_months, pr_path, tx_path, out_dir, [threshold])  # Passing [threshold] as list

#     messagebox.showinfo("Success", "Processing completed successfully.")
#     visualization_filename = f'PTOT-{start_year}-{selected_months[0]:02d}_to_{end_year}-{selected_months[-1]:02d}.tif'
#     visualize_button.config(command=lambda: visualize_data(out_dir, visualization_filename))
    
def run_processing():
    start_year = int(start_year_entry.get())
    end_year = int(end_year_entry.get())
    selected_months = [int(month) for month in months_entry.get().split(',')]  # Convert month entries to a list of integers
    threshold = int(threshold_entry.get())  # Get the threshold as an integer
    pr_path = pr_path_entry.get()
    tx_path = tx_path_entry.get()
    out_dir = out_dir_entry.get()

    years = list(range(start_year, end_year + 1))  # Generate a list of years

    if scenario_combobox.get() == 'Historical':
        process_scenario('historical', start_year, end_year, selected_months, pr_path, tx_path, out_dir, [threshold])  # Pass [threshold] as a list

    messagebox.showinfo("Success", "Processing completed successfully.")
    visualization_filename = f'PTOT-{start_year}-{selected_months[0]:02d}_to_{end_year}-{selected_months[-1]:02d}.tif'
    visualize_button.config(command=lambda: visualize_data(out_dir, visualization_filename))




def visualize_data(directory, filename):
    file_path = os.path.join(directory, filename)
    with rasterio.open(file_path) as src:
        data = src.read(1)
        plt.imshow(data, cmap='viridis')
        plt.colorbar(label='Data')
        plt.show()


# GUI Setup
window = tk.Tk()
window.title("Climate Data Processing Tool")

scenario_combobox = ttk.Combobox(window, values=["Historical", "Future"], state="readonly")
scenario_combobox.pack()
scenario_combobox.set("Historical")

# Label and Entry for Year
# year_label = tk.Label(window, text="Enter Year:")
# year_label.pack()
# year_entry = tk.Entry(window)
# year_entry.pack()

start_year_label = tk.Label(window, text="Enter Start Year:")
start_year_label.pack()
start_year_entry = tk.Entry(window)
start_year_entry.pack()

# Label and Entry for End Year
end_year_label = tk.Label(window, text="Enter End Year:")
end_year_label.pack()
end_year_entry = tk.Entry(window)
end_year_entry.pack()

# Label and Entry for Months
months_label = tk.Label(window, text="Enter Months (comma-separated):")
months_label.pack()
months_entry = tk.Entry(window)
months_entry.pack()

# # Label and Entry for Month
# month_label = tk.Label(window, text="Enter Month:")
# month_label.pack()
# month_entry = tk.Entry(window)
# month_entry.pack()

# Label and Entry for Threshold
threshold_label = tk.Label(window, text="Enter Threshold:")
threshold_label.pack()
threshold_entry = tk.Entry(window)
threshold_entry.pack()

pr_path_entry = tk.Entry(window, width=50)
pr_path_entry.pack()
pr_path_button = tk.Button(window, text="Select PR Path", command=lambda: get_directory(pr_path_entry))
pr_path_button.pack()

tx_path_entry = tk.Entry(window, width=50)
tx_path_entry.pack()
tx_path_button = tk.Button(window, text="Select TX Path", command=lambda: get_directory(tx_path_entry))
tx_path_button.pack()

out_dir_entry = tk.Entry(window, width=50)
out_dir_entry.pack()
out_dir_button = tk.Button(window, text="Select Output Directory", command=lambda: get_directory(out_dir_entry))
out_dir_button.pack()

process_button = tk.Button(window, text="Run Processing", command=run_processing)
process_button.pack()

visualize_button = tk.Button(window, text="Visualize Data")
visualize_button.pack()

window.mainloop()

