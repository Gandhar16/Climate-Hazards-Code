# import os
# import numpy as np
# import rasterio
# import pandas as pd
# from datetime import datetime
# from tkinter import Tk, Label, Entry, Button, filedialog, messagebox
# import matplotlib.pyplot as plt
# from tkinter import messagebox
# from rasterio.plot import show  



# def calc_hsm(yr, mn, tx_dir, r_ref_path, out_dir, thr=35, allyear=False):
#     # File naming
#     season = "AllYear" if allyear else "GSeason"
#     outfile = os.path.join(out_dir, f"HSM_NTx{thr}", f"{season}_HSM_NTx{thr}-{yr}-{mn}.tif")

#     if not os.path.exists(outfile):
#         os.makedirs(os.path.dirname(outfile), exist_ok=True)

#         last_day = pd.Period(f"{yr}-{mn}").days_in_month
#         dates = pd.date_range(start=f"{yr}-{mn}-01", periods=last_day)

#         fls = [os.path.join(tx_dir, f"Tmax.{date.strftime('%Y-%m-%d').replace('-', '.')}.tif") for date in dates]
#         fls = [f for f in fls if os.path.exists(f)]

#         with rasterio.open(r_ref_path) as ref_raster:
#             tmx = [rasterio.open(f).read(1) for f in fls]
#             tmx_stack = np.stack(tmx, axis=0)
#             heat_stress = np.where(tmx_stack >= thr, 1, 0)
#             hsm_sum = heat_stress.sum(axis=0)

#             profile = ref_raster.profile
#             profile.update(dtype=rasterio.uint8, count=1, compress='lzw')
#             with rasterio.open(outfile, 'w', **profile) as dst:
#                 dst.write(hsm_sum.astype(rasterio.uint8), 1)

#         # Visualize and save the map
#         visualize_and_save_map(outfile, out_dir, yr, mn)

#         messagebox.showinfo("Success", f"File saved: {outfile}")
    
# def visualize_from_tif():
#     try:
#         tif_file = filedialog.askopenfilename(title="Select GeoTIFF File", filetypes=[("GeoTIFF files", "*.tif")])
#         if tif_file:
#             out_dir = os.path.dirname(tif_file)
#             year, month = extract_date_from_filename(tif_file)
#             visualize_and_save_map(tif_file, out_dir, year, month)
#             messagebox.showinfo("Visualization", "Visualization completed successfully.")
#     except Exception as e:
#         messagebox.showerror("Error", f"An error occurred: {e}")

# def extract_date_from_filename(filename):
#     # Example filename: "HSM_NTx35_GSeason_HSM_NTx35-2024-01.tif"
#     # Extract year and month from the filename
#     parts = os.path.basename(filename).split('-')
#     year = parts[-2]
#     month = parts[-1].split('.')[0]
#     return year, month


# def visualize_and_save_map(geotiff_path, out_dir, yr_start, yr_end):
#     # Update this function to handle the new date range format
#     with rasterio.open(geotiff_path) as src:
#         fig, ax = plt.subplots(figsize=(10, 10))
#         show(src, ax=ax, title=f"HSM Map for {yr_start}-{yr_end}")
#         map_outfile = os.path.join(out_dir, f"HSM_Map_{yr_start}-{yr_end}.png")
#         plt.savefig(map_outfile)
#         plt.close(fig)
#         messagebox.showinfo("Map Visualization", f"Map image saved: {map_outfile}")


# def run_calculation():
#     try:
#         yr_start = int(year_start_entry.get())
#         yr_end = int(year_end_entry.get())
#         mn_start = int(month_start_entry.get())
#         mn_end = int(month_end_entry.get())
#         tx_dir = tx_path_entry.get()
#         r_ref_path = ref_path_entry.get()
#         out_dir = out_path_entry.get()

#         # New function to handle the calculation for the entire range
#         calc_hsm_range(yr_start, yr_end, mn_start, mn_end, tx_dir, r_ref_path, out_dir)
#     except Exception as e:
#         messagebox.showerror("Error", f"An error occurred: {e}")

# def calc_hsm_range(yr_start, yr_end, mn_start, mn_end, tx_dir, r_ref_path, out_dir, thr=35):
#     outfile = os.path.join(out_dir, f"HSM_NTx{thr}", f"HSM_NTx{thr}_{yr_start}-{yr_end}_{mn_start}-{mn_end}.tif")
    
#     if not os.path.exists(outfile):
#         os.makedirs(os.path.dirname(outfile), exist_ok=True)

#         # Create a date range for the entire period
#         start_date = f"{yr_start}-{str(mn_start).zfill(2)}-01"
#         end_month_days = pd.Period(f'{yr_end}-{mn_end}').days_in_month
#         end_date = f"{yr_end}-{str(mn_end).zfill(2)}-{end_month_days}"
#         dates = pd.date_range(start=start_date, end=end_date)

#         fls = [os.path.join(tx_dir, f"Tmax.{date.strftime('%Y-%m-%d').replace('-', '.')}.tif") for date in dates]
#         fls = [f for f in fls if os.path.exists(f)]

#         if not fls:
#             raise FileNotFoundError("No temperature files found for the specified date range.")

#         with rasterio.open(r_ref_path) as ref_raster:
#             tmx_stack = np.array([rasterio.open(f).read(1) for f in fls])
#             heat_stress = np.where(tmx_stack >= thr, 1, 0)
#             hsm_sum = heat_stress.sum(axis=0)

#             profile = ref_raster.profile
#             profile.update(dtype=rasterio.uint8, count=1, compress='lzw')
#             with rasterio.open(outfile, 'w', **profile) as dst:
#                 dst.write(hsm_sum.astype(rasterio.uint8), 1)

#         # Visualize and save the map for the entire range
#         visualize_and_save_map(outfile, out_dir, yr_start, yr_end)

#         messagebox.showinfo("Success", f"File saved: {outfile}")
#     else:
#         messagebox.showinfo("File already exists", f"The file already exists: {outfile}")


# def select_directory(entry):
#     directory = filedialog.askdirectory()
#     entry.delete(0, 'end')
#     entry.insert(0, directory)

# root = Tk()
# root.title("HSM Calculator")

# Label(root, text="Year Start:").grid(row=0, column=0)
# year_start_entry = Entry(root)
# year_start_entry.grid(row=0, column=1)

# Label(root, text="Year End:").grid(row=1, column=0)
# year_end_entry = Entry(root)
# year_end_entry.grid(row=1, column=1)

# Label(root, text="Month Start:").grid(row=2, column=0)
# month_start_entry = Entry(root)
# month_start_entry.grid(row=2, column=1)

# Label(root, text="Month End:").grid(row=3, column=0)
# month_end_entry = Entry(root)
# month_end_entry.grid(row=3, column=1)

# Label(root, text="Tmax Directory:").grid(row=4, column=0)
# tx_path_entry = Entry(root)
# tx_path_entry.grid(row=4, column=1)
# Button(root, text="Browse", command=lambda: select_directory(tx_path_entry)).grid(row=4, column=2)

# Label(root, text="Reference Path:").grid(row=5, column=0)
# ref_path_entry = Entry(root)
# ref_path_entry.grid(row=5, column=1)
# Button(root, text="Browse", command=lambda: select_directory(ref_path_entry)).grid(row=5, column=2)

# Label(root, text="Output Directory:").grid(row=6, column=0)
# out_path_entry = Entry(root)
# out_path_entry.grid(row=6, column=1)
# Button(root, text="Browse", command=lambda: select_directory(out_path_entry)).grid(row=6, column=2)

# Button(root, text="Run Calculation", command=run_calculation).grid(row=7, column=0, columnspan=3)
# Button(root, text="Visualize from TIF", command=visualize_from_tif).grid(row=8, column=0, columnspan=3)


# root.mainloop()



# Trial Code

# import os
# import numpy as np
# import rasterio
# import pandas as pd
# from datetime import datetime
# from rasterio.plot import show  

# def calc_hsm(yr, mn, tx_dir, r_ref_path, out_dir, thr=35, allyear=False):
#     # File naming
#     season = "AllYear" if allyear else "GSeason"
#     out_subdir = os.path.join(out_dir, f"HSM_NTx{thr}", season)
#     os.makedirs(out_subdir, exist_ok=True)

#     daily_outdir = os.path.join(out_subdir, "Daily_HSM")
#     os.makedirs(daily_outdir, exist_ok=True)

#     last_day = pd.Period(f"{yr}-{mn}").days_in_month
#     dates = pd.date_range(start=f"{yr}-{mn}-01", periods=last_day)

#     monthly_hsm = None

#     for date in dates:
#         day_str = date.strftime('%Y-%m-%d')
#         outfile = os.path.join(daily_outdir, f"HSM_NTx{thr}", f"{season}_HSM_NTx{thr}-{day_str}.tif")

#         if not os.path.exists(outfile):
#             fl_path = os.path.join(tx_dir, f"Tmax.{day_str.replace('-', '.')}.tif")
#             if os.path.exists(fl_path):
#                 with rasterio.open(fl_path) as src:
#                     tmx = src.read(1)
#                     heat_stress = np.where(tmx >= thr, 1, 0)

#                     profile = src.profile
#                     profile.update(dtype=rasterio.uint8, count=1, compress='lzw')
#                     with rasterio.open(outfile, 'w', **profile) as dst:
#                         dst.write(heat_stress.astype(rasterio.uint8), 1)

#                 print(f"Daily file saved: {outfile}")

#                 # Aggregate daily HSM to monthly HSM
#                 if monthly_hsm is None:
#                     monthly_hsm = heat_stress.astype(np.uint8)
#                 else:
#                     monthly_hsm += heat_stress

#             else:
#                 print(f"No temperature data found for {day_str}.")
#         else:
#             print(f"The daily file already exists: {outfile}")

#     # Save monthly HSM
#     if monthly_hsm is not None:
#         monthly_outfile = os.path.join(out_subdir, f"HSM_NTx{thr}", f"{season}_HSM_NTx{thr}-{yr}-{mn}.tif")
#         with rasterio.open(fl_path) as src:
#             profile = src.profile
#             profile.update(dtype=rasterio.uint8, count=1, compress='lzw')
#             with rasterio.open(monthly_outfile, 'w', **profile) as dst:
#                 dst.write(monthly_hsm.astype(rasterio.uint8), 1)

#         print(f"Monthly file saved: {monthly_outfile}")

# def run_calculation(yr_start, yr_end, mn_start, mn_end, tx_dir, r_ref_path, out_dir):
#     for yr in range(yr_start, yr_end + 1):
#         for mn in range(mn_start, mn_end + 1):
#             calc_hsm(yr, mn, tx_dir, r_ref_path, out_dir)

# # Example usage
# yr_start = 1983
# yr_end = 2016
# mn_start = 1
# mn_end = 1
# tx_dir = 'F:/Data/Tmax'
# r_ref_path = 'F:/Data/roi/world.tif'
# out_dir = 'F:/Data/HSM'

# run_calculation(yr_start, yr_end, mn_start, mn_end, tx_dir, r_ref_path, out_dir)



# Modified Trial with threshold
import os
import numpy as np
import rasterio
import pandas as pd
from datetime import datetime
from rasterio.plot import show  

def calc_hsm(yr, mn, tx_dir, r_ref_path, out_dir, thr=0, allyear=False):
    # File naming
    season = "AllYear" if allyear else "GSeason"
    out_subdir = os.path.join(out_dir, f"HSM_NTx{thr}", season)
    os.makedirs(out_subdir, exist_ok=True)

    daily_outdir = os.path.join(out_subdir, "Daily_HSM")
    os.makedirs(daily_outdir, exist_ok=True)

    last_day = pd.Period(f"{yr}-{mn}").days_in_month
    dates = pd.date_range(start=f"{yr}-{mn}-01", periods=last_day)

    monthly_hsm = None

    for date in dates:
        day_str = date.strftime('%Y-%m-%d')
        outfile = os.path.join(daily_outdir, f"{season}_HSM_NTx{thr}-{day_str}.tif")

        if not os.path.exists(outfile):
            fl_path = os.path.join(tx_dir, f"Tmax.{day_str.replace('-', '.')}.tif")
            if os.path.exists(fl_path):
                with rasterio.open(fl_path) as src:
                    tmx = src.read(1)
                    heat_stress = np.where(tmx >= thr, 1, 0)

                    profile = src.profile
                    profile.update(dtype=rasterio.uint8, count=1, compress='lzw')
                    with rasterio.open(outfile, 'w', **profile) as dst:
                        dst.write(heat_stress.astype(rasterio.uint8), 1)

                print(f"Daily file saved: {outfile}")

                # Aggregate daily HSM to monthly HSM
                if monthly_hsm is None:
                    monthly_hsm = heat_stress.astype(np.uint32)  # Use np.uint32 to prevent overflow
                else:
                    monthly_hsm += heat_stress.astype(np.uint32)  # Continue using np.uint32 for addition

            else:
                print(f"No temperature data found for {day_str}.")
        else:
            print(f"The daily file already exists: {outfile}")

    # Save monthly HSM
    if monthly_hsm is not None:
        monthly_outfile = os.path.join(out_subdir, f"{season}_HSM_NTx{thr}-{yr}-{mn}.tif")
        with rasterio.open(fl_path) as src:
            profile = src.profile
            profile.update(dtype=rasterio.uint8, count=1, compress='lzw')
            with rasterio.open(monthly_outfile, 'w', **profile) as dst:
                dst.write(monthly_hsm.astype(rasterio.uint8), 1)  # Convert back to uint8 for saving

        print(f"Monthly file saved: {monthly_outfile}")

def run_calculation(yr_start, yr_end, mn_start, mn_end, tx_dir, r_ref_path, out_dir):
    for yr in range(yr_start, yr_end + 1):
        for mn in range(mn_start, mn_end + 1):
            calc_hsm(yr, mn, tx_dir, r_ref_path, out_dir)

# Example usage
yr_start = 1983
yr_end = 2016
mn_start = 1
mn_end = 12
tx_dir = 'F:/Data/Tmax'
r_ref_path = 'F:/Data/roi/kenya_zambia.tif'
out_dir = 'F:/Data/HSM'


run_calculation(yr_start, yr_end, mn_start, mn_end, tx_dir, r_ref_path, out_dir)





# Modified Trial without threshold
# import os
# import numpy as np
# import rasterio
# import pandas as pd
# from datetime import datetime
# from rasterio.plot import show

# def calc_hsm(yr, mn, tx_dir, r_ref_path, out_dir, allyear=False):
#     # File naming
#     season = "AllYear" if allyear else "GSeason"
#     out_subdir = os.path.join(out_dir, "HSM_NoThr", season)
#     os.makedirs(out_subdir, exist_ok=True)

#     daily_outdir = os.path.join(out_subdir, "Daily_HSM")
#     os.makedirs(daily_outdir, exist_ok=True)

#     last_day = pd.Period(f"{yr}-{mn}").days_in_month
#     dates = pd.date_range(start=f"{yr}-{mn}-01", periods=last_day)

#     monthly_hsm = None

#     for date in dates:
#         day_str = date.strftime('%Y-%m-%d')
#         outfile = os.path.join(daily_outdir, f"{season}_HSM-{day_str}.tif")

#         if not os.path.exists(outfile):
#             fl_path = os.path.join(tx_dir, f"Tmax.{day_str.replace('-', '.')}.tif")
#             if os.path.exists(fl_path):
#                 with rasterio.open(fl_path) as src:
#                     tmx = src.read(1)

#                     # Convert to binary (1 where there is data, 0 otherwise)
#                     heat_stress = np.where(tmx >= 0, 1, 0)

#                     profile = src.profile
#                     profile.update(dtype=rasterio.uint8, count=1, compress='lzw')
#                     with rasterio.open(outfile, 'w', **profile) as dst:
#                         dst.write(heat_stress.astype(rasterio.uint8), 1)

#                 print(f"Daily file saved: {outfile}")

#                 # Aggregate daily HSM to monthly HSM
#                 if monthly_hsm is None:
#                     monthly_hsm = heat_stress.astype(np.uint8)
#                 else:
#                     monthly_hsm = monthly_hsm.astype(np.uint32) + heat_stress.astype(np.uint32)  # Ensure same data type for addition

#             else:
#                 print(f"No temperature data found for {day_str}.")
#         else:
#             print(f"The daily file already exists: {outfile}")

#     # Save monthly HSM
#     if monthly_hsm is not None:
#         monthly_outfile = os.path.join(out_subdir, f"{season}_HSM-{yr}-{mn}.tif")
#         with rasterio.open(fl_path) as src:
#             profile = src.profile
#             profile.update(dtype=rasterio.uint8, count=1, compress='lzw')
#             with rasterio.open(monthly_outfile, 'w', **profile) as dst:
#                 dst.write(monthly_hsm.astype(rasterio.uint8), 1)  # Ensure the data type is uint8 before writing

#         print(f"Monthly file saved: {monthly_outfile}")

# def run_calculation(yr_start, yr_end, mn_start, mn_end, tx_dir, r_ref_path, out_dir):
#     for yr in range(yr_start, yr_end + 1):
#         for mn in range(mn_start, mn_end + 1):
#             calc_hsm(yr, mn, tx_dir, r_ref_path, out_dir, allyear=False)

# # Example usage
# yr_start = 1983
# yr_end = 1983
# mn_start = 2
# mn_end = 2
# tx_dir = 'F:/Data/Tmax'
# r_ref_path = 'F:/Data/roi/kenya_zambia.tif'
# out_dir = 'F:/Data/HSM'

# run_calculation(yr_start, yr_end, mn_start, mn_end, tx_dir, r_ref_path, out_dir)

