# Comparing air temperature values between MERRA-2 reanalysis and WMO observations.
# The MERRA-2 air temperature values looked way high for Nuwara Eliya.
#   - How much do the Nuwara Eliya MERRA-2 temps disagree with WMO observations?
#   - Is there also notable error for Jaffna and Negombo?
#   - Is there a simple bias correction we can apply to the MERRA-2 data (e.g., using lapse rate)?


## Import libraries
import xarray as xr
import pandas as pd
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import os


# Read in MERRA-2 data
# (do this outside loop because it takes forever)

# MERRA-2 reanalysis (hourly data)
fileDir_MERRA  = "/mnt/redwood/local_drive/chanud/RawData/MERRA_1980-2021_SL/MERRA_hourly_SLV/"
fPath_MERRA  = fileDir_MERRA + 'MERRA2_tavg1_2d_slv_Nx.19800101_20211231.SUB.T2M.nc'
# Open the NetCDF file using xarray
data_MERRA_global = xr.open_dataset(fPath_MERRA)



cities = ["Negombo", "NuwaraEliya", "Jaffna"]

for city in cities:

    # Get location specific values
    # (lat/lon/elev are same as in WHATCH'EM code, run_container.pl)
    if city == "Negombo":
        lat  = 7.2008
        lon  = 79.8737
        elev = 2.
        fName_obsWMO = "Katunayake_WMO_2001-2020.csv"
    elif city == "NuwaraEliya":
        lat  = 6.9497
        lon  = 80.7891
        elev = 1868.
        fName_obsWMO = "NuwaraEliya_WMO_2001-2020.csv"
    elif city == "Jaffna":
        lat  = 9.6615
        lon  = 80.0255
        elev = 5.
        fName_obsWMO = "Jaffna_WMO_2009-2020.csv"

    ## Open the files and get temp data

    # WMO observations (daily data)
    fileDir_obs    = "/mnt/redwood/local_drive/chanud/RawData/observations/"
    fPath_obsWMO = fileDir_obs + fName_obsWMO
    # Read the CSV file into a pandas DataFrame
    data_WMO = pd.read_csv(fPath_obsWMO)
    data_WMO.set_index(pd.to_datetime(data_WMO["DATE"]), inplace=True)
    # Replace the missing value for temperatures (9999.9) with nan
    data_WMO["TEMP"].replace(9999.9, np.nan, inplace=True)
    data_WMO["MAX"].replace(9999.9, np.nan, inplace=True)
    data_WMO["MIN"].replace(9999.9, np.nan, inplace=True)
    # Get the daily average, max, and min temps in Celsius
    Tavg_WMO = (data_WMO["TEMP"] - 32) * (5/9) # convert from F to C
    Tmax_WMO = (data_WMO["MAX"] - 32) * (5/9)  # convert from F to C
    Tmin_WMO = (data_WMO["MIN"] - 32) * (5/9)  # convert from F to C
    
    # Extract MERRA-2 data for location of interest as a Pandas DataFrame
    data_MERRA = data_MERRA_global.sel(lat=lat, lon=lon, method='nearest')
    data_MERRA = data_MERRA.drop(['lat', 'lon'])
    data_MERRA = data_MERRA.to_dataframe()
    T2M_MERRA = data_MERRA["T2M"] - 273.15 # convert from K to C


    # Plot all timeseries
    if False:
        plt.figure(figsize=(10, 6))
        plt.plot(Tavg_WMO.index, Tavg_WMO.values, lw=1, ls='-', c='black', alpha=1, label='WMO avg')
        plt.plot(Tmax_WMO.index, Tmax_WMO.values, lw=1, ls='-', c='red', alpha=1, label='WMO max')
        plt.plot(Tmin_WMO.index, Tmin_WMO.values, lw=1, ls='-', c='blue', alpha=1, label='WMO min')
        plt.plot(T2M_MERRA.index, T2M_MERRA.values, lw=0.5, ls='-', c='gray', alpha=0.8, label='MERRA-2')

        plt.title(city)
        plt.xlabel('Date')
        plt.ylabel('Temperature [C]')
        if city == "Jaffna":
            plt.xlim(pd.Timestamp('2009-08-01'), pd.Timestamp('2020-12-31'))
        else:
            plt.xlim(pd.Timestamp('2001-01-01'), pd.Timestamp('2020-12-31'))
        plt.ylim(5, 40)
        plt.legend()
        plt.grid(True)

        # Save the plot
        dir_out = "/mnt/redwood/local_drive/chanud/Output/test/"
        if not os.path.isdir(dir_out):
            os.makedirs(dir_out)
        filepath_out = dir_out + "T_MERRA2_vs_obsWMO_" + city + ".png"
        plt.savefig(filepath_out, dpi=300, bbox_inches="tight")

    # Plot data as histograms
    if False:
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 12))

        # Plot histograms in the first subplot with transparency (alpha)
        colors = ['black', 'red', 'blue']
        labels = ['Tavg', 'Tmax', 'Tmin']
        axes[0].hist([Tavg_WMO, Tmax_WMO, Tmin_WMO], bins=20, edgecolor='black', color=colors, label=labels, alpha=0.7)
        axes[0].set_title('Observations (WMO)')
        axes[0].set_xlabel('Temperature [C]')
        axes[0].set_ylabel('Frequency')
        axes[0].legend()
        axes[0].set_xlim(0, 40)
                
        # Plot the fourth series in the second subplot
        axes[1].hist(T2M_MERRA, bins=20, edgecolor='black', color='gray')
        axes[1].set_title('Reanalysis (MERRA-2)')
        axes[1].set_xlabel('Temperature [C]')
        axes[1].set_ylabel('Frequency')
        axes[1].set_xlim(0, 40)

        # Adjust layout
        plt.tight_layout()

        # Save the plot
        dir_out = "/mnt/redwood/local_drive/chanud/Output/test/"
        if not os.path.isdir(dir_out):
            os.makedirs(dir_out)
        filepath_out = dir_out + "T_MERRA2_vs_obsWMO_hist_" + city + ".png"
        plt.savefig(filepath_out, dpi=300, bbox_inches="tight")



    # Bias correction
    # Do linear regression of MERRA-2 avg/max/min temps with WMO avg/max/min temps.
    # Bias correct the MERRA-2 data based on Tavg (daily average temperature).
    if False:
        # Calculate daily average, max, and min of MERRA temp
        T2Mavg_MERRA = T2M_MERRA.resample('D').mean()
        T2Mmin_MERRA = T2M_MERRA.resample('D').min()
        T2Mmax_MERRA = T2M_MERRA.resample('D').max()

        # Subset MERRA-2 and WMO data to just the common dates for which they both have non-missing data.
        common_dates_avg = T2Mavg_MERRA.index.intersection(Tavg_WMO.index)
        common_dates_avg = common_dates_avg[(~T2Mavg_MERRA.loc[common_dates_avg].isna()) & (~Tavg_WMO.loc[common_dates_avg].isna())]
        common_dates_min = T2Mmin_MERRA.index.intersection(Tmin_WMO.index)
        common_dates_min = common_dates_min[(~T2Mmin_MERRA.loc[common_dates_min].isna()) & (~Tmin_WMO.loc[common_dates_min].isna())]
        common_dates_max = T2Mmax_MERRA.index.intersection(Tmax_WMO.index)
        common_dates_max = common_dates_max[(~T2Mmax_MERRA.loc[common_dates_max].isna()) & (~Tmax_WMO.loc[common_dates_max].isna())]

        Tavg_WMO_subset = Tavg_WMO.loc[common_dates_avg]
        Tmin_WMO_subset = Tmin_WMO.loc[common_dates_min]
        Tmax_WMO_subset = Tmax_WMO.loc[common_dates_max]
        T2Mavg_MERRA_subset = T2Mavg_MERRA.loc[common_dates_avg]
        T2Mmin_MERRA_subset = T2Mmin_MERRA.loc[common_dates_min]
        T2Mmax_MERRA_subset = T2Mmax_MERRA.loc[common_dates_max]

        print("Starting linear regression...")
        # Perform linear regression
        slope_Tavg, int_Tavg, r_Tavg, _, _ = linregress(Tavg_WMO_subset.values, T2Mavg_MERRA_subset.values)
        slope_Tmin, int_Tmin, r_Tmin, _, _ = linregress(Tmin_WMO_subset.values, T2Mmin_MERRA_subset.values)
        slope_Tmax, int_Tmax, r_Tmax, _, _ = linregress(Tmax_WMO_subset.values, T2Mmax_MERRA_subset.values)

        # Function to annotate plot with regression equation and R-squared
        def annotate_regression(ax, slope, intercept, r_value):
            equation = f'T_MERRA = {slope:.2f}*T_WMO + {intercept:.2f}'
            r_squared = f'R² = {r_value**2:.4f}'
            annotation_text = f'{equation}\n{r_squared}'
            ax.annotate(annotation_text, xy=(0.100, 0.020), xycoords='axes fraction', fontsize=8, c="red")

        print("Plotting...")
        # Plot data and lines of best fit
        fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))

        axes[0].scatter(Tavg_WMO_subset, T2Mavg_MERRA_subset, label='Tavg')
        axes[0].plot(Tavg_WMO_subset, slope_Tavg * Tavg_WMO_subset + int_Tavg, color='red', label='Regression Line')
        # Add a 1-to-1 line
        min_val_4plt = min(Tavg_WMO_subset.min(), T2Mavg_MERRA_subset.min())
        max_val_4plt = max(Tavg_WMO_subset.max(), T2Mavg_MERRA_subset.max())
        axes[0].plot([min_val_4plt, max_val_4plt], [min_val_4plt, max_val_4plt], color='black', linestyle='--', label='1-to-1 Line')
        # Add annotation with regression info
        annotate_regression(axes[0], slope_Tavg, int_Tavg, r_Tavg)
        axes[0].set_title('Tavg')
        axes[0].set_xlabel('T_WMO')
        axes[0].set_ylabel('T_MERRA')
        axes[0].legend()
        # Set x and y limits to match the minimum and maximum values
        min_val_4plt = min(Tavg_WMO_subset.min(), T2Mavg_MERRA_subset.min())
        max_val_4plt = max(Tavg_WMO_subset.max(), T2Mavg_MERRA_subset.max())
        axes[0].set_xlim(min_val_4plt, max_val_4plt)
        axes[0].set_ylim(min_val_4plt, max_val_4plt)

        axes[1].scatter(Tmin_WMO_subset, T2Mmin_MERRA_subset, label='Tmin')
        axes[1].plot(Tmin_WMO_subset, slope_Tmin * Tmin_WMO_subset + int_Tmin, color='red', label='Regression Line')
        # Add a 1-to-1 line
        min_val_4plt = min(Tmin_WMO_subset.min(), T2Mmin_MERRA_subset.min())
        max_val_4plt = max(Tmin_WMO_subset.max(), T2Mmin_MERRA_subset.max())
        axes[1].plot([min_val_4plt, max_val_4plt], [min_val_4plt, max_val_4plt], color='black', linestyle='--', label='1-to-1 Line')
        # Add annotation with regression info
        annotate_regression(axes[1], slope_Tmin, int_Tmin, r_Tmin)
        axes[1].set_title('Tmin')
        axes[1].set_xlabel('T_WMO')
        axes[1].set_ylabel('T_MERRA')
        axes[1].legend()
        # Set x and y limits to match the minimum and maximum values
        axes[1].set_xlim(min_val_4plt, max_val_4plt)
        axes[1].set_ylim(min_val_4plt, max_val_4plt)

        axes[2].scatter(Tmax_WMO_subset, T2Mmax_MERRA_subset, label='Tmax')
        axes[2].plot(Tmax_WMO_subset, slope_Tmax * Tmax_WMO_subset + int_Tmax, color='red', label='Regression Line')
        # Add a 1-to-1 line
        min_val_4plt = min(Tmax_WMO_subset.min(), T2Mmax_MERRA_subset.min())
        max_val_4plt = max(Tmax_WMO_subset.max(), T2Mmax_MERRA_subset.max())
        axes[2].plot([min_val_4plt, max_val_4plt], [min_val_4plt, max_val_4plt], color='black', linestyle='--', label='1-to-1 Line')
        # Add annotation with regression info
        annotate_regression(axes[2], slope_Tmax, int_Tmax, r_Tmax)
        axes[2].set_title('Tmax')
        axes[2].set_xlabel('T_WMO')
        axes[2].set_ylabel('T_MERRA')
        axes[2].legend()
        # Set x and y limits to match the minimum and maximum values
        axes[2].set_xlim(min_val_4plt, max_val_4plt)
        axes[2].set_ylim(min_val_4plt, max_val_4plt)

        plt.tight_layout()

        print("Saving plot...")
        # Save the plot
        dir_out = "/mnt/redwood/local_drive/chanud/Output/test/"
        if not os.path.isdir(dir_out):
            os.makedirs(dir_out)
        filepath_out = dir_out + "T_MERRA2_vs_obsWMO_" + city + "_linReg" + ".png"
        plt.savefig(filepath_out, dpi=300, bbox_inches="tight")


        ##Bias correct MERRA-2 data and make same plots with this data
        print("Time for bias correction!")
        T2M_MERRA_bc = (T2M_MERRA - int_Tavg)/slope_Tavg

        # Plot timeseries
        plt.figure(figsize=(10, 6))
        plt.plot(Tavg_WMO.index, Tavg_WMO.values, lw=1, ls='-', c='black', alpha=1, label='WMO avg')
        plt.plot(Tmax_WMO.index, Tmax_WMO.values, lw=1, ls='-', c='red', alpha=1, label='WMO max')
        plt.plot(Tmin_WMO.index, Tmin_WMO.values, lw=1, ls='-', c='blue', alpha=1, label='WMO min')
        plt.plot(T2M_MERRA_bc.index, T2M_MERRA_bc.values, lw=0.5, ls='-', c='gray', alpha=0.8, label='MERRA-2 (bias corrected)')

        plt.title(city)
        plt.xlabel('Date')
        plt.ylabel('Temperature [C]')
        if city == "Jaffna":
            plt.xlim(pd.Timestamp('2009-08-01'), pd.Timestamp('2020-12-31'))
        else:
            plt.xlim(pd.Timestamp('2001-01-01'), pd.Timestamp('2020-12-31'))
        plt.ylim(5, 40)
        plt.legend()
        plt.grid(True)

        # Save the plot
        dir_out = "/mnt/redwood/local_drive/chanud/Output/test/"
        if not os.path.isdir(dir_out):
            os.makedirs(dir_out)
        filepath_out = dir_out + "T_MERRA2_bc_vs_obsWMO_" + city + ".png"
        plt.savefig(filepath_out, dpi=300, bbox_inches="tight")


        # Calculate daily average, max, and min of MERRA temp
        T2Mavg_MERRA_bc = T2M_MERRA_bc.resample('D').mean()
        T2Mmin_MERRA_bc = T2M_MERRA_bc.resample('D').min()
        T2Mmax_MERRA_bc = T2M_MERRA_bc.resample('D').max()

        # Subset to just the common dates for which they both have non-missing data.
        T2Mavg_MERRA_bc_subset = T2Mavg_MERRA_bc.loc[common_dates_avg]
        T2Mmin_MERRA_bc_subset = T2Mmin_MERRA_bc.loc[common_dates_min]
        T2Mmax_MERRA_bc_subset = T2Mmax_MERRA_bc.loc[common_dates_max]

        print("Starting linear regression...")
        # Perform linear regression
        slope_Tavg, int_Tavg, r_Tavg, _, _ = linregress(Tavg_WMO_subset.values, T2Mavg_MERRA_bc_subset.values)
        slope_Tmin, int_Tmin, r_Tmin, _, _ = linregress(Tmin_WMO_subset.values, T2Mmin_MERRA_bc_subset.values)
        slope_Tmax, int_Tmax, r_Tmax, _, _ = linregress(Tmax_WMO_subset.values, T2Mmax_MERRA_bc_subset.values)

        print("Plotting...")
        # Plot data and lines of best fit
        fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))

        axes[0].scatter(Tavg_WMO_subset, T2Mavg_MERRA_bc_subset, label='Tavg')
        axes[0].plot(Tavg_WMO_subset, slope_Tavg * Tavg_WMO_subset + int_Tavg, color='red', label='Regression Line')
        # Add a 1-to-1 line
        min_val_4plt = min(Tavg_WMO_subset.min(), T2Mavg_MERRA_bc_subset.min())
        max_val_4plt = max(Tavg_WMO_subset.max(), T2Mavg_MERRA_bc_subset.max())
        axes[0].plot([min_val_4plt, max_val_4plt], [min_val_4plt, max_val_4plt], color='black', linestyle='--', label='1-to-1 Line')
        # Add annotation with regression info
        annotate_regression(axes[0], slope_Tavg, int_Tavg, r_Tavg)
        axes[0].set_title('Tavg')
        axes[0].set_xlabel('T_WMO')
        axes[0].set_ylabel('T_MERRA')
        axes[0].legend()
        # Set x and y limits to match the minimum and maximum values
        min_val_4plt = min(Tavg_WMO_subset.min(), T2Mavg_MERRA_bc_subset.min())
        max_val_4plt = max(Tavg_WMO_subset.max(), T2Mavg_MERRA_bc_subset.max())
        axes[0].set_xlim(min_val_4plt, max_val_4plt)
        axes[0].set_ylim(min_val_4plt, max_val_4plt)

        axes[1].scatter(Tmin_WMO_subset, T2Mmin_MERRA_bc_subset, label='Tmin')
        axes[1].plot(Tmin_WMO_subset, slope_Tmin * Tmin_WMO_subset + int_Tmin, color='red', label='Regression Line')
        # Add a 1-to-1 line
        min_val_4plt = min(Tmin_WMO_subset.min(), T2Mmin_MERRA_bc_subset.min())
        max_val_4plt = max(Tmin_WMO_subset.max(), T2Mmin_MERRA_bc_subset.max())
        axes[1].plot([min_val_4plt, max_val_4plt], [min_val_4plt, max_val_4plt], color='black', linestyle='--', label='1-to-1 Line')
        # Add annotation with regression info
        annotate_regression(axes[1], slope_Tmin, int_Tmin, r_Tmin)
        axes[1].set_title('Tmin')
        axes[1].set_xlabel('T_WMO')
        axes[1].set_ylabel('T_MERRA')
        axes[1].legend()
        # Set x and y limits to match the minimum and maximum values
        axes[1].set_xlim(min_val_4plt, max_val_4plt)
        axes[1].set_ylim(min_val_4plt, max_val_4plt)

        axes[2].scatter(Tmax_WMO_subset, T2Mmax_MERRA_bc_subset, label='Tmax')
        axes[2].plot(Tmax_WMO_subset, slope_Tmax * Tmax_WMO_subset + int_Tmax, color='red', label='Regression Line')
        # Add a 1-to-1 line
        min_val_4plt = min(Tmax_WMO_subset.min(), T2Mmax_MERRA_bc_subset.min())
        max_val_4plt = max(Tmax_WMO_subset.max(), T2Mmax_MERRA_bc_subset.max())
        axes[2].plot([min_val_4plt, max_val_4plt], [min_val_4plt, max_val_4plt], color='black', linestyle='--', label='1-to-1 Line')
        # Add annotation with regression info
        annotate_regression(axes[2], slope_Tmax, int_Tmax, r_Tmax)
        axes[2].set_title('Tmax')
        axes[2].set_xlabel('T_WMO')
        axes[2].set_ylabel('T_MERRA')
        axes[2].legend()
        # Set x and y limits to match the minimum and maximum values
        axes[2].set_xlim(min_val_4plt, max_val_4plt)
        axes[2].set_ylim(min_val_4plt, max_val_4plt)

        plt.tight_layout()

        print("Saving plot...")
        # Save the plot
        dir_out = "/mnt/redwood/local_drive/chanud/Output/test/"
        if not os.path.isdir(dir_out):
            os.makedirs(dir_out)
        filepath_out = dir_out + "T_MERRA2_bc_vs_obsWMO_" + city + "_linReg" + ".png"
        plt.savefig(filepath_out, dpi=300, bbox_inches="tight")
