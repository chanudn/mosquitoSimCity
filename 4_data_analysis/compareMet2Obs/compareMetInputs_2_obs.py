# compareMetInputs_2_obs.py
#
# We have several meteorological inputs into our modeling pipeline.
#   MERRA-2: 2-m air temperature, soil temperature, low/med/high cloud fraction,
#            relative humidity (derived), vapor pressure deficit (derived)
#   IMERG:   precipitation
# These are from reanalysis/satellite products, so how well do they reflect the
# actual values of these variables? To find out, let's compare to observations.
# This script compares our met inputs to WMO observations for the variables that are
# (a) most important for our models and (b) easiest to find corresponding obs for.
# These vars are:
#   - air temperature
#   - precipitation
#   - vapor pressure deficit (derived from air temp, RH)
# We compare our inputs to obs via histograms and via scatter plots of anomalies
# (which include correlation coefficients from linear regression).
# Note: we use anomalies for the scatter plots/lin reg to reduce the impact of seasonality.


# import required packages
# ************************************************
import xarray as xr
import pandas as pd
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import math
import re
import sys
import os

# import constants
# ************************************************
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/0_constants")
from CONST_FILE_DIR import *
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/4_data_analysis")
from fns_helper import strConvert_paren2hyphen, strConvert_parenUnitsRemove


# parameters
locs = ["Negombo", "NuwaraEliya", "Jaffna"]
var_colNames = ['T_mean (C)', 'T_min (C)', 'T_max (C)', 'PRCP (mm)', 'VPD (kPa)']
var_pltNames = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in var_colNames]
start_date = '2007-01-01' # note: Jaffna starts 2009-01-01, which is accounted for later on by taking the intersection of indices
end_date   = '2020-12-31'


# Read in hourly MERRA-2 and IMERG data
#   subset to timespan of interest
#   resample from hourly to daily
dataset_T2M_MERRA_SL = xr.open_dataset(PATH_MERRA_T)   # MERRA-2 2-m air temperature
T2M_MERRA_SL = dataset_T2M_MERRA_SL['T2M']
T2M_MERRA_SL = T2M_MERRA_SL.sel(time=slice(start_date, end_date))
T2M_MERRA_SL = T2M_MERRA_SL - 273.15 # convert from K to C
T2M_mean_MERRA_SL = T2M_MERRA_SL.resample(time='D').mean()
T2M_min_MERRA_SL  = T2M_MERRA_SL.resample(time='D').min()
T2M_max_MERRA_SL  = T2M_MERRA_SL.resample(time='D').max()

dataset_T2M_MERRA_NuwaraEliya = xr.open_dataset(PATH_MERRA_T_NUWARAELIYA_BIASCORRECTED)   # MERRA-2 2-m air temperature for Nuwara Eliya, bias-corrected
T2M_MERRA_NuwaraEliya = dataset_T2M_MERRA_NuwaraEliya['T2M']
T2M_MERRA_NuwaraEliya = T2M_MERRA_NuwaraEliya.sel(time=slice(start_date, end_date))
T2M_MERRA_NuwaraEliya = T2M_MERRA_NuwaraEliya - 273.15 # convert from K to C
T2M_mean_MERRA_NuwaraEliya = T2M_MERRA_NuwaraEliya.resample(time='D').mean()
T2M_min_MERRA_NuwaraEliya  = T2M_MERRA_NuwaraEliya.resample(time='D').min()
T2M_max_MERRA_NuwaraEliya  = T2M_MERRA_NuwaraEliya.resample(time='D').max()

dataset_RH_MERRA_SL = xr.open_dataset(PATH_MERRA_RH)   # MERRA-2 2-m relative humidity (derived from T2M, QV2M, PS)
RH_MERRA_SL = dataset_RH_MERRA_SL['RH']
RH_MERRA_SL = RH_MERRA_SL.sel(time=slice(start_date, end_date))
RH_MERRA_SL = RH_MERRA_SL.resample(time='D').mean()

dataset_PRCP_IMERG_SL = xr.open_dataset(PATH_IMERG_PRECIP) # IMERG precipitation
PRCP_IMERG_SL = dataset_PRCP_IMERG_SL['precipitationCal']  # get the precipitation variable
PRCP_IMERG_SL = PRCP_IMERG_SL.sel(time=slice(start_date, end_date))
PRCP_IMERG_SL = PRCP_IMERG_SL.resample(time='D').sum()


for loc in locs:

    print(loc)

    # Get location specific values
    # (lat/lon/elev are same as in WHATCH'EM code, run_container.pl)
    if loc == "Negombo":
        lat  = 7.2008  
        lon  = 79.8737 
        elev = 2.
        PATH_WMO = PATH_WMO_NEGOMBO
        has_bc_data = False
    elif loc == "NuwaraEliya":
        lat  = 6.9497  
        lon  = 80.7891 
        elev = 1868.
        PATH_WMO = PATH_WMO_NUWARAELIYA
        has_bc_data = True
    elif loc == "Jaffna":
        lat  = 9.6615  
        lon  = 80.0255 
        elev = 5.
        PATH_WMO = PATH_WMO_JAFFNA
        has_bc_data = False

    ## Read in WMO data
    data_WMO = pd.read_csv(PATH_WMO)
    data_WMO.set_index(pd.to_datetime(data_WMO['DATE']), inplace=True)
    
    # Subset the df to only our timespan and variables of interest
    data_WMO = data_WMO.loc[start_date:end_date]
    data_WMO = data_WMO[['DEWP', 'TEMP', 'MAX', 'MIN', 'PRCP']]

    ## Replace the missing values with nan
    data_WMO['DEWP'].replace(9999.9, np.nan, inplace=True)
    data_WMO['TEMP'].replace(9999.9, np.nan, inplace=True)
    data_WMO['MAX'].replace(9999.9, np.nan, inplace=True)
    data_WMO['MIN'].replace(9999.9, np.nan, inplace=True)
    data_WMO['PRCP'].replace(99.99, np.nan, inplace=True)

    ## Convert the temps from F to C and the prcp from in. to mm
    data_WMO['DEWP'] = (data_WMO['DEWP'] - 32) * (5/9) # convert from F to C
    data_WMO['TEMP'] = (data_WMO['TEMP'] - 32) * (5/9) # convert from F to C
    data_WMO['MAX']  = (data_WMO['MAX'] - 32) * (5/9)  # convert from F to C
    data_WMO['MIN']  = (data_WMO['MIN'] - 32) * (5/9)  # convert from F to C
    data_WMO['PRCP'] = (data_WMO['PRCP'] * 25.4)       # convert from in. to mm
    data_WMO = data_WMO.rename(columns={'DEWP':'DEWP (C)', 'TEMP':'T_mean (C)', 'MAX':'T_max (C)', 'MIN':'T_min (C)', 'PRCP':'PRCP (mm)'})

    # Calculate VPD and add it as a column of the WMO df.
    # VPD is directly calculated from saturation vapor pressure (SVP) and
    # vapor pressure (VP), which in turn can be calculated from air temperature
    # (Ta) and dewpoint temperature (Td).
    #   VPD = SVP - VP  [Pa]
    # where
    #   SVP = e_s0 * exp[ (L_v/R_v)*( (1/T_0)-(1/Ta) ) ]  [Pa]
    #    VP = e_s0 * exp[ (L_v/R_v)*( (1/T_0)-(1/Td) ) ]  [Pa]
    # where
    #   e_s0 - saturation vapor pressure at T_0 [Pa]
    #   L_v  - specific enthalpy of vaporization [J/kg]
    #   R_v  - specific gas constant for water vapor [J/kg/K]
    #   T_0  - reference temperature [K] (we use 273.15 K)
    #   Ta   - air temperature [K]
    #   Td   - dewpoint temperature [K]
    
    # Define constants
    # Note: We get L_v at 25C from https://www.engineeringtoolbox.com/water-properties-d_1573.html.
    #       It's a fn of T (2.500e6 at 0C, 2.442e6 at 25C, 2.256e6 at 100C), but we use a constant value.
    e_s0 = 611       # saturation vapor pressure at T_0 [Pa]
    L_v  = 2.442e6   # Latent heat of vaporization at 25C [J/kg]
    R_v  = 461.5     # Gas constant for water vapor [J/kg/K]
    T_0  = 273.15    # reference temperature [K]
    
    # Calculate VPD
    Ta = data_WMO['T_mean (C)'] + T_0 # air temperature, converted from C to K
    Td = data_WMO['DEWP (C)'] + T_0   # dewpoint temperature, converted from C to K
    svp = e_s0*np.exp( (L_v/R_v)*((1/T_0)-(1/Ta)) )   # [Pa]
    vp  = e_s0*np.exp( (L_v/R_v)*((1/T_0)-(1/Td)) )   # [Pa]
    vpd = (svp-vp) / 1000                             # [kPa]
    data_WMO['VPD (kPa)'] = vpd

    data_WMO = data_WMO.drop(columns=['DEWP (C)'])


    # Extract MERRA-2 and IMERG data as Pandas DataFrames
    # for location of interest
    T2M_mean_MERRA = T2M_mean_MERRA_SL.sel(lat=lat, lon=lon, method='nearest')
    T2M_mean_MERRA = T2M_mean_MERRA.drop(['lat', 'lon'])
    T2M_min_MERRA  = T2M_min_MERRA_SL.sel(lat=lat, lon=lon, method='nearest')
    T2M_min_MERRA  = T2M_min_MERRA.drop(['lat', 'lon'])
    T2M_max_MERRA  = T2M_max_MERRA_SL.sel(lat=lat, lon=lon, method='nearest')
    T2M_max_MERRA  = T2M_max_MERRA.drop(['lat', 'lon'])

    T2M_mean_MERRA = pd.Series(T2M_mean_MERRA.values, index=T2M_mean_MERRA.time, name='T_mean (C)')
    T2M_min_MERRA = pd.Series(T2M_min_MERRA.values, index=T2M_min_MERRA.time, name='T_min (C)')
    T2M_max_MERRA = pd.Series(T2M_max_MERRA.values, index=T2M_max_MERRA.time, name='T_max (C)')

    RH_MERRA = RH_MERRA_SL.sel(lat=lat, lon=lon, method='nearest')
    RH_MERRA = RH_MERRA.drop(['lat', 'lon'])
    RH_MERRA = pd.Series(RH_MERRA.values, index=RH_MERRA.time, name='RH (%)')

    PRCP_IMERG = PRCP_IMERG_SL.sel(lat=lat, lon=lon, method='nearest')
    PRCP_IMERG = PRCP_IMERG.drop(['lat', 'lon'])
    PRCP_IMERG = pd.Series(PRCP_IMERG.values, index=PRCP_IMERG.time, name = 'PRCP (mm)')

    # Calculate VPD
    Ta = T2M_mean_MERRA + T_0 # air temperature, converted from C to K
    svp = e_s0*np.exp( (L_v/R_v)*((1/T_0)-(1/Ta)) )   # [Pa]
    vp  = svp * (RH_MERRA/100)                        # [Pa]
    vpd = (svp-vp) / 1000                             # [kPa]
    VPD_MERRA = vpd.rename('VPD (kPa)')

    data_MERRA_IMERG = pd.concat([T2M_mean_MERRA, T2M_min_MERRA, T2M_max_MERRA, PRCP_IMERG, VPD_MERRA], axis=1)

    # Trim datetime indices that are not present in both datasets (likely bc WMO has missing data)
    commonIndex = data_WMO.index.intersection(data_MERRA_IMERG.index)
    data_WMO = data_WMO.loc[commonIndex]
    data_MERRA_IMERG = data_MERRA_IMERG.loc[commonIndex]

    if has_bc_data: # only true for Nuwara Eliya

        # get Nuwara Eliya-specific bias-corrected data
        T2M_mean_MERRA_bc = T2M_mean_MERRA_NuwaraEliya
        T2M_min_MERRA_bc  = T2M_min_MERRA_NuwaraEliya
        T2M_max_MERRA_bc  = T2M_max_MERRA_NuwaraEliya
            
        T2M_mean_MERRA_bc = pd.Series(T2M_mean_MERRA_bc.values, index=T2M_mean_MERRA_bc.time, name='T_mean (C)')
        T2M_min_MERRA_bc = pd.Series(T2M_min_MERRA_bc.values, index=T2M_min_MERRA_bc.time, name='T_min (C)')
        T2M_max_MERRA_bc = pd.Series(T2M_max_MERRA_bc.values, index=T2M_max_MERRA_bc.time, name='T_max (C)')

        # calculate VPD from bias-corrected T2M
        Ta_bc = T2M_mean_MERRA_bc + T_0 # air temperature, converted from C to K
        svp_bc = e_s0*np.exp( (L_v/R_v)*((1/T_0)-(1/Ta_bc)) )   # [Pa]
        vp_bc  = svp_bc * (RH_MERRA/100)                        # [Pa]
        vpd_bc = (svp_bc-vp_bc) / 1000                          # [kPa]
        VPD_MERRA_bc = vpd_bc.rename('VPD (kPa)')
        
        data_MERRA_IMERG_bc = pd.concat([T2M_mean_MERRA_bc, T2M_min_MERRA_bc, T2M_max_MERRA_bc, PRCP_IMERG, VPD_MERRA_bc], axis=1)

        # Trim datetime indices that are not present in both datasets (likely bc WMO has missing data)
        commonIndex_bc = data_WMO.index.intersection(data_MERRA_IMERG_bc.index)
        data_WMO_bc = data_WMO.loc[commonIndex_bc]
        data_MERRA_IMERG_bc = data_MERRA_IMERG_bc.loc[commonIndex_bc]

    # Make histogram for each variable!
    for i_var, var_colName in enumerate(var_colNames):
        print('  ' + var_colName)
        var_pltName = var_pltNames[i_var]

        if var_colName == 'PRCP (mm)':
            label_MERRA_IMERG = 'IMERG'
        else:
            label_MERRA_IMERG = 'MERRA-2'

        fig, ax = plt.subplots() # create a new figure and axis for plot
        
        # Determine bin edges based on ranges of data
        min_val = min(data_WMO[var_colName].min(), data_MERRA_IMERG[var_colName].min())
        max_val = max(data_WMO[var_colName].max(), data_MERRA_IMERG[var_colName].max())
        num_bins = 20
        if var_colName == 'PRCP (mm)':
            tinyOffsetForLogPlots = 0.01
            bins = np.logspace(np.log10(min_val+tinyOffsetForLogPlots), np.log10(max_val+tinyOffsetForLogPlots), num_bins + 1) # bin edges
            ax.hist(data_WMO[var_colName]+tinyOffsetForLogPlots, bins=bins, alpha=0.5, label='Observations')
            ax.hist(data_MERRA_IMERG[var_colName]+tinyOffsetForLogPlots, bins=bins, alpha=0.5, label=label_MERRA_IMERG)
            ax.set_xscale('log')
        else:
            bins = np.linspace(min_val, max_val, num_bins + 1) # bin edges
            ax.hist(data_WMO[var_colName], bins=bins, alpha=0.5, label='Observations')
            ax.hist(data_MERRA_IMERG[var_colName], bins=bins, alpha=0.5, label=label_MERRA_IMERG)


        # Adding labels and legend
        plt.xlabel(var_colName)
        plt.ylabel('Frequency')
        plt.legend()

        # Save the plot
        dir_out = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/1_metInputs/'
        if not os.path.isdir(dir_out):
            os.makedirs(dir_out)
        filepath_out = dir_out + "hist_compare2obs_" + var_pltName + '_' + loc + ".png"
        plt.savefig(filepath_out, dpi=300, bbox_inches="tight")

        if has_bc_data and (var_colName != 'PRCP (mm)'): # only Nuwara Eliya has bias-corrected data, and it's all vars except precip
                
            fig, ax = plt.subplots() # create a new figure and axis for plot
        
            # Determine bin edges based on ranges of data
            min_val = min(data_WMO_bc[var_colName].min(), data_MERRA_IMERG_bc[var_colName].min())
            max_val = max(data_WMO_bc[var_colName].max(), data_MERRA_IMERG_bc[var_colName].max())
            num_bins = 20
            if var_colName == 'PRCP (mm)':
                tinyOffsetForLogPlots = 0.01
                bins = np.logspace(np.log10(min_val+tinyOffsetForLogPlots), np.log10(max_val+tinyOffsetForLogPlots), num_bins + 1) # bin edges
                ax.hist(data_WMO[var_colName]+tinyOffsetForLogPlots, bins=bins, alpha=0.5, label='Observations')
                ax.hist(data_MERRA_IMERG[var_colName]+tinyOffsetForLogPlots, bins=bins, alpha=0.5, label=label_MERRA_IMERG)
            else:
                bins = np.linspace(min_val, max_val, num_bins + 1) # bin edges
                ax.hist(data_WMO[var_colName], bins=bins, alpha=0.5, label='Observations')
                ax.hist(data_MERRA_IMERG[var_colName], bins=bins, alpha=0.5, label=label_MERRA_IMERG)

            # Adding labels and legend
            plt.xlabel(var_colName)
            plt.ylabel('Frequency')
            plt.legend()

            # Save the plot
            dir_out = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/1_metInputs/'
            if not os.path.isdir(dir_out):
                os.makedirs(dir_out)
            filepath_out = dir_out + "hist_compare2obs_" + var_pltName + '_bc' + '_' + loc + ".png"
            plt.savefig(filepath_out, dpi=300, bbox_inches="tight")