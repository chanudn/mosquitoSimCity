# createMetInput.py
#
# Ported to python from original ncl code (createMetInput.ncl).
# TO DO:
#   - add handling of doing climatology of months or seasons (rather than just full years)

# import required packages
# ************************************************
import xarray as xr
import pandas as pd
import sys
import re

# import constants
# ************************************************
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/0_constants")
from CONST_FILE_DIR import *




# ********************************************************************************************************
#   Main functions (called by other files)
# ********************************************************************************************************


# create_met_input()
# ************************************************
# Create a WHATCH'EM meteorological input text file based on hourly meteorological data.
#
# Arguments:
#   dataSrc    - string; "MERRA_IMERG" (TODO: add handling of "GEOS_S2S")
#   loc        - string; "Negombo", "NuwaraEliya", "Jaffna"
#   tspan_strt - int; must be in YYYYMMDDHH format
#   tspan_end  - int; must be in YYYYMMDDHH format
def create_met_input(dataSrc, loc, tspan_strt, tspan_end):

    # Check args
    if dataSrc not in ["MERRA_IMERG"]:
        raise ValueError(f"Invalid argument: {dataSrc}. Please provide one of: 'MERRA_IMERG'")
    if loc not in ["Negombo", "NuwaraEliya", "Jaffna"]:
        raise ValueError(f"Invalid argument: {loc}. Please provide one of: 'Negombo', 'NuwaraEliya', 'Jaffna'")
    if not re.compile(r'^\d{10}$').match(str(tspan_strt)):
        raise ValueError(f"Invalid argument: {tspan_strt}. Please provide a number in YYYYMMDDHH format.")
    if not re.compile(r'^\d{10}$').match(str(tspan_end)):
        raise ValueError(f"Invalid argument: {tspan_strt}. Please provide a number in YYYYMMDDHH format.")

    # Print input args
    print("")
    print("DATA SOURCE:  " + dataSrc)
    print("  SITE NAME:  " + loc)
    print("   TIMESPAN:  " + str(tspan_strt) + " to " + str(tspan_end))
    print("NO CLIMATOLOGY")
    print("")

    df = extract_data_to_df(dataSrc, loc, tspan_strt, tspan_end) # subset data to desired location/timespan and put in one dataframe
    df = format_df_for_txtfile(df)                               # format dataframe before saving to text file

    # Save dataframe to text file
    name_srcLocTime = dataSrc + "-" + loc + "-" + str(tspan_strt) + "_" + str(tspan_end)
    f_out = DIR_WHATCHEM_INPUT_DATA + "container_input_" + name_srcLocTime + ".txt"
    df.to_csv(f_out, sep='\t', index=False)



# create_met_input_clim()
# ************************************************
# Create WHATCH'EM meteorological input text files based on the climatology of hourly meteorological data.
# The number of generated files depends on climType (e.g., climType = "month" means 12 files, one for each month).
#
# Arguments:
#   dataSrc    - string; "MERRA_IMERG" (TODO: add handling of "GEOS_S2S")
#   loc        - string; "Negombo", "NuwaraEliya", "Jaffna"
#   tspan_strt - int; must be in YYYYMMDDHH format
#   tspan_end  - int; must be in YYYYMMDDHH format
#   climType   - string; "month" (TODO: add handling of "season" (e.g., FIM))
def create_met_input_clim(dataSrc, loc, tspan_strt, tspan_end, climType):

    # Check args
    if dataSrc not in ["MERRA_IMERG"]:
        raise ValueError(f"Invalid argument: {dataSrc}. Please provide one of: 'MERRA_IMERG'")
    if loc not in ["Negombo", "NuwaraEliya", "Jaffna"]:
        raise ValueError(f"Invalid argument: {loc}. Please provide one of: 'Negombo', 'NuwaraEliya', 'Jaffna'")
    if not re.compile(r'^\d{10}$').match(str(tspan_strt)):
        raise ValueError(f"Invalid argument: {tspan_strt}. Please provide a number in YYYYMMDDHH format.")
    if not re.compile(r'^\d{10}$').match(str(tspan_end)):
        raise ValueError(f"Invalid argument: {tspan_strt}. Please provide a number in YYYYMMDDHH format.")
    if climType not in ["month"]:
        raise ValueError(f"Invalid argument: {climType}. Please provide one of: 'month'")

    # Print input args
    print("")
    print("DATA SOURCE:  " + dataSrc)
    print("  SITE NAME:  " + loc)
    print("   TIMESPAN:  " + str(tspan_strt) + " to " + str(tspan_end))
    print("CLIMATOLOGY:  " + climType)
    print("")
    
    df = extract_data_to_df(dataSrc, loc, tspan_strt, tspan_end) # subset data to desired location/timespan and put in one dataframe
    df = calc_climatology(df)                                    # calculate climatology

    # Subset and save climatology to separate files according to desired climType
    if climType == "month":
        monthStrs = [str(x).zfill(2) for x in range(1,13)] # zero-padded (two-digit) months
        for monthStr in monthStrs:

            # Subset and format data
            df_month = df[df['MM'] == monthStr]        # subset df to month of interest
            df_month = format_df_for_txtfile(df_month) # format dataframe before saving to text file

            # Save dataframe to text file
            yyyyStr_strt = str(tspan_strt)[:4]
            yyyyStr_end  = str(tspan_end)[:4]
            name_srcLocTime = dataSrc + "-" + loc + "-" + yyyyStr_strt + "_" + yyyyStr_end + "_" + monthStr + "-clim"
            f_out = DIR_WHATCHEM_INPUT_DATA + "container_input_" + name_srcLocTime + ".txt"
            df_month.to_csv(f_out, sep='\t', index=False)







# ********************************************************************************************************
#   Helper functions (only used within this file)
# ********************************************************************************************************


# extract_data_to_df()
# ************************************************
# Subset data to desired location/timespan and put it all in one dataframe
#
# Arguments:
#   dataSrc    - string; "MERRA_IMERG" (TODO: add handling of "GEOS_S2S")
#   loc        - string; "Negombo", "NuwaraEliya", "Jaffna"
#   tspan_strt - int; must be in YYYYMMDDHH format
#   tspan_end  - int; must be in YYYYMMDDHH format
def extract_data_to_df(dataSrc, loc, tspan_strt, tspan_end):

    # Load MERRA data
    #   T..................MERRA-2 T2M
    #   RH.................MERRA-2 RH, calculated from MERRA-2 QV2M, T2M, PS
    #   LC / MC / HC.......MERRA-2 CLDLOW / CLDMID / CLDHGH
    #   TSOIL (50mm).......MERRA-2 TSOIL1 (0-0.0988m)
    if loc == "NuwaraEliya": # for Nuwara Eliya, load in the bias corrected temperature data
        file_T = xr.open_dataset(PATH_MERRA_T_NUWARAELIYA_BIASCORRECTED)
        file_TSOIL = xr.open_dataset(PATH_MERRA_TSOIL_NUWARAELIYA_BIASCORRECTED)
    else:
        file_T = xr.open_dataset(PATH_MERRA_T)
        file_TSOIL = xr.open_dataset(PATH_MERRA_TSOIL)
    file_RH    = xr.open_dataset(PATH_MERRA_RH)
    file_LC    = xr.open_dataset(PATH_MERRA_LC)
    file_MC    = xr.open_dataset(PATH_MERRA_MC)
    file_HC    = xr.open_dataset(PATH_MERRA_HC)

    # Load IMERG data
    #   PRECIP.............IMERG precipitationCal
    file_PRECIP = xr.open_dataset(PATH_IMERG_PRECIP)

    # Subset data based on location
    site_lat, site_lon = get_site_coordinates(loc)
    
    if loc == "NuwaraEliya": # the bias-corrected Nuwara Eliya temp data is already subset by location
        site_T     = file_T
        site_TSOIL = file_TSOIL
    else:
        site_T     = file_T.sel(lat=site_lat, lon=site_lon, method='nearest')
        site_TSOIL = file_TSOIL.sel(lat=site_lat, lon=site_lon, method='nearest')
    site_RH     = file_RH.sel(lat=site_lat, lon=site_lon, method='nearest')
    site_LC     = file_LC.sel(lat=site_lat, lon=site_lon, method='nearest')
    site_MC     = file_MC.sel(lat=site_lat, lon=site_lon, method='nearest')
    site_HC     = file_HC.sel(lat=site_lat, lon=site_lon, method='nearest')
    site_PRECIP = file_PRECIP.sel(lat=site_lat, lon=site_lon, method='nearest')

    # Subset data based on time span
    # Need to manually add an hour to the end time since the last index is not included
    site_time = slice(pd.to_datetime(tspan_strt, format='%Y%m%d%H'), 
                      pd.to_datetime(tspan_end, format='%Y%m%d%H') + pd.Timedelta(hours=1))

    site_T      = site_T.sel(time=site_time)
    site_RH     = site_RH.sel(time=site_time)
    site_LC     = site_LC.sel(time=site_time)
    site_MC     = site_MC.sel(time=site_time)
    site_HC     = site_HC.sel(time=site_time)
    site_TSOIL  = site_TSOIL.sel(time=site_time)
    site_PRECIP = site_PRECIP.sel(time=site_time)

    # Convert data to pandas DataFrame
    # Note: This bit of code is causing the following warning:
    #   FutureWarning: Index.ravel returning ndarray is deprecated; in a future version this will return a view on self.
    #                  values_as_series = pd.Series(values.ravel(), copy=False)
    #       I haven't figured out why this is happening, but it doesn't seem to mess anything up for now.
    df = pd.DataFrame({
        'YYYY': site_T['time'].dt.year.values,
        'MM': site_T['time'].dt.month.values,
        'DD': site_T['time'].dt.day.values,
        'HH': site_T['time'].dt.hour.values,
        'DOY': site_T['time'].dt.dayofyear.values,
        'T': site_T['T2M'].values - 273.15,            # convert from K to deg C
        'RH': site_RH['RH'].values,
        'PRECIP': site_PRECIP['precipitationCal'].values,
        'LC': site_LC['CLDLOW'].values,
        'MC': site_MC['CLDMID'].values,
        'HC': site_HC['CLDHGH'].values,
        'TSOIL': site_TSOIL['TSOIL1'].values - 273.15  # convert from K to deg C
    })

    return df





# get_site_coordinates()
# ************************************************
#
# Get location coordinates based on location name.
#
# Arguments:
#   loc - string; "Negombo", "NuwaraEliya", "Jaffna"
#                 (can handle "Colombo", "Kandy", "Trincomalee", but those are currently unused in the overall code)
def get_site_coordinates(loc):
    coordinates = {    # [deg N, deg E]
        'Negombo':     (7.2008, 79.8737),
        'Colombo':     (6.9271, 79.8612),
        'Kandy':       (7.2906, 80.6337),
        'NuwaraEliya': (6.9497, 80.7891),
        'Trincomalee': (8.5874, 81.2152),
        'Jaffna':      (9.6615, 80.0255)
    }
    return coordinates[loc]





# calc_climatology()
# ************************************************
# Calculate climatology (averaging hourly data across multiple years).
# Sets the YYYY column to a missing value of "9999".
#
# Arguments:
#   df_init  - input dataframe of hourly data
def calc_climatology(df_init):

    df = df_init.copy(deep=True)

    # Compute climatology by grouping by MMDDHH
    df['MMDDHH'] = df['MM'].astype(str).str.zfill(2) + df['DD'].astype(str).str.zfill(2) + df['HH'].astype(str).str.zfill(2)
    df = df.drop(["YYYY", "MM", "DD", "HH", "DOY"], axis=1)  # don't want these for climatology calc
    clim_df_B = df.groupby('MMDDHH').mean().reset_index()    # reset_index ensures MMDDHH remains a column of the df

    # Recreate date/time columns
    clim_df_A = pd.DataFrame({
        'YYYY': "9999", # climatology doesn't have a year; use a placeholder value with the same number of digits
        'MM': clim_df_B['MMDDHH'].str[:2],
        'DD': clim_df_B['MMDDHH'].str[2:4],
        'HH': clim_df_B['MMDDHH'].str[4:]
    })

    # Generate DOY column
    sampleLeapYr = 2000 # temp year value (a leap year) for calculating DOY
    date_str = str(sampleLeapYr) + '-' + clim_df_A['MM'].astype(str) + '-' + clim_df_A['DD'].astype(str)
    clim_df_A['DOY'] = pd.to_datetime(date_str, format='%Y-%m-%d').dt.dayofyear

    # Format climatology data into one dataframe
    clim_df = pd.concat([clim_df_A, clim_df_B], axis=1)
    clim_df = clim_df.drop("MMDDHH", axis=1) # drop unwanted columns

    return clim_df





# format_df_for_txtfile()
# ************************************************
# Format dataframe before saving to text file.
# Date columns are padded with leading zeroes (e.g., the first month is "01" rather than "1").
# Other columns are rounded to either two or four decimal places.
#
# Arguments:
#   df_init  - input dataframe of hourly data
def format_df_for_txtfile(df_init):

    df = df_init.copy(deep=True)

    cols_zfill2 = ['MM', 'DD', 'HH']
    df[cols_zfill2] = df[cols_zfill2].astype(str).apply(lambda x: x.str.zfill(2))

    cols_zfill3 = ['DOY']
    df[cols_zfill3] = df[cols_zfill3].astype(str).apply(lambda x: x.str.zfill(3))

    cols_zfill4 = ['YYYY']
    df[cols_zfill4] = df[cols_zfill4].astype(str).apply(lambda x: x.str.zfill(4))

    cols_round2 = ['T', 'RH', 'PRECIP', 'TSOIL']
    df[cols_round2] = df[cols_round2].round(2).apply(lambda x: x.apply(lambda y: '{:.2f}'.format(y)))

    cols_round4 = ['LC', 'MC', 'HC']
    df[cols_round4] = df[cols_round4].round(4).apply(lambda x: x.apply(lambda y: '{:.4f}'.format(y)))

    return df