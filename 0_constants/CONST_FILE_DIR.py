# CONST_FILE_DIR.ncl


# Code
# ****************************************
DIR_WHATCHEM_CODE = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/' # same as '~/Code/project_whatchem/whatchem-2.1_CNY/'


# WHATCHEM data and plots
# ****************************************
DIR_WHATCHEM_DATA = '/mnt/redwood/local_drive/chanud/whatchem_2.1_CNY/'
DIR_WHATCHEM_INPUT_DATA   = DIR_WHATCHEM_DATA + 'input_data/'
DIR_WHATCHEM_OUTPUT_DATA  = DIR_WHATCHEM_DATA + 'output_data/'
DIR_WHATCHEM_OUTPUT_DATA_PROCESSED = DIR_WHATCHEM_DATA + 'output_data_processed/'
DIR_WHATCHEM_OUTPUT_PLOTS = DIR_WHATCHEM_DATA + 'output_plots/'


# Downloaded raw data
# ****************************************
DIR_RAW_DATA = '/mnt/redwood/local_drive/chanud/RawData/'

# Data subset to Sri Lanka
DIR_MERRA_SL = DIR_RAW_DATA + 'MERRA_1980-2021_SL/'
DIR_IMERG_SL = DIR_RAW_DATA + 'IMERG_SL/'

# meteorological fields

# (hourly)
PATH_MERRA_T      = DIR_MERRA_SL + 'MERRA_hourly_SLV/MERRA2_tavg1_2d_slv_Nx.19800101_20211231.SUB.T2M.nc'
PATH_MERRA_T_NUWARAELIYA_BIASCORRECTED = DIR_MERRA_SL + 'MERRA_hourly_SLV/MERRA2_tavg1_2d_slv_Nx.19800101_20211231.SUB.T2M.NuwaraEliya.biasCorrected.nc'
PATH_MERRA_TSOIL  = DIR_MERRA_SL + 'MERRA_hourly_LND/MERRA2_tavg1_2d_lnd_Nx.19800101_20211231.SUB.TSOIL1.nc'
PATH_MERRA_TSOIL_NUWARAELIYA_BIASCORRECTED = DIR_MERRA_SL + 'MERRA_hourly_LND/MERRA2_tavg1_2d_lnd_Nx.19800101_20211231.SUB.TSOIL1.NuwaraEliya.biasCorrected.nc'
PATH_MERRA_RH     = DIR_MERRA_SL + 'MERRA_hourly_SLV/MERRA2_tavg1_2d_slv_Nx.19800101_20211231.SUB.RH.nc'
PATH_MERRA_LC     = DIR_MERRA_SL + 'MERRA_hourly_RAD/MERRA2_tavg1_2d_rad_Nx.19800101_20211231.SUB.CLDLOW.nc'
PATH_MERRA_MC     = DIR_MERRA_SL + 'MERRA_hourly_RAD/MERRA2_tavg1_2d_rad_Nx.19800101_20211231.SUB.CLDMID.nc'
PATH_MERRA_HC     = DIR_MERRA_SL + 'MERRA_hourly_RAD/MERRA2_tavg1_2d_rad_Nx.19800101_20211231.SUB.CLDHGH.nc'
PATH_IMERG_PRECIP = DIR_IMERG_SL + 'mergedHR_2001_2020.nc'

# WMO observations (daily)
DIR_OBS = '/mnt/redwood/local_drive/chanud/RawData/observations/'
PATH_WMO_NUWARAELIYA = DIR_OBS + 'NuwaraEliya_WMO_2001-2020.csv'
PATH_WMO_JAFFNA      = DIR_OBS + 'Jaffna_WMO_2009-2020.csv'
PATH_WMO_NEGOMBO     = DIR_OBS + 'Katunayake_WMO_2001-2020.csv'

