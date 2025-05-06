# scripts_analysis.py
#
# Assorted blocks of code for analyzing data from WHATCH'EM and the vector survival simulations.
#
# Revision history
#    2025-04-17, CNY: updated scripts to account for Focks vs Eisen survival curves



# import required packages
# ************************************************
import calendar
import sys
import re
import numpy as np
import pandas as pd



# import constants and functions
# ************************************************
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/0_constants") # directory of constants file
from CONST_FILE_DIR import DIR_WHATCHEM_OUTPUT_DATA_PROCESSED
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/4_analysis")  # directory of helper fns file
from fns_analysis import calc_tercile_crosstab_splitByMon
from fns_helper import (strConvert_paren2hyphen, strConvert_parenUnitsRemove, get_df_vectorPop_dly_merged, convert_df_vectorPop_dly2mon, 
                        create_df_withLags, load_crosstab)



# calculate and save tercile crosstabs of monthly values of variables against dengue incidence
# (each resulting df has the crosstab values for ONE variable for a given lag+outbreakYr scenario)
if False:

    # parameters related to modeling pipeline
    dataSrc = 'MERRA_IMERG'
    locs = ['Negombo', 'NuwaraEliya', 'Jaffna']
    years = list(range(2001, 2020+1))
    months = list(range(1,12+1))
    intv = 'intv00'
    initDays = list(range(1,7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000
    survCurve  = 'Eisen' # ['Focks', 'Eisen']
    var_colNames = ['TA_min (C)', 'TA (C)' ,'TA_max (C)', 'PRCP (mm)', 'TW_min_nanFill (C)', 'TW_nanFill (C)', 'TW_max_nanFill (C)', 'WH (mm)',
                   'dev_E', 'surv_E', 'dev_L', 'surv_L', 'dev_P', 'surv_P', 'pop_A(E)', 'pop_A(L)', 'pop_A(P)'] # subset of vars to actually plot
    var_pltNames = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in var_colNames]

    # parameters related to dengue data
    lagTimes = [0, 1, 2]              # in months
    den_colName_base = 'denInc'       # we will rename columns to this string + the lag time (original name is the district)
    noOutbreaks_opts = [True, False]  # whether to include data from outbreak years
    outbreakYrs = [2017, 2019]

    # parameters related to modeling pipeline + dengue data
    startDate = '2007-01-01'  # start date (first year of dengue data)
    endDate   = '2020-12-31'  # end date (last year of downloaded MERRA/IMERG data)

    # load dengue data
    fileDir_base = '/mnt/redwood/local_drive/chanud/RawData/epiReports/'
    file_dengue  = fileDir_base + 'dengueInc_monthly_2007_2022.csv'
    df_den_monthly = pd.read_csv(file_dengue, index_col="Date", parse_dates=True)

    for loc in locs:
        print(loc)
        if loc == 'Negombo':
            den_loc = 'Gampaha'
        elif loc in ['NuwaraEliya', 'Jaffna']:
            den_loc = loc
        else:
            print('Error: unknown loc.')
            exit

        # subset dengue data to location of interest and create df with lagged data (which also subsets to time of interest)
        series_den_monthly_loc = df_den_monthly[den_loc]
        df_den_monthly_loc_wLags = create_df_withLags(series_den_monthly_loc, lagTimes, den_colName_base, startDate, endDate)

        # load model pipeline data at location of interest
        # also: convert from daily to monthly, subset to time of interest
        df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve)
        df_vectorPop_monthly = convert_df_vectorPop_dly2mon(df_vectorPop_dly)
        df_vectorPop_monthly = df_vectorPop_monthly[startDate:endDate]

        # combine model pipeline data and dengue data into one df
        df_vectorPop_den_monthly = df_vectorPop_monthly.join(df_den_monthly_loc_wLags)

        if survCurve == 'Eisen':
            fileStr_survCurve = '_eisenQuadFit'
        else: # if survCurve is 'Focks'
            fileStr_survCurve = ''

        for noOutbreaks in noOutbreaks_opts:
            print(' noOutbreaks: ' + str(noOutbreaks))
            if noOutbreaks:
                fileStr_noOutbreaks = '_noOutbreaks'
                exclYrs = outbreakYrs  # outbreak years (to exclude from plot)
            else:
                fileStr_noOutbreaks = ''
                exclYrs = []

            # remove data in outbreak years we want to exclude
            df_vectorPop_den_monthly_yrFiltered = df_vectorPop_den_monthly[~df_vectorPop_den_monthly.index.year.isin(exclYrs)]
            

            for i_lag, lagTime in enumerate(lagTimes):
                print("   Lag (mo): " + str(lagTime))

                for i_var, var_colName in enumerate(var_colNames):
                    var_pltName = var_pltNames[i_var]
                    print('  Var: ' + var_pltName)

                    den_colName = den_colName_base + '-lag' + str(lagTime) # this must match col name format set in create_df_withLags()
                    den_pltName = den_colName # just reusing the column name here for simplicity
                    
                    fileName = 'crosstab_tercile_splitbyMon_' + var_pltName + '_' + den_pltName + '_' + loc + fileStr_noOutbreaks + fileStr_survCurve
                    calc_tercile_crosstab_splitByMon(df_vectorPop_den_monthly_yrFiltered, var_colName, den_colName, fileName)



# calculate and save tercile crosstab element values (e.g. High-High) of monthly values of variables against dengue incidence
# (each resulting df has the crosstab element values for ALL the variables for a given lag+outbreakYr scenario)
if False:

    # parameters related to modeling pipeline
    dataSrc = 'MERRA_IMERG'
    locs = ['Negombo', 'NuwaraEliya', 'Jaffna']
    years = list(range(2001, 2020+1))
    months = list(range(1,12+1))
    intv = 'intv00'
    initDays = list(range(1,7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000
    survCurve  = 'Eisen' # ['Focks', 'Eisen']
    var_colNames = ['TA_min (C)', 'TA (C)', 'TA_max (C)', 'PRCP (mm)', 'TW_min_nanFill (C)', 'TW_nanFill (C)', 'TW_max_nanFill (C)', 'WH (mm)',
                    'dev_E', 'surv_E', 'dev_L', 'surv_L', 'dev_P', 'surv_P', 'pop_A(E)', 'pop_A(L)', 'pop_A(P)'] # subset of vars to actually plot
    var_pltNames = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in var_colNames]

    # parameters related to dengue data
    lagTimes = [0, 1, 2]              # in months
    den_colName_base = 'denInc'       # we will rename columns to this string + the lag time (original name is the district)
    noOutbreaks_opts = [True, False]  # whether to include data from outbreak years
    outbreakYrs = [2017, 2019]

    # parameters related to modeling pipeline + dengue data
    startDate = '2007-01-01'  # start date (first year of dengue data)
    endDate   = '2020-12-31'  # end date (last year of downloaded MERRA/IMERG data)

    # parameters related to the crosstab element of interest
    # format is 'terc1-terc2' where terc1 is tercile of the variable, terc2 is tercile of dengue
    crosstab_elts = ['HI-HI', 'MED-HI', 'LO-HI', 'HI-LO', 'MED-LO', 'LO-LO']
    map_tercAbbr_to_tercName = {'HI': 'High', 'MED': 'Medium', 'LO': 'Low'} # map to names of rows/columns in crosstab csv file

    # load dengue data
    fileDir_base = '/mnt/redwood/local_drive/chanud/RawData/epiReports/'
    file_dengue  = fileDir_base + 'dengueInc_monthly_2007_2022.csv'
    df_den_monthly = pd.read_csv(file_dengue, index_col="Date", parse_dates=True)

    for loc in locs:
        print(loc)
        if loc == 'Negombo':
            den_loc = 'Gampaha'
        elif loc in ['NuwaraEliya', 'Jaffna']:
            den_loc = loc
        else:
            print('Error: unknown loc.')
            exit

        # subset dengue data to location of interest and create df with lagged data (which also subsets to time of interest)
        series_den_monthly_loc = df_den_monthly[den_loc]
        df_den_monthly_loc_wLags = create_df_withLags(series_den_monthly_loc, lagTimes, den_colName_base, startDate, endDate)

        # load model pipeline data at location of interest
        # also: convert from daily to monthly, subset to time of interest
        df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve)
        df_vectorPop_monthly = convert_df_vectorPop_dly2mon(df_vectorPop_dly)
        df_vectorPop_monthly = df_vectorPop_monthly[startDate:endDate]

        # combine model pipeline data and dengue data into one df
        df_vectorPop_den_monthly = df_vectorPop_monthly.join(df_den_monthly_loc_wLags)

        if survCurve == 'Eisen':
            fileStr_survCurve = '_eisenQuadFit'
        else: # if survCurve is 'Focks'
            fileStr_survCurve = ''

        for noOutbreaks in noOutbreaks_opts:
            print(' noOutbreaks: ' + str(noOutbreaks))
            if noOutbreaks:
                fileStr_noOutbreaks = '_noOutbreaks'
                exclYrs = outbreakYrs  # outbreak years (to exclude from plot)
            else:
                fileStr_noOutbreaks = ''
                exclYrs = []

            # remove data in outbreak years we want to exclude
            df_vectorPop_den_monthly_yrFiltered = df_vectorPop_den_monthly[~df_vectorPop_den_monthly.index.year.isin(exclYrs)]
            
            for i_lag, lagTime in enumerate(lagTimes):
                print("   Lag (mo): " + str(lagTime))

                for crosstab_elt in crosstab_elts:
                    print("    Crosstab element: " + crosstab_elt)

                    tercileStr_modelVar, tercileStr_den = crosstab_elt.split('-')
                    dir_crosstabs = 'crosstabs_' + tercileStr_modelVar + '_' + tercileStr_den + '/'
                    tercileName_modelVar = map_tercAbbr_to_tercName.get(tercileStr_modelVar) # get corresponding row name in crosstabs csv file
                    tercileName_den      = map_tercAbbr_to_tercName.get(tercileStr_den)      # get corresponding col name in crosstabs csv file
                    colName_4_df = tercileName_modelVar + '-' + tercileName_den + ' Value' # (e.g., 'High-High Value')

                    crosstab_elt_vals = {}      # dict to store tercile-tercile values for each variable
                    crosstab_elt_vals_norm = {} # dict to store tercile-tercile values for each variable (normalized to a %)

                    for i_var, var_colName in enumerate(var_colNames):
                        var_pltName = var_pltNames[i_var]
                        print('  Var: ' + var_pltName)

                        den_colName = den_colName_base + '-lag' + str(lagTime) # this must match col name format set in create_df_withLags()
                        den_pltName = den_colName # just reusing the column name here for simplicity
                        
                        fileName = 'crosstab_tercile_splitbyMon_' + var_pltName + '_' + den_pltName + '_' + loc + fileStr_noOutbreaks + fileStr_survCurve
                        filePath = DIR_WHATCHEM_OUTPUT_DATA_PROCESSED + 'crosstabs/' + fileName + '.csv'
                        filePath_norm = DIR_WHATCHEM_OUTPUT_DATA_PROCESSED + 'crosstabs/' + fileName + '_norm.csv'
                        crosstab = load_crosstab(filePath)
                        crosstab_norm = load_crosstab(filePath_norm)

                        crosstab_elt_vals[var_colName]      = crosstab.loc[tercileName_modelVar, tercileName_den] if tercileName_modelVar in crosstab.index and tercileName_den in crosstab.columns else np.nan
                        crosstab_elt_vals_norm[var_colName] = crosstab_norm.loc[tercileName_modelVar, tercileName_den] if tercileName_modelVar in crosstab_norm.index and tercileName_den in crosstab_norm.columns else np.nan
                
                    # Convert the dictionary to a DataFrame
                    df_crosstab_elt = pd.DataFrame(list(crosstab_elt_vals.items()), columns=['Variable', colName_4_df])
                    df_crosstab_elt_norm = pd.DataFrame(list(crosstab_elt_vals_norm.items()), columns=['Variable', colName_4_df])
        
                    fileName = 'crosstab_tercile_splitbyMon_' + tercileStr_modelVar + '-' + tercileStr_den + '_' + den_pltName + '_' + loc + fileStr_noOutbreaks + fileStr_survCurve
                    filePath = DIR_WHATCHEM_OUTPUT_DATA_PROCESSED + dir_crosstabs + fileName + '.csv'
                    filePath_norm = DIR_WHATCHEM_OUTPUT_DATA_PROCESSED + dir_crosstabs + fileName + '_norm.csv'
                    df_crosstab_elt.to_csv(filePath, index=False)
                    df_crosstab_elt_norm.to_csv(filePath_norm, index=False)

