# scripts_plotting.py
#
# Assorted blocks of code for plotting data from WHATCH'EM and the vector survival simulations.



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
from CONST_FILE_DIR import DIR_WHATCHEM_OUTPUT_DATA
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting")       # directory of plotting files
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/4_data_analysis")  # directory of helper fns file
from fns_plotting import (plot_varVsVar_byYear, plot_waterVarvsTime_monthSummary, plot_varVsTime_climEndmembers_oneMon,
                          plot_interannualVar_allMon_tseries, plot_interannualVar_allMon_boxplot, plot_interannualVar_allMon_boxplot_2var, plot_scatter_colorMon, plot_correlogram_allVars,
                          plot_bar_tercileXY_splitByMon, plot_bar_binaryX1X2_tercileY_splitByMon)
from fns_helper import (WHATCHEM_specs_to_fpath, strConvert_paren2hyphen, strConvert_parenUnitsRemove, combine_dfs_overInitDays, 
                          get_df_vectorPop_dly_merged, convert_df_vectorPop_dly2mon, create_df_withLags)






## Plot water height and temperature evolution for each month of the timespan
#  Note: plot_waterVarvsTime_monthSummary() can generate many different plots 
#        (e.g., time evolution for each year; time evolution for each year and mean of that)
#        See the code for plot_waterVarvsTime_monthSummary() to see what
#        plotting is actually enabled.
if False:
    dataSrc   = "MERRA_IMERG"
    locs      = ["Negombo", "Jaffna", "NuwaraEliya"]
    yrStrt    = 2001
    yrStop    = 2020
    years     = range(yrStrt, yrStop+1) # note: first number is included, last number is not
    months    = range(1,13)             # note: first number is included, last number is not
    intv      = "intv00"
    waterVars = ["wh", "tw"]

    print("Plotting...")
    for loc in locs:
        for waterVar in waterVars:
            print("Location: " + loc + "  Variable: " + waterVar)
            plot_waterVarvsTime_monthSummary(dataSrc, loc, years, months, intv, waterVar)






if False: # if plotting
    
        
    # load in previously created dfs
    # then plot interannual variability
    locs = ["Negombo", "Jaffna", "NuwaraEliya"]
    monthNums = np.arange(1,13)
    climVar_colNames = ["TA (C)", "PRCP (mm)", "VPD (kPa)"]
    contVar_colNames = ["TW (C)", "WH (mm)"]
    devVar_colNames  = ["dev_E", "dev_L", "dev_P"]
    survVar_colNames = ["surv_E", "surv_L", "surv_P"]
    popVar_total_colNames  = ['pop_E', 'pop_L', 'pop_P', 'pop_A']
    popVar_subset_colNames = ['pop_L(E)', 'pop_L(L)', 'pop_P(E)', 'pop_P(L)', 'pop_P(P)', 
                            'pop_A(E)', 'pop_A(L)', 'pop_A(P)']

    for loc in locs:
        print(loc)
        df_dev_hrly = pd.read_csv("df_dev_hrly_" + loc + ".txt", sep='\t', na_values="NaN", index_col=0, parse_dates=True)
        df_dev_surv_dly = pd.read_csv("df_dev_surv_dly_" + loc + ".txt", sep='\t', na_values="NaN", index_col=0, parse_dates=True)

        df_dev_hrly_clim = pd.read_csv("df_dev_hrly_clim_" + loc + ".txt", sep='\t', na_values="NaN", index_col=["Month", "Day", "Hour"], parse_dates=True)
        df_dev_surv_dly_clim = pd.read_csv("df_dev_surv_dly_clim_" + loc + ".txt", sep='\t', na_values="NaN", index_col=["Month", "Day"], parse_dates=True)
                
        if False: # if plotting
            for var_colName in climVar_colNames:
                print(var_colName)
                pltTitle_allMon = var_colName + "_" + loc + "_allMon"
                fileName = "TEST_interannualVar_" + strConvert_paren2hyphen(var_colName) + "_" + loc + "_allMon"
                plot_interannualVar_allMon(df_dev_surv_dly, df_dev_surv_dly_clim, var_colName, pltTitle_allMon, fileName)
                #for monthNum in monthNums:
                    #print(str(monthNum))
                    #pltTitle_oneMon = var_colName + "_" + loc + "_" + str(monthNum).zfill(2)
                    #fileName = "TEST_interannualVar_" + strConvert_paren2hyphen(var_colName) + "_" + loc + "_" + str(monthNum).zfill(2)
                    #plot_interannualVar_oneMon(df_dev_surv_dly, df_dev_surv_dly_clim, var_colName, pltTitle_oneMon, monthNum, fileName)

        if False: # if plotting
            for var_colName in contVar_colNames:
                print(var_colName)
                pltTitle_allMon = var_colName + "_" + loc + "_allMon"
                fileName = "TEST_interannualVar_" + strConvert_paren2hyphen(var_colName) + "_" + loc + "_allMon"
                plot_interannualVar_allMon(df_dev_surv_dly, df_dev_surv_dly_clim, var_colName, pltTitle_allMon, fileName)
                #for monthNum in monthNums:
                    #print(str(monthNum))
                    #pltTitle_oneMon = var_colName + "_" + loc + "_" + str(monthNum).zfill(2)
                    #fileName = "TEST_interannualVar_" + strConvert_paren2hyphen(var_colName) + "_" + loc + "_" + str(monthNum).zfill(2)
                    #plot_interannualVar_oneMon(df_dev_surv_dly, df_dev_surv_dly_clim, var_colName, pltTitle_oneMon, monthNum, fileName)

        if False: # if plotting
            for var_colName in devVar_colNames:
                print(var_colName)
                pltTitle_allMon = var_colName + "_" + loc + "_allMon"
                fileName = "TEST_interannualVar_" + strConvert_paren2hyphen(var_colName) + "_" + loc + "_allMon"
                plot_interannualVar_allMon(df_dev_surv_dly, df_dev_surv_dly_clim, var_colName, pltTitle_allMon, fileName)
                #for monthNum in monthNums:
                    #print(str(monthNum))
                    #pltTitle_oneMon = var_colName + "_" + loc + "_" + str(monthNum).zfill(2)
                    #fileName = "TEST_interannualVar_" + strConvert_paren2hyphen(var_colName) + "_" + loc + "_" + str(monthNum).zfill(2)
                    #plot_interannualVar_oneMon(df_dev_surv_dly, df_dev_surv_dly_clim, var_colName, pltTitle_oneMon, monthNum, fileName)

        if False: # if plotting
            for var_colName in survVar_colNames:
                print(var_colName)
                pltTitle_allMon = var_colName + "_" + loc + "_allMon"
                fileName = "TEST_interannualVar_" + strConvert_paren2hyphen(var_colName) + "_" + loc + "_allMon"
                plot_interannualVar_allMon(df_dev_surv_dly, df_dev_surv_dly_clim, var_colName, pltTitle_allMon, fileName)
                #for monthNum in monthNums:
                    #print(str(monthNum))
                    #pltTitle_oneMon = var_colName + "_" + loc + "_" + str(monthNum).zfill(2)
                    #fileName = "TEST_interannualVar_" + strConvert_paren2hyphen(var_colName) + "_" + loc + "_" + str(monthNum).zfill(2)
                    #plot_interannualVar_oneMon(df_dev_surv_dly, df_dev_surv_dly_clim, var_colName, pltTitle_oneMon, monthNum, fileName)






# Just do boxplots for initDay = 1
if False: # if plotting

    dataSrc = 'MERRA_IMERG'
    locs    = ['Negombo', 'Jaffna', 'NuwaraEliya']
    years   = list(range(2001, 2020+1))
    months  = list(range(1, 12+1))
    intv    = 'intv00'
    initDay    = 1
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000
    for loc in locs:
        
        df_list = []
        # Compile all data across months and years into one df
        for year in years:
            for month in months:
                _, subDir, fileStr_allSpecs = WHATCHEM_specs_to_fpath(dataSrc, loc, year, month, intv)
                popStr = 'day' + str(initDay) + '_' + str(initEggs) + 'e' + str(initLarvae) + 'l' + str(initPupae) + 'p'
                fpath = DIR_WHATCHEM_OUTPUT_DATA + subDir + "df_vectorPop_dly_" + fileStr_allSpecs + '-' + popStr + ".txt"
                df_vectorPop_dly = pd.read_csv(fpath, sep='\t', na_values="NaN", index_col=0, parse_dates=True, comment='#')
                df_list.append(df_vectorPop_dly)
        df_vectorPop_dly_all = pd.concat(df_list)
            
        # Compute df of mean across all years
        df_vectorPop_dly_all['MMDD'] = df_vectorPop_dly_all.index.strftime('%m%d')
        df_vectorPop_dly_mean = df_vectorPop_dly_all.groupby('MMDD').mean()
        df_vectorPop_dly_mean.index = pd.to_datetime('2000' + df_vectorPop_dly_mean.index, format='%Y%m%d')
        df_vectorPop_dly_mean = df_vectorPop_dly_mean.rename_axis("Date")

        # Do plotting
        if False:
            for var_colName in popVar_total_colNames:
                print(var_colName)
                pltTitle_allMon = var_colName + '_' + loc + '_allMon'
                fName = 'interannualVar_' + strConvert_paren2hyphen(var_colName) + '_' + loc + '_initDay' + str(initDay) + '_allMon_tseries'
                plot_interannualVar_allMon_tseries(df_vectorPop_dly_all, df_vectorPop_dly_mean, var_colName, pltTitle_allMon, fName)
        if True:
            var_colName = 'pop_A'
            pltTitle_allMon = var_colName + "_" + loc
            fName = 'interannualVar_' + strConvert_paren2hyphen(var_colName) + '_' + loc + '_initDay' + str(initDay) + '_allMon_boxplot'
            plot_interannualVar_allMon_boxplot(df_vectorPop_dly_all, var_colName, pltTitle_allMon, fName)
            var_colName = 'pop_A(E)'
            pltTitle_allMon = var_colName + "_" + loc
            fName = 'interannualVar_' + strConvert_paren2hyphen(var_colName) + '_' + loc + '_initDay' + str(initDay) + '_allMon_boxplot'
            plot_interannualVar_allMon_boxplot(df_vectorPop_dly_all, var_colName, pltTitle_allMon, fName)
            var_colName = 'pop_A(L)'
            pltTitle_allMon = var_colName + "_" + loc
            fName = 'interannualVar_' + strConvert_paren2hyphen(var_colName) + '_' + loc + '_initDay' + str(initDay) + '_allMon_boxplot'
            plot_interannualVar_allMon_boxplot(df_vectorPop_dly_all, var_colName, pltTitle_allMon, fName)
            var_colName = 'pop_A(P)'
            pltTitle_allMon = var_colName + "_" + loc
            fName = 'interannualVar_' + strConvert_paren2hyphen(var_colName) + '_' + loc + '_initDay' + str(initDay) + '_allMon_boxplot'
            plot_interannualVar_allMon_boxplot(df_vectorPop_dly_all, var_colName, pltTitle_allMon, fName)






## TEST THIS!
# make timeseries plots showing interannual variability of daily data
if False: # if plotting

    # parameters related to modeling pipeline
    dataSrc = 'MERRA_IMERG'
    locs    = ['Negombo']#['Negombo', 'Jaffna', 'NuwaraEliya']
    years   = list(range(2001, 2020+1))
    months  = list(range(1, 12+1))
    intv    = 'intv00'
    initDays   = list(range(1, 7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000
    var_colNames = ['TA_max (C)', 'TA (C)', 'PRCP (mm)', 'VPD (kPa)', 'TW_nanFill (C)', 'WH (mm)', 
                    'dev_E', 'surv_E', 'dev_L', 'surv_L', 'dev_P', 'surv_P', 
                    'pop_A(E)', 'pop_A(L)', 'pop_A(P)'] # subset of vars to actually plot
    var_pltNames = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in var_colNames]

    for loc in locs:
        print(loc)

        # load model pipeline data at location of interest
        # also: convert from daily to monthly
        df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae)
        
        # Compute df of mean across all years
        df_vectorPop_dly['MMDD'] = df_vectorPop_dly.index.strftime('%m%d')
        df_vectorPop_dly_mean = df_vectorPop_dly.groupby('MMDD').mean()
        df_vectorPop_dly_mean.index = pd.to_datetime('2000' + df_vectorPop_dly_mean.index, format='%Y%m%d') # use dummy year 2000 to make datetime index
        df_vectorPop_dly_mean = df_vectorPop_dly_mean.rename_axis("Date")

        for i_var, var_colName in enumerate(var_colNames):
            var_pltName = var_pltNames[i_var]
            print('  Var: ' + var_pltName)

            fName = 'tseries_' + var_pltName + '_' + loc
            plot_interannualVar_allMon_tseries(df_vectorPop_dly, df_vectorPop_dly_mean, var_colName, var_pltName, fName)






# make boxplots showing interannual variability of monthly climate/modeled data
if False: # if plotting

    # parameters related to modeling pipeline
    dataSrc = 'MERRA_IMERG'
    locs    = ['Negombo']#['Negombo', 'Jaffna', 'NuwaraEliya']
    years   = list(range(2001, 2020+1))
    months  = list(range(1, 12+1))
    intv    = 'intv00'
    initDays   = list(range(1, 7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000
    var_colNames = ['WH (mm)', 'surv_P']#['TA_max (C)', 'TA (C)', 'PRCP (mm)', 'VPD (kPa)', 'TW_nanFill (C)', 'WH (mm)', 
                   #'dev_E', 'surv_E', 'dev_L', 'surv_L', 'dev_P', 'surv_P', 
                   #'pop_A(E)', 'pop_A(L)', 'pop_A(P)',
                   #'pop_A(E)_frac', 'pop_A(L)_frac', 'pop_A(P)_frac'] # subset of vars to actually plot
    var_pltNames = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in var_colNames]

    # parameters related to dengue data
    dengueDataYrs_noOutbreaks_opts = [True]#[True, False] # whether to plot (dengue data years w/o outbreak years (2007-20 w/o 2017, 2019)) or (all years)
    dengueDataYrs = list(range(2007, 2020+1))
    outbreakYrs = [2017, 2019]

    for loc in locs:
        print(loc)

        # load model pipeline data at location of interest
        # also: convert from daily to monthly
        df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae)
        df_vectorPop_monthly = convert_df_vectorPop_dly2mon(df_vectorPop_dly)
        
        for i_var, var_colName in enumerate(var_colNames):
            var_pltName = var_pltNames[i_var]
            print('  Var: ' + var_pltName)

            # if the pop variable is a fractional one, we have to create it manually
            # (TODO: at some point we need to edit/rerun the modeling code so that these vars are natively in the output files)
            if var_colName in ['pop_A(E)_frac', 'pop_A(L)_frac', 'pop_A(P)_frac']:
                popTotal_E = len(initDays) * initEggs # this code only works for eggs AND larvae and pupae as long as initEggs = initLarvae = initPupae
                df_vectorPop_monthly[var_colName] = df_vectorPop_monthly[var_colName.replace('_frac', '')] / popTotal_E
            
            for dengueDataYrs_noOutbreaks in dengueDataYrs_noOutbreaks_opts:
                print(' dengueDataYrs_noOutbreaks: ' + str(dengueDataYrs_noOutbreaks))
                if dengueDataYrs_noOutbreaks:
                    fileStr_dengueDataYrs_noOutbreaks = '_' + str(dengueDataYrs[0]) + '-' + str(dengueDataYrs[-1]) + '_noOutbreaks'
                    exclYrs = list(set([year for year in years if year not in dengueDataYrs] + outbreakYrs))  # years to exclude from plot
                else:
                    fileStr_dengueDataYrs_noOutbreaks = '_' + str(years[0]) + '-' + str(years[-1])
                    exclYrs = []

                # remove data in outbreak years we want to exclude
                df_vectorPop_monthly_yrFiltered = df_vectorPop_monthly[~df_vectorPop_monthly.index.year.isin(exclYrs)]

                fName = 'boxplot_' + var_pltName + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                plot_interannualVar_allMon_boxplot(df_vectorPop_monthly_yrFiltered, var_colName, var_pltName, fName)






# make boxplots showing interannual variability of monthly max/min temperatures on the same plot
if False: # if plotting

    # parameters related to modeling pipeline
    dataSrc = 'MERRA_IMERG'
    locs    = ['Negombo', 'Jaffna', 'NuwaraEliya']
    years   = list(range(2001, 2020+1))
    months  = list(range(1, 12+1))
    intv    = 'intv00'
    initDays   = list(range(1, 7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000
    var1_colNames = ['TA_max (C)', 'TW_max_nanFill (C)'] # subset of vars to actually plot
    var2_colNames = ['TA_min (C)', 'TW_min_nanFill (C)'] # subset of vars to actually plot
    var1_pltNames = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in var1_colNames]
    var2_pltNames = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in var2_colNames]

    # parameters related to dengue data
    dengueDataYrs_noOutbreaks_opts = [True, False] # whether to plot (dengue data years w/o outbreak years (2007-20 w/o 2017, 2019)) or (all years)
    dengueDataYrs = list(range(2007, 2020+1))
    outbreakYrs = [2017, 2019]

    for loc in locs:
        print(loc)

        # load model pipeline data at location of interest
        # also: convert from daily to monthly
        df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae)
        df_vectorPop_monthly = convert_df_vectorPop_dly2mon(df_vectorPop_dly)

        for i_var, var1_colName in enumerate(var1_colNames):
            var1_pltName = var1_pltNames[i_var]
            var2_colName = var2_colNames[i_var]
            var2_pltName = var2_pltNames[i_var]

            print('  Vars: ' + var1_pltName + ', ' + var2_pltName)

            for dengueDataYrs_noOutbreaks in dengueDataYrs_noOutbreaks_opts:
                print(' dengueDataYrs_noOutbreaks: ' + str(dengueDataYrs_noOutbreaks))
                if dengueDataYrs_noOutbreaks:
                    fileStr_dengueDataYrs_noOutbreaks = '_' + str(dengueDataYrs[0]) + '-' + str(dengueDataYrs[-1]) + '_noOutbreaks'
                    exclYrs = list(set([year for year in years if year not in dengueDataYrs] + outbreakYrs))  # years to exclude from plot
                else:
                    fileStr_dengueDataYrs_noOutbreaks = '_' + str(years[0]) + '-' + str(years[-1])
                    exclYrs = []

                # remove data in outbreak years we want to exclude
                df_vectorPop_monthly_yrFiltered = df_vectorPop_monthly[~df_vectorPop_monthly.index.year.isin(exclYrs)]

                fName = 'boxplot_' + var1_pltName + '-' + var2_pltName + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                plot_interannualVar_allMon_boxplot_2var(df_vectorPop_monthly_yrFiltered, var1_colName, var1_pltName, var2_colName, var2_pltName, fName)






# make boxplots showing interannual variability of monthly dengue data (case count and incidence per 100k)
if False: # if plotting

    # parameters related to modeling pipeline
    dataSrc = 'MERRA_IMERG'
    locs    = ['Negombo', 'Jaffna', 'NuwaraEliya']
    years   = list(range(2007, 2020+1))   # note: not 2001 to 2020!
    months  = list(range(1, 12+1))
    intv    = 'intv00'
    initDays   = list(range(1, 7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000
    
    # parameters related to dengue data
    dataTypes = ['dengueData', 'dengueInc']
    noOutbreaks_opts = [True, False] # whether to include data from outbreak years
    outbreakYrs = [2017, 2019]

    fileDir_base = '/mnt/redwood/local_drive/chanud/RawData/epiReports/'
    for dataType in dataTypes:
        print(dataType)
        var_colName = dataType
        var_pltName = dataType
        # load dengue data
        if dataType == 'dengueData':
            file_dengue  = fileDir_base + 'dengueData_monthly_2007_2022.csv'
        elif dataType == ' dengueInc':
            file_dengue  = fileDir_base + 'dengueInc_monthly_2007_2022.csv'
        df_den_monthly = pd.read_csv(file_dengue, index_col="Date", parse_dates=True)
        df_den_monthly = df_den_monthly[df_den_monthly.index.year.isin(years)] # cut extraneous years
        
        for loc in locs:
            print(loc)
            if loc == 'Negombo':
                den_loc = 'Gampaha'
            elif loc in ['NuwaraEliya', 'Jaffna']:
                den_loc = loc
            else:
                print('Error: unknown loc.')
                exit

            # subset dengue data to location of interest
            series_den_monthly_loc = df_den_monthly[den_loc]
            series_den_monthly_loc.name = var_colName
            df_den_monthly_loc = series_den_monthly_loc.to_frame()

            for noOutbreaks in noOutbreaks_opts:
                print(' noOutbreaks: ' + str(noOutbreaks))
                if noOutbreaks:
                    fileStr_noOutbreaks = '_' + str(years[0]) + '-' + str(years[-1]) + '_noOutbreaks'
                    exclYrs = outbreakYrs  # outbreak years (to exclude from plot)
                else:
                    fileStr_noOutbreaks = '_' + str(years[0]) + '-' + str(years[-1])
                    exclYrs = []

                # remove data in outbreak years we want to exclude
                df_den_monthly_loc_yrFiltered = df_den_monthly_loc[~df_den_monthly_loc.index.year.isin(exclYrs)]

                fName = 'boxplot_' + var_pltName + '_' + loc + fileStr_noOutbreaks
                plot_interannualVar_allMon_boxplot(df_den_monthly_loc_yrFiltered, var_colName, var_pltName, fName)





# do scatter plot of vars against adult pops
if False: # if plotting

    dataSrc = 'MERRA_IMERG'
    locs    = ['Jaffna'] # ['Negombo', 'Jaffna', 'NuwaraEliya']
    years   = list(range(2001, 2020+1))
    months  = [12]#list(range(1, 12+1))
    intv    = 'intv00'
    initDays   = list(range(1, 7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000
    
    df_list = []
    # Compile all data across years into one df
    loc = locs[0]
    month = months[0]
    for year in years:
        df_allInitDays = combine_dfs_overInitDays(dataSrc, loc, year, month, intv, initDays)
        df_list.append(df_allInitDays, axis=1)
    df_vectorPop_dly_all = pd.concat(df_list)

    df_vectorPop_monthly_all = convert_df_vectorPop_dly2mon(df_vectorPop_dly_all)
    df_vectorPop_monthly_all = df_vectorPop_monthly_all[df_vectorPop_monthly_all.index.month == month]

    varX_colName = 'surv_E'
    varX_pltName = strConvert_paren2hyphen(varX_colName)
    varY_colName = 'pop_A(E)'
    varY_pltName = strConvert_paren2hyphen(varY_colName)
    pltTitle = ''
    fileName = 'interannualVar_monthly_' + varY_pltName + '-' + varX_pltName + '_' + loc + '_month' + str(month).zfill(2)
    plot_varVsVar_byYear(df_vectorPop_monthly_all, varX_colName, varX_pltName, varY_colName, varY_pltName, pltTitle, fileName)






# do scatter plot of vars against dengue cases
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
    var_colNames = ['pop_A(E)', 'pop_A(L)', 'pop_A(P)'] # subset of vars to actually plot
    var_pltNames = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in var_colNames]
    
    # parameters related to dengue data
    lagTimes = [0, 1, 2] # in months
    den_colName_base = 'denCases' # we will rename columns to this string + the lag time (original name is the district)
    noOutbreaks_opts = [True, False] # whether to include data from outbreak years
    outbreakYrs = [2017, 2019]

    # parameters related to modeling pipeline + dengue data
    startDate = '2007-01-01'  # start date (first year of dengue data)
    endDate   = '2020-12-31'  # end date (last year of downloaded MERRA/IMERG data)

    # load dengue data
    fileDir_base = '/mnt/redwood/local_drive/chanud/RawData/epiReports/'
    file_dengue  = fileDir_base + 'dengueData_monthly_2007_2022.csv'
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
        df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae)
        df_vectorPop_monthly = convert_df_vectorPop_dly2mon(df_vectorPop_dly)
        df_vectorPop_monthly = df_vectorPop_monthly[startDate:endDate]

        # combine model pipeline data and dengue data into one df
        df_vectorPop_den_monthly = df_vectorPop_monthly.join(df_den_monthly_loc_wLags)

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

            for i_var, var_colName in enumerate(var_colNames):
                print('  Var: ' + var_colName)
                var_pltName = var_pltNames[i_var]

                for i_lag, lagTime in enumerate(lagTimes):
                    print('   Lag (mo): ' + str(lagTime))

                    den_colName = den_colName_base + '-lag' + str(lagTime) # this must match col name format set in create_df_withLags()
                    den_pltName = den_colName # just reusing the column name here for simplicity

                    fileName = 'scatter_' + var_pltName + '_' + den_pltName + '_' + loc + fileStr_noOutbreaks
                    plot_scatter_colorMon(df_vectorPop_den_monthly_yrFiltered, var_colName, var_pltName, den_colName, den_pltName, fileName)

                    for month in months:
                        print('  ' + str(month))
                        
                        df_vectorPop_den_monthly_yrFiltered_oneMon = df_vectorPop_den_monthly_yrFiltered[df_vectorPop_den_monthly_yrFiltered.index.month == month]
#
                        fileName = 'scatter_' + var_pltName + '_' + den_colName + '_' + loc + '_month' + str(month).zfill(2) + fileStr_noOutbreaks
                        plot_scatter_colorMon(df_vectorPop_den_monthly_yrFiltered_oneMon, var_colName, var_pltName, den_colName, den_pltName, fileName)






# do timeseries plots of var, highlighting the years corresponding to max/min climate endmember values
# note: unlike other scripts in this file, here we don't do any removal of outbreak years. I don't think it'd
#       really make sense to do that bc then we might lose out on an endmember year?
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
    pltVar_colNames = ['TA_max (C)', 'TA (C)', 'PRCP (mm)', 'VPD (kPa)', 'TW_nanFill (C)', 'WH (mm)', 
                       'dev_E', 'surv_E', 'dev_L', 'surv_L', 'dev_P', 'surv_P', 
                       'pop_A(E)', 'pop_A(L)', 'pop_A(P)'] # subset of variables to actually plot
    pltVar_pltNames = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in pltVar_colNames]
    endmemberVar_colNames = ['TA (C)', 'PRCP (mm)'] # variables to use as endmembers (highlight timeseries for which these are min/max)
    endmemberVar_pltNames = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in endmemberVar_colNames]
    
    # define a custom color palette for the high/low values of each endmember var
    # from colorbrewer (4-color, qualitative, colorblind-friendly), but made a bit darker (and with a separate color for endmember combos)
    dict_endmemberColors = {'Hi-TA': '#23901c',
                            'Lo-TA': '#a2cf7a',
                            'Hi-PRCP': '#0f68a4',
                            'Lo-PRCP': '#96bed3',
                            'Lo-TA, Hi-PRCP': '#224455',
                            'Hi-TA, Lo-PRCP': '#225544'}

    for loc in locs:
        print(loc)
        
        # load model pipeline data at location of interest
        # also: convert from daily to monthly (will use this dataset to find min/max monthly values)
        df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae)
        df_vectorPop_monthly = convert_df_vectorPop_dly2mon(df_vectorPop_dly)

        for month in months:
            print(' ' + str(month))
            
            # subset data to month of interest
            df_vectorPop_dly_oneMon = df_vectorPop_dly[df_vectorPop_dly.index.month == month]
            df_vectorPop_monthly_oneMon = df_vectorPop_monthly[df_vectorPop_monthly.index.month == month]

            # For each climate endmember var, find years with min/max monthly value and store these in dicts
            dict_endmemberYrs = {}
            for i_endmember, endmemberVar_colName in enumerate(endmemberVar_colNames):
                maxYr = df_vectorPop_monthly_oneMon[endmemberVar_colName].idxmax().year
                minYr = df_vectorPop_monthly_oneMon[endmemberVar_colName].idxmin().year
                
                # update dict, appending variable name to existing value if the key (year) already exists
                # (e.g., if 2022 is the max year for both 'TA' and 'PRCP', the dict entry will be 2022: 'Hi-TA, Hi-PRCP')
                endmemberVar_pltName = endmemberVar_pltNames[i_endmember]
                if maxYr in dict_endmemberYrs:
                    dict_endmemberYrs[maxYr] += ', ' + 'Hi-' + endmemberVar_pltName
                else:
                    dict_endmemberYrs[maxYr] = 'Hi-' + endmemberVar_pltName

                if minYr in dict_endmemberYrs:
                    dict_endmemberYrs[minYr] += ', ' + 'Lo-' + endmemberVar_pltName
                else:
                    dict_endmemberYrs[minYr] = 'Lo-' + endmemberVar_pltName
            
            print(dict_endmemberYrs)

            # For each variable we want to plot, do the plot!
            for i_pltVar, pltVar_colName in enumerate(pltVar_colNames):
                pltVar_pltName = pltVar_pltNames[i_pltVar]
                print('  Plotted variable: ' + pltVar_pltName)
                print('  Climate endmembers: ' + ', '.join(endmemberVar_pltNames))

                fileName = 'tseries_' + pltVar_pltName + '_climEndmembers_' + loc + '_month' + str(month).zfill(2)
                plot_varVsTime_climEndmembers_oneMon(df_vectorPop_dly_oneMon, dict_endmemberYrs, dict_endmemberColors, pltVar_colName, pltVar_pltName, fileName)






# do correlogram of monthly values EXCEPT dengue cases
# Note: the timespan here is the entire MERRA/IMERG record 
# that we downloaded: 2001-2020.
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

    for loc in locs:
        print('  ' + loc)

        # load model pipeline data at location of interest
        # also: convert from daily to monthly and drop unwanted cols
        df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae)
        df_vectorPop_monthly = convert_df_vectorPop_dly2mon(df_vectorPop_dly)
        cols2drop = ['pop_E', 'pop_L', 'pop_P', 'pop_A', 'pop_L(E)', 'pop_L(L)', 'pop_P(E)', 'pop_P(L)', 'pop_P(P)',
                     'surv_temp_E','surv_temp_L','surv_temp_P','surv_desicc_E','surv_desicc_L','surv_desicc_P',
                     'WH_del (mm)', 'TW (C)', 'TW_min (C)', 'TW_max (C)']
        df_vectorPop_monthly = df_vectorPop_monthly.drop(cols2drop, axis=1)

        # plot correlogram for all data (across months and years)
        fileName = 'correlogram_' + loc
        plot_correlogram_allVars(df_vectorPop_monthly, fileName)

        for month in months:
            print('   ' + str(month))

            # subset data to month of interest
            df_vectorPop_monthly_oneMon = df_vectorPop_monthly[df_vectorPop_monthly.index.month == month]

            # plot correlogram for one calendar month
            fileName = 'correlogram_' + loc + '_month' + str(month).zfill(2)
            plot_correlogram_allVars(df_vectorPop_monthly_oneMon, fileName)






# do correlogram of monthly values INCLUDING dengue cases
# Note: the timespan here is the overlap between our MERRA/IMERG record 
# and the dengue case data: 2007-2020.
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

    # parameters related to dengue data
    lagTimes = [0, 1, 2] # in months
    den_colName_base = 'denCases' # we will rename columns to this string + the lag time (original name is the district)
    noOutbreaks_opts = [True, False] # whether to include data from outbreak years
    outbreakYrs = [2017, 2019]

    # parameters related to modeling pipeline + dengue data
    startDate = '2007-01-01'  # start date (first year of dengue data)
    endDate   = '2020-12-31'  # end date (last year of downloaded MERRA/IMERG data)

    # load dengue data
    fileDir_base = '/mnt/redwood/local_drive/chanud/RawData/epiReports/'
    file_dengue  = fileDir_base + 'dengueData_monthly_2007_2022.csv'
    df_den_monthly = pd.read_csv(file_dengue, index_col="Date", parse_dates=True)

    for loc in locs:
        print('  ' + loc)
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
        # also: convert from daily to monthly, subset to time of interest, and drop unwanted cols
        df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae)
        df_vectorPop_monthly = convert_df_vectorPop_dly2mon(df_vectorPop_dly)
        df_vectorPop_monthly = df_vectorPop_monthly[startDate:endDate]
        cols2drop = ['pop_E', 'pop_L', 'pop_P', 'pop_A', 'pop_L(E)', 'pop_L(L)', 'pop_P(E)', 'pop_P(L)', 'pop_P(P)',
                     'surv_temp_E','surv_temp_L','surv_temp_P','surv_desicc_E','surv_desicc_L','surv_desicc_P',
                     'WH_del (mm)', 'TW (C)', 'TW_min (C)', 'TW_max (C)']
        df_vectorPop_monthly = df_vectorPop_monthly.drop(cols2drop, axis=1)

        # combine model pipeline data and dengue data into one df
        df_vectorPop_den_monthly = df_vectorPop_monthly.join(df_den_monthly_loc_wLags)

        for noOutbreaks in noOutbreaks_opts:
            if noOutbreaks:
                fileStr_noOutbreaks = '_noOutbreaks'
                exclYrs = outbreakYrs  # outbreak years (to exclude from plot)
            else:
                fileStr_noOutbreaks = ''
                exclYrs = []

            # remove data in outbreak years we want to exclude
            df_vectorPop_den_monthly_yrFiltered = df_vectorPop_den_monthly[~df_vectorPop_den_monthly.index.year.isin(exclYrs)]
            
            # plot correlogram for all data (across months and years)
            fileName = 'correlogram_withDen_' + loc + fileStr_noOutbreaks
            plot_correlogram_allVars(df_vectorPop_den_monthly_yrFiltered, fileName)

            for month in months:
                print('   ' + str(month))

                # subset data to month of interest
                df_vectorPop_den_monthly_yrFiltered_oneMon = df_vectorPop_den_monthly_yrFiltered[df_vectorPop_den_monthly_yrFiltered.index.month == month]

                # plot correlogram for one calendar month
                fileName = 'correlogram_withDen_' + loc + '_month' + str(month).zfill(2) + fileStr_noOutbreaks
                plot_correlogram_allVars(df_vectorPop_den_monthly_yrFiltered_oneMon, fileName)






# do quantile stacked bar plot of monthly values of temp/prcp against adult pops
if True:

    # parameters related to modeling pipeline
    dataSrc = 'MERRA_IMERG'
    locs = ['Jaffna']#['Negombo', 'NuwaraEliya', 'Jaffna']
    years = list(range(2001, 2020+1))
    months = list(range(1,12+1))
    intv = 'intv00'
    initDays = list(range(1,7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000
    varX1_colName = 'TA (C)'
    varX1_pltName = strConvert_paren2hyphen(strConvert_parenUnitsRemove(varX1_colName))
    varX2_colName = 'PRCP (mm)'
    varX2_pltName = strConvert_paren2hyphen(strConvert_parenUnitsRemove(varX2_colName))
    varY_colNames = ['pop_A(L)']# ['pop_A(E)', 'pop_A(L)', 'pop_A(P)'] # subset of vars to actually plot
    varY_pltNames = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in varY_colNames]

    # parameters related to dengue data
    noOutbreaks_opts = [False]#[True, False] # whether to include data from outbreak years
    outbreakYrs = [2017, 2019]


    for loc in locs:
        print(loc)

        # load model pipeline data at location of interest
        # also: convert from daily to monthly
        df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae)
        df_vectorPop_monthly = convert_df_vectorPop_dly2mon(df_vectorPop_dly)

        for noOutbreaks in noOutbreaks_opts:
            print(' noOutbreaks: ' + str(noOutbreaks))
            if noOutbreaks:
                fileStr_noOutbreaks = '_noOutbreaks'
                exclYrs = outbreakYrs  # outbreak years (to exclude from plot)
            else:
                fileStr_noOutbreaks = ''
                exclYrs = []

            # remove data in outbreak years we want to exclude
            df_vectorPop_monthly_yrFiltered = df_vectorPop_monthly[~df_vectorPop_monthly.index.year.isin(exclYrs)]
            
            for i_varY, varY_colName in enumerate(varY_colNames):
                varY_pltName = varY_pltNames[i_varY]
                print('  Var: ' + varY_pltName)

                fileName = 'bar_binaryTercile_' + varX1_pltName + '-' + varX2_pltName + '_' + varY_pltName + '_' + loc + fileStr_noOutbreaks
                plot_bar_binaryX1X2_tercileY_splitByMon(df_vectorPop_monthly_yrFiltered, varX1_colName, varX1_pltName, varX2_colName, varX2_pltName, varY_colName, varY_pltName, 'allMon', fileName)





# do tercile stacked bar plot of monthly values of variables against dengue incidence
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
    var_colNames = ['TA_min (C)', 'TA (C)', 'TA_max (C)', 'PRCP (mm)', 'TW_min_nanFill (C)', 'TW_nanFill (C)', 'TW_max_nanFill (C)', 'WH (mm)',
                   'dev_E', 'surv_E', 'dev_L', 'surv_L', 'dev_P', 'surv_P', 'pop_A(E)', 'pop_A(L)', 'pop_A(P)'] # subset of vars to actually plot
    var_pltNames = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in var_colNames]

    # parameters related to dengue data
    lagTimes = [0, 1, 2]             # in months
    den_colName_base = 'denInc'    # we will rename columns to this string + the lag time (original name is the district)
    noOutbreaks_opts = [True, False] # whether to include data from outbreak years
    outbreakYrs = [2017, 2019]

    # other parameters
    seasons = ['allMon']#['allMon', 'NEM', 'FIM', 'SWM', 'SIM']

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
        df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae)
        df_vectorPop_monthly = convert_df_vectorPop_dly2mon(df_vectorPop_dly)
        df_vectorPop_monthly = df_vectorPop_monthly[startDate:endDate]

        # combine model pipeline data and dengue data into one df
        df_vectorPop_den_monthly = df_vectorPop_monthly.join(df_den_monthly_loc_wLags)

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
            

            for i_var, var_colName in enumerate(var_colNames):
                var_pltName = var_pltNames[i_var]
                print('  Var: ' + var_pltName)

                for i_lag, lagTime in enumerate(lagTimes):
                    print("   Lag (mo): " + str(lagTime))

                    den_colName = den_colName_base + '-lag' + str(lagTime) # this must match col name format set in create_df_withLags()
                    den_pltName = den_colName # just reusing the column name here for simplicity
                    
                    for season in seasons:
                        print('     Season: ' + season)
                        if season == 'allMon':
                            fileStr_season = ''
                        else:
                            fileStr_season = '_' + season

                        fileName = 'bar_tercile_splitbyMon_' + var_pltName + '_' + den_pltName + '_' + loc + fileStr_season + fileStr_noOutbreaks
                        plot_bar_tercileXY_splitByMon(df_vectorPop_den_monthly_yrFiltered, var_colName, var_pltName, den_colName, den_pltName, season, fileName)


