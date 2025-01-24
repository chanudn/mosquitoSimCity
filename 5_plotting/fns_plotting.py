#*******************************************************
# helperFns_WHATCHEM_plotting.py
#*******************************************************
# Revision history
#    ????-??-??, CNY: redoing plotting code in Python
#    2024-10-19, CNY: reworked this while reorganizing code
#
#
# Functions:
#   plot_waterVarvsTime_monthSummary()
#   plot_varVsVar()
#   plot_varVsVar_byYear()
#   plot_varVsTime_climEndmembers_oneMon()
#   plot_interannualVar_oneMon()
#   plot_interannualVar_allMon_OLD()
#   plot_interannualVar_allMon_tseries()
#   plot_interannualVar_allMon_boxplot()
#   plot_interannualVar_allMon_boxplot_2var()
#   plot_varVsDengue_scatter()
#   plot_correlogram_allVars()
#   plot_tercilebar_allMon()
#   MORE!


# import required packages
# ************************************************
import matplotlib.pyplot as plt
from matplotlib import dates
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import FuncFormatter, FixedLocator, NullFormatter
import pandas as pd
import numpy as np
import math
from datetime import date
import calendar
import os
import sys
import re
import seaborn as sns
from scipy.stats import pearsonr, chisquare, chi2_contingency, kendalltau
from cycler import cycler


# import constants and custom functions
# ************************************************
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/0_constants") # directory of constants file
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/4_analysis") # directory of helper functions file
from CONST_FILE_DIR import DIR_WHATCHEM_OUTPUT_DATA, DIR_WHATCHEM_OUTPUT_PLOTS
from fns_helper import get_tbounds_yyyymmddhh, read_WHATCHEM_output, categorizeSeries_LowMedHigh_byMon, categorizeSeries_LowHigh_byMon




# Plotting functions
# ************************************************

# Plot WHATCH'EM water height vs time for a given location, years, and month.
#   TO DO: add in doClim into args
#   years  - vector of all years to plot
#   months - vector of all months to consider (represented as integers)
#   intv   - intv00, intvF0, intv0E, or intvFE (0 is no intv, F is fill, E is empty)
#   waterVar - variable to plot (options: "wh", "tw")
def plot_waterVarvsTime_monthSummary(dataSrc, loc, years, months, intv, waterVar):

  # change default color cycle
  plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.tab20.colors)

  # get parameters related to variable of interest
  if waterVar == "wh":
    waterVar_str = "wh"
    waterVar_dfColName = "WH (mm)"
    waterVar_ylabel = "Water height [mm]"
    heightThreshold_min = 15              # [mm]; height below which model can't compute temp (not min for immatures' survival)
    heightThreshold_max = 368             # [mm]; top of container
    thresholds = [heightThreshold_min, heightThreshold_max] # these will be plotted as horiz lines
    waterVar_yMin   = 0                   # [mm]; min water height for plot
    waterVar_yMax   = heightThreshold_max # [mm]; max water height for plot
    waterVar_sd_yMax = 100                # [mm]; max standard deviation of water height for plot
  elif waterVar == "tw":
    waterVar_str = "tw"
    waterVar_dfColName = "TW (C)"
    waterVar_ylabel = "Water temperature [" + chr(176) + "C]" 
    tempThreshold_max = 40               # [deg C]; when Focks et al. survival factor ~ 0.5 for both eggs and larvae/pupae
    thresholds = [tempThreshold_max]     # these will be plotted as horiz lines
    waterVar_yMin   = 10                 # [deg C]; min temperature for plot
    waterVar_yMax   = 50                 # [deg C]; max temperature for plot
    waterVar_sd_yMax = 6                 # [deg C]; max standard deviation of temperature for plot
  else:
    raise Exception("Invalid waterVar value in plot_waterVarvsTime_monthSummary()")

  # get remaining parameters
  yrStrt = min(years)
  yrEnd  = max(years)
  container = "Bucket"
  refl      = "Gray"
  shade     = "HalfShade"
  contSpecs = container + "-" + refl + "-" + shade
  if dataSrc == "MERRA_IMERG":
    clouds = "Observed"
  elif dataSrc == "GEOS_S2S":
    clouds = "PartlyCloudy"    # "Clear", "Overcast", "PartlyCloudy"


  # read in WHATCH'EM output files as dataframes
  df_in_list   = []   # holds dfs of waterVar timeseries for each calendar month
  df_mean_list = []   # holds dfs of waterVar timeseries mean for each calendar month
  df_sd_list   = []   # holds dfs of waterVar timeseries stdev for each calendar month
  for month in months:
    df_in_list_monthN = [] # holds dfs of waterVar timeseries for a given month
    for year in years:

      tStrt, tEnd = get_tbounds_yyyymmddhh(year, month)
      tspan       = str(tStrt) + "_" + str(tEnd)                 # "yyyymmddhh_yyyymmddhh"
      subDir_in   = dataSrc + '-' + loc + "-" + tspan + "/"
      dir_in      = DIR_WHATCHEM_OUTPUT_DATA + subDir_in

      fileStr_allSpecs  = dataSrc + "-" + loc + "-" +  tspan + "-" + intv + "-" + contSpecs
      if clouds != "Observed": # if clouds are specified (rather than from data), the filename incl the 'clouds' string
        fileStr_allSpecs = fileStr_allSpecs + "-" + clouds
      filename_in_noExt = "container_output_" + fileStr_allSpecs
      filepath_in       = dir_in + filename_in_noExt + ".txt"

      df_in = read_WHATCHEM_output(filepath_in)
      df_in_list_monthN.append(df_in)

    df_in_list.append(df_in_list_monthN)  # store list (of individual timeseries) within list (in which one entry is one calendar month)

    # calculate mean and stdev of model output (across years)
    df_concat = pd.concat(df_in_list_monthN)
    df_concat.index = [x.replace(year=2000) for x in df_concat.index] # substitute in dummy year 2000 (chosen bc it's a leap year)
    df_concat_grp = df_concat.groupby([df_concat.index])
    df_mean_list.append(df_concat_grp.mean())
    df_sd_list.append(df_concat_grp.std())


  # start plotting
  if True:
    for month in months:
      
      # define timespans based on years and month
      tspan_out = str(yrStrt) + "_" + str(yrEnd) + "-" + str(month).zfill(2)         # e.g., "2000_2020-01"
      tspan_4plt = str(yrStrt) + "-" + str(yrEnd) + " " + calendar.month_abbr[month] # e.g., "2000-2020 Jan"
    
      # define output filepath
      subDir_out = dataSrc + '-' + loc + "-" + tspan_out + "/"
      dir_out = DIR_WHATCHEM_OUTPUT_PLOTS + subDir_out
      fileStr_out_allSpecs = dataSrc + "-" + loc + "-" +  tspan_out + "-" + intv + "-" + contSpecs
      if clouds != "Observed": # if clouds are specified (rather than from data), the filename incl the 'clouds' string
        fileStr_out_allSpecs = fileStr_out_allSpecs + "-" + clouds


      ## PLOT 1
      # plot all individual timeseries for a given calendar month
      if False:

        filename_out_noExt = waterVar_str + "-time_allRuns_" + fileStr_out_allSpecs
        filepath_out = dir_out + filename_out_noExt + ".png"

        df_in_list_monthN = df_in_list[month-1]

        fig, ax = plt.subplots()

        # plot data
        for threshold in thresholds:
          ax.axhline(y = threshold, color='r', linestyle = '--') # horizontal
        for df in df_in_list_monthN:
          df.index = [x.replace(year=2000) for x in df.index] # substitute in dummy year 2000 (chosen bc it's a leap year)
          ax.plot(df[waterVar_dfColName], c='k', lw=0.5)

        ax.set_title(loc + ", " + tspan_4plt)
        ax.set_xlabel("Day of Month")
        ax.set_ylabel(waterVar_ylabel)
        ax.set_ylim(waterVar_yMin, waterVar_yMax)

        # label x-axis with day number, and only at the beginning of every week (i.e., every 7 days)
        ax.xaxis.set_major_locator(dates.DayLocator(interval=7))
        ax.xaxis.set_major_formatter(dates.DateFormatter('%d'))

        if not os.path.isdir(dir_out):
          os.makedirs(dir_out)
        plt.savefig(filepath_out, dpi=300, bbox_inches="tight")
        plt.close()


      ## PLOT 2
      # plot mean timeseries along with individual timeseries for a given calendar month
      if True:

        filename_out_noExt = waterVar_str + "-time_meanWithAllRuns_" + fileStr_out_allSpecs
        filepath_out = dir_out + filename_out_noExt + ".png"

        df_in_list_monthN = df_in_list[month-1]
        df_mean = df_mean_list[month-1]
        waterVar_mean = df_mean[waterVar_dfColName]

        fig, ax = plt.subplots()

        # plot data
        for threshold in thresholds:
          ax.axhline(y = threshold, color='r', linestyle = '--') # horizontal
        for df in df_in_list_monthN:
          df.index = [x.replace(year=2000) for x in df.index] # substitute in dummy year 2000 (chosen bc it's a leap year)
          ax.plot(df[waterVar_dfColName], c='k', lw=0.5, alpha=0.3)
        ax.plot(waterVar_mean)

        ax.set_title(loc + ", " + tspan_4plt)
        ax.set_xlabel("Day of Month")
        ax.set_ylabel(waterVar_ylabel)
        ax.set_ylim(waterVar_yMin, waterVar_yMax)

        # label x-axis with day number, and only at the beginning of every week (i.e., every 7 days)
        ax.xaxis.set_major_locator(dates.DayLocator(interval=7))
        ax.xaxis.set_major_formatter(dates.DateFormatter('%d'))

        if not os.path.isdir(dir_out):
          os.makedirs(dir_out)
        plt.savefig(filepath_out, dpi=300, bbox_inches="tight")
        plt.close()


      ## PLOT 3
      # plot mean timeseries with 2*stdev as shading for a given calendar month
      if False:

        filename_out_noExt = waterVar_str + "-time_meanWith2sd_" + fileStr_out_allSpecs
        filepath_out = dir_out + filename_out_noExt + ".png"

        df_mean = df_mean_list[month-1]
        waterVar_mean = df_mean[waterVar_dfColName]
        df_sd   = df_sd_list[month-1]
        waterVar_sd   = df_sd[waterVar_dfColName]

        fig, ax = plt.subplots()

        # plot data
        for threshold in thresholds:
          ax.axhline(y = threshold, color='r', linestyle = '--') # horizontal
        ax.fill_between(x=df_sd.index, y1=waterVar_mean-(2*waterVar_sd), y2=waterVar_mean+(2*waterVar_sd), alpha=0.4)
        ax.plot(waterVar_mean)

        ax.set_title(loc + ", " + tspan_4plt)
        ax.set_xlabel("Day of Month")
        ax.set_ylabel(waterVar_ylabel)
        ax.set_ylim(waterVar_yMin, waterVar_yMax)

        # label x-axis with day number, and only at the beginning of every week (i.e., every 7 days)
        ax.xaxis.set_major_locator(dates.DayLocator(interval=7))
        ax.xaxis.set_major_formatter(dates.DateFormatter('%d'))

        if not os.path.isdir(dir_out):
          os.makedirs(dir_out)
        plt.savefig(filepath_out, dpi=300, bbox_inches="tight")
        plt.close()


    # redefine values for next plot, which will combine data from all months
    # (compared to previous plots that were month-specific)

    # define timespans based on years
    tspan_out = str(yrStrt) + "_" + str(yrEnd)
    tspan_4plt = str(yrStrt) + "-" + str(yrEnd)
  
    # define output filepath
    subDir_out = dataSrc + '-' + loc + "-" + tspan_out + "/"
    dir_out = DIR_WHATCHEM_OUTPUT_PLOTS + subDir_out
    fileStr_out_allSpecs = dataSrc + "-" + loc + "-" +  tspan_out + "-" + intv + "-" + contSpecs
    if clouds != "Observed": # if clouds are specified (rather than from data), the filename incl the 'clouds' string
      fileStr_out_allSpecs = fileStr_out_allSpecs + "-" + clouds


    ## PLOT 4
    # plot mean timeseries for all calendar months
    if False:

      filename_out_noExt = waterVar_str + "-time_mean_" + fileStr_out_allSpecs
      filepath_out = dir_out + filename_out_noExt + ".png"

      fig, ax = plt.subplots()

      # plot data
      for month in months:
        df = df_mean_list[month-1]
        df.index = [x.replace(month=1) for x in df.index] # substitute in dummy month 1 (chosen bc January has max possible number of days: 31)
        ax.plot(df[waterVar_dfColName], label=calendar.month_abbr[month])

      ax.set_title(loc + ", " + tspan_4plt)
      ax.set_xlabel("Day of Month")
      ax.set_ylabel(waterVar_ylabel)
      ax.set_ylim(waterVar_yMin, waterVar_yMax)
      ax.legend(bbox_to_anchor=(1.0, 1.0), loc="upper left")

      # label x-axis with day number, and only at the beginning of every week (i.e., every 7 days)
      ax.xaxis.set_major_locator(dates.DayLocator(interval=7))
      ax.xaxis.set_major_formatter(dates.DateFormatter('%d'))

      if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
      plt.savefig(filepath_out, dpi=300, bbox_inches="tight")
      plt.close()


    ## PLOT 5
    # plot stdev timeseries for all calendar months
    if False:

      filename_out_noExt = waterVar_str + "-time_stdev_" + fileStr_out_allSpecs
      filepath_out = dir_out + filename_out_noExt + ".png"

      fig, ax = plt.subplots()

      # plot data
      for month in months:
        df = df_sd_list[month-1]
        df.index = [x.replace(month=1) for x in df.index] # substitute in dummy month 1 (chosen bc January has max possible number of days: 31)
        ax.plot(df[waterVar_dfColName], label=calendar.month_abbr[month])

      ax.set_title(loc + ", " + tspan_4plt)
      ax.set_xlabel("Day of Month")
      ax.set_ylabel(waterVar_ylabel + " (STDEV)")
      ax.set_ylim(0, waterVar_sd_yMax)
      ax.legend(bbox_to_anchor=(1.0, 1.0), loc="upper left")

      # label x-axis with day number, and only at the beginning of every week (i.e., every 7 days)
      ax.xaxis.set_major_locator(dates.DayLocator(interval=7))
      ax.xaxis.set_major_formatter(dates.DateFormatter('%d'))

      if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
      plt.savefig(filepath_out, dpi=300, bbox_inches="tight")
      plt.close()




# Plot one variable against another as a scatter plot.
def plot_varVsVar(df_in, varX_colName, varX_pltName, varY_colName, varY_pltName, pltTitle, fileName, addJitter, isClim):
    dir_out = "/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/simulateVectorSurvival/" # TEMPORARY!

    fig, ax = plt.subplots() # create a new figure and axis for plot

    df = df_in.copy(deep=True) # create a temporary copy of the df so that we don't alter the original

    # Extract month from index (which we will use to color the points)
    if isClim: # if climatology, df has a MultiIndex
        df['Month'] = df.index.get_level_values('Month')
    else: # if not climatology, df has a datetime index
        df['Month'] = df.index.month
    
    # Define a colormap for months
    cmap = ListedColormap(plt.get_cmap("tab20").colors[:12])

    # Create a scatter plot
    xVar = df[varX_colName]
    yVar = df[varY_colName]
    if addJitter:
        xVar = xVar + np.random.normal(scale=0.1, size=len(xVar))
    scatter = ax.scatter(xVar, yVar, c=df['Month'], cmap=cmap, s=6)

    # Add labels and a colorbar
    ax.set_xlabel(varX_pltName)
    ax.set_ylabel(varY_pltName)
    ax.set_title(pltTitle)
    tick_pos = np.arange(1.5, 12.5, (12.5-1.5)/12) # set tick positions at the center of each color interval
    month_labels = [calendar.month_abbr[month_num] for month_num in np.arange(1,13)]  # ["Jan", "Feb", ...]
    cbar = plt.colorbar(scatter, ax=ax, label='Month')
    cbar.set_ticks(ticks = tick_pos, labels = month_labels)
    cbar.ax.invert_yaxis() # Reverse the direction of the colorbar (to have "Jan" on top)

    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)




# Plot one variable against another as a scatter plot.
def plot_varVsVar_byYear(df_in, varX_colName, varX_pltName, varY_colName, varY_pltName, pltTitle, fileName):
    dir_out = "/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/simulateVectorSurvival/" # TEMPORARY!

    fig, ax = plt.subplots() # create a new figure and axis for plot

    df = df_in.copy(deep=True) # create a temporary copy of the df so that we don't alter the original

    df['Year'] = df.index.year
    
    # Define a colormap for years
    numYrs = 20
    cmap = ListedColormap(plt.get_cmap('viridis')(np.linspace(0, 1, numYrs)))

    # Create a scatter plot
    xVar = df[varX_colName]
    yVar = df[varY_colName]
    scatter = ax.scatter(xVar, yVar, c=df['Year'], cmap=cmap, s=6)

    # Add labels and a colorbar
    ax.set_xlabel(varX_pltName)
    ax.set_ylabel(varY_pltName)
    ax.set_title(pltTitle)
    tick_pos = np.arange(1.5, 12.5, (12.5-1.5)/12) # set tick positions at the center of each color interval
    month_labels = [calendar.month_abbr[month_num] for month_num in np.arange(1,13)]  # ["Jan", "Feb", ...]
    cbar = plt.colorbar(scatter, ax=ax, label='Month')
    cbar.set_ticks(ticks = tick_pos, labels = month_labels)
    cbar.ax.invert_yaxis() # Reverse the direction of the colorbar (to have "Jan" on top)

    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)




# Plot one variable for a given month, showing each individual year of data as well as their mean
def plot_interannualVar_oneMon(df_dly_INIT, df_dly_clim_INIT, var_colName, var_pltName, monthNum, fileName):
    dir_out = "/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/" # TEMPORARY!

    fig, ax = plt.subplots() # create a new figure and axis for plot

    df_dly = df_dly_INIT.copy(deep=True)           # create a temporary copy of the df so that we don't alter the original
    df_dly_clim = df_dly_clim_INIT.copy(deep=True) # create a temporary copy of the df so that we don't alter the original

    # Extract the specific month we're interested in (e.g., January, represented by 1)
    df_dly_target_month = df_dly[df_dly.index.month == monthNum].copy(deep=True)
    df_dly_target_month["Day"] = df_dly_target_month.index.day-1
    df_dly_clim_target_month = df_dly_clim[df_dly_clim.index.get_level_values('Month') == monthNum]

    # Plot the time series for each year separately
    for year, group in df_dly_target_month.groupby(df_dly_target_month.index.year):
        group.plot(x="Day", y=var_colName, ax=ax, c='k', lw=0.5, alpha=0.3)

    # Plot the average time series over the twenty years
    df_dly_clim_target_month[var_colName].plot(color='black', linewidth=2)

    # label x-axis with day number, and only at the beginning of every week (i.e., every 7 days)
    ax.xaxis.set_minor_locator(dates.DayLocator())
    ax.xaxis.set_major_locator(dates.DayLocator(interval=7))
    ax.xaxis.set_major_formatter(dates.DateFormatter('%d'))
    ax.set_xlabel('Day of Month')
    
    # Set the lower limit of the y-axis to zero
    ax.set_ylim(bottom=0)
    ax.set_ylabel(var_pltName)

    ax.set_title("Interannual variability of " + var_pltName + " (" + calendar.month_abbr[monthNum] + ")")


    ax.legend().set_visible(False)

    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)




# Plot one variable for each month of the year, showing each individual year of data as well as their mean
# df_indiv - a dataframe of daily data for all months across all years (index name: "Date"; index format: YYYY-MM-DD)
# df_mean  - a dataframe containing the mean of df_indiv across years (index name: "Date"; index format: YYYY-MM-DD where YYYY is a dummy leap year (e.g., 2000))
def plot_interannualVar_allMon_tseries(df_indiv, df_mean, var_colName, var_pltName, fileName):
    dir_out = "/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/" # TEMPORARY!

    fig, ax = plt.subplots() # create a new figure and axis for plot

    df_indiv_cpy = df_indiv.copy(deep=True)  # create a temporary copy of the df so that we don't alter the original
    df_mean_cpy  = df_mean.copy(deep=True)   # create a temporary copy of the df so that we don't alter the original

    # Extract the day of the year
    df_indiv_cpy["Day of Year"] = df_indiv_cpy.index.dayofyear
    df_mean_cpy["Day of Year"]  = np.arange(1,366+1)

    # Define a colormap for months
    cmap = ListedColormap(plt.get_cmap("tab20").colors[:12])

    # Plot the time series.
    # For temperature, only plot min/max temps.
    if var_colName == "TA (C)" or var_colName == "TW (C)":

        if var_colName == "TA (C)":
            varMax_colName = "TA_max (C)"
            varMin_colName = "TA_min (C)"
        elif var_colName == "TW (C)":
            varMax_colName = "TW_max (C)"
            varMin_colName = "TW_min (C)"

        # Plot the time series for each month
        for month, group_mon in df_indiv_cpy.groupby(df_indiv_cpy.index.month): # for a given calendar month (e.g., January)
            print("Month: " + str(month))
            for year, group_monYr in group_mon.groupby(group_mon.index.year): # for a specific month (e.g., January 2001)
                group_monYr.plot(x="Day of Year", y=varMax_colName, ax=ax, lw=0.6, alpha=0.20, color=cmap(month - 1))
                group_monYr.plot(x="Day of Year", y=varMin_colName, ax=ax, lw=0.6, alpha=0.20, color=cmap(month - 1))

        # Plot the average time series over the twenty years
        for month, group in df_mean_cpy.groupby(df_mean_cpy.index.month):
            group.plot(x="Day of Year", y=varMax_colName, ax=ax, lw=2, color=cmap(month - 1))
            group.plot(x="Day of Year", y=varMin_colName, ax=ax, lw=2, color=cmap(month - 1))

    else: # if any variable other than temperature
            
        # Plot the time series for each month
        for month, group_mon in df_indiv_cpy.groupby(df_indiv_cpy.index.month): # for a given calendar month (e.g., January)
            print("Month: " + str(month))
            for year, group_monYr in group_mon.groupby(group_mon.index.year): # for a specific month (e.g., January 2001)
                group_monYr.plot(x="Day of Year", y=var_colName, ax=ax, lw=0.6, alpha=0.20, color=cmap(month - 1))

        # Plot the average time series over the twenty years
        for month, group in df_mean_cpy.groupby(df_mean_cpy.index.month):
            group.plot(x="Day of Year", y=var_colName, ax=ax, lw=2, color=cmap(month - 1))


    # label x-axis with day number, and only at the beginning of every week (i.e., every 7 days)
    ax.xaxis.set_major_locator(dates.MonthLocator())
    ax.xaxis.set_minor_locator(dates.MonthLocator(bymonthday=15))
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.xaxis.set_minor_formatter(dates.DateFormatter('%b'))
    ax.set_xlim(left=-1, right=367)
    ax.set_xlabel('Day of Year', fontsize=16)

    # format y-axis
    if var_colName == "TW (C)" or var_colName == "TA (C)":
        ax.set_ylim(bottom=5, top=50)
    else:
        ax.set_ylim(bottom=0)
    ax.set_ylabel(var_pltName, fontsize=16)

    ax.tick_params(axis='y', which='both', labelsize=16)
    ax.tick_params(axis='x', which='both', labelsize=13)

    ax.set_title("Interannual variability of " + var_pltName)

    ax.legend().set_visible(False)

    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)




# Plot one variable for each month of the year, one boxplot for each month
# Arguments:
#   df_in - df of all data (presumably monthly)
def plot_interannualVar_allMon_boxplot(df_in, var_colName, var_pltName, fileName):
  dir_out = "/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/" # TEMPORARY!
  
  fig, ax = plt.subplots() # create a new figure and axis for plot
  
  df = df_in.copy(deep=True)  # create a temporary copy of the df so that we don't alter the original

  # Define hex color codes for each month 
  # Colors are based on colorbrewer sequential, single hue palettes with 9 colors (with some adjustments)
  month_colors = ['#9e3138', '#faaf9c', '#aa3400', '#fdae6b',            # NEM partial (JF), FIM (MA)
                  '#00441b', '#006d2c', '#41ab5d', '#74c476', '#a1d99b', # SWM (MJJAS)
                  '#08306b', '#c6dbef', '#800817']                       # SIM (ON), NEM partial (D)

  # Extract month from the index to use as labels
  month_labels = df.index.strftime('%b').unique()

  df['M'] = df.index.month

  # Create box plots for each month
  boxplot_dict = df.boxplot(column=var_colName, by='M', grid=False, ax=ax, 
                            boxprops=dict(color='k'), medianprops=dict(color='k'),
                            whiskerprops=dict(color='k'), flierprops=dict(color='k'),
                            return_type='both', patch_artist=True)

  # Set colors for each box plot
  for row_key, (ax,row) in boxplot_dict.iteritems():
      ax.set_xlabel('')
      for i,box in enumerate(row['boxes']):
          box.set_facecolor(month_colors[i])

  # Set x and y tick labels 
  ax.set_xticks(np.arange(len(month_labels)))   # set correct tick positions
  ax.set_xticklabels(month_labels, fontsize=14) # set x-axis labels using the month labels
  ax.tick_params(axis='y', labelsize=14)        # set y-axis tick label font size

  # Set y-axis limits based on the variable
  if var_colName in ['surv_E', 'surv_L', 'surv_P', ]:
    ax.set_ylim(0, 1)
  if var_colName in ['PRCP (mm)', 'WH (mm)', 'pop_A(E)_frac', 'pop_A(L)_frac', 'pop_A(P)_frac']:
     ax.set_ylim(bottom=0)

  # Set plot title and labels
  ax.set_title('Interannual variability of ' + var_pltName)
  ax.set_xlabel('Month', fontsize=16)
  ax.set_ylabel(var_pltName, fontsize=16)
      
  # if outbreak years are in the dataset, mark those years of data
  if {2017, 2019}.issubset(df.index.year.unique()):
    # Mark datapoints for outbreak years 2017 and 2019
    for i, (month, df_month) in enumerate(df.groupby(df.index.month)):
      outliers_2017 = df_month[df_month.index.year == 2017][var_colName]
      ax.scatter([i + 1] * len(outliers_2017), outliers_2017, marker='*', color='gold', zorder=5)
      outliers_2019 = df_month[df_month.index.year == 2019][var_colName]
      ax.scatter([i + 1] * len(outliers_2019), outliers_2019, marker='*', color='darkgoldenrod', zorder=5)
    # Add legend
    star_2017 = plt.Line2D([0], [0], marker='*', color='gold', linestyle='', label='2017')
    star_2019 = plt.Line2D([0], [0], marker='*', color='darkgoldenrod', linestyle='', label='2019')
    ax.legend(handles=[star_2017, star_2019])

  if not os.path.isdir(dir_out):
    os.makedirs(dir_out)
  plt.savefig(dir_out + fileName + ".png", dpi=300, bbox_inches="tight")
  plt.close(fig)




# Plot two variables for each month of the year, one boxplot for each month
# This is mainly intended for showing boxplots of max/min temperatures on the same plot.
# Arguments:
#   df_in - df of all data (presumably monthly)
def plot_interannualVar_allMon_boxplot_2var(df_in, var1_colName, var1_pltName, var2_colName, var2_pltName, fileName):
  dir_out = "/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/" # TEMPORARY!
  
  fig, ax = plt.subplots() # create a new figure and axis for plot
  
  df = df_in.copy(deep=True)  # create a temporary copy of the df so that we don't alter the original

  # Define hex color codes for each month 
  # Colors are based on colorbrewer sequential, single hue palettes with 9 colors (with some adjustments)
  month_colors = ['#9e3138', '#faaf9c', '#aa3400', '#fdae6b',           # NEM partial (JF), FIM (MA)
                  '#00441b', '#006d2c', '#41ab5d', '#74c476', '#a1d99b', # SWM (MJJAS)
                  '#08306b', '#c6dbef', '#800817']                       # SIM (ON), NEM partial (D)

  # Extract month from the index to use as labels
  month_labels = df.index.strftime('%b').unique()

  df['M'] = df.index.month

  # Prepare data for both variables in a combined format
  data_var1 = [df[df['M'] == month][var1_colName].values for month in range(1, 13)]
  data_var2 = [df[df['M'] == month][var2_colName].values for month in range(1, 13)]

  # Plot for the first variable
  boxplot_dict1 = ax.boxplot(data_var1, positions=range(1, 13), widths=0.4,
                              patch_artist=True, medianprops=dict(color='k'),
                              whiskerprops=dict(color='k'), flierprops=dict(color='k'))

  # Plot for the second variable
  boxplot_dict2 = ax.boxplot(data_var2, positions=range(1, 13), widths=0.4,
                              patch_artist=True, medianprops=dict(color='k'),
                              whiskerprops=dict(color='k'), flierprops=dict(color='k'))

  # Set colors for each box plot based on the month
  for i in range(12):
      boxplot_dict1['boxes'][i].set_facecolor(month_colors[i])
      boxplot_dict2['boxes'][i].set_facecolor(month_colors[i])

  # Set x-axis labels using the month labels
  ax.set_xticks(range(1, 13))
  ax.set_xticklabels(month_labels)

  # Set plot title and labels
  ax.set_title('Interannual variability of ' + var1_pltName + ' and ' + var2_pltName)
  ax.set_xlabel('Month', fontsize=16)
  ax.set_ylabel(var1_pltName + ' and ' + var2_pltName, fontsize=16)
  
  # if outbreak years are in the dataset, mark those years of data
  if {2017, 2019}.issubset(df.index.year.unique()):
    # Mark datapoints for outbreak years 2017 and 2019
    for i, (month, df_month) in enumerate(df.groupby(df.index.month)):
      for var_colName in [var1_colName, var2_colName]:
          outliers_2017 = df_month[df_month.index.year == 2017][var_colName]
          ax.scatter([i + 1] * len(outliers_2017), outliers_2017, marker='*', color='gold', zorder=5)
          outliers_2019 = df_month[df_month.index.year == 2019][var_colName]
          ax.scatter([i + 1] * len(outliers_2019), outliers_2019, marker='*', color='darkgoldenrod', zorder=5)
    # Add legend
    star_2017 = plt.Line2D([0], [0], marker='*', color='gold', linestyle='', label='2017')
    star_2019 = plt.Line2D([0], [0], marker='*', color='darkgoldenrod', linestyle='', label='2019')
    ax.legend(handles=[star_2017, star_2019])

  if not os.path.isdir(dir_out):
    os.makedirs(dir_out)
  plt.savefig(dir_out + fileName + ".png", dpi=300, bbox_inches="tight")
  plt.close(fig)




# Plot timeseries of one variable from the passed df.
# Highlight specific timeseries corresponding to the passed endmember years (i.e., years in which endmember vars are min/max).
# Arguments:
#   df_in - df of all data
#   dict_endmemberYrs    - dict indicating which years are endmembers (and should be highlighted in the plot)
#   dict_endmemberColors - dict indicating the plot colors to use for each endmember
def plot_varVsTime_climEndmembers_oneMon(df_in, dict_endmemberYrs, dict_endmemberColors, var_colName, var_pltName, fileName):
    dir_out = "/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/" # TEMPORARY!
 
    df = df_in.copy(deep=True)  # create a temporary copy of the df so that we don't alter the original

    fig, ax = plt.subplots() # create a new figure and axis for plot

    # first plot all the non-endmember years so that they're on lower layers
    df_nonendmemberYrs = df[~df.index.year.isin(dict_endmemberYrs.keys())]
    if var_colName in df_nonendmemberYrs.columns:
        grouped = df_nonendmemberYrs.groupby(df_nonendmemberYrs.index.year)
        for year, group in grouped:
            day_of_month = group.index.day
            ax.plot(day_of_month, group[var_colName], lw=1, c='black', alpha=0.10)
    else:
        print(f"Warning: Column '{var_colName}' not found in non-endmember years data")

    # next plot the endmember years so that they're on upper layers
    for endmemberYr, endmemberVars in dict_endmemberYrs.items():
        df_endmemberYr = df[df.index.year == endmemberYr]
        if var_colName in df.columns:
            day_of_month = df_endmemberYr.index.day
            color = dict_endmemberColors[endmemberVars]
            endmemberLabel = str(endmemberYr) + ' ('+ endmemberVars + ')'
            ax.plot(day_of_month, df_endmemberYr[var_colName], lw=3, label=endmemberLabel, color=color)
        else:
            print(f"Warning: Column '{var_colName}' not found in endmember years data")

    # Add labels
    ax.set_xlabel('Day of Month')
    ax.set_ylabel(var_colName) # using colName instead of pltName here bc colName has units
    ax.legend()

    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)




# Make scatter plot based on two columns of the passed df.
# The data is colored based on calendar month.
def plot_scatter_colorMon(df_in, varX_colName, varX_pltName, varY_colName, varY_pltName, fileName):
    dir_out = "/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/" # TEMPORARY!

    fig, ax = plt.subplots() # create a new figure and axis for plot

    df = df_in.copy(deep=True) # create a temporary copy of the df so that we don't alter the original

    # Create a discrete color map
    numMon = 12
    colors = plt.cm.twilight_shifted(np.linspace(0, 1, numMon))
    cmap = ListedColormap(colors)
    bounds = np.arange(1, numMon + 2) - 0.5
    norm = BoundaryNorm(bounds, cmap.N)

    # Create scatter plot
    month_colors = [cmap(norm(month)) for month in df.index.month]
    ax.scatter(df[varX_colName], df[varY_colName], c=month_colors, s=10)

    # Add labels
    ax.set_xlabel(varX_pltName)
    ax.set_ylabel(varY_pltName)

    # Create discrete colorbar (only if there's more than one unique calendar month in the data)
    if len(df.index.month.unique()) > 1:
        cbar = plt.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, ticks=np.arange(1, numMon + 1))
        cbar.set_label('Month')
        cbar.set_ticks(np.arange(1, numMon + 1))
        cbar.set_ticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])

    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)




# Plot correlogram for all variables in the passed df.
# Statistically significant (p < 0.05) correlations are denoted with '*'.
# Arguments:
#   df - df with all columns for correlogram
def plot_correlogram_allVars(df_in, fileName):
    dir_out = "/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/" # TEMPORARY!

    fig, ax = plt.subplots() # create a new figure and axis for plot

    df = df_in.copy(deep=True)  # create a temporary copy of the df so that we don't alter the original

    # calculate correlation matrix
    corrMatrix = df.corr()

    # Calculate and store p values
    p_values = pd.DataFrame(np.zeros((df.shape[1], df.shape[1])), columns=df.columns, index=df.columns)
    for col1 in df.columns:
        for col2 in df.columns:
            if col1 == col2:
                p_values.loc[col1, col2] = np.nan  # p-value for correlation with itself is not relevant
            else:
                _, p_value = pearsonr(df[col1], df[col2])
                p_values.loc[col1, col2] = p_value

    # Create a mask for the upper triangle (since correlation matrix is symmetric)
    mask = np.triu(np.ones_like(corrMatrix, dtype=bool))

    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(corrMatrix, mask=mask, annot=False, cmap='coolwarm', linewidths=0.5, vmin=-1, vmax=1)

    # Add annotations for significant correlations on the lower triangle
    for i in range(corrMatrix.shape[0]):
        for j in range(i):
            p_value = p_values.iloc[i, j]
            if p_value < 0.05:  # Common significance level
                plt.text(j + 0.5, i + 0.65, '*', fontsize=10, color='black', ha='center', va='center')

    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)




# Plot stacked bar plot of terciles, showing the percentage of each tercile of varY that corresponds
# to each tercile of varX.
# season - months from which to include data (options: 'allMon', 'NEM', 'FIM', 'SWM', 'SIM')
def plot_bar_tercileXY_splitByMon(df_in, varX_colName, varX_pltName, varY_colName, varY_pltName, season, fileName):
    dir_out = "/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/" # TEMPORARY!

    # set custom colors for stacked bar plot
    tercileColors = ['#edf8b1', '#7fcdbb', '#2c7fb8']

    fig, ax = plt.subplots() # create a new figure and axis for plotx

    df = df_in.copy(deep=True)  # create a temporary copy of the df so that we don't alter the original

    # subset to season of interest
    if season == 'allMon':
        months2incl = list(range(1,12+1))
    elif season == 'NEM':
        months2incl = [12, 1, 2]
    elif season == 'FIM':
        months2incl = [3, 4]
    elif season == 'SWM':
        months2incl = [5, 6, 7, 8, 9]
    elif season == 'SIM':
        months2incl = [10, 11]
    else:
       print(f'Error: unknown season argument: {season}')
       return
    df = df[df.index.month.isin(months2incl)]
    
    # Categorize data into Low, Medium, High, on a monthly basis (and then recombine all monthly results into one series)
    varX_terciles = categorizeSeries_LowMedHigh_byMon(df[varX_colName])    
    varY_terciles = categorizeSeries_LowMedHigh_byMon(df[varY_colName])

    # If categorization failed for either variable, we're not going to make a plot.
    if varX_terciles.isnull().all():
       print(f'Error: categorization failed for variable {varX_colName}. Skipping plot generation.')
       return
    if varY_terciles.isnull().all():
       print(f'Error: categorization failed for variable {varY_colName}. Skipping plot generation.')
       return

    # Add tercile columns to the DataFrame
    varXtercile_colName = varX_colName + '_tercile'
    varYtercile_colName = varY_colName + '_tercile'
    df[varXtercile_colName] = varX_terciles
    df[varYtercile_colName] = varY_terciles

    # Step 2: Create a crosstab to calculate the counts of varY terciles within each varX tercile
    crosstab = pd.crosstab(df[varXtercile_colName], df[varYtercile_colName])

    # Step 3: Normalize the crosstab to get percentages
    crosstab_normalized = crosstab.div(crosstab.sum(axis=1), axis=0) * 100

    # Step 4: Perform chi-square test of independence
    chi2_stat, chi2_p_val, _, _ = chi2_contingency(crosstab)
    print(" ")
    print("Chi-square test")
    print("*******************************")
    print("  Chi-Square Statistic:", chi2_stat)
    print("  p-value:", chi2_p_val)

    # Step 5: Calculate Cramer's V (strength of association; 0 to 1)
    cramers_n = crosstab.sum().sum()       # Total number of observations
    cramers_min_dim = min(crosstab.shape)  # The minimum dimension of the contingency table
    cramers_v_val = np.sqrt(chi2_stat / (cramers_n * (cramers_min_dim - 1)))
    print(" ")
    print("Cramer's V test")
    print("*******************************")
    print("  Cramer's V value: ", cramers_v_val)

    # Step 6: Calculate Kendall's tau-b (direction of association; -1 to 1)
    kendall_tau_b, kendall_p_val = kendalltau(df[varXtercile_colName], df[varYtercile_colName])
    print(" ")
    print("Kendall's tau-b test")
    print("*******************************")
    print("  Kendall's tau-b:", kendall_tau_b)
    print("  p-value:", kendall_p_val)

    # Step 7: Plot the stacked bar chart
    bars = crosstab_normalized.plot(kind='bar', stacked=True, color=tercileColors, ax=ax)
    # Adding labels with percents/counts on each section of the stacked bars
    for i, rect in enumerate(bars.containers):
        for j, bar in enumerate(rect):
            # Calculate the percents/counts for this section
            quantile_percent = round(crosstab_normalized.iloc[j, i])
            quantile_count = crosstab.iloc[j, i]

            if quantile_count > 0: # only do labeling if this category is not empty
                # Get x and y placement of the label based on rectangle location
                x_value = bar.get_x() + bar.get_width() / 2
                y_value = bar.get_y() + bar.get_height() / 2

                # Format and add label
                label = f'{quantile_percent}%\n({quantile_count})'
                ax.text(x_value, y_value, label, ha='center', va='center', color='black', fontsize=10)

    # Add text indicating statistical test results outside the plot
    signif_level = 0.05
    stat_results_text = (r"$\bf{Chi-square\ test}$" + "\n"
                        f"(dependence)\n"
                        f"  Chi-square statistic: {chi2_stat:.2f}\n"
                        f"  p-value: " + (r"$\bf{" if chi2_p_val < signif_level else "") +
                        f"{chi2_p_val:.2e}" + (r"*}$" if chi2_p_val < signif_level else "") + "\n"
                        f"\n"
                        r"$\bf{Cramer's\ V\ test}$" + "\n"
                        f"(assoc. strength; 0 to 1)\n"
                        f"  Cramer's V value: {cramers_v_val:.2f}\n"
                        f"\n"
                        r"$\bf{Kendall's\ tau-b\ test}$" + "\n"
                        f"(assoc. direction; -1 to 1)\n"
                        f"  Kendall's tau-b: {kendall_tau_b:.2f}\n"
                        f"  p-value: " + (r"$\bf{" if kendall_p_val < signif_level else "") +
                        f"{kendall_p_val:.2e}" + (r"*}$" if kendall_p_val < signif_level else ""))
    ax.text(1.025, 0.50, stat_results_text, transform=ax.transAxes, fontsize=8, 
            verticalalignment='center', horizontalalignment='left', bbox=dict(facecolor='white', alpha=0.5))

    # Adding labels and title
    ax.set_xlabel(varX_pltName + ' Tercile')
    ax.set_ylabel('Percentage of ' + varY_pltName + ' Terciles')
    
    # Adjust legend position to avoid overlapping with bars
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))

    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)
   



# Plot stacked bar plot of quantiles, showing the percentage of each binary quantile of varY that corresponds
# to each combined binary categories of varX1 and varX2 (e.g., one combined binary category is low varX1 + varX2).
# Note: we don't do the Kendall tau-b test here since the combined varX categories are unordered.
# Note: this fn is mainly intended for looking at how temp and prcp (X1 and X2) are associated with other model variables.
def plot_bar_binaryX1X2_tercileY_splitByMon(df_in, varX1_colName, varX1_pltName, varX2_colName, varX2_pltName, varY_colName, varY_pltName, season, fileName):
    dir_out = "/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/" # TEMPORARY!

    # set custom colors for stacked bar plot
    # from colorbrewer, sequential, 3-class YlGnBu
    varYColors = ['#edf8b1', '#7fcdbb', '#2c7fb8']

    fig, ax = plt.subplots() # create a new figure and axis for plot

    df = df_in.copy(deep=True)  # create a temporary copy of the df so that we don't alter the original

    # subset to season of interest
    if season == 'allMon':
        months2incl = list(range(1,12+1))
    elif season == 'NEM':
        months2incl = [12, 1, 2]
    elif season == 'FIM':
        months2incl = [3, 4]
    elif season == 'SWM':
        months2incl = [5, 6, 7, 8, 9]
    elif season == 'SIM':
        months2incl = [10, 11]
    else:
       print(f'Error: unknown season argument: {season}')
       return
    df = df[df.index.month.isin(months2incl)]

    # Step 1: Create binary categories for varX1 and varX2
    varX1_binary = categorizeSeries_LowHigh_byMon(df[varX1_colName])
    varX2_binary = categorizeSeries_LowHigh_byMon(df[varX2_colName])

    # If categorization failed for either variable, we're not going to make a plot.
    if varX1_binary.isnull().all():
       print(f'Error: categorization failed for variable {varX1_colName}. Skipping plot generation.')
       return
    if varX2_binary.isnull().all():
       print(f'Error: categorization failed for variable {varX2_colName}. Skipping plot generation.')
       return

    # Combine varX1 and varX2 categories
    label_varX1_lo = 'Low '  + varX1_pltName
    label_varX1_hi = 'High ' + varX1_pltName
    label_varX2_lo = 'Low '  + varX2_pltName
    label_varX2_hi = 'High ' + varX2_pltName
    varX1_binary = varX1_binary.cat.rename_categories({'Low': label_varX1_lo, 'High': label_varX1_hi}) # rename categories to differentiate the variables
    varX2_binary = varX2_binary.cat.rename_categories({'Low': label_varX2_lo, 'High': label_varX2_hi}) # rename categories to differentiate the variables
    varX_combined = varX1_binary.astype(str) + ' + ' + varX2_binary.astype(str) # combined categories
    varX_combined_labels = [label_varX1_lo + ' + ' + label_varX2_lo, label_varX1_lo + ' + ' + label_varX2_hi, 
                            label_varX1_hi + ' + ' + label_varX2_lo, label_varX1_hi + ' + ' + label_varX2_hi]
    varX_combined = pd.Categorical(varX_combined, categories=varX_combined_labels, ordered=True)

    # Step 2: Create tercile categories for varY
    varY_tercile = categorizeSeries_LowMedHigh_byMon(df[varY_colName])

    # If categorization failed for either variable, we're not going to make a plot.
    if varY_tercile.isnull().all():
      print(f'Error: categorization failed for variable {varY_colName}. Skipping plot generation.')
      return

    # Add combined varX binary and varY tercile columns to the DataFrame
    varX_combined_colName = varX1_colName + '_' + varX2_colName + '_combined'
    varY_tercile_colName = varY_colName + '_tercile'
    df[varX_combined_colName] = varX_combined
    df[varY_tercile_colName] = varY_tercile

    # Step 3: Create a crosstab to calculate the counts of varY tercile within each combined varX binary
    crosstab = pd.crosstab(df[varX_combined_colName], df[varY_tercile_colName])

    # Step 4: Normalize the crosstab to get percentages
    crosstab_normalized = crosstab.div(crosstab.sum(axis=1), axis=0) * 100

    # Step 5: Perform chi-square test of independence
    chi2_stat, chi2_p_val, _, _ = chi2_contingency(crosstab)
    print(" ")
    print("Chi-square test")
    print("*******************************")
    print("  Chi-Square Statistic:", chi2_stat)
    print("  p-value:", chi2_p_val)

    # Step 6: Calculate Cramer's V (strength of association; 0 to 1)
    cramers_n = crosstab.sum().sum()       # Total number of observations
    cramers_min_dim = min(crosstab.shape)  # The minimum dimension of the contingency table
    cramers_v_val = np.sqrt(chi2_stat / (cramers_n * (cramers_min_dim - 1)))
    print(" ")
    print("Cramer's V test")
    print("*******************************")
    print("  Cramer's V value: ", cramers_v_val)

    # Step 7: Plot the stacked bar chart
    bars = crosstab_normalized.plot(kind='bar', stacked=True, color=varYColors, ax=ax)
    # Adding labels with percents/counts on each section of the stacked bars
    for i, rect in enumerate(bars.containers):
        for j, bar in enumerate(rect):
            # Calculate the percents/counts for this section
            quantile_percent = round(crosstab_normalized.iloc[j, i])
            quantile_count = crosstab.iloc[j, i]

            if quantile_count > 0: # only do labeling if this category is not empty
                # Get x and y placement of the label based on rectangle location
                x_value = bar.get_x() + bar.get_width() / 2
                y_value = bar.get_y() + bar.get_height() / 2

                # Format and add label
                label = f'{quantile_percent}%\n({quantile_count})'
                ax.text(x_value, y_value, label, ha='center', va='center', color='black', fontsize=10)

    # Add text indicating statistical test results outside the plot
    signif_level = 0.05
    stat_results_text = (r"$\bf{Chi-square\ test}$" + "\n"
                        f"(dependence)\n"
                        f"  Chi-square statistic: {chi2_stat:.2f}\n"
                        f"  p-value: " + (r"$\bf{" if chi2_p_val < signif_level else "") +
                        f"{chi2_p_val:.2e}" + (r"*}$" if chi2_p_val < signif_level else "") + "\n"
                        f"\n"
                        r"$\bf{Cramer's\ V\ test}$" + "\n"
                        f"(assoc. strength; 0 to 1)\n"
                        f"  Cramer's V value: {cramers_v_val:.2f}")
    ax.text(1.025, 0.50, stat_results_text, transform=ax.transAxes, fontsize=8, 
            verticalalignment='center', horizontalalignment='left', bbox=dict(facecolor='white', alpha=0.5))

    # Adding labels and title
    ax.set_xlabel(varX1_pltName + ' + ' + varX2_pltName + ' Categories')
    ax.set_ylabel('Percentage of ' + varY_pltName + ' Categories')
    
    # Adjust legend position to avoid overlapping with bars
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))

    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)


