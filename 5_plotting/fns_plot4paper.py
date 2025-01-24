#*******************************************************
# helperFns_WHATCHEM_plotting.py
#*******************************************************
# Revision history
#    redoing plotting code in Python
#
# Plotting functions:
#
# Non-plotting functions:


# import required packages
# ************************************************
import matplotlib.pyplot as plt
from matplotlib import dates
from matplotlib.colors import ListedColormap, BoundaryNorm, to_rgba, to_hex
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import FuncFormatter, FixedLocator, NullFormatter, FormatStrFormatter
import matplotlib.colors as mcolors
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
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
sys.path.append('/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/0_constants') # directory of constants file
sys.path.append('/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/4_data_analysis') # directory of helper fns file
from CONST_FILE_DIR import DIR_WHATCHEM_OUTPUT_DATA_PROCESSED
from fns_helper import (get_df_vectorPop_dly_merged, convert_df_vectorPop_dly2mon, 
                        categorizeSeries_LowHigh, categorizeSeries_LowHigh_byMon, categorizeSeries_LowMedHigh_byMon, 
                        strConvert_paren2hyphen, strConvert_parenUnitsRemove,
                        create_df_withLags)



# Plotting functions
# ************************************************



### DENGUE PLOTS #################################

def plot_tseries_denInc(df_in, yearStart_data, yearEnd_data, outbreakYears, fileName):
    dir_out = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/'  # TEMPORARY!

    # Create figure
    fig, ax1 = plt.subplots(figsize=(16, 12))

    # Filter for the specified locations.
    locs = ['SRILANKA', 'Gampaha', 'Jaffna', 'NuwaraEliya']
    df_SL  = df_in[[locs[0]]].copy(deep=True)
    df_Neg = df_in[[locs[1]]].copy(deep=True)
    df_Jaf = df_in[[locs[2]]].copy(deep=True)
    df_NuE = df_in[[locs[3]]].copy(deep=True)

    for df in [df_SL, df_Neg, df_Jaf, df_NuE]:
        df['Month'] = df.index.month
        df['Year'] = df.index.year

    # Compute averages across years (excluding outbreak years)
    df_SL_noOutbreakYrs  = df_SL[~df_SL['Year'].isin(outbreakYears)]
    df_Neg_noOutbreakYrs = df_Neg[~df_Neg['Year'].isin(outbreakYears)]
    df_Jaf_noOutbreakYrs = df_Jaf[~df_Jaf['Year'].isin(outbreakYears)]
    df_NuE_noOutbreakYrs = df_NuE[~df_NuE['Year'].isin(outbreakYears)]
    df_SL_noOutbreakYrs_avg   = df_SL_noOutbreakYrs.groupby('Month')[locs[0]].mean()
    df_Neg_noOutbreakYrs_avg = df_Neg_noOutbreakYrs.groupby('Month')[locs[1]].mean()
    df_Jaf_noOutbreakYrs_avg = df_Jaf_noOutbreakYrs.groupby('Month')[locs[2]].mean()
    df_NuE_noOutbreakYrs_avg = df_NuE_noOutbreakYrs.groupby('Month')[locs[3]].mean()

    # Compute standard deviations for all of Sri Lanka (excluding outbreak years)
    df_SL_noOutbreakYrs_std = df_SL_noOutbreakYrs.groupby('Month')[locs[0]].std()

    # Define plot line parameters for main plot
    color_main       = '#222222'
    color_cities     = '#0000dd'
    lw_main          = 10

    # Define plot line parameters for inset plot
    color_inset_2017       = '#E60000'
    color_inset_2019       = '#660000'
    lw_inset_major         = 5
    lw_inset_2017_2019     = 2

    # Set y-axis bounds based on loc
    ymin = 0
    ymax_SL, ystep_SL = 51, 10
    ymax_SL_inset, ystep_SL_inset = 160, 50

    ### Make main plot
    # Plot average across (non-outbreak) years as well as each (non-outbreak) year.
    ax1.plot(df_SL_noOutbreakYrs_avg.index, df_SL_noOutbreakYrs_avg.values, label=locs[0], color=color_main, lw=lw_main)

    # Plot the shaded area for ±2 standard deviation
    ax1.fill_between(
        df_SL_noOutbreakYrs_avg.index,
        df_SL_noOutbreakYrs_avg.values - 2*df_SL_noOutbreakYrs_std.values,  # Lower bound
        df_SL_noOutbreakYrs_avg.values + 2*df_SL_noOutbreakYrs_std.values,  # Upper bound
        color='gray', alpha=0.2, lw=0
    )

    ls_locs = ['--',':','-.']
    for i, df in enumerate([df_Neg_noOutbreakYrs_avg, df_Jaf_noOutbreakYrs_avg, df_NuE_noOutbreakYrs_avg]):
        locNum = i+1
        loc = locs[locNum]
        ax1.plot(df.index, df.values, label=loc, color=color_cities, lw=6, alpha=0.9, ls=ls_locs[i])


    # Add inset plots (with all years, including outbreak years).
    ax1_inset = inset_axes(ax1, width='100%', height='100%', bbox_to_anchor=(0.62, 0.70, 0.27, 0.26), bbox_transform=ax1.transAxes)
    ax1_inset.plot(df_SL_noOutbreakYrs_avg.index, df_SL_noOutbreakYrs_avg.values, label=locs[0], 
                   color=color_main, lw=lw_inset_major)
    for year in outbreakYears:
        if year == 2017:
            color_SL     = color_inset_2017
        elif year == 2019:
            color_SL     = color_inset_2019
        lw           = lw_inset_2017_2019
        zorder       = 4
        ax1_inset.plot(df_SL[df_SL['Year'] == year]['Month'], df_SL[df_SL['Year'] == year][locs[0]],
                    color=color_SL, lw=lw, zorder=zorder)
    ax1_inset.fill_between(
        df_SL_noOutbreakYrs_avg.index,
        df_SL_noOutbreakYrs_avg.values - 2*df_SL_noOutbreakYrs_std.values,  # Lower bound
        df_SL_noOutbreakYrs_avg.values + 2*df_SL_noOutbreakYrs_std.values,  # Upper bound
        color='gray', alpha=0.3, lw=0
    )

    # Format all y-axes
    ax1.set_ylim(ymin, ymax_SL)
    ax1.set_yticks(np.arange(ymin, ymax_SL+0.01, ystep_SL))
    ax1_inset.set_ylim(ymin, ymax_SL_inset)
    ax1_inset.set_yticks(np.arange(ymin, ymax_SL_inset+0.01, ystep_SL_inset))
    for ax in [ax1]:
        ax.set_ylabel('Dengue incidence (per 100,000)', fontsize=20)
    ax1_inset.set_ylabel('Den. inc. (per 100k)', fontsize=18)
    for ax in [ax1, ax1_inset]:
        ax.tick_params(axis='y', labelsize=18)  # set y-axis tick label size
        ax.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, left=True, right=False) # adjust minor ticks for y-axis
        ax.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, left=True, right=False) # adjust major ticks for y-axis

    # Format all x-axes (months)
    ax1.set_xticks(np.arange(1, 12+1))
    ax1.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'], fontsize=20)
    ax1_inset.set_xticks(np.arange(1, 12+1), minor=True)
    ax1_inset.set_xlim(left=0.75, right=12.25)   # set x-axis limits so that it begins close to January (1) and ends close to December (12)
    ax1_inset.set_xticks([2, 5, 8, 11])
    ax1_inset.set_xticklabels(['Feb', 'May', 'Aug', 'Nov'], fontsize=18)

    # Save the figure
    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)




##################################################
### TSERIES PLOTS ################################

# Plot mean timeseries of variable for each climate category (based on combinations of low/high temp/precip).
def plot_tseries_climCategories(df_in, month, var_colName, var_pltName, varTemp_colName, varTemp_pltName, varPrcp_colName, varPrcp_pltName, fileName):
    dir_out = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/' # TEMPORARY!

    fig, ax = plt.subplots() # create a new figure and axis for plot

    df_dly = df_in.copy(deep=True)  # create a temporary copy of the df so that we don't alter the original

    df_monthly = convert_df_vectorPop_dly2mon(df_dly) # convert to monthly

    df_dly = df_dly[df_dly.index.month == month] # subset to month of interest
    df_monthly = df_monthly[df_monthly.index.month == month] # subset to month of interest

    # Create binary categories for temp and prcp variables
    varTemp_binary = categorizeSeries_LowHigh(df_monthly[varTemp_colName])
    varPrcp_binary = categorizeSeries_LowHigh(df_monthly[varPrcp_colName])

    # If categorization failed for either variable, we're not going to make a plot.
    if varTemp_binary.isnull().all():
       print(f'Error: categorization failed for variable {varTemp_colName}. Skipping plot generation.')
       return
    if varPrcp_binary.isnull().all():
       print(f'Error: categorization failed for variable {varPrcp_colName}. Skipping plot generation.')
       return

    # Combine varTemp and varPrcp categories
    label_varTemp_lo = 'Low '  + varTemp_pltName
    label_varTemp_hi = 'High ' + varTemp_pltName
    label_varPrcp_lo = 'Low '  + varPrcp_pltName
    label_varPrcp_hi = 'High ' + varPrcp_pltName
    varTemp_binary = varTemp_binary.cat.rename_categories({'Low': label_varTemp_lo, 'High': label_varTemp_hi}) # rename categories to differentiate the variables
    varPrcp_binary = varPrcp_binary.cat.rename_categories({'Low': label_varPrcp_lo, 'High': label_varPrcp_hi}) # rename categories to differentiate the variables
    varClimCombined = varTemp_binary.astype(str) + ' + ' + varPrcp_binary.astype(str) # combined categories
    varClimCombined_labels = [label_varTemp_lo + ' + ' + label_varPrcp_lo, label_varTemp_lo + ' + ' + label_varPrcp_hi, 
                              label_varTemp_hi + ' + ' + label_varPrcp_lo, label_varTemp_hi + ' + ' + label_varPrcp_hi]
    varClimCombined = pd.Categorical(varClimCombined, categories=varClimCombined_labels, ordered=True)

    # Add combined varClim binary to the DataFrame
    varClimCombined_colName = varTemp_colName + '_' + varPrcp_colName + '_combined'
    df_monthly[varClimCombined_colName] = varClimCombined

    # Define color mapping for each climate category
    color_mapping = {
        varClimCombined_labels[0]: '#529ef7',   # Low Temp + Low Prcp
        varClimCombined_labels[1]: '#7737de',   # Low Temp + High Prcp
        varClimCombined_labels[2]: '#d95f02',   # High Temp + Low Prcp
        varClimCombined_labels[3]: '#1b9e77',   # High Temp + High Prcp
    }

    ls_mapping = {
        varClimCombined_labels[0]: ':',   # Low Temp + Low Prcp
        varClimCombined_labels[1]: '-.',  # Low Temp + High Prcp
        varClimCombined_labels[2]: '--',  # High Temp + Low Prcp
        varClimCombined_labels[3]: '-',   # High Temp + High Prcp
    }

    # Subset daily data to var of interest
    df_dly_var = df_dly.copy(deep=True)
    df_dly_var = df_dly_var[[var_colName]]

    # Map climate categories to daily data
    df_dly_var[varClimCombined_colName] = df_dly_var.index.to_period('M').map(
                                           lambda x: df_monthly.loc[x.to_timestamp(), varClimCombined_colName]
                                           if x.to_timestamp() in df_monthly.index else None
                                           )

    # Compute daily averages and stdevs for each climate category
    df_dly_var['Day'] = df_dly_var.index.day
    df_dly_var_climCatStats = df_dly_var.groupby(['Day', varClimCombined_colName]).agg(['mean', 'std'])
    df_dly_var_climCatAvgs = df_dly_var_climCatStats[var_colName]['mean'].unstack(varClimCombined_colName)
    df_dly_var_climCatStds = df_dly_var_climCatStats[var_colName]['std'].unstack(varClimCombined_colName)

    # Plot the timeseries
    for climCat in varClimCombined_labels:
        if climCat in df_dly_var_climCatAvgs:

            # Plot the mean timeseries
            ax.plot(df_dly_var_climCatAvgs.index, df_dly_var_climCatAvgs[climCat], 
                    label=climCat, color=color_mapping.get(climCat, 'black'), ls=ls_mapping.get(climCat, '-'), lw=3)
            
            # Add shading for standard deviation
            ax.fill_between(df_dly_var_climCatStds.index, 
                            df_dly_var_climCatAvgs[climCat] - df_dly_var_climCatStds[climCat], 
                            df_dly_var_climCatAvgs[climCat] + df_dly_var_climCatStds[climCat], 
                            color=color_mapping.get(climCat, 'black'), 
                            alpha=0.1)

    # Adding labels and title
    ax.set_xlabel('Day of Month')
    ax.set_ylabel(var_colName)
    
    # Adjust legend position to avoid overlapping with bars
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))

    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


# Function to plot timeseries of all variables for each climate category.
def plot_tseries_climCategories_allVars(df_in, loc, month, varTemp_colName, varTemp_pltName, varPrcp_colName, varPrcp_pltName, fileName):
    dir_out = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/'  # TEMPORARY!

    df_dly = df_in.copy(deep=True)  # create a temporary copy of the df so that we don't alter the original
    df_monthly = convert_df_vectorPop_dly2mon(df_dly) # convert to monthly

    df_dly = df_dly[df_dly.index.month == month] # subset to month of interest
    df_monthly = df_monthly[df_monthly.index.month == month] # subset to month of interest

    # Create binary categories for temp and prcp variables
    varTemp_binary = categorizeSeries_LowHigh(df_monthly[varTemp_colName])
    varPrcp_binary = categorizeSeries_LowHigh(df_monthly[varPrcp_colName])

    # If categorization failed for either variable, we're not going to make a plot.
    if varTemp_binary.isnull().all():
       print(f'Error: categorization failed for variable {varTemp_colName}. Skipping plot generation.')
       return
    if varPrcp_binary.isnull().all():
       print(f'Error: categorization failed for variable {varPrcp_colName}. Skipping plot generation.')
       return

    # Combine varTemp and varPrcp categories
    label_varTemp_lo = 'Low '  + varTemp_pltName
    label_varTemp_hi = 'High ' + varTemp_pltName
    label_varPrcp_lo = 'Low '  + varPrcp_pltName
    label_varPrcp_hi = 'High ' + varPrcp_pltName
    varTemp_binary = varTemp_binary.cat.rename_categories({'Low': label_varTemp_lo, 'High': label_varTemp_hi}) # rename categories to differentiate the variables
    varPrcp_binary = varPrcp_binary.cat.rename_categories({'Low': label_varPrcp_lo, 'High': label_varPrcp_hi}) # rename categories to differentiate the variables
    varClimCombined = varTemp_binary.astype(str) + ' + ' + varPrcp_binary.astype(str) # combined categories
    varClimCombined_labels = [label_varTemp_lo + ' + ' + label_varPrcp_lo, label_varTemp_lo + ' + ' + label_varPrcp_hi, 
                              label_varTemp_hi + ' + ' + label_varPrcp_lo, label_varTemp_hi + ' + ' + label_varPrcp_hi]
    varClimCombined = pd.Categorical(varClimCombined, categories=varClimCombined_labels, ordered=True)

    # Add combined varClim binary to the DataFrame
    varClimCombined_colName = varTemp_colName + '_' + varPrcp_colName + '_combined'
    df_monthly[varClimCombined_colName] = varClimCombined
    
    # Define colors, linewidths, and plotting order for each climate category
    # NOTE: The ordering of these matters for the legend!
    climCat_params = {
        varClimCombined_labels[1]: {'color': '#7737de', 'ls': '-.', 'zorder': 8},  # Low Temp + High Prcp
        varClimCombined_labels[0]: {'color': '#529ef7', 'ls': ':',  'zorder': 6},  # Low Temp + Low Prcp
        varClimCombined_labels[3]: {'color': '#1b9e77', 'ls': '-',  'zorder': 9},  # High Temp + High Prcp
        varClimCombined_labels[2]: {'color': '#d95f02', 'ls': '--', 'zorder': 7}   # High Temp + Low Prcp
    }

    # Define colors for variables
    color_met  = '#48c335' # color for meteorology variables
    color_cont = '#2677ff' # color for container water dynamics variables
    color_vectBio_L = '#b9216e'  # color for vector biology (larvae) variables
    color_vectPop_L  = '#1b9e77' # color for vector pop (larvae) variables

    # Define variable groups
    #  some values are defined based on location
    #  y_params are in the format [ymin, ymax, ymin_tick, ymax_tick, ystep_tick]
    if loc in ['Jaffna', 'Negombo']:
        y_params_T = [17, 43, 20, 40, 5]
    elif loc == 'NuwaraEliya':
        y_params_T = [7, 33, 10, 30, 5]

    met_vars = {
        'PRCP': {'var_pltName': 'Precipitation (mm)', 'var_list': ['PRCP (mm)'], 'ax_color': color_met, 'y_params': [0, 45, 0, 45, 10], 'y_str_fmt': '%.0f'},
        'TA':   {'var_pltName': 'Air temperature (°C)', 'var_list': ['TA (C)', 'TA_min (C)', 'TA_max (C)'], 'ax_color': color_met, 'y_params': y_params_T, 'y_str_fmt': '%.0f'}
    }
    cont_vars = {
        'WH':   {'var_pltName': 'Water height (mm)', 'var_list': ['WH (mm)'], 'ax_color': color_cont, 'y_params': [0, 120, 0, 120, 30], 'y_str_fmt': '%.0f'},
        'TW':   {'var_pltName': 'Water temperature (°C)', 'var_list': ['TW (C)', 'TW_min (C)', 'TW_max (C)'], 'ax_color': color_cont, 'y_params': y_params_T, 'y_str_fmt': '%.0f'}
    }
    dev_vars = {
        'dev_E':   {'var_pltName': 'Dev. rate (day$^{-1}$)', 'var_list': ['dev_E'], 'ax_color': color_vectBio_L, 'y_params': [0, 0.65, 0, 0.65, 0.2], 'y_str_fmt': '%.1f'},
        'dev_L':   {'var_pltName': 'Dev. rate (day$^{-1}$)', 'var_list': ['dev_L'], 'ax_color': color_vectBio_L, 'y_params': [0, 0.65, 0, 0.65, 0.2], 'y_str_fmt': '%.1f'},
        'dev_P':   {'var_pltName': 'Dev. rate (day$^{-1}$)', 'var_list': ['dev_P'], 'ax_color': color_vectBio_L, 'y_params': [0, 0.65, 0, 0.65, 0.2], 'y_str_fmt': '%.1f'}
    }
    surv_vars = {
        'surv_E':  {'var_pltName': 'Survival rate (day$^{-1}$)', 'var_list': ['surv_E'], 'ax_color': color_vectBio_L, 'y_params': [0, 1.10, 0, 1.10, 0.2], 'y_str_fmt': '%.1f'},
        'surv_L':  {'var_pltName': 'Survival rate (day$^{-1}$)', 'var_list': ['surv_L'], 'ax_color': color_vectBio_L, 'y_params': [0, 1.10, 0, 1.10, 0.2], 'y_str_fmt': '%.1f'},
        'surv_P':  {'var_pltName': 'Survival rate (day$^{-1}$)', 'var_list': ['surv_P'], 'ax_color': color_vectBio_L, 'y_params': [0, 1.10, 0, 1.10, 0.2], 'y_str_fmt': '%.1f'}
    }
    pop_vars = {
        'pop_A(E)_frac': {'var_pltName': 'Adult population (norm.)', 'var_list': ['pop_A(E)_frac'], 'ax_color': color_vectPop_L, 'y_params': [0, 0.90, 0, 0.90, 0.2], 'y_str_fmt': '%.1f'},
        'pop_A(L)_frac': {'var_pltName': 'Adult population (norm.)', 'var_list': ['pop_A(L)_frac'], 'ax_color': color_vectPop_L, 'y_params': [0, 0.90, 0, 0.90, 0.2], 'y_str_fmt': '%.1f'},
        'pop_A(P)_frac': {'var_pltName': 'Adult population (norm.)', 'var_list': ['pop_A(P)_frac'], 'ax_color': color_vectPop_L, 'y_params': [0, 0.90, 0, 0.90, 0.2], 'y_str_fmt': '%.1f'}
    }


    # Create the figure and gridspec layout
    fig = plt.figure(figsize=(30, 36))
    gs_outer = GridSpec(3, 1, figure=fig, height_ratios=[2, 2, 1], hspace=0.32)  # hspace adds space for the section titles
    gs_top    = GridSpecFromSubplotSpec(2, 2, subplot_spec=gs_outer[0, 0], hspace=0.0, wspace=0.03)
    gs_center = GridSpecFromSubplotSpec(2, 3, subplot_spec=gs_outer[1, 0], hspace=0.0, wspace=0.00)
    gs_bot    = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs_outer[2, 0], hspace=0.0, wspace=0.00)

    # Top left plots
    ax_prcp = fig.add_subplot(gs_top[0, 0])
    ax_ta = fig.add_subplot(gs_top[1, 0])
    plot_climCategory_timeseries(ax_prcp, df_dly, df_monthly, met_vars['PRCP'], varClimCombined_colName, varClimCombined_labels, climCat_params)
    plot_climCategory_timeseries(ax_ta, df_dly, df_monthly, met_vars['TA'], varClimCombined_colName, varClimCombined_labels, climCat_params)

    # Top right plots
    ax_wh = fig.add_subplot(gs_top[0, 1])
    ax_tw = fig.add_subplot(gs_top[1, 1])
    plot_climCategory_timeseries(ax_wh, df_dly, df_monthly, cont_vars['WH'], varClimCombined_colName, varClimCombined_labels, climCat_params)
    plot_climCategory_timeseries(ax_tw, df_dly, df_monthly, cont_vars['TW'], varClimCombined_colName, varClimCombined_labels, climCat_params)

    # Center plots
    ax_dev_E   = fig.add_subplot(gs_center[0, 0])
    ax_dev_L   = fig.add_subplot(gs_center[0, 1])
    ax_dev_P   = fig.add_subplot(gs_center[0, 2])
    plot_climCategory_timeseries(ax_dev_E, df_dly, df_monthly, dev_vars['dev_E'], varClimCombined_colName, varClimCombined_labels, climCat_params)
    plot_climCategory_timeseries(ax_dev_L, df_dly, df_monthly, dev_vars['dev_L'], varClimCombined_colName, varClimCombined_labels, climCat_params)
    plot_climCategory_timeseries(ax_dev_P, df_dly, df_monthly, dev_vars['dev_P'], varClimCombined_colName, varClimCombined_labels, climCat_params)

    ax_surv_E  = fig.add_subplot(gs_center[1, 0])
    ax_surv_L  = fig.add_subplot(gs_center[1, 1])
    ax_surv_P  = fig.add_subplot(gs_center[1, 2])
    plot_climCategory_timeseries(ax_surv_E, df_dly, df_monthly, surv_vars['surv_E'], varClimCombined_colName, varClimCombined_labels, climCat_params)
    plot_climCategory_timeseries(ax_surv_L, df_dly, df_monthly, surv_vars['surv_L'], varClimCombined_colName, varClimCombined_labels, climCat_params)
    plot_climCategory_timeseries(ax_surv_P, df_dly, df_monthly, surv_vars['surv_P'], varClimCombined_colName, varClimCombined_labels, climCat_params)
    
    # Bottom plots
    ax_pop_A_E = fig.add_subplot(gs_bot[0, 0])
    ax_pop_A_L = fig.add_subplot(gs_bot[0, 1])
    ax_pop_A_P = fig.add_subplot(gs_bot[0, 2])
    plot_climCategory_timeseries(ax_pop_A_E, df_dly, df_monthly, pop_vars['pop_A(E)_frac'], varClimCombined_colName, varClimCombined_labels, climCat_params)
    plot_climCategory_timeseries(ax_pop_A_L, df_dly, df_monthly, pop_vars['pop_A(L)_frac'], varClimCombined_colName, varClimCombined_labels, climCat_params)
    plot_climCategory_timeseries(ax_pop_A_P, df_dly, df_monthly, pop_vars['pop_A(P)_frac'], varClimCombined_colName, varClimCombined_labels, climCat_params)

    # For interior plot borders, overwrite color with black and reduce thickness
    for ax in [ax_prcp, ax_wh, ax_dev_E, ax_dev_L, ax_dev_P]: # modify lower border
        ax.spines['bottom'].set_color('#666666')
        ax.spines['bottom'].set_linewidth(1)
    for ax in [ax_ta, ax_tw, ax_surv_E, ax_surv_L, ax_surv_P]: # modify upper border
        ax.spines['top'].set_color('#666666')
        ax.spines['top'].set_linewidth(1)
    for ax in [ax_dev_E, ax_surv_E, ax_pop_A_E, ax_dev_L, ax_surv_L, ax_pop_A_L]: # modify right border
        ax.spines['right'].set_color('#666666')
        ax.spines['right'].set_linewidth(1)
    for ax in [ax_dev_L, ax_surv_L, ax_pop_A_L, ax_dev_P, ax_surv_P, ax_pop_A_P]: # modify left border
        ax.spines['left'].set_color('#666666')
        ax.spines['left'].set_linewidth(1)
        

    # Remove x-axis labels for all but the bottommost plots in each section
    for ax in [ax_prcp, ax_wh, ax_dev_E, ax_dev_L, ax_dev_P]:
        ax.set_xlabel('')

    # Move y-axis labels to right side for container variable plots
    for ax in [ax_wh, ax_tw]:
        ax.tick_params(axis='y', labelleft=False, labelright=True) # set tick labels only on right side
        ax.yaxis.set_label_position('right') # set axis label on right side
        ax.yaxis.label.set_rotation(270)
        ax.yaxis.labelpad = 30 # move a smidge to the right to avoid overlap with tick labels

    # Remove y-axis labels for larvae and pupae plots
    for ax in [ax_dev_L, ax_dev_P, ax_surv_L, ax_surv_P, ax_pop_A_L, ax_pop_A_P]:
        ax.set_yticklabels([])
        ax.set_ylabel('')

    # Make labels for meterology and container sections
    for ax in [ax_prcp, ax_wh]:

        if ax == ax_prcp:
            labelStr = 'METEOROLOGY'
            rect_color = color_met
        elif ax == ax_wh:
            labelStr = 'CONTAINER HABITAT'
            rect_color = color_cont

        # make rectangle
        bbox = ax.get_position()                   # returns bounding box of axis (x0, y0, width, height) in figure coordinates
        rect_height = bbox.height*0.15             # height of rectangle (x% of the axis height)
        rect_pos_y  = bbox.y1 + bbox.height*0.03   # set y-position above the top of the axes
        rect_pos_x_offset = bbox.width*0.005
        rect_pos_x  = bbox.x0 + rect_pos_x_offset
        rect_width  = bbox.width - 2*rect_pos_x_offset
        rect = patches.Rectangle(xy=(rect_pos_x, rect_pos_y), width=rect_width, height=rect_height, 
                                 linewidth=0, edgecolor='none', facecolor=rect_color, alpha=0.8, zorder=5, 
                                 transform=ax.figure.transFigure)  # Use figure coordinates
        ax.figure.patches.append(rect)

        # make label
        stage_label_x = bbox.x0 + bbox.width*0.5          # Scale x to normalized coords
        stage_label_y = rect_pos_y + (rect_height / 2)    # Center the text in the rectangle (in normalized y)
        ax.figure.text(x=stage_label_x, y=stage_label_y, s=labelStr, 
                        color='black', ha='center', va='center', fontsize=24, zorder=6, fontweight='bold', 
                        transform=ax.figure.transFigure)  # Use figure coordinates for text

    # Make labels for vector biology and vector population sections
    for ax_left, ax_right in [(ax_dev_E, ax_dev_P), (ax_pop_A_E, ax_pop_A_P)]:
        if ax_left == ax_dev_E:
            labelStr = 'VECTOR BIOLOGY'
            rect_color = color_vectBio_L
        elif ax_left == ax_pop_A_E:
            labelStr = 'VECTOR POPULATION'
            rect_color = color_vectPop_L

        # make rectangle
        bbox_left  = ax_left.get_position()  # returns bounding box of axis (x0, y0, width, height) in figure coordinates
        bbox_right = ax_right.get_position()  # returns bounding box of axis (x0, y0, width, height) in figure coordinates
        rect_height = bbox_left.height*0.15             # height of rectangle (x% of the axis height)
        rect_pos_y  = bbox_left.y1 + bbox_left.height*0.10   # set y-position above the top of the axes
        rect_pos_x_offset = bbox_left.width*0.005
        rect_pos_x  = bbox_left.x0 + rect_pos_x_offset
        bboxes_width = (bbox_right.x0 - bbox_left.x0) + bbox_right.width
        rect_width  = bboxes_width - 2*rect_pos_x_offset
        rect = patches.Rectangle(xy=(rect_pos_x, rect_pos_y), width=rect_width, height=rect_height, 
                                 linewidth=0, edgecolor='none', facecolor=rect_color, alpha=0.8, zorder=5, 
                                 transform=ax.figure.transFigure)  # Use figure coordinates
        ax.figure.patches.append(rect)

        # make label
        stage_label_x = bbox_left.x0 + bboxes_width*0.5   # Scale x to normalized coords
        stage_label_y = rect_pos_y + (rect_height / 2)    # Center the text in the rectangle (in normalized y)
        ax.figure.text(x=stage_label_x, y=stage_label_y, s=labelStr, 
                        color='black', ha='center', va='center', fontsize=24, zorder=6, fontweight='bold', 
                        transform=ax.figure.transFigure)  # Use figure coordinates for text

    # Make labels for eggs, larvae, and pupae columns
    for ax in [ax_dev_E, ax_dev_L, ax_dev_P, ax_pop_A_E, ax_pop_A_L, ax_pop_A_P]:

        bbox = ax.get_position() # returns bounding box of axis (x0, y0, width, height) in figure coordinates
        stage_label_x = bbox.x0 + bbox.width*0.5   # Scale x to normalized coords
        stage_label_y = bbox.y1 + bbox.height*0.04

        if ax in [ax_dev_E, ax_pop_A_E]:
            labelStr = 'EGGS'
        elif ax in [ax_dev_L, ax_pop_A_L]:
            labelStr = 'LARVAE'
        elif ax in [ax_dev_P, ax_pop_A_P]:
            labelStr = 'PUPAE'
        
        ax.figure.text(x=stage_label_x, y=stage_label_y, s=labelStr, 
                        color='black', ha='center', va='center', fontsize=22, zorder=6, fontweight='bold', 
                        transform=ax.figure.transFigure)  # Use figure coordinates for text

    # Add legend within PRCP plot
    bbox_prcp = ax_prcp.get_position() # we'll use this to position the legend
    legend_x_pos = bbox_prcp.x0 + 0.7*(bbox_prcp.x1-bbox_prcp.x0) + 0.010
    legend_y_pos = bbox_prcp.y0 + 0.85*(bbox_prcp.y1-bbox_prcp.y0)
    legend_labels = ['','','','']
    legend_handles = [plt.Line2D([0], [0], color=params['color'], lw=7, ls=params['ls'], label=label) 
                    for label, params in climCat_params.items()]
    legend = fig.legend(labels=legend_labels, handles=legend_handles, loc='upper center', bbox_to_anchor=(legend_x_pos, legend_y_pos), 
                        ncol=2, handlelength=12, labelspacing=1.5, frameon=False) # labelspacing affects vertical spacing

    bbox_legend = legend.get_window_extent().transformed(fig.transFigure.inverted())
    legend_width  = bbox_legend.x1 - bbox_legend.x0
    legend_height = bbox_legend.y1 - bbox_legend.y0

    # Add rectangle outlining legend
    legend_rect = patches.Rectangle((bbox_legend.x0-0.004, bbox_legend.y0-0.0005), legend_width+0.005, legend_height+0.002, transform=fig.transFigure,
                                    facecolor='none', edgecolor='black', linewidth=1, zorder=5)
    fig.add_artist(legend_rect)

    # Add lines within the legend rectangle
    legend_line_vert  = plt.Line2D([legend_rect.get_x() + 0.5*legend_rect.get_width() - 0.002, legend_rect.get_x() + 0.5*legend_rect.get_width() - 0.002], 
                                [legend_rect.get_y(), legend_rect.get_y() + legend_rect.get_height()], 
                                color='#666666', lw=0.5, zorder=4)
    legend_line_horiz = plt.Line2D([legend_rect.get_x(), legend_rect.get_x() + legend_rect.get_width()], 
                                [legend_rect.get_y() + 0.5*legend_rect.get_height(), legend_rect.get_y() + 0.5*legend_rect.get_height()], 
                                color='#666666', lw=0.5, zorder=4)
    fig.add_artist(legend_line_vert)
    fig.add_artist(legend_line_horiz)
    
    # Add labels to the legend
    legend_main_title = 'CLIMATE CATEGORIES'
    legend_col_title  = 'Mean air temp.'
    legend_row_title  = 'Prcp.'
    legend_rowCol_subtitles = ['Low', 'High']

    fig.text(legend_rect.get_x() + 0.5*legend_rect.get_width(), bbox_legend.y1 + 0.008, legend_main_title, ha='center', va='center', fontweight='bold', fontsize=21) # main title
    fig.text(bbox_legend.x0 + 0.50*legend_width, bbox_legend.y0 - 0.017, legend_col_title, ha='center', va='center', fontsize=22)  # col title
    fig.text(bbox_legend.x0 - 0.041, bbox_legend.y0 + 0.50*legend_height, legend_row_title, ha='right', va='center', fontsize=22)  # row title
    for i, subtitle in enumerate(legend_rowCol_subtitles):
        fig.text(bbox_legend.x0 + (0.22 + 0.52*i)*legend_width, bbox_legend.y0 - 0.007, subtitle, ha='center', va='center', fontsize=21)  # col subtitles
        fig.text(bbox_legend.x0 - 0.011, bbox_legend.y0 + (0.25 + 0.50*i)*legend_height, subtitle, ha='right', va='center', fontsize=21)  # row subtitles

    # Adjust layout and save figure
    #plt.tight_layout()
    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


# Helper function to plot timeseries for a given variable list and climate categories.
def plot_climCategory_timeseries(ax, df_dly, df_monthly, var_dict, varClimCombined_colName, varClimCombined_labels, climCat_params):

    var_pltName = var_dict['var_pltName']
    var_list    = var_dict['var_list']
    ax_color    = var_dict['ax_color']
    ymin, ymax, ymin_tick, ymax_tick, ystep_tick = var_dict['y_params']
    y_str_fmt = var_dict['y_str_fmt']
    for var_colName in var_list:

        # Subset daily data to var of interest
        df_dly_var = df_dly.copy(deep=True)
        df_dly_var = df_dly_var[[var_colName]]
            
        # Map climate categories to daily data
        df_dly_var[varClimCombined_colName] = df_dly_var.index.to_period('M').map(
                                                lambda x: df_monthly.loc[x.to_timestamp(), varClimCombined_colName]
                                                if x.to_timestamp() in df_monthly.index else None
                                                )

        # Compute daily averages and std deviations
        df_dly_var['Day'] = df_dly_var.index.day
        df_dly_var_climCatStats = df_dly_var.groupby(['Day', varClimCombined_colName]).agg(['mean', 'std'])
        df_dly_var_climCatAvgs = df_dly_var_climCatStats[var_colName]['mean'].unstack(varClimCombined_colName)
        df_dly_var_climCatStds = df_dly_var_climCatStats[var_colName]['std'].unstack(varClimCombined_colName)

        # Plot timeseries and shading for std deviation
        for climCat in varClimCombined_labels:
            if climCat in df_dly_var_climCatAvgs:

                climCat_color  = climCat_params[climCat]['color']
                climCat_ls     = climCat_params[climCat]['ls']
                climCat_zorder = climCat_params[climCat]['zorder']
                # Plot the mean timeseries
                ax.plot(df_dly_var_climCatAvgs.index, df_dly_var_climCatAvgs[climCat], 
                        label=climCat, color=climCat_color, ls=climCat_ls, lw=5, alpha=0.85, zorder=climCat_zorder)
                
                # Add shading for standard deviation
                ax.fill_between(df_dly_var_climCatStds.index, 
                                df_dly_var_climCatAvgs[climCat] - df_dly_var_climCatStds[climCat], 
                                df_dly_var_climCatAvgs[climCat] + df_dly_var_climCatStds[climCat], 
                                color=climCat_color, alpha=0.1, lw=0)

    # Enable minor ticks (but only for y-axis)
    ax.minorticks_on()
    ax.tick_params(axis='x', which='minor', top=False, bottom=False)

    # Edit x axis
    numDaysInMon = len(df_dly.index.day.unique())
    ax.set_xlabel('Day of Month', fontsize=28)
    xtick_locs = [1, 8, 15, 22, 29]
    ax.set_xticks(xtick_locs)
    ax.set_xticklabels(xtick_locs, fontsize=24)
    ax.tick_params(axis='x', which='major', length=8, color='gray', direction='in', width=1.5, top=True, bottom=True) # adjust major ticks for x-axis
    ax.set_xlim(left=0.5, right=numDaysInMon + 0.5)   # set x-axis limits

    # Edit y axis
    ax.set_ylabel(var_pltName, fontsize=26, color=ax_color)
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(np.round(np.arange(ymin_tick, ymax_tick + 0.01, ystep_tick),1))
    ax.set_yticklabels(ax.get_yticks(), fontsize=24)  # set y-axis tick label size
    ax.tick_params(axis='y', which='both', labelcolor=ax_color) # Set tick label color for y-axis
    ax.tick_params(axis='y', which='minor', length=6, color='gray', direction='in', width=1.0, left=True, right=True) # adjust minor ticks for y-axis
    ax.tick_params(axis='y', which='major', length=10, color='gray', direction='in', width=1.5, left=True, right=True) # adjust major ticks for y-axis
    ax.yaxis.set_major_formatter(FormatStrFormatter(y_str_fmt))

    # Set border color and thickness
    for spine in ax.spines.values():
        spine.set_color(ax_color)
        spine.set_linewidth(3)






##################################################
### CLIMCAT TERCILE PLOTS ########################


# adultPop_colName indicates which adult population variable to plot: 'pop_A(E)', 'pop_A(L)', 'pop_A(P)'
def plot_bar_climCategories_adultPop(df_in, varTemp_colName, varTemp_pltName, varPrcp_colName, varPrcp_pltName, adultPop_colName, fileName):
    dir_out = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/'  # TEMPORARY!

    df_dly = df_in.copy(deep=True)  # create a temporary copy of the df so that we don't alter the original
    df_monthly = convert_df_vectorPop_dly2mon(df_dly) # convert to monthly

    # Step 1: Create binary categories for temp and prcp variables
    varTemp_binary = categorizeSeries_LowHigh_byMon(df_monthly[varTemp_colName])
    varPrcp_binary = categorizeSeries_LowHigh_byMon(df_monthly[varPrcp_colName])

    # If categorization failed for either variable, we're not going to make a plot.
    if varTemp_binary.isnull().all():
       print(f'Error: categorization failed for variable {varTemp_colName}. Skipping plot generation.')
       return
    if varPrcp_binary.isnull().all():
       print(f'Error: categorization failed for variable {varPrcp_colName}. Skipping plot generation.')
       return

    # Combine varTemp and varPrcp categories
    label_varTemp_lo = 'Low '  + varTemp_pltName
    label_varTemp_hi = 'High ' + varTemp_pltName
    label_varPrcp_lo = 'Low '  + varPrcp_pltName
    label_varPrcp_hi = 'High ' + varPrcp_pltName
    varTemp_binary = varTemp_binary.cat.rename_categories({'Low': label_varTemp_lo, 'High': label_varTemp_hi}) # rename categories to differentiate the variables
    varPrcp_binary = varPrcp_binary.cat.rename_categories({'Low': label_varPrcp_lo, 'High': label_varPrcp_hi}) # rename categories to differentiate the variables
    varClimCombined = varTemp_binary.astype(str) + ' + ' + varPrcp_binary.astype(str) # combined categories
    varClimCombined_labels = [label_varTemp_lo + ' + ' + label_varPrcp_lo, label_varTemp_lo + ' + ' + label_varPrcp_hi, 
                              label_varTemp_hi + ' + ' + label_varPrcp_lo, label_varTemp_hi + ' + ' + label_varPrcp_hi]
    varClimCombined = pd.Categorical(varClimCombined, categories=varClimCombined_labels, ordered=True)

    # Add combined varClim binary to the DataFrame
    varClimCombined_colName = varTemp_colName + '_' + varPrcp_colName + '_combined'
    df_monthly[varClimCombined_colName] = varClimCombined

    # Step 2: Create tercile categories for adult pop variables
    # If categorization fails, try alternate categorization where all zeros are binned together.
    # If even that fails, we don't make the plot.
    varPopAx_tercile = categorizeSeries_LowMedHigh_byMon(df_monthly[adultPop_colName])
    varPopAx_tercile_colName = adultPop_colName + '_tercile'
    df_monthly[varPopAx_tercile_colName] = varPopAx_tercile # add adultPop tercile to the DataFrame

    # Step 3: Create a crosstab to calculate the counts of adult pop tercile within each climate category
    crosstab = pd.crosstab(df_monthly[varClimCombined_colName], df_monthly[varPopAx_tercile_colName])
    
    # Step 4: Normalize the crosstab to get percentages
    crosstab_norm = crosstab.div(crosstab.sum(axis=1), axis=0) * 100

    # Step 4: Perform chi-square test of independence
    chi2_stat, chi2_p_val, _, _ = chi2_contingency(crosstab)
    print(" ")
    print("Chi-square test")
    print("*******************************")
    print("  Chi-Square Statistic:", chi2_stat)
    print("  p-value:", chi2_p_val)

    # Extract climate categories
    categories = crosstab_norm.index.tolist()

    # Adjusts color brightness (color = hex color code or RGBA tuple; factor = mult. factor, <1 for darker, >1 for lighter)
    def adj_bright(color, factor):
        rgba = to_rgba(color)
        adjusted = tuple(min(1, max(0, channel * factor)) for channel in rgba[:3]) + (rgba[3],)
        return to_hex(adjusted)

    # Colors for categories
    # make the second color (light blue) slightly darker to be more readable
    catColors = [adj_bright('#529ef7', 0.9), '#7737de', '#d95f02', '#1b9e77'] # Low Temp + Low Prcp, Low Temp + High Prcp, High Temp + Low Prcp,  High Temp + High Prcp

    # Colors for terciles
    tercColors = ['#edf8b1', '#7fcdbb', '#2c7fb8']

    # Initialize plot
    fig, ax = plt.subplots(figsize=(8, 8))

    bar_width = 0.80
    bars = crosstab_norm.plot(kind='bar', stacked=True, width=bar_width, color=tercColors, ax=ax, zorder=2)

    # Adding labels with percents/counts on each section of the stacked bars
    for i, rect in enumerate(bars.containers):
        for j, bar in enumerate(rect):
            # Calculate the percents/counts for this section
            quantile_percent = round(crosstab_norm.iloc[j, i])
            quantile_count = crosstab.iloc[j, i]

            if quantile_count > 0: # only do labeling if this category is not empty
                # Get x and y placement of the label based on rectangle location
                x_value = bar.get_x() + bar.get_width() / 2
                y_value = bar.get_y() + bar.get_height() / 2

                # Format and add label
                label = f'{quantile_percent}%\n({quantile_count})'
                ax.text(x_value, y_value, label, ha='center', va='center', fontsize=13, fontweight='bold', color='black', zorder=4)

    # Add text indicating statistical test results outside the plot
    signif_level = 0.05
    stat_results_text = (r"$\bf{\chi^2\ test:}$" + " " +
                         r'$\chi^2$' + f" = {chi2_stat:.1f}, " + 
                         "p = " + (r"$\bf{" if chi2_p_val < signif_level else "") + f"{chi2_p_val:.2e}" + (r"*}$" if chi2_p_val < signif_level else ""))
    ax.text(0.5, 1.01, stat_results_text, transform=ax.transAxes, fontsize=14, va='bottom', ha='center')

    # Formatting
    category_labels = [label.replace(' + ', '\n') for label in categories]
    category_labels_4plt = []
    for item in category_labels:
        # Replace "TA" with "TEMP"
        item = item.replace("TA", "TEMP")
        # Split the string by newline (\n) and handle each part separately
        parts = item.split('\n')
        parts[0] = parts[0].lower() if parts[0].startswith("Low") else parts[0].upper()
        parts[1] = parts[1].lower() if parts[1].startswith("Low") else parts[1].upper()
        # Join the parts back together and append to the new list
        category_labels_4plt.append('\n'.join(parts))

    ax.set_xticklabels(category_labels_4plt, ha='center', fontsize=15, fontweight='bold', rotation=0)
    for tick_label, color in zip(ax.get_xticklabels(), catColors): # Apply colors to the x-axis tick labels
        tick_label.set_color(color)
    ax.tick_params(axis='y', labelsize=14)
    ax.set_xlabel('') #ax.set_xlabel('Climate Categories', fontsize=16)
    ax.set_ylabel('Frequency of association with\nsimulated adult population tercile (%)', fontsize=16)
    #ax.set_title(loc, fontsize=16)
    for spine in ax.spines.values():
        spine.set_zorder(5)  # Set zorder for all spines to 5 (to be above bars)
    ax.grid(axis='y', linestyle='--', alpha=0.7, zorder=0)
    ax.set_ylim(0, 105)
    ax.legend(handles=[], frameon=False) # remove legend

    # Adjust layout and save figure
    plt.tight_layout()
    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)




def plot_bar_climCategories_adultPops_all(df_in, varTemp_colName, varTemp_pltName, varPrcp_colName, varPrcp_pltName, fileName):
    dir_out = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/'  # TEMPORARY!

    df_dly = df_in.copy(deep=True)  # create a temporary copy of the df so that we don't alter the original
    df_monthly = convert_df_vectorPop_dly2mon(df_dly) # convert to monthly

    # Step 1: Create binary categories for temp and prcp variables
    varTemp_binary = categorizeSeries_LowHigh_byMon(df_monthly[varTemp_colName])
    varPrcp_binary = categorizeSeries_LowHigh_byMon(df_monthly[varPrcp_colName])

    # If categorization failed for either variable, we're not going to make a plot.
    if varTemp_binary.isnull().all():
       print(f'Error: categorization failed for variable {varTemp_colName}. Skipping plot generation.')
       return
    if varPrcp_binary.isnull().all():
       print(f'Error: categorization failed for variable {varPrcp_colName}. Skipping plot generation.')
       return

    # Combine varTemp and varPrcp categories
    label_varTemp_lo = 'Low '  + varTemp_pltName
    label_varTemp_hi = 'High ' + varTemp_pltName
    label_varPrcp_lo = 'Low '  + varPrcp_pltName
    label_varPrcp_hi = 'High ' + varPrcp_pltName
    varTemp_binary = varTemp_binary.cat.rename_categories({'Low': label_varTemp_lo, 'High': label_varTemp_hi}) # rename categories to differentiate the variables
    varPrcp_binary = varPrcp_binary.cat.rename_categories({'Low': label_varPrcp_lo, 'High': label_varPrcp_hi}) # rename categories to differentiate the variables
    varClimCombined = varTemp_binary.astype(str) + ' + ' + varPrcp_binary.astype(str) # combined categories
    varClimCombined_labels = [label_varTemp_lo + ' + ' + label_varPrcp_lo, label_varTemp_lo + ' + ' + label_varPrcp_hi, 
                              label_varTemp_hi + ' + ' + label_varPrcp_lo, label_varTemp_hi + ' + ' + label_varPrcp_hi]
    varClimCombined = pd.Categorical(varClimCombined, categories=varClimCombined_labels, ordered=True)

    # Add combined varClim binary to the DataFrame
    varClimCombined_colName = varTemp_colName + '_' + varPrcp_colName + '_combined'
    df_monthly[varClimCombined_colName] = varClimCombined

    # Step 2: Create tercile categories for adult pop variables
    # If categorization fails, try alternate categorization where all zeros are binned together.
    # If even that fails, we don't make the plot.
    adultPop_colNames = ['pop_A(E)', 'pop_A(L)', 'pop_A(P)']
    crosstabs_dict = {}
    crosstabs_norm_dict = {}
    for adultPop_colName in adultPop_colNames:
        varPopAx_tercile = categorizeSeries_LowMedHigh_byMon(df_monthly[adultPop_colName])
        varPopAx_tercile_colName = adultPop_colName + '_tercile'
        df_monthly[varPopAx_tercile_colName] = varPopAx_tercile # add adultPop tercile to the DataFrame

        # Step 3: Create a crosstab to calculate the counts of adult pop tercile within each climate category
        crosstabs_dict[adultPop_colName] = pd.crosstab(df_monthly[varClimCombined_colName], df_monthly[varPopAx_tercile_colName])
        
        # Step 4: Normalize the crosstab to get percentages
        crosstabs_norm_dict[adultPop_colName] = crosstabs_dict[adultPop_colName].div(crosstabs_dict[adultPop_colName].sum(axis=1), axis=0) * 100


    # Extract climate categories and tercile labels
    categories = crosstabs_norm_dict['pop_A(L)'].index.tolist()
    terciles = ['Low', 'Medium', 'High']

    # Adjusts color brightness (color = hex color code or RGBA tuple; factor = mult. factor, <1 for darker, >1 for lighter)
    def adj_bright(color, factor):
        rgba = to_rgba(color)
        adjusted = tuple(min(1, max(0, channel * factor)) for channel in rgba[:3]) + (rgba[3],)
        return to_hex(adjusted)

    # Colors for categories
    # make the second color (light blue) slightly darker to be more readable
    catColors = [adj_bright('#529ef7', 0.9), '#7737de', '#d95f02', '#1b9e77'] # Low Temp + Low Prcp, Low Temp + High Prcp, High Temp + Low Prcp,  High Temp + High Prcp

    # Colors for terciles (pop_A(E) is slightly lighter, pop_A(P) is slightly darker)
    tercColors_L = ['#edf8b1', '#7fcdbb', '#2c7fb8']
    tercColors_E = [adj_bright(tercColors_L[0], 1.1), adj_bright(tercColors_L[1], 1.2), adj_bright(tercColors_L[2], 1.3)]
    tercColors_P = [adj_bright(tercColors_L[0], 0.85), adj_bright(tercColors_L[1], 0.80), adj_bright(tercColors_L[2], 0.70)]

    # Initialize plot
    fig, ax = plt.subplots(figsize=(12, 8))

    # Define bar positions and width
    n_datasets = len(crosstabs_norm_dict)  # Number of datasets (pop_A(E), pop_A(L), pop_A(P))
    n_categories = len(categories)
    bar_width = 0.25
    spacing_factor = 1.1 # space between bars in a given climate category
    x = np.arange(n_categories)

    # Plot each dataset as a group of stacked bars
    for i, (dataset, values) in enumerate(crosstabs_norm_dict.items()):

        # Assign colors based on the dataset
        if dataset == 'pop_A(E)':
            tercColors = tercColors_E
        elif dataset == 'pop_A(L)':
            tercColors = tercColors_L
        elif dataset == 'pop_A(P)':
            tercColors = tercColors_P

        # Get data for each category
        tercile_values = np.array(list(values.values)).T  # Shape: (3, n_categories)
        
        # Compute offsets for the bar groups
        bar_positions = x + ((i + 1 - (n_datasets-1) / 2) * bar_width * spacing_factor)
        
        # Plot stacked bars for each tercile
        bottom = np.zeros(n_categories)
        for j, tercile in enumerate(terciles):
            bars = ax.bar(bar_positions, tercile_values[j], bar_width, bottom=bottom, color=tercColors[j], label=f"{tercile}" if i == 0 else "", zorder=2)

            # Add annotations to each segment
            for bar, value in zip(bars, tercile_values[j]):
                if value > 0:  # Only annotate non-zero segments
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,  # X position: center of the bar
                        bar.get_y() + bar.get_height() / 2,  # Y position: middle of the segment
                        f"{value:.0f}%",  # Text: value as a percentage
                        ha='center', va='center', fontsize=11, fontweight='bold', color='black', zorder=4
                    )

            bottom += tercile_values[j]

    # Formatting
    x_tick_pos = x + ((1 + 1 - (n_datasets - 1) / 2) * bar_width * spacing_factor)
    ax.set_xticks(x_tick_pos)
    category_labels = [label.replace(' + ', '\n') for label in categories]
    ax.set_xticklabels(category_labels, ha='center', fontsize=14, fontweight='bold')
    for tick_label, color in zip(ax.get_xticklabels(), catColors): # Apply colors to the x-axis tick labels
        tick_label.set_color(color)
    ax.tick_params(axis='y', labelsize=14)
    ax.set_xlabel('Climate Categories', fontsize=16)
    ax.set_ylabel('Frequency of association with\nsimulated adult population tercile (%)', fontsize=16)
    ax.set_title(loc, fontsize=16)
    for spine in ax.spines.values():
        spine.set_zorder(5)  # Set zorder for all spines to 5 (to be above bars)
    ax.grid(axis='y', linestyle='--', alpha=0.7, zorder=0)
    ax.set_ylim(0, 105)

    # Adjust layout and save figure
    plt.tight_layout()
    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    if False:
        # Create the figure and gridspec layout
        fig = plt.figure(figsize=(30, 36))
        gs_outer = GridSpec(3, 1, figure=fig, height_ratios=[2, 2, 1], hspace=0.32)  # hspace adds space for the section titles
        gs_top    = GridSpecFromSubplotSpec(2, 2, subplot_spec=gs_outer[0, 0], hspace=0.0, wspace=0.03)
        gs_center = GridSpecFromSubplotSpec(2, 3, subplot_spec=gs_outer[1, 0], hspace=0.0, wspace=0.00)
        gs_bot    = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs_outer[2, 0], hspace=0.0, wspace=0.00)

        # Top left plots
        ax_prcp = fig.add_subplot(gs_top[0, 0])
        ax_ta = fig.add_subplot(gs_top[1, 0])
        plot_climCategory_timeseries(ax_prcp, df_dly, df_monthly, met_vars['PRCP'], varClimCombined_colName, varClimCombined_labels, climCat_params)
        plot_climCategory_timeseries(ax_ta, df_dly, df_monthly, met_vars['TA'], varClimCombined_colName, varClimCombined_labels, climCat_params)

        # Top right plots
        ax_wh = fig.add_subplot(gs_top[0, 1])
        ax_tw = fig.add_subplot(gs_top[1, 1])
        plot_climCategory_timeseries(ax_wh, df_dly, df_monthly, cont_vars['WH'], varClimCombined_colName, varClimCombined_labels, climCat_params)
        plot_climCategory_timeseries(ax_tw, df_dly, df_monthly, cont_vars['TW'], varClimCombined_colName, varClimCombined_labels, climCat_params)

        # Center plots
        ax_dev_E   = fig.add_subplot(gs_center[0, 0])
        ax_dev_L   = fig.add_subplot(gs_center[0, 1])
        ax_dev_P   = fig.add_subplot(gs_center[0, 2])
        plot_climCategory_timeseries(ax_dev_E, df_dly, df_monthly, dev_vars['dev_E'], varClimCombined_colName, varClimCombined_labels, climCat_params)
        plot_climCategory_timeseries(ax_dev_L, df_dly, df_monthly, dev_vars['dev_L'], varClimCombined_colName, varClimCombined_labels, climCat_params)
        plot_climCategory_timeseries(ax_dev_P, df_dly, df_monthly, dev_vars['dev_P'], varClimCombined_colName, varClimCombined_labels, climCat_params)

        ax_surv_E  = fig.add_subplot(gs_center[1, 0])
        ax_surv_L  = fig.add_subplot(gs_center[1, 1])
        ax_surv_P  = fig.add_subplot(gs_center[1, 2])
        plot_climCategory_timeseries(ax_surv_E, df_dly, df_monthly, surv_vars['surv_E'], varClimCombined_colName, varClimCombined_labels, climCat_params)
        plot_climCategory_timeseries(ax_surv_L, df_dly, df_monthly, surv_vars['surv_L'], varClimCombined_colName, varClimCombined_labels, climCat_params)
        plot_climCategory_timeseries(ax_surv_P, df_dly, df_monthly, surv_vars['surv_P'], varClimCombined_colName, varClimCombined_labels, climCat_params)
        
        # Bottom plots
        ax_pop_A_E = fig.add_subplot(gs_bot[0, 0])
        ax_pop_A_L = fig.add_subplot(gs_bot[0, 1])
        ax_pop_A_P = fig.add_subplot(gs_bot[0, 2])
        plot_climCategory_timeseries(ax_pop_A_E, df_dly, df_monthly, pop_vars['pop_A(E)_frac'], varClimCombined_colName, varClimCombined_labels, climCat_params)
        plot_climCategory_timeseries(ax_pop_A_L, df_dly, df_monthly, pop_vars['pop_A(L)_frac'], varClimCombined_colName, varClimCombined_labels, climCat_params)
        plot_climCategory_timeseries(ax_pop_A_P, df_dly, df_monthly, pop_vars['pop_A(P)_frac'], varClimCombined_colName, varClimCombined_labels, climCat_params)

        # Adjust layout and save figure
        #plt.tight_layout()
        if not os.path.isdir(dir_out):
            os.makedirs(dir_out)
        plt.savefig(dir_out + fileName + '.png', dpi=300, bbox_inches='tight')
        plt.close(fig)






##################################################
### VIOLIN PLOTS #################################

def plot_violinplt_prcp_WH(df_in, loc, fileName):
    dir_out = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/'  # TEMPORARY!

    fig, ax = plt.subplots(figsize=(16, 6))  # Create a new figure and axis for the plot

    df = df_in.copy(deep=True)  # Create a temporary copy of the df

    # Define column names
    prcp_colName = 'PRCP (mm)'
    WH_colName = 'WH (mm)'

    # Set y-axis bounds
    ymin_prcp, ymax_prcp, ystep_prcp = 0, 640, 200
    ymin_WH, ymax_WH, ystep_WH       = 0, ymax_prcp/2, ystep_prcp/2

    # Extract month from the index to use as labels
    df['Month'] = df.index.strftime('%b')
    df['Month'] = pd.Categorical(df['Month'], categories=df['Month'].unique(), ordered=True)  # Convert to categorical

    # Rescale water height to twice the value (since we'll scale down the right-side axis by half)
    df[WH_colName] =df[WH_colName]*2

    # Prep melted dataframes for split violin plots    
    df_melted = df.melt(id_vars='Month', value_vars=[prcp_colName, WH_colName], var_name='Variable', value_name='Value')
    df_melted['Type'] = df_melted['Variable'].apply(lambda x: 'PRCP' if 'PRCP' in x else 'WH')

    # Define violin plot color mapping
    color_met  = '#74b728' # color for meteorology variables
    color_cont = '#554db8' # color for container water dynamics variables
    hue_colors = {'PRCP': color_met, 'WH': color_cont}

    # Plot the split violin plots
    violin_width = 0.85
    sns.violinplot(x='Month', y='Value', hue='Type', data=df_melted, palette=hue_colors, split=False, inner=None, ax=ax, 
                   density_norm='width', cut=0, gap=0.2, width=violin_width, legend=False, linewidth=0)

    # Set transparency for the violin plots
    for violin in ax.collections:
        violin.set_alpha(0.5)

    # Plot specific data points with offsets
    offset = violin_width / 4  # Adjust this value as needed
    ax.scatter(df['Month'].cat.codes - offset, df[prcp_colName], color='black', label='prcp', alpha=0.7, zorder=3, marker='.', s=16) # prcp points slightly to the left
    ax.scatter(df['Month'].cat.codes + offset, df[WH_colName], color='black', label='WH', alpha=0.7, zorder=3, marker='.', s=16)     # WH points slightly to the right

    # Edit axes
    ax.set_xlabel('') # remove x axis title
    ax.set_ylabel('Precipitation (mm)', fontsize=16, color=color_met) # edit y axis title
    ax.minorticks_on() # enable minor tick marks
    ax.set_xticks(range(len(df['Month'].unique())))  # set positions for x-axis ticks
    ax.set_xticklabels(df['Month'].cat.categories, fontsize=14)  # edit x-axis tick labels
    ax.set_yticks(np.arange(ymin_prcp, ymax_prcp + 0.01, ystep_prcp))  # ensure y-ticks are set explicitly
    ax.set_yticklabels(ax.get_yticks(), fontsize=14)  # set y-axis tick labels
    ax.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, left=True, right=False) # set minor ticks only for the y-axis
    ax.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, left=True, right=False) # adjust major ticks for y-axis
    ax.tick_params(axis='y', labelcolor=color_met)  # Set tick label color for y-axis
    ax.tick_params(axis='x', which='both', bottom=False, top=False)  # remove x axis ticks
    ax.set_ylim(bottom=ymin_prcp, top=ymax_prcp)   # set y-axis limits
    ax.set_xlim(left=-0.5, right=11.5)   # set x-axis limits so that it begins at January (index 0) and ends at December (index 11)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f')) # set y-axis tick label formatting

    # Create a secondary axis on the right
    ax_right = ax.twinx()
    ax_right.set_ylabel('Water height (mm)', fontsize=16, color=color_cont) # edit y axis title
    ax_right.set_yticks(np.arange(ymin_WH, ymax_WH + 0.01, ystep_WH))  # ensure y-ticks are set explicitly
    ax_right.set_yticklabels(ax_right.get_yticks(), fontsize=14)  # set y-axis tick labels
    ax_right.minorticks_on() # enable minor tick marks
    ax_right.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, left=False, right=True) # set minor ticks only for the y-axis
    ax_right.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, left=False, right=True) # adjust major ticks for y-axis
    ax_right.tick_params(axis='y', labelcolor=color_cont)  # Set tick label color for y-axis
    ax_right.set_ylim(bottom=ymin_WH, top=ymax_WH)   # set y-axis limits
    ax_right.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    # Add horizontal rectangles to mark seasons
    season_rects = [
        {'months': (0, 1),   'label': 'NEM', 'color': '#9e3138'}, # Jan-Feb
        {'months': (2, 3),   'label': 'FIM', 'color': '#d47136'}, # Mar-Apr
        {'months': (4, 8),   'label': 'SWM', 'color': '#41ab5d'}, # May-Sep
        {'months': (9, 10),  'label': 'SIM', 'color': '#6786ad'}, # Oct-Nov
        {'months': (11, 11), 'label': 'NEM', 'color': '#9e3138'}  # Dec
    ]
    
    rect_height = 0.05 # fraction of plot height
    y_pos = 1.0        # place at top of plot
    for season in season_rects:
        # Add the rectangle
        rect = patches.Rectangle(xy=(season['months'][0] / 12, y_pos - rect_height),  # Scale x to normalized coords
                                 width=(season['months'][1] - season['months'][0] + 1) / 12,  # Adjust width for months (normalized)
                                 height=rect_height, linewidth=1, edgecolor='none', facecolor=season['color'], alpha=1, zorder=5, 
                                 transform=ax.transAxes)  # transform needed for text positioning/sizing
        ax.add_patch(rect)
        
        # Add label for the season
        season_label_x = (season['months'][1] + season['months'][0] + 1) / 24  # Scale x to normalized coords
        season_label_y = y_pos - (rect_height / 2)  # Center the text in the rectangle (in normalized y)
        ax.text(x=season_label_x, y=season_label_y, s=season['label'], 
                color='black', ha='center', va='center', fontsize=10, zorder=6, fontweight='bold', 
                transform=ax.transAxes)  # transform needed for text positioning/sizing

    # Add vertical lines between months
    num_months = len(df['Month'].unique())
    for i in range(1, num_months):
        ax.axvline(x=i-0.5, color='grey', linestyle='--', lw=0.8, alpha=0.7, zorder=0)  # Add a vertical line between month ticks

    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)




def plot_violinplt_TA_TW(df_in, loc, fileName):
    dir_out = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/'  # TEMPORARY!

    fig, ax = plt.subplots(figsize=(16, 6))  # Create a new figure and axis for the plot

    df = df_in.copy(deep=True)  # Create a temporary copy of the df

    # Define column names
    TA_max_colName = 'TA_max (C)'
    TA_min_colName = 'TA_min (C)'
    TA_mean_colName = 'TA (C)'
    TW_max_colName = 'TW_max (C)'
    TW_min_colName = 'TW_min (C)'
    TW_mean_colName = 'TW (C)'

    # Get y-axis bounds based on loc
    if loc == 'Negombo':
        ymin, ymax, ystep = 18, 44, 5
        ymin_tick, ymax_tick = 20, 40
    elif loc == 'Jaffna':
        ymin, ymax, ystep = 20, 49, 5
        ymin_tick, ymax_tick = 25, 45
    elif loc == 'NuwaraEliya':
        ymin, ymax, ystep = 4, 36, 5
        ymin_tick, ymax_tick = 5, 30

    # Extract month from the index to use as labels
    df['Month'] = df.index.strftime('%b')
    df['Month'] = pd.Categorical(df['Month'], categories=df['Month'].unique(), ordered=True)  # Convert to categorical

    # Prep melted dataframes for split violin plots
    df_melted_max = df.melt(id_vars='Month', value_vars=[TA_max_colName, TW_max_colName], var_name='Variable', value_name='Value')
    df_melted_max['Type'] = df_melted_max['Variable'].apply(lambda x: 'TA' if 'TA' in x else 'TW')
    
    df_melted_min = df.melt(id_vars='Month', value_vars=[TA_min_colName, TW_min_colName], var_name='Variable', value_name='Value')
    df_melted_min['Type'] = df_melted_min['Variable'].apply(lambda x: 'TA' if 'TA' in x else 'TW')
    
    df_melted_mean = df.melt(id_vars='Month', value_vars=[TA_mean_colName, TW_mean_colName], var_name='Variable', value_name='Value')
    df_melted_mean['Type'] = df_melted_mean['Variable'].apply(lambda x: 'TA' if 'TA' in x else 'TW')

    # Define violin plot color mapping
    color_met  = '#74b728' # color for meteorology variables
    color_cont = '#554db8' # color for container water dynamics variables
    hue_colors = {'TA': color_met, 'TW': color_cont}

    # Plot the split violin plots
    violin_width = 0.85
    sns.violinplot(x='Month', y='Value', hue='Type', data=df_melted_max, palette=hue_colors, split=False, inner=None, ax=ax, 
                   density_norm='width', cut=0, gap=0.2, width=violin_width, legend=False, linewidth=0)
    sns.violinplot(x='Month', y='Value', hue='Type', data=df_melted_min, palette=hue_colors, split=False, inner=None, ax=ax, 
                   density_norm='width', cut=0, gap=0.2, width=violin_width, legend=False, linewidth=0)
    sns.violinplot(x='Month', y='Value', hue='Type', data=df_melted_mean, palette=hue_colors, split=False, inner=None, ax=ax, 
                   density_norm='width', cut=0, gap=0.2, width=violin_width, legend=False, linewidth=0)

    # Set transparency for the violin plots
    for violin in ax.collections:
        violin.set_alpha(0.5)

    # Plot specific data points with offsets
    offset = violin_width / 4  # Adjust this value as needed
    ax.scatter(df['Month'].cat.codes - offset, df[TA_max_colName], color='red', label='TA Max', alpha=0.7, zorder=3, marker='.', s=16)      # TA points slightly to the left
    ax.scatter(df['Month'].cat.codes + offset, df[TW_max_colName], color='red', label='TW Max', alpha=0.7, zorder=3, marker='.', s=16)      # TW points slightly to the right
    ax.scatter(df['Month'].cat.codes - offset, df[TA_min_colName], color='blue', label='TA Min', alpha=0.7, zorder=3, marker='.', s=16)     # TA points slightly to the left
    ax.scatter(df['Month'].cat.codes + offset, df[TW_min_colName], color='blue', label='TW Min', alpha=0.7, zorder=3, marker='.', s=16)     # TW points slightly to the right
    ax.scatter(df['Month'].cat.codes - offset, df[TA_mean_colName], color='black', label='TA Mean', alpha=0.7, zorder=3, marker='.', s=16)  # TA points slightly to the left
    ax.scatter(df['Month'].cat.codes + offset, df[TW_mean_colName], color='black', label='TW Mean', alpha=0.7, zorder=3, marker='.', s=16)  # TW points slightly to the right

    # Edit axes
    ax.set_xlabel('') # remove x axis title
    ax.set_ylabel('Air temperature (°C)', fontsize=16, color=color_met) # edit y axis title
    ax.minorticks_on() # enable minor tick marks
    ax.set_xticks(range(len(df['Month'].unique())))  # set positions for x-axis ticks
    ax.set_xticklabels(df['Month'].cat.categories, fontsize=14)  # edit x-axis tick labels
    ax.set_yticks(np.arange(ymin_tick, ymax_tick + 0.01, ystep))  # ensure y-ticks are set explicitly
    ax.set_yticklabels(ax.get_yticks(), fontsize=14)  # set y-axis tick labels
    ax.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, left=True, right=False) # set minor ticks only for the y-axis
    ax.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, left=True, right=False) # adjust major ticks for y-axis
    ax.tick_params(axis='y', labelcolor=color_met)  # Set tick label color for y-axis
    ax.tick_params(axis='x', which='both', bottom=False, top=False)  # remove x axis ticks
    ax.set_ylim(bottom=ymin, top=ymax)   # set y-axis limits
    ax.set_xlim(left=-0.5, right=11.5)   # set x-axis limits so that it begins at January (index 0) and ends at December (index 11)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f')) # set y-axis tick label formatting

    # Create a secondary axis on the right
    ax_right = ax.twinx()
    ax_right.set_ylabel('Water temperature (°C)', fontsize=16, color=color_cont) # edit y axis title
    ax_right.set_yticks(np.arange(ymin_tick, ymax_tick + 0.01, ystep))  # ensure y-ticks are set explicitly
    ax_right.set_yticklabels(ax_right.get_yticks(), fontsize=14)  # set y-axis tick labels
    ax_right.minorticks_on() # enable minor tick marks
    ax_right.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, left=False, right=True) # set minor ticks only for the y-axis
    ax_right.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, left=False, right=True) # adjust major ticks for y-axis
    ax_right.tick_params(axis='y', labelcolor=color_cont)  # Set tick label color for y-axis
    ax_right.set_ylim(bottom=ymin, top=ymax)   # set y-axis limits
    ax_right.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    # Add horizontal rectangles to mark seasons
    season_rects = [
        {'months': (0, 1),   'label': 'NEM', 'color': '#9e3138'}, # Jan-Feb
        {'months': (2, 3),   'label': 'FIM', 'color': '#d47136'}, # Mar-Apr
        {'months': (4, 8),   'label': 'SWM', 'color': '#41ab5d'}, # May-Sep
        {'months': (9, 10),  'label': 'SIM', 'color': '#6786ad'}, # Oct-Nov
        {'months': (11, 11), 'label': 'NEM', 'color': '#9e3138'}  # Dec
    ]
    
    rect_height = 0.05 # fraction of plot height
    y_pos = 1.0        # place at top of plot
    for season in season_rects:
        # Add the rectangle
        rect = patches.Rectangle(xy=(season['months'][0] / 12, y_pos - rect_height),  # Scale x to normalized coords
                                 width=(season['months'][1] - season['months'][0] + 1) / 12,  # Adjust width for months (normalized)
                                 height=rect_height, linewidth=1, edgecolor='none', facecolor=season['color'], alpha=1, zorder=5, 
                                 transform=ax.transAxes)  # transform needed for text positioning/sizing
        ax.add_patch(rect)
        
        # Add label for the season
        season_label_x = (season['months'][1] + season['months'][0] + 1) / 24  # Scale x to normalized coords
        season_label_y = y_pos - (rect_height / 2)  # Center the text in the rectangle (in normalized y)
        ax.text(x=season_label_x, y=season_label_y, s=season['label'], 
                color='black', ha='center', va='center', fontsize=10, zorder=6, fontweight='bold', 
                transform=ax.transAxes)  # transform needed for text positioning/sizing

    # Add vertical lines between months
    num_months = len(df['Month'].unique())
    for i in range(1, num_months):
        ax.axvline(x=i-0.5, color='grey', linestyle='--', lw=0.8, alpha=0.7, zorder=0)  # Add a vertical line between month ticks

    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)




def plot_violinplt_dev(df_in, loc, fileName):
    dir_out = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/'  # TEMPORARY!

    fig, ax = plt.subplots(figsize=(16, 6))  # Create a new figure and axis for the plot

    df = df_in.copy(deep=True)  # Create a temporary copy of the df

    # Define column names
    dev_E_colName = 'dev_E'
    dev_L_colName = 'dev_L'
    dev_P_colName = 'dev_P'

    # Get y-axis bounds based on loc
    if loc == 'Negombo':
        ymin, ymax, ystep = 0, 0.75, 0.2
    elif loc == 'Jaffna':
        ymin, ymax, ystep = 0, 1.07, 0.2
    elif loc == 'NuwaraEliya':
        ymin, ymax, ystep = 0, 0.37, 0.1

    # Extract month from the index to use as labels
    df['Month'] = df.index.strftime('%b')
    df['Month'] = pd.Categorical(df['Month'], categories=df['Month'].unique(), ordered=True)  # Convert to categorical

    # Prep melted dataframes for split violin plots    
    df_melted = df.melt(id_vars='Month', value_vars=[dev_E_colName, dev_L_colName, dev_P_colName], var_name='Variable', value_name='Value')
    df_melted['Type'] = df_melted['Variable'].apply(lambda x: 'dev_E' if 'dev_E' in x else ('dev_L' if 'dev_L' in x else 'dev_P'))

    # Define violin plot color mapping
    color_vectBio_E = '#f394c5' # color for vector biology (eggs) variables
    color_vectBio_L = '#b9216e' # color for vector biology (larvae) variables
    color_vectBio_P = '#450c29' # color for vector biology (pupae) variables
    hue_colors = {'dev_E': color_vectBio_E, 'dev_L': color_vectBio_L, 'dev_P': color_vectBio_P}

    # Plot the split violin plots
    violin_width = 0.95
    sns.violinplot(x='Month', y='Value', hue='Type', data=df_melted, palette=hue_colors, split=False, inner=None, ax=ax, 
                   density_norm='width', cut=0, gap=0.2, width=violin_width, legend=False, linewidth=0)

    # Set transparency for the violin plots
    for violin in ax.collections:
        violin.set_alpha(0.5)

    # Plot specific data points with offsets
    offset = violin_width / 3  # Adjust this value as needed
    ax.scatter(df['Month'].cat.codes - offset, df[dev_E_colName], color='black', label='dev_E', alpha=0.7, zorder=3, marker='.', s=16) # dev_E points slightly to the left
    ax.scatter(df['Month'].cat.codes,          df[dev_L_colName], color='black', label='dev_L', alpha=0.7, zorder=3, marker='.', s=16) # dev_L points in center
    ax.scatter(df['Month'].cat.codes + offset, df[dev_P_colName], color='black', label='dev_P', alpha=0.7, zorder=3, marker='.', s=16) # dev_E points slightly to the right


    # Edit axes
    ax.set_xlabel('') # remove x axis title
    ax.set_ylabel(r'Development rate (day$^{-1}$)', fontsize=16) # edit y axis title
    ax.minorticks_on() # enable minor tick marks
    ax.set_xticks(range(len(df['Month'].unique())))  # set positions for x-axis ticks
    ax.set_xticklabels(df['Month'].cat.categories, fontsize=14)  # edit x-axis tick labels
    ax.set_yticks(np.arange(ymin, ymax + 0.01, ystep))  # ensure y-ticks are set explicitly
    ax.set_yticklabels(ax.get_yticks(), fontsize=14)  # set y-axis tick labels
    ax.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, left=True, right=True) # set minor ticks only for the y-axis
    ax.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, left=True, right=True) # adjust major ticks for y-axis
    ax.tick_params(axis='x', which='both', bottom=False, top=False)  # remove x axis ticks
    ax.set_ylim(bottom=ymin, top=ymax)   # set y-axis limits
    ax.set_xlim(left=-0.5, right=11.5)   # set x-axis limits so that it begins at January (index 0) and ends at December (index 11)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f')) # set y-axis tick label formatting

    # Add horizontal rectangles to mark seasons
    season_rects = [
        {'months': (0, 1),   'label': 'NEM', 'color': '#9e3138'}, # Jan-Feb
        {'months': (2, 3),   'label': 'FIM', 'color': '#d47136'}, # Mar-Apr
        {'months': (4, 8),   'label': 'SWM', 'color': '#41ab5d'}, # May-Sep
        {'months': (9, 10),  'label': 'SIM', 'color': '#6786ad'}, # Oct-Nov
        {'months': (11, 11), 'label': 'NEM', 'color': '#9e3138'}  # Dec
    ]
    
    rect_height = 0.05 # fraction of plot height
    y_pos = 1.0        # place at top of plot
    for season in season_rects:
        # Add the rectangle
        rect = patches.Rectangle(xy=(season['months'][0] / 12, y_pos - rect_height),  # Scale x to normalized coords
                                 width=(season['months'][1] - season['months'][0] + 1) / 12,  # Adjust width for months (normalized)
                                 height=rect_height, linewidth=1, edgecolor='none', facecolor=season['color'], alpha=1, zorder=5, 
                                 transform=ax.transAxes)  # transform needed for text positioning/sizing
        ax.add_patch(rect)
        
        # Add label for the season
        season_label_x = (season['months'][1] + season['months'][0] + 1) / 24  # Scale x to normalized coords
        season_label_y = y_pos - (rect_height / 2)  # Center the text in the rectangle (in normalized y)
        ax.text(x=season_label_x, y=season_label_y, s=season['label'], 
                color='black', ha='center', va='center', fontsize=10, zorder=6, fontweight='bold', 
                transform=ax.transAxes)  # transform needed for text positioning/sizing

    # Add vertical lines between months
    num_months = len(df['Month'].unique())
    for i in range(1, num_months):
        ax.axvline(x=i-0.5, color='grey', linestyle='--', lw=0.8, alpha=0.7, zorder=0)  # Add a vertical line between month ticks

    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)




def plot_violinplt_surv(df_in, loc, fileName):
    dir_out = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/'  # TEMPORARY!

    fig, ax = plt.subplots(figsize=(16, 6))  # Create a new figure and axis for the plot

    df = df_in.copy(deep=True)  # Create a temporary copy of the df

    # Define column names
    surv_E_colName = 'surv_E'
    surv_L_colName = 'surv_L'
    surv_P_colName = 'surv_P'

    # Set y-axis bounds
    ymin, ymax, ystep = 0, 1.07, 0.2

    # Extract month from the index to use as labels
    df['Month'] = df.index.strftime('%b')
    df['Month'] = pd.Categorical(df['Month'], categories=df['Month'].unique(), ordered=True)  # Convert to categorical

    # Prep melted dataframes for split violin plots    
    df_melted = df.melt(id_vars='Month', value_vars=[surv_E_colName, surv_L_colName, surv_P_colName], var_name='Variable', value_name='Value')
    df_melted['Type'] = df_melted['Variable'].apply(lambda x: 'surv_E' if 'surv_E' in x else ('surv_L' if 'surv_L' in x else 'surv_P'))

    # Define violin plot color mapping
    color_vectBio_E = '#f394c5' # color for vector biology (eggs) variables
    color_vectBio_L = '#b9216e' # color for vector biology (larvae) variables
    color_vectBio_P = '#450c29' # color for vector biology (pupae) variables
    hue_colors = {'surv_E': color_vectBio_E, 'surv_L': color_vectBio_L, 'surv_P': color_vectBio_P}

    # Plot the split violin plots
    violin_width = 0.90
    sns.violinplot(x='Month', y='Value', hue='Type', data=df_melted, palette=hue_colors, split=False, inner=None, ax=ax, 
                   density_norm='width', cut=0, gap=0.2, width=violin_width, legend=False, linewidth=0)

    # Set transparency for the violin plots
    for violin in ax.collections:
        violin.set_alpha(0.5)

    # Plot specific data points with offsets
    offset = violin_width / 3  # Adjust this value as needed
    ax.scatter(df['Month'].cat.codes - offset, df[surv_E_colName], color='black', label='dev_E', alpha=0.7, zorder=3, marker='.', s=16) # dev_E points slightly to the left
    ax.scatter(df['Month'].cat.codes,          df[surv_L_colName], color='black', label='dev_L', alpha=0.7, zorder=3, marker='.', s=16) # dev_L points in center
    ax.scatter(df['Month'].cat.codes + offset, df[surv_P_colName], color='black', label='dev_P', alpha=0.7, zorder=3, marker='.', s=16) # dev_E points slightly to the right


    # Edit axes
    ax.set_xlabel('') # remove x axis title
    ax.set_ylabel(r'Survival rate (day$^{-1}$)', fontsize=16) # edit y axis title
    ax.minorticks_on() # enable minor tick marks
    ax.set_xticks(range(len(df['Month'].unique())))  # set positions for x-axis ticks
    ax.set_xticklabels(df['Month'].cat.categories, fontsize=14)  # edit x-axis tick labels
    ax.set_yticks(np.arange(ymin, ymax + 0.01, ystep))  # ensure y-ticks are set explicitly
    ax.set_yticklabels(ax.get_yticks(), fontsize=14)  # set y-axis tick labels
    ax.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, left=True, right=True) # set minor ticks only for the y-axis
    ax.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, left=True, right=True) # adjust major ticks for y-axis
    ax.tick_params(axis='x', which='both', bottom=False, top=False)  # remove x axis ticks
    ax.set_ylim(bottom=ymin, top=ymax)   # set y-axis limits
    ax.set_xlim(left=-0.5, right=11.5)   # set x-axis limits so that it begins at January (index 0) and ends at December (index 11)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f')) # set y-axis tick label formatting

    # Add horizontal rectangles to mark seasons
    season_rects = [
        {'months': (0, 1),   'label': 'NEM', 'color': '#9e3138'}, # Jan-Feb
        {'months': (2, 3),   'label': 'FIM', 'color': '#d47136'}, # Mar-Apr
        {'months': (4, 8),   'label': 'SWM', 'color': '#41ab5d'}, # May-Sep
        {'months': (9, 10),  'label': 'SIM', 'color': '#6786ad'}, # Oct-Nov
        {'months': (11, 11), 'label': 'NEM', 'color': '#9e3138'}  # Dec
    ]
    
    rect_height = 0.05 # fraction of plot height
    y_pos = 1.0        # place at top of plot
    for season in season_rects:
        # Add the rectangle
        rect = patches.Rectangle(xy=(season['months'][0] / 12, y_pos - rect_height),  # Scale x to normalized coords
                                 width=(season['months'][1] - season['months'][0] + 1) / 12,  # Adjust width for months (normalized)
                                 height=rect_height, linewidth=1, edgecolor='none', facecolor=season['color'], alpha=1, zorder=5, 
                                 transform=ax.transAxes)  # transform needed for text positioning/sizing
        ax.add_patch(rect)
        
        # Add label for the season
        season_label_x = (season['months'][1] + season['months'][0] + 1) / 24  # Scale x to normalized coords
        season_label_y = y_pos - (rect_height / 2)  # Center the text in the rectangle (in normalized y)
        ax.text(x=season_label_x, y=season_label_y, s=season['label'], 
                color='black', ha='center', va='center', fontsize=10, zorder=6, fontweight='bold', 
                transform=ax.transAxes)  # transform needed for text positioning/sizing

    # Add vertical lines between months
    num_months = len(df['Month'].unique())
    for i in range(1, num_months):
        ax.axvline(x=i-0.5, color='grey', linestyle='--', lw=0.8, alpha=0.7, zorder=0)  # Add a vertical line between month ticks

    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)




def plot_violinplt_dev_surv(df_in, loc, fileName):
    dir_out = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/'  # TEMPORARY!

    # Create figure with two vertically stacked subplots, sharing the x-axis
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 9), sharex=True, gridspec_kw={'height_ratios': [1, 2], 'hspace': 0})

    df = df_in.copy(deep=True)

    # Define column names
    dev_E_colName = 'dev_E'
    dev_L_colName = 'dev_L'
    dev_P_colName = 'dev_P'
    surv_E_colName = 'surv_E'
    surv_L_colName = 'surv_L'
    surv_P_colName = 'surv_P'

    # Set y-axis bounds based on loc
    if loc == 'Negombo':
        ymin_dev, ymax_dev, ystep_dev = 0, 0.79, 0.2
        ymin_surv, ymax_surv, ystep_surv = 0, 1.02, 0.2
    elif loc == 'Jaffna':
        ymin_dev, ymax_dev, ystep_dev = 0, 1.12, 0.2
        ymin_surv, ymax_surv, ystep_surv = 0, 1.02, 0.2
    elif loc == 'NuwaraEliya':
        ymin_dev, ymax_dev, ystep_dev = 0, 0.37, 0.1
        ymin_surv, ymax_surv, ystep_surv = 0, 1.02, 0.2

    # Extract month from the index to use as labels
    df['Month'] = df.index.strftime('%b')
    df['Month'] = pd.Categorical(df['Month'], categories=df['Month'].unique(), ordered=True)

    # Melt the DataFrame for dev and surv plots
    df_melted_dev = df.melt(id_vars='Month', value_vars=[dev_E_colName, dev_L_colName, dev_P_colName], var_name='Variable', value_name='Value')
    df_melted_dev['Type'] = df_melted_dev['Variable'].apply(lambda x: 'dev_E' if 'dev_E' in x else ('dev_L' if 'dev_L' in x else 'dev_P'))
    df_melted_surv = df.melt(id_vars='Month', value_vars=[surv_E_colName, surv_L_colName, surv_P_colName], var_name='Variable', value_name='Value')
    df_melted_surv['Type'] = df_melted_surv['Variable'].apply(lambda x: 'surv_E' if 'surv_E' in x else ('surv_L' if 'surv_L' in x else 'surv_P'))

    # Define color mappings
    color_vectBio_E = '#f394c5' # color for vector biology (eggs) variables
    color_vectBio_L = '#b9216e' # color for vector biology (larvae) variables
    color_vectBio_P = '#450c29' # color for vector biology (pupae) variables
    hue_colors_surv = {'surv_E': color_vectBio_E, 'surv_L': color_vectBio_L, 'surv_P': color_vectBio_P}
    hue_colors_dev  = {'dev_E': color_vectBio_E, 'dev_L': color_vectBio_L, 'dev_P': color_vectBio_P}

    violin_width = 0.90
    offset = violin_width / 3

    # Development rate plot (top - ax1)
    sns.violinplot(x='Month', y='Value', hue='Type', data=df_melted_dev, palette=hue_colors_dev, split=False, inner=None, ax=ax1, 
                   density_norm='width', cut=0, gap=0.2, width=violin_width, legend=False, linewidth=0)
    
    # Set transparency for the violin plots
    for violin in ax1.collections:
        violin.set_alpha(0.5)
        
    # Plot specific data points with offsets
    ax1.scatter(df['Month'].cat.codes - offset, df[dev_E_colName], color='black', alpha=0.7, zorder=3, marker='.', s=16)
    ax1.scatter(df['Month'].cat.codes,          df[dev_L_colName], color='black', alpha=0.7, zorder=3, marker='.', s=16)
    ax1.scatter(df['Month'].cat.codes + offset, df[dev_P_colName], color='black', alpha=0.7, zorder=3, marker='.', s=16)

    # Survival rate plot (bottom - ax2)
    sns.violinplot(x='Month', y='Value', hue='Type', data=df_melted_surv, palette=hue_colors_surv, split=False, inner=None, ax=ax2, 
                   density_norm='width', cut=0, gap=0.2, width=violin_width, legend=False, linewidth=0)
    
    # Set transparency for the violin plots
    for violin in ax2.collections:
        violin.set_alpha(0.5)
    
    # Plot specific data points with offsets
    ax2.scatter(df['Month'].cat.codes - offset, df[surv_E_colName], color='black', alpha=0.7, zorder=3, marker='.', s=16)
    ax2.scatter(df['Month'].cat.codes,          df[surv_L_colName], color='black', alpha=0.7, zorder=3, marker='.', s=16)
    ax2.scatter(df['Month'].cat.codes + offset, df[surv_P_colName], color='black', alpha=0.7, zorder=3, marker='.', s=16)


    # Set y-axis for development rate plot on the right side
    ax1.set_ylabel(r'Dev. rate (day$^{-1}$)', fontsize=16)
    ax1.yaxis.set_label_position('right')
    ax1.yaxis.tick_right()
    ax1.set_ylim(ymin_dev, ymax_dev)
    ax1.set_yticks(np.round(np.arange(ymin_dev, ymax_dev + 0.01, ystep_dev),1))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    # Set y-axis for survival rate plot on the left side
    ax2.set_ylabel(r'Survival rate (day$^{-1}$)', fontsize=16)
    ax2.set_ylim(ymin_surv, ymax_surv)
    ax2.set_yticks(np.round(np.arange(ymin_surv, ymax_surv + 0.01, ystep_surv),1))
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    # Formatting x-axis (months)
    for ax in [ax1, ax2]:
        ax.set_xlabel('') # remove x axis title
        ax.minorticks_on() # enable minor tick marks
        ax.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, left=True, right=True) # set minor ticks only for the y-axis
        ax.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, left=True, right=True) # adjust major ticks for y-axis
        ax.tick_params(axis='x', which='both', bottom=False, top=False)  # remove x axis ticks
        ax.set_xlim(left=-0.5, right=11.5)   # set x-axis limits so that it begins at January (index 0) and ends at December (index 11)
        ax.set_yticklabels(ax.get_yticks(), fontsize=14)  # set y-axis tick labels
    ax2.set_xticks(range(len(df['Month'].unique())))
    ax2.set_xticklabels(df['Month'].cat.categories, fontsize=14)

    # Add season patches to both plots
    season_rects = [
        {'months': (0, 1),   'label': 'NEM', 'color': '#9e3138'},
        {'months': (2, 3),   'label': 'FIM', 'color': '#d47136'},
        {'months': (4, 8),   'label': 'SWM', 'color': '#41ab5d'},
        {'months': (9, 10),  'label': 'SIM', 'color': '#6786ad'},
        {'months': (11, 11), 'label': 'NEM', 'color': '#9e3138'}
    ]
    
    rect_height = 0.05*2
    y_pos = 1.0
    for season in season_rects:
        rect = patches.Rectangle(xy=(season['months'][0] / 12, y_pos - rect_height), 
                                    width=(season['months'][1] - season['months'][0] + 1) / 12, 
                                    height=rect_height, linewidth=1, edgecolor='none', facecolor=season['color'], alpha=1, zorder=5, 
                                    transform=ax1.transAxes)
        ax1.add_patch(rect)
        season_label_x = (season['months'][1] + season['months'][0] + 1) / 24
        season_label_y = y_pos - (rect_height / 2)
        ax1.text(x=season_label_x, y=season_label_y, s=season['label'], color='black', ha='center', va='center', fontsize=10, 
                 zorder=6, fontweight='bold', transform=ax1.transAxes)

    # Add vertical lines between months
    num_months = len(df['Month'].unique())
    for i in range(1, num_months):
        ax1.axvline(x=i-0.5, color='grey', linestyle='--', lw=0.8, alpha=0.7, zorder=0)  # Add a vertical line between month ticks
        ax2.axvline(x=i-0.5, color='grey', linestyle='--', lw=0.8, alpha=0.7, zorder=0)  # Add a vertical line between month ticks

    # Save the figure
    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)




def plot_violinplt_popFrac(df_in, loc, fileName):
    dir_out = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/'  # TEMPORARY!

    fig, ax = plt.subplots(figsize=(16, 6))  # Create a new figure and axis for the plot

    df = df_in.copy(deep=True)  # Create a temporary copy of the df

    # Define column names
    popFrac_A_E_colName = 'pop_A(E)_frac'
    popFrac_A_L_colName = 'pop_A(L)_frac'
    popFrac_A_P_colName = 'pop_A(P)_frac'

    # Define y-axis bounds
    ymin, ymax, ystep = 0, 0.96, 0.2

    # Extract month from the index to use as labels
    df['Month'] = df.index.strftime('%b')
    df['Month'] = pd.Categorical(df['Month'], categories=df['Month'].unique(), ordered=True)  # Convert to categorical

    # Prep melted dataframes for split violin plots    
    df_melted = df.melt(id_vars='Month', value_vars=[popFrac_A_E_colName, popFrac_A_L_colName, popFrac_A_P_colName], var_name='Variable', value_name='Value')
    df_melted['Type'] = df_melted['Variable'].apply(lambda x: 'pop_A(E)_frac' if 'pop_A(E)_frac' in x else ('pop_A(L)_frac' if 'pop_A(L)_frac' in x else 'pop_A(P)_frac'))

    # Define violin plot color mapping
    color_vectPop_E  = '#76c5ad' # color for vector pop (eggs) variables
    color_vectPop_L  = '#1b9e77' # color for vector pop (larvae) variables
    color_vectPop_P  = '#0b3f30' # color for vector pop (pupae) variables
    hue_colors = {'pop_A(E)_frac': color_vectPop_E, 'pop_A(L)_frac': color_vectPop_L, 'pop_A(P)_frac': color_vectPop_P}

    # Plot the split violin plots
    violin_width = 0.90
    sns.violinplot(x='Month', y='Value', hue='Type', data=df_melted, palette=hue_colors, split=False, inner=None, ax=ax, 
                   density_norm='width', cut=0, gap=0.2, width=violin_width, legend=False, linewidth=0)

    # Set transparency for the violin plots
    for violin in ax.collections:
        violin.set_alpha(0.5)

    # Plot specific data points with offsets
    offset = violin_width / 3  # Adjust this value as needed
    ax.scatter(df['Month'].cat.codes - offset, df[popFrac_A_E_colName], color='black', label='pop_A_E_frac', alpha=0.7, zorder=3, marker='.', s=16) # points slightly to left
    ax.scatter(df['Month'].cat.codes,          df[popFrac_A_L_colName], color='black', label='pop_A_L_frac', alpha=0.7, zorder=3, marker='.', s=16) # points centered
    ax.scatter(df['Month'].cat.codes + offset, df[popFrac_A_P_colName], color='black', label='pop_A_P_frac', alpha=0.7, zorder=3, marker='.', s=16) # points slightly to right

    # Edit axes
    ax.set_xlabel('') # remove x axis title
    ax.set_ylabel('Norm. adult population', fontsize=16) # edit y axis title
    ax.minorticks_on() # enable minor tick marks
    ax.set_xticks(range(len(df['Month'].unique())))  # set positions for x-axis ticks
    ax.set_xticklabels(df['Month'].cat.categories, fontsize=14)  # edit x-axis tick labels
    ax.set_yticks(np.arange(ymin, ymax + 0.01, ystep))  # ensure y-ticks are set explicitly
    ax.set_yticklabels(ax.get_yticks(), fontsize=14)  # set y-axis tick labels
    ax.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, left=True, right=False) # set minor ticks only for the y-axis
    ax.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, left=True, right=False) # adjust major ticks for y-axis
    ax.tick_params(axis='y', labelcolor='black')  # Set tick label color for y-axis
    ax.tick_params(axis='x', which='both', bottom=False, top=False)  # remove x axis ticks
    ax.set_ylim(bottom=ymin, top=ymax)   # set y-axis limits
    ax.set_xlim(left=-0.5, right=11.5)   # set x-axis limits so that it begins at January (index 0) and ends at December (index 11)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f')) # set y-axis tick label formatting

    # Add horizontal rectangles to mark seasons
    season_rects = [
        {'months': (0, 1),   'label': 'NEM', 'color': '#9e3138'}, # Jan-Feb
        {'months': (2, 3),   'label': 'FIM', 'color': '#d47136'}, # Mar-Apr
        {'months': (4, 8),   'label': 'SWM', 'color': '#41ab5d'}, # May-Sep
        {'months': (9, 10),  'label': 'SIM', 'color': '#6786ad'}, # Oct-Nov
        {'months': (11, 11), 'label': 'NEM', 'color': '#9e3138'}  # Dec
    ]
    
    rect_height = 0.05 # fraction of plot height
    y_pos = 1.0        # place at top of plot
    for season in season_rects:
        # Add the rectangle
        rect = patches.Rectangle(xy=(season['months'][0] / 12, y_pos - rect_height),  # Scale x to normalized coords
                                 width=(season['months'][1] - season['months'][0] + 1) / 12,  # Adjust width for months (normalized)
                                 height=rect_height, linewidth=1, edgecolor='none', facecolor=season['color'], alpha=1, zorder=5, 
                                 transform=ax.transAxes)  # transform needed for text positioning/sizing
        ax.add_patch(rect)
        
        # Add label for the season
        season_label_x = (season['months'][1] + season['months'][0] + 1) / 24  # Scale x to normalized coords
        season_label_y = y_pos - (rect_height / 2)  # Center the text in the rectangle (in normalized y)
        ax.text(x=season_label_x, y=season_label_y, s=season['label'], 
                color='black', ha='center', va='center', fontsize=10, zorder=6, fontweight='bold', 
                transform=ax.transAxes)  # transform needed for text positioning/sizing

    # Add vertical lines between months
    num_months = len(df['Month'].unique())
    for i in range(1, num_months):
        ax.axvline(x=i-0.5, color='grey', linestyle='--', lw=0.8, alpha=0.7, zorder=0)  # Add a vertical line between month ticks

    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)




def plot_violinplt_denInc(df_in, loc, fileName):
    dir_out = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/'  # TEMPORARY!

    fig, ax = plt.subplots(figsize=(16, 6))  # Create a new figure and axis for the plot

    df = df_in.copy(deep=True)  # Create a temporary copy of the df

    # Define column names
    denInc_colName = 'dengueInc'

    # Get y-axis bounds based on loc
    if loc == 'Negombo':
        ymin, ymax, ystep = 0, 49, 10
    elif loc == 'Jaffna':
        ymin, ymax, ystep = 0, 179, 50
    elif loc == 'NuwaraEliya':
        ymin, ymax, ystep = 0, 13, 5

    # Extract month from the index to use as labels
    df['Month'] = df.index.strftime('%b')
    df['Month'] = pd.Categorical(df['Month'], categories=df['Month'].unique(), ordered=True)  # Convert to categorical

    # Define violin plot color
    color_den = '#d95f02' # color for dengue variables

    # Plot the split violin plots
    violin_width = 0.70
    sns.violinplot(x='Month', y=denInc_colName, data=df, color=color_den, inner=None, ax=ax, 
                   density_norm='width', cut=0, gap=0.2, width=violin_width, legend=False, linewidth=0)

    # Set transparency for the violin plots
    for violin in ax.collections:
        violin.set_alpha(0.5)

    # Plot specific data points with offsets
    ax.scatter(df['Month'].cat.codes, df[denInc_colName], color='black', label='dengueInc', alpha=0.7, zorder=3, marker='.', s=16)     # points to right

    # Edit axes
    ax.set_xlabel('') # remove x axis title
    ax.set_ylabel('Dengue incidence (per 100,000)', fontsize=16, color='black') # edit y axis title
    ax.minorticks_on() # enable minor tick marks
    ax.set_xticks(range(len(df['Month'].unique())))  # set positions for x-axis ticks
    ax.set_xticklabels(df['Month'].cat.categories, fontsize=14)  # edit x-axis tick labels
    ax.set_yticks(np.arange(ymin, ymax + 0.01, ystep))  # ensure y-ticks are set explicitly
    ax.set_yticklabels(ax.get_yticks(), fontsize=14)  # set y-axis tick labels
    ax.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, left=True, right=False) # set minor ticks only for the y-axis
    ax.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, left=True, right=False) # adjust major ticks for y-axis
    ax.tick_params(axis='y', labelcolor='black')  # Set tick label color for y-axis
    ax.tick_params(axis='x', which='both', bottom=False, top=False)  # remove x axis ticks
    ax.set_ylim(bottom=ymin, top=ymax)   # set y-axis limits
    ax.set_xlim(left=-0.5, right=11.5)   # set x-axis limits so that it begins at January (index 0) and ends at December (index 11)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f')) # set y-axis tick label formatting

    # Add horizontal rectangles to mark seasons
    season_rects = [
        {'months': (0, 1),   'label': 'NEM', 'color': '#9e3138'}, # Jan-Feb
        {'months': (2, 3),   'label': 'FIM', 'color': '#d47136'}, # Mar-Apr
        {'months': (4, 8),   'label': 'SWM', 'color': '#41ab5d'}, # May-Sep
        {'months': (9, 10),  'label': 'SIM', 'color': '#6786ad'}, # Oct-Nov
        {'months': (11, 11), 'label': 'NEM', 'color': '#9e3138'}  # Dec
    ]
    
    rect_height = 0.05 # fraction of plot height
    y_pos = 1.0        # place at top of plot
    for season in season_rects:
        # Add the rectangle
        rect = patches.Rectangle(xy=(season['months'][0] / 12, y_pos - rect_height),  # Scale x to normalized coords
                                 width=(season['months'][1] - season['months'][0] + 1) / 12,  # Adjust width for months (normalized)
                                 height=rect_height, linewidth=1, edgecolor='none', facecolor=season['color'], alpha=1, zorder=5, 
                                 transform=ax.transAxes)  # transform needed for text positioning/sizing
        ax.add_patch(rect)
        
        # Add label for the season
        season_label_x = (season['months'][1] + season['months'][0] + 1) / 24  # Scale x to normalized coords
        season_label_y = y_pos - (rect_height / 2)  # Center the text in the rectangle (in normalized y)
        ax.text(x=season_label_x, y=season_label_y, s=season['label'], 
                color='black', ha='center', va='center', fontsize=10, zorder=6, fontweight='bold', 
                transform=ax.transAxes)  # transform needed for text positioning/sizing

    # Add vertical lines between months
    num_months = len(df['Month'].unique())
    for i in range(1, num_months):
        ax.axvline(x=i-0.5, color='grey', linestyle='--', lw=0.8, alpha=0.7, zorder=0)  # Add a vertical line between month ticks

    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)




def plot_violinplt_popFrac_denInc(df_in, loc, fileName):
    dir_out = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/'  # TEMPORARY!

    fig, ax = plt.subplots(figsize=(16, 6))  # Create a new figure and axis for the plot

    df = df_in.copy(deep=True)  # Create a temporary copy of the df

    # Define column names
    popFrac_A_E_colName = 'pop_A(E)_frac'
    popFrac_A_L_colName = 'pop_A(L)_frac'
    popFrac_A_P_colName = 'pop_A(P)_frac'
    denInc_colName = 'dengueInc'

    # Get y-axis bounds based on loc
    ymin_popFrac, ymax_popFrac, ystep_popFrac = 0, 0.96, 0.2
    if loc == 'Negombo':
        ymin_denInc, ymax_denInc, ystep_denInc = 0, 49, 10
    elif loc == 'Jaffna':
        ymin_denInc, ymax_denInc, ystep_denInc = 0, 179, 50
    elif loc == 'NuwaraEliya':
        ymin_denInc, ymax_denInc, ystep_denInc = 0, 13, 5

    # Extract month from the index to use as labels
    df['Month'] = df.index.strftime('%b')
    df['Month'] = pd.Categorical(df['Month'], categories=df['Month'].unique(), ordered=True)  # Convert to categorical

    # Rescale dengue incidence based on the scaling of the left y-axis vs the right y-axis (we plot dengue incidence with the
    # left axis, but the values are meant to be read off the right axis)
    ax_right_scaling = (ymax_popFrac-ymin_popFrac)/(ymax_denInc-ymin_denInc)
    df[denInc_colName] =df[denInc_colName]*ax_right_scaling

    # Prep melted dataframes for split violin plots    
    df_melted = df.melt(id_vars='Month', value_vars=[popFrac_A_E_colName, popFrac_A_L_colName, popFrac_A_P_colName, denInc_colName], var_name='Variable', value_name='Value')
    df_melted['Type'] = df_melted['Variable'].apply(lambda x: 'pop_A(E)_frac' if 'pop_A(E)_frac' in x else 
                                                               ('pop_A(L)_frac' if 'pop_A(L)_frac' in x else 
                                                                ('pop_A(P)_frac' if 'pop_A(P)_frac' in x else 'dengueInc')))

    # Define violin plot color mapping
    color_vectPop_E  = '#76c5ad' # color for vector pop (eggs) variables
    color_vectPop_L  = '#1b9e77' # color for vector pop (larvae) variables
    color_vectPop_P  = '#0b3f30' # color for vector pop (pupae) variables
    color_denInc     = '#d95f02' # color for dengue variables
    hue_colors = {'pop_A(E)_frac': color_vectPop_E, 'pop_A(L)_frac': color_vectPop_L, 'pop_A(P)_frac': color_vectPop_P, 'dengueInc': color_denInc}

    # Plot the split violin plots
    violin_width = 0.90
    sns.violinplot(x='Month', y='Value', hue='Type', data=df_melted, palette=hue_colors, split=False, inner=None, ax=ax, 
                   density_norm='width', cut=0, gap=0.2, width=violin_width, legend=False, linewidth=0)

    # Set transparency for the violin plots
    for violin in ax.collections:
        violin.set_alpha(0.5)

    # Plot specific data points with offsets
    offset = violin_width / 8  # Adjust this value as needed
    ax.scatter(df['Month'].cat.codes - 3*offset, df[popFrac_A_E_colName], color='black', label='pop_A_E_frac', alpha=0.7, zorder=3, marker='.', s=16) # points to the left
    ax.scatter(df['Month'].cat.codes - 1*offset, df[popFrac_A_L_colName], color='black', label='pop_A_L_frac', alpha=0.7, zorder=3, marker='.', s=16) # points slightly to left
    ax.scatter(df['Month'].cat.codes + 1*offset, df[popFrac_A_P_colName], color='black', label='pop_A_P_frac', alpha=0.7, zorder=3, marker='.', s=16) # points slightly to right
    ax.scatter(df['Month'].cat.codes + 3*offset, df[denInc_colName], color='black', label='dengueInc', alpha=0.7, zorder=3, marker='.', s=16)     # points to right

    # Edit axes
    ax.set_xlabel('') # remove x axis title
    ax.set_ylabel('Norm. adult population', fontsize=16, color=color_vectPop_L) # edit y axis title
    ax.minorticks_on() # enable minor tick marks
    ax.set_xticks(range(len(df['Month'].unique())))  # set positions for x-axis ticks
    ax.set_xticklabels(df['Month'].cat.categories, fontsize=14)  # edit x-axis tick labels
    ax.set_yticks(np.arange(ymin_popFrac, ymax_popFrac + 0.01, ystep_popFrac))  # ensure y-ticks are set explicitly
    ax.set_yticklabels(ax.get_yticks(), fontsize=14)  # set y-axis tick labels
    ax.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, left=True, right=False) # set minor ticks only for the y-axis
    ax.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, left=True, right=False) # adjust major ticks for y-axis
    ax.tick_params(axis='y', labelcolor=color_vectPop_L)  # Set tick label color for y-axis
    ax.tick_params(axis='x', which='both', bottom=False, top=False)  # remove x axis ticks
    ax.set_ylim(bottom=ymin_popFrac, top=ymax_popFrac)   # set y-axis limits
    ax.set_xlim(left=-0.5, right=11.5)   # set x-axis limits so that it begins at January (index 0) and ends at December (index 11)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f')) # set y-axis tick label formatting

    # Create a secondary axis on the right
    ax_right = ax.twinx()
    ax_right.set_ylabel('Dengue incidence (per 100,000)', fontsize=16, color=color_denInc) # edit y axis title
    ax_right.set_yticks(np.arange(ymin_denInc, ymax_denInc + 0.01, step=ystep_denInc))  # ensure y-ticks are set explicitly
    ax_right.set_yticklabels(ax_right.get_yticks(), fontsize=14)  # set y-axis tick labels
    ax_right.minorticks_on() # enable minor tick marks
    ax_right.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, left=False, right=True) # set minor ticks only for the y-axis
    ax_right.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, left=False, right=True) # adjust major ticks for y-axis
    ax_right.tick_params(axis='y', which='both', labelcolor=color_denInc)     # Set tick label color for y-axis
    ax_right.set_ylim(bottom=ymin_denInc, top=ymax_denInc)   # set y-axis limits
    ax_right.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    # Add horizontal rectangles to mark seasons
    season_rects = [
        {'months': (0, 1),   'label': 'NEM', 'color': '#9e3138'}, # Jan-Feb
        {'months': (2, 3),   'label': 'FIM', 'color': '#d47136'}, # Mar-Apr
        {'months': (4, 8),   'label': 'SWM', 'color': '#41ab5d'}, # May-Sep
        {'months': (9, 10),  'label': 'SIM', 'color': '#6786ad'}, # Oct-Nov
        {'months': (11, 11), 'label': 'NEM', 'color': '#9e3138'}  # Dec
    ]
    
    rect_height = 0.05 # fraction of plot height
    y_pos = 1.0        # place at top of plot
    for season in season_rects:
        # Add the rectangle
        rect = patches.Rectangle(xy=(season['months'][0] / 12, y_pos - rect_height),  # Scale x to normalized coords
                                 width=(season['months'][1] - season['months'][0] + 1) / 12,  # Adjust width for months (normalized)
                                 height=rect_height, linewidth=1, edgecolor='none', facecolor=season['color'], alpha=1, zorder=5, 
                                 transform=ax.transAxes)  # transform needed for text positioning/sizing
        ax.add_patch(rect)
        
        # Add label for the season
        season_label_x = (season['months'][1] + season['months'][0] + 1) / 24  # Scale x to normalized coords
        season_label_y = y_pos - (rect_height / 2)  # Center the text in the rectangle (in normalized y)
        ax.text(x=season_label_x, y=season_label_y, s=season['label'], 
                color='black', ha='center', va='center', fontsize=10, zorder=6, fontweight='bold', 
                transform=ax.transAxes)  # transform needed for text positioning/sizing

    # Add vertical lines between months
    num_months = len(df['Month'].unique())
    for i in range(1, num_months):
        ax.axvline(x=i-0.5, color='grey', linestyle='--', lw=0.8, alpha=0.7, zorder=0)  # Add a vertical line between month ticks

    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)




def plot_violinplt_allVars(df_in, loc, fileName):
    dir_out = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/'  # TEMPORARY!

    # Create figure with multiple vertically stacked subplots, sharing the x-axis
    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(20, 24), sharex=True, gridspec_kw={'height_ratios': [1, 2, 1, 2, 2], 'hspace': 0})

    df = df_in.copy(deep=True)

    # Define column names
    prcp_colName = 'PRCP (mm)'
    WH_colName = 'WH (mm)'
    TA_max_colName = 'TA_max (C)'
    TA_min_colName = 'TA_min (C)'
    TA_mean_colName = 'TA (C)'
    TW_max_colName = 'TW_max (C)'
    TW_min_colName = 'TW_min (C)'
    TW_mean_colName = 'TW (C)'
    dev_E_colName = 'dev_E'
    dev_L_colName = 'dev_L'
    dev_P_colName = 'dev_P'
    surv_E_colName = 'surv_E'
    surv_L_colName = 'surv_L'
    surv_P_colName = 'surv_P'
    popFrac_A_E_colName = 'pop_A(E)_frac'
    popFrac_A_L_colName = 'pop_A(L)_frac'
    popFrac_A_P_colName = 'pop_A(P)_frac'
    denInc_colName = 'dengueInc'

    # Set y-axis bounds based on loc
    ymin_prcp, ymax_prcp, ystep_prcp = 0, 640, 200
    ymin_WH, ymax_WH, ystep_WH       = 0, ymax_prcp/2, ystep_prcp/2
    ymin_popFrac, ymax_popFrac, ystep_popFrac = 0, 0.90, 0.2
    ymin_surv, ymax_surv, ystep_surv = 0, 1.05, 0.2
    if loc == 'Negombo':
        ymin_temp, ymax_temp, ystep_temp = 18, 44, 5
        ymin_temp_tick, ymax_temp_tick = 20, 40
        ymin_dev, ymax_dev, ystep_dev = 0, 0.79, 0.2
        ymin_denInc, ymax_denInc, ystep_denInc = 0, 49, 10
    elif loc == 'Jaffna':
        ymin_temp, ymax_temp, ystep_temp = 20, 49, 5
        ymin_temp_tick, ymax_temp_tick = 25, 45
        ymin_dev, ymax_dev, ystep_dev = 0, 1.12, 0.2
        ymin_denInc, ymax_denInc, ystep_denInc = 0, 179, 50
    elif loc == 'NuwaraEliya':
        ymin_temp, ymax_temp, ystep_temp = 4, 36, 5
        ymin_temp_tick, ymax_temp_tick = 5, 30
        ymin_dev, ymax_dev, ystep_dev = 0, 0.37, 0.1
        ymin_denInc, ymax_denInc, ystep_denInc = 0, 13, 5

    # Extract month from the index to use as labels
    df['Month'] = df.index.strftime('%b')
    df['Month'] = pd.Categorical(df['Month'], categories=df['Month'].unique(), ordered=True)

    # Rescale variables that we plot on the left axis, but are meant to be read off the right axis
    df[WH_colName] =df[WH_colName]*((ymax_prcp-ymin_prcp)/(ymax_WH-ymin_WH))
    df[denInc_colName] =df[denInc_colName]*((ymax_popFrac-ymin_popFrac)/(ymax_denInc-ymin_denInc))

    # Melt the DataFrames
    df_melted_prcp_WH = df.melt(id_vars='Month', value_vars=[prcp_colName, WH_colName], var_name='Variable', value_name='Value')
    df_melted_prcp_WH['Type'] = df_melted_prcp_WH['Variable'].apply(lambda x: 'PRCP' if 'PRCP' in x else 'WH')
    df_melted_max = df.melt(id_vars='Month', value_vars=[TA_max_colName, TW_max_colName], var_name='Variable', value_name='Value')
    df_melted_max['Type'] = df_melted_max['Variable'].apply(lambda x: 'TA' if 'TA' in x else 'TW')
    df_melted_min = df.melt(id_vars='Month', value_vars=[TA_min_colName, TW_min_colName], var_name='Variable', value_name='Value')
    df_melted_min['Type'] = df_melted_min['Variable'].apply(lambda x: 'TA' if 'TA' in x else 'TW')
    df_melted_mean = df.melt(id_vars='Month', value_vars=[TA_mean_colName, TW_mean_colName], var_name='Variable', value_name='Value')
    df_melted_mean['Type'] = df_melted_mean['Variable'].apply(lambda x: 'TA' if 'TA' in x else 'TW')
    df_melted_dev = df.melt(id_vars='Month', value_vars=[dev_E_colName, dev_L_colName, dev_P_colName], var_name='Variable', value_name='Value')
    df_melted_dev['Type'] = df_melted_dev['Variable'].apply(lambda x: 'dev_E' if 'dev_E' in x else ('dev_L' if 'dev_L' in x else 'dev_P'))
    df_melted_surv = df.melt(id_vars='Month', value_vars=[surv_E_colName, surv_L_colName, surv_P_colName], var_name='Variable', value_name='Value')
    df_melted_surv['Type'] = df_melted_surv['Variable'].apply(lambda x: 'surv_E' if 'surv_E' in x else ('surv_L' if 'surv_L' in x else 'surv_P'))
    df_melted_popFrac_denInc = df.melt(id_vars='Month', value_vars=[popFrac_A_E_colName, popFrac_A_L_colName, popFrac_A_P_colName, denInc_colName], var_name='Variable', value_name='Value')
    df_melted_popFrac_denInc['Type'] = df_melted_popFrac_denInc['Variable'].apply(lambda x: 'pop_A(E)_frac' if 'pop_A(E)_frac' in x else 
                                                                                  ('pop_A(L)_frac' if 'pop_A(L)_frac' in x else 
                                                                                   ('pop_A(P)_frac' if 'pop_A(P)_frac' in x else 'dengueInc')))

    # Define color mappings
    color_met  = '#48c335' # color for meteorology variables
    color_cont = '#2677ff' # color for container water dynamics variables
    hue_colors_prcp_WH = {'PRCP': color_met, 'WH': color_cont}
    hue_colors_TA_TW = {'TA': color_met, 'TW': color_cont}
    color_vectBio_E = '#f394c5' # color for vector biology (eggs) variables
    color_vectBio_L = '#b9216e' # color for vector biology (larvae) variables
    color_vectBio_P = '#450c29' # color for vector biology (pupae) variables
    hue_colors_surv = {'surv_E': color_vectBio_E, 'surv_L': color_vectBio_L, 'surv_P': color_vectBio_P}
    hue_colors_dev  = {'dev_E': color_vectBio_E, 'dev_L': color_vectBio_L, 'dev_P': color_vectBio_P}
    color_vectPop_E  = '#76c5ad' # color for vector pop (eggs) variables
    color_vectPop_L  = '#1b9e77' # color for vector pop (larvae) variables
    color_vectPop_P  = '#0b3f30' # color for vector pop (pupae) variables
    color_denInc     = '#d95f02' # color for dengue variables
    hue_colors_popFrac_denInc = {'pop_A(E)_frac': color_vectPop_E, 'pop_A(L)_frac': color_vectPop_L, 'pop_A(P)_frac': color_vectPop_P, 'dengueInc': color_denInc}

    # prcp/WH plot (ax1)
    violin_width = 0.85
    offset = violin_width / 4
    sns.violinplot(x='Month', y='Value', hue='Type', data=df_melted_prcp_WH, palette=hue_colors_prcp_WH, split=False, inner=None, ax=ax1, 
                   density_norm='width', cut=0, gap=0.2, width=violin_width, legend=False, linewidth=0)
    ax1.scatter(df['Month'].cat.codes - offset, df[prcp_colName], color='black', label='prcp', alpha=0.7, zorder=3, marker='.', s=16) # prcp points slightly to the left
    ax1.scatter(df['Month'].cat.codes + offset, df[WH_colName], color='black', label='WH', alpha=0.7, zorder=3, marker='.', s=16)     # WH points slightly to the right

    ax1.set_ylabel('Precipitation (mm)', fontsize=16, color=color_met)
    ax1.set_ylim(ymin_prcp, ymax_prcp)
    ax1.set_yticks(np.arange(ymin_prcp, ymax_prcp + 0.01, ystep_prcp))
    ax1.set_yticklabels(ax1.get_yticks(), fontsize=14)  # set y-axis tick label size
    ax1.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, left=True, right=False) # adjust minor ticks for y-axis
    ax1.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, left=True, right=False) # adjust major ticks for y-axis
    ax1.tick_params(axis='y', which='both', labelcolor=color_met)  # Set tick label color for y-axis
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    ax1_right = ax1.twinx()
    ax1_right.set_ylabel('Water height (mm)', fontsize=16, color=color_cont, rotation=270, labelpad=20) # edit y axis title
    ax1_right.set_ylim(ymin_WH, ymax_WH)
    ax1_right.set_yticks(np.arange(ymin_WH, ymax_WH + 0.01, step=ystep_WH))
    ax1_right.set_yticklabels(ax1_right.get_yticks(), fontsize=14)  # set y-axis tick label size
    ax1_right.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, left=False, right=True) # # adjust minor ticks for y-axis
    ax1_right.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, left=False, right=True) # adjust major ticks for y-axis
    ax1_right.tick_params(axis='y', which='both', labelcolor=color_cont)     # Set tick label color for y-axis
    ax1_right.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    # TA/TW plot (ax2)
    violin_width = 0.85
    offset = violin_width / 4
    sns.violinplot(x='Month', y='Value', hue='Type', data=df_melted_max, palette=hue_colors_TA_TW, split=False, inner=None, ax=ax2, 
                   density_norm='width', cut=0, gap=0.2, width=violin_width, legend=False, linewidth=0)
    sns.violinplot(x='Month', y='Value', hue='Type', data=df_melted_min, palette=hue_colors_TA_TW, split=False, inner=None, ax=ax2, 
                   density_norm='width', cut=0, gap=0.2, width=violin_width, legend=False, linewidth=0)
    sns.violinplot(x='Month', y='Value', hue='Type', data=df_melted_mean, palette=hue_colors_TA_TW, split=False, inner=None, ax=ax2, 
                   density_norm='width', cut=0, gap=0.2, width=violin_width, legend=False, linewidth=0)
    ax2.scatter(df['Month'].cat.codes - offset, df[TA_max_colName], color='#900000', label='TA Max', alpha=0.7, zorder=3, marker='.', s=16)  # TA points slightly to the left
    ax2.scatter(df['Month'].cat.codes + offset, df[TW_max_colName], color='#900000', label='TW Max', alpha=0.7, zorder=3, marker='.', s=16)  # TW points slightly to the right
    ax2.scatter(df['Month'].cat.codes - offset, df[TA_min_colName], color='blue', label='TA Min', alpha=0.7, zorder=3, marker='.', s=16)     # TA points slightly to the left
    ax2.scatter(df['Month'].cat.codes + offset, df[TW_min_colName], color='blue', label='TW Min', alpha=0.7, zorder=3, marker='.', s=16)     # TW points slightly to the right
    ax2.scatter(df['Month'].cat.codes - offset, df[TA_mean_colName], color='black', label='TA Mean', alpha=0.7, zorder=3, marker='.', s=16)  # TA points slightly to the left
    ax2.scatter(df['Month'].cat.codes + offset, df[TW_mean_colName], color='black', label='TW Mean', alpha=0.7, zorder=3, marker='.', s=16)  # TW points slightly to the right

    ax2.set_ylabel('Air temperature (°C)', fontsize=16, color=color_met)
    ax2.set_ylim(ymin_temp, ymax_temp)
    ax2.set_yticks(np.arange(ymin_temp_tick, ymax_temp_tick + 0.01, ystep_temp))
    ax2.set_yticklabels(ax2.get_yticks(), fontsize=14)  # set y-axis tick label size
    ax2.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, left=True, right=False) # adjust minor ticks for y-axis
    ax2.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, left=True, right=False) # adjust major ticks for y-axis
    ax2.tick_params(axis='y', which='both', labelcolor=color_met)  # Set tick label color for y-axis
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    ax2_right = ax2.twinx()
    ax2_right.set_ylabel('Water temperature (°C)', fontsize=16, color=color_cont, rotation=270, labelpad=20) # edit y axis title
    ax2_right.set_ylim(ymin_temp, ymax_temp)
    ax2_right.set_yticks(np.arange(ymin_temp_tick, ymax_temp_tick + 0.01, step=ystep_temp))
    ax2_right.set_yticklabels(ax2_right.get_yticks(), fontsize=14)  # set y-axis tick label size
    ax2_right.set_yticklabels(ax2_right.get_yticks(), fontsize=14)  # set y-axis tick labels
    ax2_right.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, left=False, right=True) # # adjust minor ticks for y-axis
    ax2_right.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, left=False, right=True) # adjust major ticks for y-axis
    ax2_right.tick_params(axis='y', which='both', labelcolor=color_cont)     # Set tick label color for y-axis
    ax2_right.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    # Dev rate plot (ax3)
    violin_width = 0.90
    offset = violin_width / 3
    sns.violinplot(x='Month', y='Value', hue='Type', data=df_melted_dev, palette=hue_colors_dev, split=False, inner=None, ax=ax3, 
                   density_norm='width', cut=0, gap=0.2, width=violin_width, legend=False, linewidth=0)
    ax3.scatter(df['Month'].cat.codes - offset, df[dev_E_colName], color='black', alpha=0.7, zorder=3, marker='.', s=16)
    ax3.scatter(df['Month'].cat.codes,          df[dev_L_colName], color='black', alpha=0.7, zorder=3, marker='.', s=16)
    ax3.scatter(df['Month'].cat.codes + offset, df[dev_P_colName], color='black', alpha=0.7, zorder=3, marker='.', s=16)

    ax3.set_ylabel(r'Dev. rate (day$^{-1}$)', fontsize=16, color=color_vectBio_L)
    ax3.set_ylim(ymin_dev, ymax_dev)
    ax3.set_yticks(np.round(np.arange(ymin_dev, ymax_dev + 0.01, ystep_dev),1))
    ax3.set_yticklabels(ax3.get_yticks(), fontsize=14)  # set y-axis tick label size
    ax3.tick_params(axis='y', which='both', labelcolor=color_vectBio_L)     # Set tick label color for y-axis
    ax3.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, left=True, right=True) # adjust minor ticks for y-axis
    ax3.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, left=True, right=True) # adjust major ticks for y-axis
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    # Survival rate plot (ax4)
    violin_width = 0.90
    offset = violin_width / 3
    sns.violinplot(x='Month', y='Value', hue='Type', data=df_melted_surv, palette=hue_colors_surv, split=False, inner=None, ax=ax4, 
                   density_norm='width', cut=0, gap=0.2, width=violin_width, legend=False, linewidth=0)
    ax4.scatter(df['Month'].cat.codes - offset, df[surv_E_colName], color='black', alpha=0.7, zorder=3, marker='.', s=16)
    ax4.scatter(df['Month'].cat.codes,          df[surv_L_colName], color='black', alpha=0.7, zorder=3, marker='.', s=16)
    ax4.scatter(df['Month'].cat.codes + offset, df[surv_P_colName], color='black', alpha=0.7, zorder=3, marker='.', s=16)

    ax4.set_ylabel(r'Survival rate (day$^{-1}$)', fontsize=16, color=color_vectBio_L)
    ax4.set_ylim(ymin_surv, ymax_surv)
    ax4.set_yticks(np.round(np.arange(ymin_surv, ymax_surv + 0.01, ystep_surv),1))
    ax4.set_yticklabels(ax4.get_yticks(), fontsize=14)  # set y-axis tick label size
    ax4.tick_params(axis='y', which='both', labelcolor=color_vectBio_L)     # Set tick label color for y-axis
    ax4.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, left=True, right=True) # adjust minor ticks for y-axis
    ax4.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, left=True, right=True) # adjust major ticks for y-axis
    ax4.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    # popFrac/denInc plot (ax5)
    violin_width = 0.90
    offset = violin_width / 8  # Adjust this value as needed
    sns.violinplot(x='Month', y='Value', hue='Type', data=df_melted_popFrac_denInc, palette=hue_colors_popFrac_denInc, split=False, inner=None, ax=ax5, 
                   density_norm='width', cut=0, gap=0.2, width=violin_width, legend=False, linewidth=0)
    ax5.scatter(df['Month'].cat.codes - 3*offset, df[popFrac_A_E_colName], color='black', label='pop_A_E_frac', alpha=0.7, zorder=3, marker='.', s=16) # points to the left
    ax5.scatter(df['Month'].cat.codes - 1*offset, df[popFrac_A_L_colName], color='black', label='pop_A_L_frac', alpha=0.7, zorder=3, marker='.', s=16) # points slightly to left
    ax5.scatter(df['Month'].cat.codes + 1*offset, df[popFrac_A_P_colName], color='black', label='pop_A_P_frac', alpha=0.7, zorder=3, marker='.', s=16) # points slightly to right
    ax5.scatter(df['Month'].cat.codes + 3*offset, df[denInc_colName], color='black', label='dengueInc', alpha=0.7, zorder=3, marker='.', s=16)         # points to right

    ax5.set_ylabel('Adult population (norm.)', fontsize=16, color=color_vectPop_L)
    ax5.set_ylim(ymin_popFrac, ymax_popFrac)
    ax5.set_yticks(np.round(np.arange(ymin_popFrac, ymax_popFrac + 0.01, ystep_popFrac),1))
    ax5.set_yticklabels(ax5.get_yticks(), fontsize=14)  # set y-axis tick label size
    ax5.tick_params(axis='y', which='both', labelcolor=color_vectPop_L)     # Set tick label color for y-axis
    ax5.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, left=True, right=False) # adjust minor ticks for y-axis
    ax5.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, left=True, right=False) # adjust major ticks for y-axis
    ax5.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    ax5_right = ax5.twinx()
    ax5_right.set_ylabel('Dengue incidence (per 100,000)', fontsize=16, color=color_denInc, rotation=270, labelpad=20) # edit y axis title
    ax5_right.set_ylim(ymin_denInc, ymax_denInc)   # set y-axis limits
    ax5_right.set_yticks(np.arange(ymin_denInc, ymax_denInc + 0.01, step=ystep_denInc))
    ax5_right.set_yticklabels(ax5_right.get_yticks(), fontsize=14)  # set y-axis tick labels
    ax5_right.tick_params(axis='y', which='both', labelcolor=color_denInc)     # Set tick label color for y-axis
    ax5_right.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, left=False, right=True) # # adjust minor ticks for y-axis
    ax5_right.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, left=False, right=True) # adjust major ticks for y-axis
    ax5_right.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    # Set transparency of all violin plots
    for ax in [ax1, ax2, ax3, ax4, ax5]:
        for violin in ax.collections:
            violin.set_alpha(0.8)

    # Format all x-axes (months)
    ax1.xaxis.set_ticks_position('top')  # Move x-axis ticks to the top
    ax1.xaxis.set_label_position('top')   # Move x-axis label to the top
    for ax in [ax1, ax1_right, ax2, ax2_right, ax3, ax4, ax5, ax5_right]:
        ax.set_xlabel('') # remove x axis title
        ax.minorticks_on() # enable minor tick marks
        ax.tick_params(axis='x', which='both', bottom=False, top=False)  # remove x axis ticks
        ax.set_xlim(left=-0.5, right=11.5)   # set x-axis limits so that it begins at January (index 0) and ends at December (index 11)
    ax1.set_xticks(range(len(df['Month'].unique())))
    ax1.set_xticklabels(df['Month'].cat.categories, fontsize=18)
    ax5.set_xticks(range(len(df['Month'].unique())))
    ax5.set_xticklabels(df['Month'].cat.categories, fontsize=18)

    # Add season patches to top plot
    season_rects = [
        {'months': (0, 1),   'label': 'NEM', 'color': '#9e3138'},
        {'months': (2, 3),   'label': 'FIM', 'color': '#d47136'},
        {'months': (4, 8),   'label': 'SWM', 'color': '#41ab5d'},
        {'months': (9, 10),  'label': 'SIM', 'color': '#6786ad'},
        {'months': (11, 11), 'label': 'NEM', 'color': '#9e3138'}
    ]
    
    bbox = ax1.get_position()                 # returns bounding box of axis (x0, y0, width, height) in figure coordinates
    rect_height = bbox.height*0.20            # height of rectangles (x% of the axis height)
    rect_pos_y = bbox.y1 + bbox.height*0.20   # set y-position above the top of the axes
    for season in season_rects:
        rect_pos_x = bbox.x0 + bbox.width*(season['months'][0] / 12)
        rect_width = bbox.width*((season['months'][1] - season['months'][0] + 1) / 12)  # adjust width for months (normalized)
        rect = patches.Rectangle(xy=(rect_pos_x, rect_pos_y), width=rect_width, height=rect_height, 
                                 linewidth=0, edgecolor='none', facecolor=season['color'], alpha=0.8, zorder=5, 
                                 transform=ax1.figure.transFigure)  # Use figure coordinates
        ax1.figure.patches.append(rect)  # Add the rectangle to the figure

        season_label_x = bbox.x0 + bbox.width*((season['months'][1] + season['months'][0] + 1) / 24)  # Scale x to normalized coords
        season_label_y = rect_pos_y + (rect_height / 2)                                               # Center the text in the rectangle (in normalized y)
        ax1.figure.text(x=season_label_x, y=season_label_y, s=season['label'], 
                        color='black', ha='center', va='center', fontsize=16, zorder=6, fontweight='bold', 
                        transform=ax1.figure.transFigure)  # Use figure coordinates for text

    # Add vertical lines between months
    num_months = len(df['Month'].unique())
    for ax in [ax1, ax2, ax3, ax4, ax5]:
        for i in range(1, num_months):
            ax.axvline(x=i-0.5, color='grey', linestyle='--', lw=0.8, alpha=0.7, zorder=0)  # Add a vertical line between month ticks

    # Save the figure
    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)



##################################################
### TERCILE-TERCILE PLOTS ########################



# Plot stacked bar plot of terciles, showing the percentage of each tercile of dengue incidence that corresponds
# to each tercile of the model variable.
def plot_bar_varX_denInc(df_in, varX_colName, varX_pltName, denInc_colName, denInc_pltName, fileName):
    dir_out = "/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/" # TEMPORARY!

    df_monthly = df_in.copy(deep=True)  # create a temporary copy of the df so that we don't alter the original

    # Step 1: Categorize data into Low, Medium, High, on a monthly basis (and then recombine all monthly results into one series)
    varX_terciles   = categorizeSeries_LowMedHigh_byMon(df_monthly[varX_colName])    
    denInc_terciles = categorizeSeries_LowMedHigh_byMon(df_monthly[denInc_colName])

    # If categorization failed for either variable, we're not going to make a plot.
    if varX_terciles.isnull().all():
       print(f'Error: categorization failed for variable {varX_colName}. Skipping plot generation.')
       return
    if denInc_terciles.isnull().all():
       print(f'Error: categorization failed for variable {denInc_colName}. Skipping plot generation.')
       return

    # Add tercile columns to the DataFrame
    varX_tercile_colName   = varX_colName + '_tercile'
    denInc_tercile_colName = denInc_colName  + '_tercile'
    df_monthly[varX_tercile_colName]   = varX_terciles
    df_monthly[denInc_tercile_colName] = denInc_terciles

    # Step 2: Create a crosstab to calculate the counts of denInc terciles within each varX tercile
    crosstab = pd.crosstab(df_monthly[varX_tercile_colName], df_monthly[denInc_tercile_colName])

    # Step 3: Normalize the crosstab to get percentages
    crosstab_norm = crosstab.div(crosstab.sum(axis=1), axis=0) * 100

    # Step 4: Perform chi-square test of independence
    chi2_stat, chi2_p_val, _, _ = chi2_contingency(crosstab)
    print(" ")
    print("Chi-square test")
    print("*******************************")
    print("  Chi-Square Statistic:", chi2_stat)
    print("  p-value:", chi2_p_val)

    # Step 5: Calculate Kendall's tau-b (direction of association; -1 to 1)
    kendall_tau_b, kendall_p_val = kendalltau(df_monthly[varX_tercile_colName], df_monthly[denInc_tercile_colName])
    print(" ")
    print("Kendall's tau-b test")
    print("*******************************")
    print("  Kendall's tau-b:", kendall_tau_b)
    print("  p-value:", kendall_p_val)
    print(" ")


    # Step 6: Plot the stacked bar chart
    # Initialize plot
    fig, ax = plt.subplots(figsize=(8, 8))
    
    tercColors = ['#edf8b1', '#7fcdbb', '#2c7fb8'] # colors for terciles
    bar_width = 0.80
    bars = crosstab_norm.plot(kind='bar', stacked=True, width=bar_width, color=tercColors, ax=ax, zorder=2)

    # Adding labels with percents/counts on each section of the stacked bars
    for i, rect in enumerate(bars.containers):
        for j, bar in enumerate(rect):
            # Calculate the percents/counts for this section
            quantile_percent = round(crosstab_norm.iloc[j, i])
            quantile_count = crosstab.iloc[j, i]

            if quantile_count > 0: # only do labeling if this category is not empty
                # Get x and y placement of the label based on rectangle location
                x_value = bar.get_x() + bar.get_width() / 2
                y_value = bar.get_y() + bar.get_height() / 2

                # Format and add label
                label = f'{quantile_percent}%\n({quantile_count})'
                ax.text(x_value, y_value, label, ha='center', va='center', fontsize=13, fontweight='bold', color='black', zorder=4)

    # Add text indicating statistical test results outside the plot
    signif_level = 0.05
    stat_results_text = (r"$\bf{\chi^2\ test:}$" + " " +
                         r'$\chi^2$' + f" = {chi2_stat:.1f}, " + 
                         "p = " + (r"$\bf{" if chi2_p_val < signif_level else "") + f"{chi2_p_val:.2e}" + (r"*}$" if chi2_p_val < signif_level else "") + "\n" +
                         r"$\bf{Kendall's\ \tau_B\ test:}$" + " " +
                         r'$\tau_B$' + f" = {kendall_tau_b:.2f}, " +
                         "p = " + (r"$\bf{" if kendall_p_val < signif_level else "") + f"{kendall_p_val:.2e}" + (r"*}$" if kendall_p_val < signif_level else ""))
    ax.text(0.5, 1.01, stat_results_text, transform=ax.transAxes, fontsize=14, va='bottom', ha='center')

    # Adding labels and title
    ax.set_xticklabels(ax.get_xticklabels(), ha='center', fontsize=15, fontweight='bold', rotation=0)
    ax.tick_params(axis='y', labelsize=14)
    ax.set_xlabel('')
    ax.set_ylabel('Frequency of association with\ndengue incidence tercile (%)', fontsize=16)
    for spine in ax.spines.values():
        spine.set_zorder(5)  # Set zorder for all spines to 5 (to be above bars)
    ax.grid(axis='y', linestyle='--', alpha=0.7, zorder=0)
    ax.set_ylim(0, 105)
    ax.legend(handles=[], frameon=False) # remove legend

    # Adjust layout and save figure
    plt.tight_layout()
    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)



##################################################
### CROSSTAB SCATTER PLOTS #######################
    


# plot with variables as the x-axis categories (with multiple vars in same x-axis category)
def plot_scatter_crosstabs(df_in, loc, tercileName_modelVar, tercileName_den, fileName):
    dir_out = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/5_plotting/'  # TEMPORARY!

    fig, ax = plt.subplots(figsize=(27, 9))  # Create a new figure and axis for the plot

    # Define y-axis limits (tercile-tercile assoc. %)
    ymin = 0
    ymax = 100

    # Define color mapping
    color_met  = '#48c335'       # color for meteorology variables
    color_cont = '#2677ff'       # color for container water dynamics variables
    color_vectBio_E = '#f394c5'  # color for vector biology (eggs) variables
    color_vectBio_L = '#b9216e'  # color for vector biology (larvae) variables
    color_vectBio_P = '#450c29'  # color for vector biology (pupae) variables
    color_vectPop_E  = '#76c5ad' # color for vector pop (eggs) variables
    color_vectPop_L  = '#1b9e77' # color for vector pop (larvae) variables
    color_vectPop_P  = '#0b3f30' # color for vector pop (pupae) variables

    color_map = {'TA_min (C)': color_met,
                 'TA (C)': color_met,
                 'TA_max (C)': color_met,
                 'PRCP (mm)': color_met,
                 'TW_min_nanFill (C)': color_cont,
                 'TW_nanFill (C)': color_cont,
                 'TW_max_nanFill (C)': color_cont,
                 'WH (mm)': color_cont,
                 'dev_E': color_vectBio_E,
                 'surv_E': color_vectBio_E,
                 'dev_L': color_vectBio_L,
                 'surv_L': color_vectBio_L,
                 'dev_P': color_vectBio_P,
                 'surv_P': color_vectBio_P,
                 'pop_A(E)': color_vectPop_E,
                 'pop_A(L)': color_vectPop_L,
                 'pop_A(P)': color_vectPop_P
                }

    # Define datapoint offset mapping
    varOffset_left  = -1.0
    varOffset_mid   =  0.0
    varOffset_right =  1.0
    varOffset_map = {
        'TA_min (C)':         varOffset_left,
        'TA (C)':             varOffset_mid,
        'TA_max (C)':         varOffset_right,
        'PRCP (mm)':          varOffset_mid,
        'TW_min_nanFill (C)': varOffset_left,
        'TW_nanFill (C)':     varOffset_mid,
        'TW_max_nanFill (C)': varOffset_right,
        'WH (mm)':            varOffset_mid,
        'dev_E':              varOffset_left,
        'dev_L':              varOffset_mid,
        'dev_P':              varOffset_right,
        'surv_E':             varOffset_left,
        'surv_L':             varOffset_mid,
        'surv_P':             varOffset_right,
        'pop_A(E)':           varOffset_left,
        'pop_A(L)':           varOffset_mid,
        'pop_A(P)':           varOffset_right
    }

    # Define how we group the variables
    groupWidth_wide   = 3.2
    groupWidth_narrow = 1.2
    varGroup_params = {
        'TA_group'  : {'groupLabel': 'Air temperature',   'groupWidth': groupWidth_wide, 'groupVars': ['TA_min (C)', 'TA (C)', 'TA_max (C)']},
        'PRCP_group': {'groupLabel': 'Precip.',           'groupWidth': groupWidth_narrow, 'groupVars': ['PRCP (mm)']},
        'TW_group'  : {'groupLabel': 'Water temperature', 'groupWidth': groupWidth_wide, 'groupVars': ['TW_min_nanFill (C)', 'TW_nanFill (C)', 'TW_max_nanFill (C)']},
        'WH_group':   {'groupLabel': 'Water height',      'groupWidth': groupWidth_narrow, 'groupVars': ['WH (mm)']},
        'dev_group':  {'groupLabel': 'Development rate',  'groupWidth': groupWidth_wide, 'groupVars': ['dev_E', 'dev_L', 'dev_P']},
        'surv_group': {'groupLabel': 'Survival rate',     'groupWidth': groupWidth_wide, 'groupVars': ['surv_E', 'surv_L', 'surv_P']},
        'pop_group':  {'groupLabel': 'Adult population',  'groupWidth': groupWidth_wide, 'groupVars': ['pop_A(E)', 'pop_A(L)', 'pop_A(P)']}
    }

    # Custom marker style map based on column names
    fill_noOutbreaks = False
    fill_inclOutbreaks = True
    txtColor_noOutbreaks = 'black'
    txtColor_inclOutbreaks = 'white'
    offset_lag0 = -0.25
    offset_lag1 = 0
    offset_lag2 = 0.25
    marker_map = {
        'lag0_noOutbreaks': {'marker': 'D', 'fill': fill_noOutbreaks,   'number': '0', 'txtColor': txtColor_noOutbreaks,   'offset': offset_lag0},
        'lag0':             {'marker': 'D', 'fill': fill_inclOutbreaks, 'number': '0', 'txtColor': txtColor_inclOutbreaks, 'offset': offset_lag0},
        'lag1_noOutbreaks': {'marker': 'D', 'fill': fill_noOutbreaks,   'number': '1', 'txtColor': txtColor_noOutbreaks,   'offset': offset_lag1},
        'lag1':             {'marker': 'D', 'fill': fill_inclOutbreaks, 'number': '1', 'txtColor': txtColor_inclOutbreaks, 'offset': offset_lag1},
        'lag2_noOutbreaks': {'marker': 'D', 'fill': fill_noOutbreaks,   'number': '2', 'txtColor': txtColor_noOutbreaks,   'offset': offset_lag2},
        'lag2':             {'marker': 'D', 'fill': fill_inclOutbreaks, 'number': '2', 'txtColor': txtColor_inclOutbreaks, 'offset': offset_lag2}
    }


    # Dynamically set the parameter 'groupPosX' based on the starting position (0) and each groupWidth
    previous_group_width = None
    previous_group_posX = 0      # initialize the first group's position as 0
    for group_key, group_params in varGroup_params.items(): # iterate over the groups and calculate 'groupPosX' based on 'groupWidth'
        if previous_group_width is None:  # First group
            group_params['groupPosX'] = 0
        else: # Calculate the new groupPosX
            group_params['groupPosX'] = previous_group_posX + (previous_group_width + group_params['groupWidth']) / 2
        # Update previous group values for the next iteration
        previous_group_width = group_params['groupWidth']
        previous_group_posX = group_params['groupPosX']

    # Function to get groupLabel, groupPosX, and groupWidth for a given variable
    def get_group_params(varName, varGroup_params):
        for group, params in varGroup_params.items():
            if varName in params['groupVars']:
                return params['groupLabel'], float(params['groupPosX']), float(params['groupWidth'])
        return None, None, None  # If the variable is not found


    # Plot each variable
    x_positions_unique = []
    x_groupWidths_unique = []
    x_labels_unique = []
    seen_labels = set()
    for i, variable in enumerate(df_in.index):

        x_label, groupPosX, groupWidth = get_group_params(variable, varGroup_params)

        # the Nuwara Eliya pop_A(E) data is mostly zeroes and doesn't map nicely to terciles, 
        # so we don't plot that data (and draw a gray box instead)
        if loc == 'NuwaraEliya' and variable == 'pop_A(E)': 
            grayBox_xPos_left = groupPosX-(0.5*groupWidth)
            grayBox_xPos_right = grayBox_xPos_left + (1/3)*groupWidth
            ax.fill_betweenx(y=[ymin, ymax], x1=grayBox_xPos_left, x2=grayBox_xPos_right, color='#888888', alpha=0.3)

        else:
            for j, column in enumerate(df_in.columns):
                
                # Get marker style based on column name
                marker_info = marker_map[column]
                fill_color = color_map[variable] if marker_info['fill'] else 'none'  # Use fill or no fill
                x_position = groupPosX + marker_info['offset'] + varOffset_map[variable]

                # Plot the marker
                ax.scatter(x_position, df_in[column].iloc[i], 
                        edgecolor=color_map[variable], facecolor=fill_color, linewidth=2,
                        marker=marker_info['marker'], s=140, label=column if i == 0 else '')

                # Add the corresponding number inside the marker
                if variable in ['TA_min (C)', 'TW_min_nanFill (C)']:
                    if column in ['lag0_noOutbreaks', 'lag1_noOutbreaks', 'lag2_noOutbreaks']:
                        txtColor = '#0b43bb'
                    elif column in ['lag0', 'lag1', 'lag2']:
                        txtColor = '#d6e4ff'
                elif variable in ['TA_max (C)', 'TW_max_nanFill (C)']:
                    if column in ['lag0_noOutbreaks', 'lag1_noOutbreaks', 'lag2_noOutbreaks']:
                        txtColor = '#900000'
                    elif column in ['lag0', 'lag1', 'lag2']:
                        txtColor = '#ffbaba'
                else:
                    txtColor = marker_info['txtColor']
                ax.text(x_position, df_in[column].iloc[i], marker_info['number'], 
                        fontsize=8, fontweight='bold', color=txtColor, ha='center', va='center')
                
        # Only add unique labels and positions
        if x_label not in seen_labels:
            x_positions_unique.append(groupPosX)
            x_groupWidths_unique.append(groupWidth)
            x_labels_unique.append(x_label)
            seen_labels.add(x_label)

    # Add vertical lines separating each category
    for i, pos in enumerate(x_positions_unique[:-1]):  # Omit the last position, no need for a line after the last category
        vertLine_xPos = pos + (x_groupWidths_unique[i]/2)
        ax.axvline(x=vertLine_xPos, color='gray', linewidth=1)
    # Add a horizontal line at y=33.33 (indicating expected Hi-Hi correspondence if the association was random)
    ax.axhline(y=100/3, color='gray', linestyle='--', linewidth=0.5)

    # Set axis ticks and labels
    ax.tick_params(axis='x', which='both', bottom=False, top=False)  # remove x axis ticks
    ax.set_xticks(x_positions_unique)  # only unique x-positions
    ax.set_xticklabels(x_labels_unique, fontsize=15)
    ax.tick_params(axis='y', labelsize=16)  # Set y-axis tick label font size to 12

    # Color x axis tick labels based on variable category
    colors_xTickLabel = [color_met, color_met, color_cont, color_cont, color_vectBio_L, color_vectBio_L, color_vectPop_L]  # color order must match plotting order!
    for tick_label, color in zip(ax.get_xticklabels(), colors_xTickLabel): # apply colors to each x-axis tick label
        tick_label.set_color(color)

    # Edit axes
    ax.set_xlabel('')
    ax.set_ylabel('Frequency of ' + tercileName_modelVar.lower() + ' tercile association with\n' + tercileName_den.lower() + ' dengue incidence (%)', fontsize=16)
    ax.set_xlim(-0.5*(varGroup_params['TA_group']['groupWidth']), max(x_positions_unique) + 0.5*(varGroup_params['pop_group']['groupWidth']))  # Add padding to the first and last positions to avoid off-centering
    ax.set_ylim(bottom=ymin, top=ymax)

    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)



##################################################








# Scripts that should be spun off into a separate file!!!!
###########################################################



# make timeseries plots of dengue data
if False:

    # Read in dengue data.
    fileDir_base = '/mnt/redwood/local_drive/chanud/RawData/epiReports/'
    yearStart_file = 2007
    yearEnd_file = 2022
    timespanStr_file = str(yearStart_file) + '_' + str(yearEnd_file)
    file_denInc_monthly  = fileDir_base + 'dengueInc_monthly_' + timespanStr_file + '.csv'
    df_denInc_monthly_allData  = pd.read_csv(file_denInc_monthly, index_col='Date', parse_dates=True)
    
    # Subset to time range of interest.
    yearStart_data = 2007
    yearEnd_data = 2020
    outbreakYears = [2017, 2019]
    timespanStr_data = str(yearStart_data) + '_' + str(yearEnd_data)
    df_denInc_monthly = df_denInc_monthly_allData[str(yearStart_data)+'-01-01' : str(yearEnd_data)+'-12-31'] # subset to years of interest

    fName = 'tseries_' + 'denInc_monthly_' + timespanStr_data
    plot_tseries_denInc(df_denInc_monthly, yearStart_data, yearEnd_data, outbreakYears, fName)




# make timeseries plots of model output by climate category (lo/hi temp x lo/hi prcp)
if False:

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
    varTemp_colName = 'TA (C)'
    varTemp_pltName = strConvert_paren2hyphen(strConvert_parenUnitsRemove(varTemp_colName))
    varPrcp_colName = 'PRCP (mm)'
    varPrcp_pltName = strConvert_paren2hyphen(strConvert_parenUnitsRemove(varPrcp_colName))
    var_colNames = ['TA (C)', 'TA_min (C)', 'TA_max (C)', 'PRCP (mm)', 'TW_nanFill (C)', 'TW_min_nanFill (C)', 'TW_max_nanFill (C)', 
                    'WH (mm)', 'dev_E', 'dev_L', 'dev_P', 'surv_E', 'surv_L', 'surv_P', 'pop_A(E)', 'pop_A(L)', 'pop_A(P)'] # subset of vars to actually plot
    var_pltNames = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in var_colNames]

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
            df_vectorPop_dly_yrFiltered = df_vectorPop_dly[~df_vectorPop_dly.index.year.isin(exclYrs)]
            
            for i_var, var_colName in enumerate(var_colNames):
                var_pltName = var_pltNames[i_var]
                print('  Var: ' + var_pltName)

                for month in [1]:#months:
                    fileName = 'tseries_climCat_' + varTemp_pltName + '-' + varPrcp_pltName + '_' + var_pltName + '_' + loc + fileStr_noOutbreaks + '_mon' + str(month).zfill(2)
                    plot_tseries_climCategories(df_vectorPop_dly_yrFiltered, month, var_colName, var_pltName, varTemp_colName, varTemp_pltName, varPrcp_colName, varPrcp_pltName, fileName)




# make BIG (all vars) timeseries plot of model output by climate category (lo/hi temp x lo/hi prcp)
if False:
    
    # parameters related to modeling pipeline
    dataSrc = 'MERRA_IMERG'
    locs = ['Negombo', 'NuwaraEliya']#['Negombo', 'NuwaraEliya', 'Jaffna']
    years = list(range(2001, 2020+1))
    months = list(range(1,12+1))
    intv = 'intv00'
    initDays = list(range(1,7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000
    varTemp_colName = 'TA (C)'
    varTemp_pltName = strConvert_paren2hyphen(strConvert_parenUnitsRemove(varTemp_colName))
    varPrcp_colName = 'PRCP (mm)'
    varPrcp_pltName = strConvert_paren2hyphen(strConvert_parenUnitsRemove(varPrcp_colName))

    # parameters related to dengue data
    noOutbreaks_opts = [False]#[True, False] # whether to exclude data from outbreak years
    outbreakYrs = [2017, 2019]


    for loc in locs:
        print(loc)

        # load model pipeline data at location of interest
        # also: convert from daily to monthly
        df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae)

        for noOutbreaks in noOutbreaks_opts:
            print(' noOutbreaks: ' + str(noOutbreaks))
            if noOutbreaks:
                fileStr_noOutbreaks = '_noOutbreaks'
                exclYrs = outbreakYrs  # outbreak years (to exclude from plot)
            else:
                fileStr_noOutbreaks = ''
                exclYrs = []

            # remove data in outbreak years we want to exclude
            df_vectorPop_dly_yrFiltered = df_vectorPop_dly[~df_vectorPop_dly.index.year.isin(exclYrs)]

            # manually create the popFrac columns and add them to vectorPop df
            popTotal_E = len(initDays) * initEggs # this code only works for eggs AND larvae and pupae as long as initEggs = initLarvae = initPupae
            df_vectorPop_dly_yrFiltered_withFrac = df_vectorPop_dly_yrFiltered.copy(deep=True)
            df_vectorPop_dly_yrFiltered_withFrac['pop_A(E)_frac'] = df_vectorPop_dly_yrFiltered['pop_A(E)'] / popTotal_E
            df_vectorPop_dly_yrFiltered_withFrac['pop_A(L)_frac'] = df_vectorPop_dly_yrFiltered['pop_A(L)'] / popTotal_E
            df_vectorPop_dly_yrFiltered_withFrac['pop_A(P)_frac'] = df_vectorPop_dly_yrFiltered['pop_A(P)'] / popTotal_E
        
            for month in [7]:#months:
                fileName = 'tseries_climCat_' + varTemp_pltName + '-' + varPrcp_pltName + '_' + 'allVars' + '_' + loc + fileStr_noOutbreaks + '_mon' + str(month).zfill(2)
                plot_tseries_climCategories_allVars(df_vectorPop_dly_yrFiltered_withFrac, loc, month, varTemp_colName, varTemp_pltName, varPrcp_colName, varPrcp_pltName, fileName)




# make tercile plots for climate categories vs adult pop
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
    varTemp_colName = 'TA (C)'
    varTemp_pltName = strConvert_paren2hyphen(strConvert_parenUnitsRemove(varTemp_colName))
    varPrcp_colName = 'PRCP (mm)'
    varPrcp_pltName = strConvert_paren2hyphen(strConvert_parenUnitsRemove(varPrcp_colName))
    adultPop_colNames = ['pop_A(E)', 'pop_A(L)', 'pop_A(P)']
    adultPop_pltNames = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in adultPop_colNames]

    # parameters related to dengue data
    noOutbreaks_opts = [False]#[True, False] # whether to exclude data from outbreak years
    outbreakYrs = [2017, 2019]

    for loc in locs:
        print(loc)

        # load model pipeline data at location of interest
        # also: convert from daily to monthly
        df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae)

        for noOutbreaks in noOutbreaks_opts:
            print(' noOutbreaks: ' + str(noOutbreaks))
            if noOutbreaks:
                fileStr_noOutbreaks = '_noOutbreaks'
                exclYrs = outbreakYrs  # outbreak years (to exclude from plot)
            else:
                fileStr_noOutbreaks = ''
                exclYrs = []

            # remove data in outbreak years we want to exclude
            df_vectorPop_dly_yrFiltered = df_vectorPop_dly[~df_vectorPop_dly.index.year.isin(exclYrs)]

            # make tercile bar plots with all adult pop variables (E, L, P)
            fileName = 'barTerciles_climCat_' + varTemp_pltName + '-' + varPrcp_pltName + '_' + 'popA' + '_' + loc + fileStr_noOutbreaks
            #plot_bar_climCategories_adultPops(df_vectorPop_dly_yrFiltered, varTemp_colName, varTemp_pltName, varPrcp_colName, varPrcp_pltName, fileName)
            
            # make separate tercile bar plots for each adult pop variable (E, L, P)
            for i, adultPop_colName in enumerate(adultPop_colNames):
                adultPop_pltName = adultPop_pltNames[i]
                fileName = 'barTerciles_climCat_' + varTemp_pltName + '-' + varPrcp_pltName + '_' + adultPop_pltName + '_' + loc + fileStr_noOutbreaks
                plot_bar_climCategories_adultPop(df_vectorPop_dly_yrFiltered, varTemp_colName, varTemp_pltName, varPrcp_colName, varPrcp_pltName, adultPop_colName, fileName)




# make violin plots
if False:
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

    # parameters related to dengue data
    dengueDataYrs_noOutbreaks_opts = [True]#[True, False] # whether to plot (dengue data years w/o outbreak years (2007-20 w/o 2017, 2019)) or (all years)
    dengueDataYrs = list(range(2007, 2020+1))
    outbreakYrs = [2017, 2019]


    # load dengue incidence data
    file_dengue  = '/mnt/redwood/local_drive/chanud/RawData/epiReports/' + 'dengueInc_monthly_2007_2022.csv'
    df_den_monthly_allLocs = pd.read_csv(file_dengue, index_col='Date', parse_dates=True)
    df_den_monthly_allLocs = df_den_monthly_allLocs[df_den_monthly_allLocs.index.year.isin(years)] # cut extraneous years

    for loc in locs:
        print(loc)
        
        # subset dengue data to location of interest
        if loc == 'Negombo':
            den_loc = 'Gampaha'
        elif loc in ['NuwaraEliya', 'Jaffna']:
            den_loc = loc
        else:
            print('Error: unknown loc.')
            exit
        series_den_monthly = df_den_monthly_allLocs[den_loc]
        series_den_monthly.name = 'dengueInc'
        df_den_monthly = series_den_monthly.to_frame()

        # load model pipeline data at location of interest
        # also: convert from daily to monthly
        df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae)
        df_vectorPop_monthly = convert_df_vectorPop_dly2mon(df_vectorPop_dly)

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
            df_den_monthly_yrFiltered = df_den_monthly[~df_den_monthly.index.year.isin(exclYrs)]

            fName = 'violinplot_' + 'TA-TW' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
            #plot_violinplt_TA_TW(df_vectorPop_monthly_yrFiltered, loc, fName)
            
            fName = 'violinplot_' + 'prcp-WH' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
            #plot_violinplt_prcp_WH(df_vectorPop_monthly_yrFiltered, loc, fName)

            fName = 'violinplot_' + 'dev' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
            #plot_violinplt_dev(df_vectorPop_monthly_yrFiltered, loc, fName)

            fName = 'violinplot_' + 'surv' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
            #plot_violinplt_surv(df_vectorPop_monthly_yrFiltered, loc, fName)

            fName = 'violinplot_' + 'dev-surv' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
            #plot_violinplt_dev_surv(df_vectorPop_monthly_yrFiltered, loc, fName)


            # manually create the popFrac columns and add them to vectorPop df
            popTotal_E = len(initDays) * initEggs # this code only works for eggs AND larvae and pupae as long as initEggs = initLarvae = initPupae
            df_vectorPop_monthly_yrFiltered_withFrac = df_vectorPop_monthly_yrFiltered.copy(deep=True)
            df_vectorPop_monthly_yrFiltered_withFrac['pop_A(E)_frac'] = df_vectorPop_monthly_yrFiltered['pop_A(E)'] / popTotal_E
            df_vectorPop_monthly_yrFiltered_withFrac['pop_A(L)_frac'] = df_vectorPop_monthly_yrFiltered['pop_A(L)'] / popTotal_E
            df_vectorPop_monthly_yrFiltered_withFrac['pop_A(P)_frac'] = df_vectorPop_monthly_yrFiltered['pop_A(P)'] / popTotal_E

            fName = 'violinplot_' + 'popFrac' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
            #plot_violinplt_popFrac(df_vectorPop_monthly_yrFiltered_withFrac, loc, fName)

            if dengueDataYrs_noOutbreaks_opts: # only do this if we're restricting years to no outbreak dengue data years (2007-2020 minus 2017 and 2019)

                # create new df combining vectorPop and denInc dfs
                df_vectorPop_denInc_monthly_yrFiltered = df_vectorPop_monthly_yrFiltered_withFrac.copy(deep=True)
                df_vectorPop_denInc_monthly_yrFiltered['dengueInc'] = df_den_monthly_yrFiltered['dengueInc']

                fName = 'violinplot_' + 'denInc' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                #plot_violinplt_denInc(df_vectorPop_denInc_monthly_yrFiltered, loc, fName)

                fName = 'violinplot_' + 'popFrac-denInc' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                #plot_violinplt_popFrac_denInc(df_vectorPop_denInc_monthly_yrFiltered, loc, fName)

                fName = 'violinplot_' + 'allVars' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                plot_violinplt_allVars(df_vectorPop_denInc_monthly_yrFiltered, loc, fName)




# make scatter plots of crosstab hi-hi values
if False:

    # parameters related to modeling pipeline
    dataSrc = 'MERRA_IMERG'
    locs    = ['NuwaraEliya']#['Negombo', 'Jaffna', 'NuwaraEliya']
    years   = list(range(2001, 2020+1))
    months  = list(range(1, 12+1))
    intv    = 'intv00'
    initDays   = list(range(1, 7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000

    # parameters related to dengue data
    lagTimes = [0, 1, 2]              # in months
    den_colName_base = 'denInc'
    noOutbreaks_opts = [True, False]  # whether to include data from outbreak years
    outbreakYrs = [2017, 2019]
    
    # parameters related to modeling pipeline + dengue data
    startDate = '2007-01-01'  # start date (first year of dengue data)
    endDate   = '2020-12-31'  # end date (last year of downloaded MERRA/IMERG data)
    
    # parameters related to the crosstab element of interest
    # format is 'terc1-terc2' where terc1 is tercile of the variable, terc2 is tercile of dengue
    crosstab_elts = ['HI-HI']#['HI-HI', 'MED-HI', 'LO-HI', 'HI-LO', 'MED-LO', 'LO-LO']
    map_tercAbbr_to_tercName = {'HI': 'High', 'MED': 'Medium', 'LO': 'Low'} # map to names of rows/columns in crosstab csv file

    for crosstab_elt in crosstab_elts:
        print('Crosstab element: ' + crosstab_elt)

        tercileStr_modelVar, tercileStr_den = crosstab_elt.split('-')
        dir_crosstabs = 'crosstabs_' + tercileStr_modelVar + '_' + tercileStr_den + '/'
        tercileName_modelVar = map_tercAbbr_to_tercName.get(tercileStr_modelVar) # get corresponding row name in crosstabs csv file
        tercileName_den      = map_tercAbbr_to_tercName.get(tercileStr_den)      # get corresponding col name in crosstabs csv file
        colName_4_df = tercileName_modelVar + '-' + tercileName_den + ' Value' # (e.g., 'High-High Value')

        for loc in locs:
            print(loc)
            
            df_list = []
            df_name_list = []
            for i_lag, lagTime in enumerate(lagTimes):
                print(' Lag (mo): ' + str(lagTime))

                for noOutbreaks in noOutbreaks_opts:
                    print('   noOutbreaks: ' + str(noOutbreaks))
                    if noOutbreaks:
                        fileStr_noOutbreaks = '_noOutbreaks'
                        exclYrs = outbreakYrs  # outbreak years (to exclude from plot)
                    else:
                        fileStr_noOutbreaks = ''
                        exclYrs = []
                    
                    den_colName = den_colName_base + '-lag' + str(lagTime) # this must match col name format set in create_df_withLags()
                    den_pltName = den_colName # just reusing the column name here for simplicity
                        
                    # read in crosstab hi-hi values
                    fileName_in = 'crosstab_tercile_splitbyMon_' + tercileStr_modelVar + '-' + tercileStr_den + '_' + den_pltName + '_' + loc + fileStr_noOutbreaks
                    filePath_norm = DIR_WHATCHEM_OUTPUT_DATA_PROCESSED + dir_crosstabs + fileName_in + '_norm.csv'
                    df_crosstab_elt_norm_oneScenario = pd.read_csv(filePath_norm, index_col=0, header=0)

                    # rename based on lag+outbreakYr scenario and add to list
                    df_name = 'lag' + str(lagTime) + fileStr_noOutbreaks
                    df_crosstab_elt_norm_oneScenario.rename(columns={colName_4_df: df_name}, inplace=True)
                    df_list.append(df_crosstab_elt_norm_oneScenario)

            df_crosstab_elt_norm = pd.concat(df_list, axis=1)

            fileName_out = 'crosstab_tercile_splitbyMon_' + tercileStr_modelVar + '-' + tercileStr_den + '_' + den_colName_base + '_' + loc + '_norm'
            plot_scatter_crosstabs(df_crosstab_elt_norm, loc, tercileName_modelVar, tercileName_den, fileName_out)




# make tercile-tercile plots of model vars against dengue incidence
if True:
    
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
    den_colName_base = 'denInc'      # we will rename columns to this string + the lag time (original name is the district)
    noOutbreaks_opts = [True, False] # whether to include data from outbreak years
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

                    fileName = 'barTerciles_' + var_pltName + '_' + den_pltName + '_' + loc + fileStr_noOutbreaks
                    plot_bar_varX_denInc(df_vectorPop_den_monthly_yrFiltered, var_colName, var_pltName, den_colName, den_pltName, fileName)


