#*******************************************************
# fns_analysis.py
#*******************************************************
#
# The functions in this file process the data from the
# modeling pipeline. In some cases the resulting data
# is saved to file. In other cases they're returned by
# the function (e.g., as a precursor to plotting the data).
#
#
# Revision history
# 2024-10-19, CNY: file created
#
# Functions:
# 



# import required packages
# ************************************************
import pandas as pd
import sys


# import constants and custom functions
# ************************************************
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/4_analysis") # directory of helper functions file
from fns_helper import categorizeSeries_LowMedHigh_byMon




def calc_tercile_crosstab_splitByMon(df_in, varX_colName, varY_colName, fileName):

   filedir_out = '/home/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/4_data_analysis/'
   df = df_in.copy(deep=True)  # create a temporary copy of the df so that we don't alter the original
   
   # Categorize data into Low, Medium, High, on a monthly basis (and then recombine all monthly results into one series)
   varX_terciles = categorizeSeries_LowMedHigh_byMon(df[varX_colName])    
   varY_terciles = categorizeSeries_LowMedHigh_byMon(df[varY_colName])

   # If categorization failed for either variable, we're not going to do the calculation.
   if varX_terciles.isnull().all():
      print(f'Error: categorization failed for variable {varX_colName}. Skipping calculation.')
      return
   if varY_terciles.isnull().all():
      print(f'Error: categorization failed for variable {varY_colName}. Skipping calculation.')
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

   crosstab.to_csv(filedir_out + fileName + '.csv')
   crosstab_normalized.to_csv(filedir_out + fileName + '_norm.csv')

