# scripts_createMetInput.py
#
# Assorted blocks of code for creating meteorological input files for WHATCH'EM.
# Essentially just runs the create_met_input() function from createMetInput.py
# with different arguments.


# import required packages
# ************************************************
import calendar
import sys

# import constants
# ************************************************
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY/metInputs")
from createMetInput import create_met_input, create_met_input_clim




## Create met inputs for each month of our time range of interest (2001 to 2020)
if False:
    for year in range(2001, 2021):
        for month in range(1, 13):
            print(f"Year: {year}  Month: {month}")

            # Define start and end times in YYYYMMDDHH format
            days_in_month = calendar.monthrange(year, month)[1]
            t_start = year * 1000000 + month * 10000 + 100
            t_end = year * 1000000 + month * 10000 + days_in_month * 100 + 23

            # Next, run the code to create the met input files
            create_met_input("MERRA_IMERG", "Jaffna", t_start, t_end)
            create_met_input("MERRA_IMERG", "Negombo", t_start, t_end)
            create_met_input("MERRA_IMERG", "NuwaraEliya", t_start, t_end)



## Create met inputs for the climatology of our time range of interest (2001 to 2020)
if False:
    # Define start and end times in YYYYMMDDHH format
    t_start  = 2001010100
    t_end    = 2020123123
    climType = "month"

    # Next, run the code to create the met input files
    create_met_input_clim("MERRA_IMERG", "Jaffna", t_start, t_end, climType)
    create_met_input_clim("MERRA_IMERG", "Negombo", t_start, t_end, climType)
    create_met_input_clim("MERRA_IMERG", "NuwaraEliya", t_start, t_end, climType)

