# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 15:28:05 2024

@author: kjloo
"""
# This file calcualtes the diffuse skylight ratio from ERA5 climate data downloaded from the Copernicus climate store (https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=download)
# and processes them to be compatible for upload into Google Earth Engine.
# We used the hourly reanalysis data at single levels, taking the "Mean surface direct short-wave radiation flux" and the "Mean surface downward short-wave radiation flux" and selected one month of one year
# The naming conventions used in the ERA5 datasets changed during the course of this project

# Import packages
import os
import xarray as xr
import rioxarray as rio
import numpy as np

#%% Input parameters
# File name prefix - the files downloaded from Copernicus climate data store need to be saved with a file name prefix and two digits to indicate the month (e.g. "01" for January and "12" for December)
# and four digits to indicate the year (e.g. "prefix_01_2012.nc" for the file containing data for January of 2012). The prefix needs to be specified here
file_prefix = "ssw_radiation_hourly"

# Variable names
# Mean surface direct short-wave radiation flux
direct_rad = "avg_sdirswrf" #Note that this could also be 'msdrswrf' depending on when the data was downloaded due to changes in the naming conventions

# Mean surface downward short-wave radiation flux
all_rad = "avg_sdswrf" #Note that this could also be'msdwswrf' depending on when the data was downloaded due to changes in the naming conventions

# Month and year of the dataset
month = 1 # months are numerical between 1 and 12
yr = 2013

# Specifiy directory where the input files are located
file_source_location = "C:/"

# Specify directory to save output files
file_save_location = "C:/"

#%% Define function to calculate diffuse skylight ratio from files downloaded from the ERA5 through Copernicus climate data store

def dsr_calc(file_prefix,direct_rad,all_rad,month,yr,file_source_location,file_save_location):
    os.chdir(file_source_location)
    if month < 10:
        month_str = str('0'+str(month))
    else:
        month_str = str(month)

    file_name_1 = file_prefix + "_" + month_str + "_" + str(yr) + ".nc"
    
    print('Opening datasets for '+str(yr)+'-'+month_str)
    
    # in the ERA 5, the mean ratesfor hourly data are taken over the hour preceeding the time stamp on the data. This means that for the zero hour, this data is actually from the 
    # last hour of the day before. To correct for this, we discard the 0 hour data and tack on the 0 hour data from the next month, shifting the hourly data by one hour
    if direct_rad == "avg_sdirswrf":
        nc_file_1 =  xr.open_dataset(file_name_1).isel(valid_time = slice(1,None))
    elif direct_rad == "msdrswrf":
        nc_file_1 =  xr.open_dataset(file_name_1).isel(time = slice(1,None))
    else:
        print("Invalid variable name")
        return
    
    if (month+1) < 10:
        file_name_2 = file_prefix + "_0" + str(month + 1) + "_" + str(yr) + ".nc"
    elif (month+1) == 13:
        file_name_2 = file_prefix + "_01_" + str(yr + 1) + ".nc"
    else:
        file_name_2 = file_prefix + "_" + str(month + 1) + "_" + str(yr) + ".nc"
    
    nc_file_2 = xr.open_dataset(file_name_2).isel(time = 0)
    
    if direct_rad == "avg_sdirswrf":
        nc_file = xr.concat([nc_file_1,nc_file_2], dim = 'valid_time')
    elif direct_rad == "msdrswrf":
        nc_file = xr.concat([nc_file_1,nc_file_2], dim = 'time')
    
    nc_file_1.close()
    nc_file_2.close()
    
    del(nc_file_1) #since these are big files, these variables are deleted to save memory
    del(nc_file_2)
    
    # Downward radiation all
    if direct_rad == "avg_sdirswrf":
        ssr_down = nc_file.avg_sdswrf
    elif direct_rad == "msdrswrf":
        ssr_down = nc_file.msdwswrf
    
    # Mask out invalid data
    ssr_down = ssr_down.where(ssr_down>=0,np.nan)
    
    # Define geospatial parameters
    ssr_down = ssr_down.rio.set_spatial_dims(x_dim='longitude', y_dim='latitude')
    ssr_down.rio.write_crs('EPSG:4326',inplace=True)
    
    # Downward radiation direct
    if direct_rad == "avg_sdirswrf":
        ssr_dir = nc_file.avg_sdirswrf
    elif direct_rad == "msdrswrf":
        ssr_dir = nc_file.msdrswrf
   
    # Define geospatial parameters
    ssr_dir = ssr_dir.rio.set_spatial_dims(x_dim='longitude', y_dim='latitude')
    ssr_dir.rio.write_crs('EPSG:4326',inplace=True)
    
    del(nc_file)
    
    # Calculate diffuse skylight ratios
    print('Calculating...')
    dsr = (ssr_down - ssr_dir) / ssr_down
    
    # Mask out invalid data
    dsr = dsr.where(dsr>0,np.nan)
    dsr = dsr.where(dsr<1,np.nan)
    
    del(ssr_down)
    del(ssr_dir)
    
    os.chdir(file_save_location)
    print('Saving files')
    
    # Save daily rasters with hourly data as bands
    for n in range(int(len(dsr[:,0,0])/24)):
        img = dsr[n*24:(n+1)*24,:,:]
        img_daily = []
        for h in range(24):
            img_h = img[h,:,:]
            img_daily.append(img_h)
        img_n = xr.concat(img_daily,dim='band')

        if (n+1) < 10:
            img_n.rio.to_raster(str('dsr_'+str(yr)+'_'+month_str+'_0'+str(n+1)+'.tiff'))
        else:
            img_n.rio.to_raster(str('dsr_'+str(yr)+'_'+month_str+'_'+str(n+1)+'.tiff'))
    
    del(dsr)
    del(img_daily)
    del(img_n)

#%% Run for one month of one year

dsr_calc(file_prefix,direct_rad,all_rad,month,yr,file_source_location,file_save_location)

#%% Run for all months of one year

for m in range(1,13):
    dsr_calc(file_prefix,direct_rad,all_rad,m,yr,file_source_location,file_save_location)

