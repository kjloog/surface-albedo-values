# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 08:58:36 2024

@author: kjloo
"""
# Import necessary packages
import os
import pandas as pd

#%% Define function to combine geographic chunks into global datasets at a monthly level
# Takes the land cover type (lc_in), the year (yr) and the month (m) associated with the data output file as inputs (numeric inputs)
# Also requires the location of the source files (the chunks output from dailymean_blueksy_albedo.py. These files will be saved into 
# the Google Drive that is associated with the users Google Earth Engine account. This Google Drive either needs to be synced with your computer and 
# setup as directory on your computer, or the files need to be moved into an accessible directory). A separate directory to save the processed files also needs
# to be specified (inputs as strings with the path to the directory). If a suffix (sfx) was used in dailymean_bluesky_albedo.py, that also needs to be specified here
 
def chunks_to_months(lc_in,yr,m,source_file_loc,save_file_loc,sfx):
    os.chdir(source_file_loc) 
    alphas_list = ['a','b','c','d','e','f','g','h','i','j'] # list of the file name suffixes for the geographic file chunks
    band_names = ['01-Jan','02-Feb','03-Mar','04-Apr','05-May','06-Jun','07-Jul','08-Aug','09-Sep','10-Oct','11-Nov','12-Dec']
    
    # convert numerical land cover types into strings
    lc_num_list = list(range(1,15)) + [16,11.27,11.5,12.7,12.8,14.25,14.35,'forest']
    lc_str_list = list(range(1,15)) + [16,'11-27','11-50','12-7','12-8','14-25','14-35','forest']
    lc_dict = dict(zip(lc_num_list,lc_str_list))
    
    # error handling for invalid land cover types
    try:
        lc = lc_dict[lc_in]
    except:
        print('Invalid land cover type')
        return
    
    df_list = [] #creates an empty list to store the geographic chunks
    for alpha in alphas_list: #iterates over all the chunks
        file_name = str('albedo_bls_daily_lc-'+str(lc)+'_25kms_'+str(yr)+'-'+band_names[m-1][:2]+sfx+'_'+alpha+'-10.csv') #defines the file name of the chunk
        try:
            df_list.append(pd.read_csv(file_name)) # attempts to open the chunk
        except FileNotFoundError:
            print('lc-'+str(lc)+' '+str(m) + " failed. Missing " + alpha) #if the chunk does not exist or otherwise cannot be opened, the function terminates and prints a message to the user identifying the problematic file
            return

    df = pd.concat(df_list) #combine the chunks into one file
        
    os.chdir(save_file_loc)
    df.to_csv(str('albedo_bls_daily_lc-'+str(lc)+'_25kms_'+str(yr)+'-'+band_names[m-1][:2]+sfx+'.csv')) #save the combined file
    print('lc-'+str(lc)+' '+str(m) + " successful") #prints a message to let the user know the file was successfully processed
    
#%% Run function to aggregate one month of data for one land cover type
# Input parameters
lc = 11 #land cover type (see list of land cover types in readme)
year = 2016
month = 6
source_file_location = 'H:/My Drive' #source file location if Google Drive is mapped as a directory
save_file_location = 'C:/Users/Public/Documents'
sfx = ""

# Run function
chunks_to_months(lc,year,month,source_file_location,save_file_location,sfx)

#%% Run function to aggregate chunks into monthly data files for all months and years (loops function across months and years)
# Input parameters
lc = 11 #land cover type (see list of land cover types in readme)
source_file_location = 'H:/My Drive' #source file location if Google Drive is mapped as a directory
save_file_location = 'C:/Users/Public/Documents'
sfx = ""

for year in range(2011,2021):
    for month in range(13):
        chunks_to_months(lc,year,month,source_file_location,save_file_location,sfx)
        print("Year: " + str(year) + " Month: " + str(month) )

#%% Define function to combine monthly data into yearly datasets
# Takes the land cover type (lc_in) and the year (yr) inputs (numeric inputs) (Note that you need to have combined the geographic chunks into global monthly files first)
# Also requires the location of the source files (the files output from the previous function (chunks_to_months()). A separate directory to save the processed files also needs
# to be specified (inputs as strings with the path to the directory). If a suffix (sfx) was used in dailymean_bluesky_albedo.py, that also needs to be specified here
# The csv file containing the monthly weighting factors is required and needs to be in the same location as the other source files

def months_to_year(lc_in,yr,source_file_loc,sfx,save_file_loc):
    # convert numerical land cover types into strings
    lc_num_list = list(range(1,15)) + [16,11.27,11.5,12.7,12.8,14.25,14.35,'forest']
    lc_str_list = list(range(1,15)) + [16,'11-27','11-50','12-7','12-8','14-25','14-35','forest']
    lc_dict = dict(zip(lc_num_list,lc_str_list))
    
    # error handling for invalid land cover types
    try:
        lc = lc_dict[lc_in]
    except:
        print('Invalid land cover type')
        return
    
    #open January dataset
    os.chdir(source_file_loc)
    
    file_name = str('albedo_bls_daily_lc-'+str(lc)+'_25kms_'+str(yr)+'-01'+sfx+'.csv')    
    df_1 = pd.read_csv(file_name).set_index(['lat','lon']) #this file contains both the means and the variances
    
    #load the weighting factors used to calculate yearly averages
    if len(df_1) == 69937:
        df_wt = pd.read_csv('monthly_weighting_factors.csv').set_index(['lat','lon'])
    elif len(df_1) == 6256091:
        df_wt = pd.read_csv('monthly_weighting_factors_005.csv').set_index(['lat','lon'])
    
    #create separate dataframes for the means and the variances
    band_names = ['01-Jan','02-Feb','03-Mar','04-Apr','05-May','06-Jun','07-Jul','08-Aug','09-Sep','10-Oct','11-Nov','12-Dec']
    df_means = df_1[['sum_mean']].rename(columns = {'sum_mean':band_names[0]})
    df_var = df_1[['sum_variance']].rename(columns = {'sum_variance':band_names[0]})
    
    #Join the weighting factors
    df_means = df_means.join(df_wt)
    df_var = df_var.join(df_wt)
    
    df_means[band_names[0]+'-alb-wt'] = df_means[band_names[0]] * df_means[band_names[0]+'-wt']
    df_var[band_names[0]+'-alb-wt'] = df_var[band_names[0]] * df_var[band_names[0]+'-wt']
    
    # Open and join other monthly datasets
    os.chdir(source_file_loc)
    for m in range(2,13):
        file_name = str('albedo_bls_daily_lc-'+str(lc)+'_25kms_'+str(yr)+'-'+band_names[m-1][:2]+sfx+'.csv')
        
        df_m = pd.read_csv(file_name).set_index(['lat','lon'])
        
        df_means_m = df_m[['sum_mean']].rename(columns = {'sum_mean':band_names[m-1]})
        df_var_m = df_m[['sum_variance']].rename(columns = {'sum_variance':band_names[m-1]})
        
        df_means = df_means.join(df_means_m)
        df_var = df_var.join(df_var_m)
        
        df_means[band_names[m-1]+'-alb-wt'] = df_means[band_names[m-1]] * df_means[band_names[m-1]+'-wt']
        df_var[band_names[m-1]+'-alb-wt'] = df_var[band_names[m-1]] * df_var[band_names[m-1]+'-wt']
    
    # Calculate yearly averages
    df_means['00-Yr'] = df_means[[band_name +'-alb-wt' for band_name in band_names]].sum(axis=1,skipna=False)
    
    df_var['00-Yr'] = df_var[[band_name +'-alb-wt' for band_name in band_names]].sum(axis=1,skipna=False)
    
    df_means = df_means[band_names + ['00-Yr']]
    df_var = df_var[band_names + ['00-Yr']]
    # save files
    os.chdir(save_file_loc)
    means_file_name = str('albedo_bls_daily_lc-'+str(lc)+'_'+str(yr)+'.csv')
    df_means.to_csv(means_file_name)
    
    var_file_name = str('albedo_bls_daily_var_lc-'+str(lc)+'_'+str(yr)+'.csv')
    df_var.to_csv(var_file_name)

#%% Run function to aggregate one year of one land cover type
# Input parameters
lc = 11 #land cover type (see list of land cover types in readme)
year = 2016
source_file_location = 'H:/My Drive' #source file location if Google Drive is mapped as a directory
save_file_location = 'C:/Users/Public/Documents'
sfx = ""

# Run function
months_to_year(lc,year,source_file_location,sfx,save_file_location)

#%% oRun function for one land cover type for all years
# Input parameters
lc = 10
source_file_location = 'H:/My Drive' #source file location if Google Drive is mapped as a directory
save_file_location = 'C:/Users/Public/Documents'
sfx = ''

# Run function
for year in range(2011, 2021):
    months_to_year(lc,year,source_file_location,sfx,save_file_location)
    print(str(year))

#%% Define function to combine yearly data into a climatology
# This function takes the yearly files created with the prvious function and averages them into a climatology.
# It takes as inputs the start and end years (start_year and end_year) of the climatology preiod, and land cover type (lc_in),
# the directory containing the yearly files (source_file_loc) and a directery where to save the output files (save_file_loc)
# The csv file containing the monthly weighting factors is required and needs to be in the same location as the other source files

def years_to_clima(start_year,end_year,lc,source_file_loc,save_file_loc):
    os.chdir(source_file_loc)
    
    file_list_means = []
    file_list_vars = []
    
    # Combine all yearly files (means and variances)
    for yr in range(start_year,end_year+1):
        means_file_name = str('albedo_bls_daily_lc-'+str(lc)+'_'+str(yr)+'.csv')
        var_file_name = str('albedo_bls_daily_var_lc-'+str(lc)+'_'+str(yr)+'.csv')
        
        df_means = pd.read_csv(means_file_name).set_index(['lat','lon'])
        df_var = pd.read_csv(var_file_name).set_index(['lat','lon'])
        
        file_list_means.append(df_means)
        file_list_vars.append(df_var)
    
    # Take means across the months across all years
    df_concat_means = pd.concat(file_list_means)
    df_rows_means = df_concat_means.groupby(df_concat_means.index).mean()
    df_rows_means.index = pd.MultiIndex.from_tuples(df_rows_means.index,names=['lat','lon'])
    
    df_concat_vars = pd.concat(file_list_vars)
    df_rows_vars = df_concat_vars.groupby(df_concat_vars.index).mean()
    df_rows_vars.index = pd.MultiIndex.from_tuples(df_rows_vars.index,names=['lat','lon'])
    
    # Retrieve weighting factors
    if len(df_means) == 69937:
        df_wt = pd.read_csv('monthly_weighting_factors.csv').set_index(['lat','lon'])
    elif len(df_means) == 6256091:
        df_wt = pd.read_csv('monthly_weighting_factors_005.csv').set_index(['lat','lon'])    
    
    # Join weighting factors dataset
    df_rows_means = df_rows_means.join(df_wt)
    df_rows_vars = df_rows_vars.join(df_wt)
    
    band_names = ['01-Jan','02-Feb','03-Mar','04-Apr','05-May','06-Jun','07-Jul','08-Aug','09-Sep','10-Oct','11-Nov','12-Dec']
    
    # Multiply weights by monthly values
    for band_name in band_names:
        df_rows_means[band_name+'-alb-wt'] = df_rows_means[band_name] * df_rows_means[band_name+'-wt']
        df_rows_vars[band_name+'-alb-wt'] = df_rows_vars[band_name] * df_rows_vars[band_name+'-wt']
    
    # Combine weighted values into yearly average
    df_rows_means['00-Yr'] = df_rows_means[[band_name+'-alb-wt' for band_name in band_names]].sum(axis=1,skipna=False)
    df_rows_vars['00-Yr'] = df_rows_vars[[band_name+'-alb-wt' for band_name in band_names]].sum(axis=1,skipna=False)
    
    df_rows_means = df_rows_means[band_names + ['00-Yr']]
    df_rows_vars = df_rows_vars[band_names + ['00-Yr']]
    
    # Save files
    os.chdir(save_file_loc)
    means_file_name = str('albedo_bls_daily_lc-'+str(lc)+'_'+str(start_year)+'-'+str(end_year)+'.csv')
    df_rows_means.to_csv(means_file_name)
    
    var_file_name = str('albedo_bls_daily_var_lc-'+str(lc)+'_'+str(start_year)+'-'+str(end_year)+'.csv')
    df_rows_vars.to_csv(var_file_name)

#%% Run function for one land cover type
# Input parameters
start_year = 2011
end_year = 2020
source_file_loc = 'C:/Users/Public/Documents'
save_file_loc = 'C:/Users/Public/Documents'
lc = 9

# Run function
years_to_clima(start_year, end_year, lc, source_file_loc, save_file_loc)

#%% Run function for all land cover types
start_year = 2011
end_year = 2020
source_file_loc = 'C:/Users/Public/Documents'
save_file_loc = 'C:/Users/Public/Documents'

lc_list = list(range(1,15)) + [16, 11.27,11.50,12.7,12.8,14.25,14.35,'forest']

for lc in lc_list:
    years_to_clima(start_year, end_year, lc, source_file_loc, save_file_loc)
    print(lc)
