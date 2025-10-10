# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 14:41:54 2023

@author: kjloo
"""


# Due to the large processing requirements, this code can only run one year at a time for one land cover 
# type which is also split into geographical chunks. You then need to run "file_aggregation" to aggregate
# the chunks into global datasets, and to aggregate the yearly data into climatologies

# Input parameters
# Land cover type by numerical code (see readme)
lc = 7 # put land cover type here

# Year
yr = 2016

# Output scale. The results can be output at either a 0.5 degree scale (input 0.5) or a 0.05 degree scale (input 0.05)
output_scale = 0.5

#%% Import required packages
import ee
import math

#%% Activate google earth engine (note that you need an active GEE account to run this code)
ee.Authenticate()
ee.Initialize()

#%% Load datasets from GEE

# Surface albedo data (MODIS MCD43A1 v6.1 )
albedo_BRDF = ee.ImageCollection("MODIS/061/MCD43A1")
albedo_scale = albedo_BRDF.first().projection().nominalScale()

# Land cover data (MODIS MCD12Q1 v6.1)
lc_data = ee.ImageCollection('MODIS/061/MCD12Q1')
lc_quality = lc_data.select('QC')
lc_scale = lc_data.first().projection().nominalScale()

# Output geometries
#load grid at 0.5 degrees (~50 km)
albedo_geom_05 = ee.FeatureCollection('projects/kloog-project-2/assets/albedo_geom_05_v1')
#load grid at 0.05 degrees (~5 km)
albedo_geom_005 = ee.FeatureCollection('projects/kloog-project-2/assets/albedo_geom_005')

#%% Land cover change mask
# This function takes an image from the MODIS MCD12Q1 land cover dataset as an input and creates a mask that filters out 
# any data points where the land cover type is not the same in either the previous or following year

def change_mask(image):
    # Identify the start and end dates of the input image
    start_date = image.date()
    end_date = start_date.advance(1, 'year')

    # Identify start and end dates of the image from one year previous
    start_date_before = start_date.advance(-1, 'year')
    end_date_before = start_date
    
    # Identify start and end dates of the image from one year after
    start_date_after = end_date
    end_date_after = start_date_after.advance(1, 'year')
    
    # Select images from the year befreo and the year after based on the dates determined
    image_before = lc_data.select('LC_Type1').filterDate(start_date_before, end_date_before).first()
    image_after = lc_data.select('LC_Type1').filterDate(start_date_after, end_date_after).first()

    # Compare current image to the images from other years
    lc_change_mask_before = image.eq(image_before) # returns 1 where values are euqal in both years and 0 where unequal
    lc_change_mask_after = image.eq(image_after)
    
    # Combine before and after into one mask
    lc_change_mask = lc_change_mask_before.where(lc_change_mask_after.eq(0), 0) 
    return lc_change_mask

#%% Land cover spatial homogeneity mask
# This function takes an image from the MODIS MCD12Q1 land cover dataset as an input and creates a mask that filters out 
# any data points where the land cover type is not the same in the surrounding cells

def homog_mask(image):
     # Define the base land cover types
     lc_orig = image
     # Calculates the minimum of all the land cover types (expressed as integers) in the surrounding cells
     lc_min = image.reduceNeighborhood(reducer=ee.Reducer.min(),
                                       kernel=ee.Kernel.square(radius=1,
                                                               units='pixels'),
                                       inputWeight='mask')
     # Calculates the maximum of all the land cover types (expressed as integers) in the surrounding cells
     lc_max = image.reduceNeighborhood(reducer=ee.Reducer.max(),
                                       kernel=ee.Kernel.square(radius=1,
                                                               units='pixels'),
                                       inputWeight='mask')
     # If the min, max, and base land cover types are all the same value, this indicates that all surrounding cells are the same.
     # Where either the min or max value differs from the base land cover type, these values are masked out
     lc_homog_mask = lc_orig.where(lc_orig.eq(lc_min), 1)
     lc_homog_mask = lc_homog_mask.where(lc_orig.neq(lc_min), 0)
     lc_homog_mask = lc_homog_mask.where(lc_orig.neq(lc_max), 0)
     return lc_homog_mask

#%% Land cover data quality mask
# This function creates a mask based on the MODIS MCD12Q1 data quality flag

def qual_mask(image):
    lc_qual_x = image
    lc_qual_mask_x = ee.Image.constant(1).where(lc_qual_x.eq(1), 0).where(lc_qual_x.eq(2), 0).where(lc_qual_x.eq(3), 0).where(lc_qual_x.eq(4), 0).where(lc_qual_x.eq(5), 0).where(lc_qual_x.eq(7), 0)
    return lc_qual_mask_x

#%% Function to calculate daily-mean blue-sky albedo calculation and export results to file
# This is a series of nested functions that use global parameters. Although this is generally
# considered poor coding practice this was necessary to make use of the GEE mapping function which
# allows us to map one function across all the images of an image collection or all elements of a list. 
# One constraint of the mapping functionality is that only functions with a single input can be mapped.
# To work around this, we have created a nested function where upper levels provide necessary input 
# parameters to further nested functions


# Define hour angles
# Black-sky albedo is calculated every hour but on the half hour to correspond with the diffuse sylight data
# taken from the ERA5 dataset. A function to calculate hourly surface albedo values iterates over this list
hours = ee.List.sequence(start=0.5, end=23.5, step=1)

# Function to calculate the daily-mean blue-sy albedo values for a given land cover type for one month of the year
# yr is the year, d is sampling distance (in km), lc_in is the land cover type, m is the month of the year
# tile_scale is used to specify the tile scaling for GEE and sfx is a suffix that can be added to the output file
def albedo_bls_export(yr, d, lc_in, geom, m, tile_scale, sfx):
    # Identify dates at beginning and end of year to select corresponding land cover data from the MODIS MCD12Q1
    # land cover dataset which is available on a yearly basis
    start_date = ee.Date.fromYMD(yr, 1, 1)
    end_date = start_date.advance(1, 'year')
    
    #Band selection for the MODIS MCD12Q1. Depending on the land cover type selected in the input parameters different 
    #bands are required to retrieve the data
    if lc_in == 12.7:
        band = 'LC_Type5'
        lc = 7
    elif lc_in == 12.8:
        band = 'LC_Type5'
        lc = 8
    elif lc_in == 14.25:
        band = 'LC_Prop2'
        lc = 25
    elif lc_in == 14.35:
        band = 'LC_Prop2'
        lc = 35
    else:
        band = 'LC_Type1'
        lc = lc_in

    def albedo_bls_daily(image):
        # This function iterates over the MODIS MCD43A1 dataset (images are daily) and calculates the daily-mean blue-sky albedo
        # Extract starting date of the image to select images from the same date across datasets
        start_time = ee.Date(image.get('system:time_start'))
        end_time = start_time.advance(1, 'day')

        # BRDF model parameters from the MODIS MCD43A1 dataset
        BRDF_iso = image.select('BRDF_Albedo_Parameters_shortwave_iso').multiply(0.001)
        BRDF_vol = image.select('BRDF_Albedo_Parameters_shortwave_vol').multiply(0.001)
        BRDF_geo = image.select('BRDF_Albedo_Parameters_shortwave_geo').multiply(0.001)

        # Calculation of white-sky albedo from the BRDF model parameters
        albedo_ws_x = BRDF_iso.add(BRDF_vol.multiply(0.189184)).add(BRDF_geo.multiply(-1.377622))

        # Retrieve the hourly diffuse skylight ratios for the same date (images are daily and contain 24 bands, one for each hour of the day)
        dsr_x = ee.ImageCollection('projects/kloog-project-2/assets/dsr_hourly/dsr_'+str(yr)).filterDate(start_time, end_time).first()
        dsr_inv_x = dsr_x.subtract(1).multiply(-1)

        # Define additional parameters needed to calculate the blue sky albedo
        # Longitude (in degrees)
        lon = image.pixelCoordinates(projection='EPSG:4326').select('x')
        # Latitude (in radians)
        lat = image.pixelCoordinates(projection='EPSG:4326').select('y').multiply(math.pi).divide(180)
        #Day of year
        DOY = start_time.getRelative('day', 'year')
        # Solar declination angle (based on the day of year)
        SDA = (DOY.add(284)).multiply(360).divide(365.24).multiply(math.pi).divide(180).sin().multiply(23.44).multiply(math.pi).divide(180)# in radians
        

        # Function to calculate the cosine of the solar zenith angle (SZA) from UTC time (input is a constant)
        def cos_SZA_calc(hour):
            # Convert input UTC time to image
            UTC_time_x = ee.Image.constant(hour)
            # Calculate solar time for each location based on the longitute and UTC time
            solar_time_x = UTC_time_x.add(lon.divide(15))
            solar_time_x = solar_time_x.where(solar_time_x.lt(0), solar_time_x.add(24))
            # Calculate hour angles from solar time
            hour_angles_x = (solar_time_x.subtract(ee.Number(12))).multiply(15).divide(180).multiply(math.pi)

            # Calculate cos(SZA) from hour angles
            cos_SZA_x = lat.sin().multiply(SDA.sin()).add(lat.cos().multiply(SDA.cos()).multiply(hour_angles_x.cos()))
            # Mask where cos(SZA) is less than 0
            mask = cos_SZA_x.mask().where(cos_SZA_x.lt(0), 0)
            cos_SZA_x = cos_SZA_x.updateMask(mask)
            return cos_SZA_x
        
        # Map function across the hours list to get the global solar zenith angles for all hours of the day
        cos_SZA_hourly = ee.ImageCollection.fromImages(hours.map(cos_SZA_calc))

        # Function to calculate the instantaneous black-sky albedo based on the cosine of the SZA
        def BKS_calc(cos_sza_x):
            # Define solar zenith angle to input into BRDF polynomial function
            SZA = cos_sza_x.acos()
            # Calculate the three terms of the BRDF polynomial (iso, vol, and geo)
            BRDF_term_1 = BRDF_iso
            BRDF_term_2 = BRDF_vol.multiply(ee.Image.constant(-0.007574).add(SZA.multiply(SZA).multiply(-0.070987)).add(SZA.multiply(SZA).multiply(SZA).multiply(0.307588)))
            BRDF_term_3 = BRDF_geo.multiply(ee.Image.constant(-1.284909).add(SZA.multiply(SZA).multiply(-0.166314)).add(SZA.multiply(SZA).multiply(SZA).multiply(0.041840)))

            # Sum the terms of the BRDF polynomial to get the instantaneous black-sky albedo
            bks_alb_x = BRDF_term_1.add(BRDF_term_2).add(BRDF_term_3)

            # Mask out any anomolous data (greater than 1 or less than 0)
            albedo_bks_mask = bks_alb_x.where(bks_alb_x.lte(0), 0)  # invalid
            albedo_bks_mask = albedo_bks_mask.where(bks_alb_x.gt(0), 1)  # valid 
            albedo_bks_mask = albedo_bks_mask.where(bks_alb_x.gte(1), 0)  # invalid
            bks_alb_x = bks_alb_x.updateMask(albedo_bks_mask)
            return bks_alb_x
        
        # Map function across all solar zenith angles to get global hourly black sky albedo values and
        # convert to image with same band names as the image containing the diffuse skylight ratios
        bands_hourly = ['b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10', 'b11', 'b12', 
                        'b13', 'b14', 'b15', 'b16', 'b17', 'b18', 'b19', 'b20', 'b21', 'b22', 'b23', 'b24']
        
        albedo_bks_x = cos_SZA_hourly.map(BKS_calc).toBands().rename(bands_hourly)
        
        # Mask out and anomolous data (greater than 1 or less than 0)
        albedo_ws_mask = albedo_ws_x.mask().where(albedo_ws_x.lte(0), 0)  # invalid
        albedo_ws_mask = albedo_ws_mask.where(albedo_ws_x.gte(1), 0)  # invalid
        albedo_ws_x = albedo_ws_x.updateMask(albedo_ws_mask)
        
        # Mask out high solar zenith angles
        SZA_noon = lat.sin().multiply(SDA.sin()).add(lat.cos().multiply(SDA.cos())).acos().multiply(180).divide(math.pi)
        SZA_mask = SZA_noon.mask().where(SZA_noon.gt(75),0)
        albedo_bks_x = albedo_bks_x.updateMask(SZA_mask)

        # Calculate blue-sky albedo values
        # Hourly blue-sky values calculate =s from the daily white-sky albedo and hourly black-sky values and hourly DSR values
        albedo_BLS_x = ((albedo_bks_x.multiply(dsr_inv_x)).add(albedo_ws_x.multiply(dsr_x)))
        
        # Determine weighting factors to calcualte daily-average blue sky values (weighted by the incoming solar radiation)
        # Total daily solar insolation
        tot_daily_solar = cos_SZA_hourly.toBands().updateMask(dsr_x.mask()).reduce(reducer=ee.Reducer.sum()).multiply(1361)
        # Hourly solar insolation divided by total daily solar to give fractional weighting for each hour of the day
        wtg_factor = cos_SZA_hourly.toBands().updateMask(dsr_x.mask()).rename(bands_hourly).multiply(1361).divide(tot_daily_solar)
        
        # Daily-mean blue-sky albedo values
        albedo_BLS_mean = albedo_BLS_x.multiply(wtg_factor).reduce(ee.Reducer.sum().unweighted())

        return albedo_BLS_mean

    # Calculate monthly average albedo values by mapping the daily-mean blue-sky albedo function across one month of data
    albedo_bls_m = albedo_BRDF.filterDate(start_date, end_date).filter(ee.Filter.calendarRange(m, m, 'month')).map(albedo_bls_daily).mean()

    # Isolate per land cover type by masking out all other land cover types
    lc_igbp_x = lc_data.select(band).filterDate(start_date, end_date).first()
    
    if lc == 'forest':
        lc_type_mask_x = lc_igbp_x.where(lc_igbp_x.neq(1), 0).where(lc_igbp_x.eq(1), 1).where(lc_igbp_x.eq(2), 1).where(lc_igbp_x.eq(3), 1).where(lc_igbp_x.eq(4), 1).where(lc_igbp_x.eq(5), 1)
        albedo_bls_m = albedo_bls_m.updateMask(lc_type_mask_x)
        
    elif lc == 11.27 or lc == 11.50:
        lc_type_mask_x = lc_igbp_x.where(lc_igbp_x.neq(11), 0).where(lc_igbp_x.eq(11), 1)
        albedo_bls_m = albedo_bls_m.updateMask(lc_type_mask_x)
        if lc == 11.27:
            lc_prop_x = lc_data.select('LC_Prop3').filterDate(start_date, end_date).first()
            lc_prop_type_mask_x = lc_prop_x.where(lc_igbp_x.neq(27), 0).where(lc_prop_x.eq(27), 1)
            albedo_bls_m = albedo_bls_m.updateMask(lc_prop_type_mask_x)
            
        if lc == 11.50:
            lc_prop_x = lc_data.select('LC_Prop3').filterDate(start_date, end_date).first()
            lc_prop_type_mask_x = lc_prop_x.where(lc_igbp_x.neq(50), 0).where(lc_prop_x.eq(50), 1)
            albedo_bls_m = albedo_bls_m.updateMask(lc_prop_type_mask_x)
            
    else:
        lc_type_mask_x = lc_igbp_x.where(lc_igbp_x.neq(lc), 0).where(lc_igbp_x.eq(lc), 1)
        albedo_bls_m = albedo_bls_m.updateMask(lc_type_mask_x)
        
    # Calculate and land cover masks
    # Calculate data quality mask by mapping masking function over the MODIS MCD12Q1 quality dataset corresponding to the same timeframe
    lc_qual_mask_x = lc_quality.filterDate(start_date, end_date).map(qual_mask).first()
    # Calculate land cover change mask by mapping masking function over the MODIS MCD12Q1 land cover dataset corresponding to the same timeframe
    lc_change_mask_x = lc_data.select('LC_Type1').filterDate(start_date, end_date).map(change_mask).first()
    # Calculate land cover spatial homogeneity mask by mapping masking function over the MODIS MCD12Q1 land cover dataset corresponding to the same timeframe
    lc_homog_mask_x = lc_data.select('LC_Type1').filterDate(start_date, end_date).map(homog_mask).first()
    
    # Apply land cover masks
    albedo_bls_m = albedo_bls_m.updateMask(lc_qual_mask_x).updateMask(lc_change_mask_x).updateMask(lc_homog_mask_x)

    # Apply spatial averaging, extending values by the distnace specified by d
    albedo_bls_m = albedo_bls_m.reduceNeighborhood(reducer=ee.Reducer.mean().combine(reducer2=ee.Reducer.variance(), sharedInputs=True),
                                                   kernel=ee.Kernel.circle(radius=d*1000,
                                                                           units='meters'),
                                                   skipMasked=False).reproject('EPSG:4326', scale=lc_scale).reduceRegions(collection=geom,
                                                                                                                          reducer=ee.Reducer.mean(),
                                                                                                                          scale=lc_scale,
                                                                                                                          tileScale=tile_scale)
    # Strings for file naming
    if band == 'LC_Type5':
        if lc == 7:
            lc = '12-7'
        if lc == 8:
            lc = '12-8'
    if band == 'LC_Prop2':
        if lc == 25:
            lc = '14-25'
        if lc == 35:
            lc = '14-35'
    if lc == 11.27:
        lc = '11-27'
    if lc == 11.50:
        lc = '11-50'

    band_names = ['01-Jan', '02-Feb', '03-Mar', '04-Apr', '05-May',
                  '06-Jun', '07-Jul', '08-Aug', '09-Sep', '10-Oct', '11-Nov', '12-Dec']

    # Export results
    ee.batch.Export.table.toDrive(collection=albedo_bls_m,
                                  description=str('albedo_bls_daily_lc-'+str(lc)+'_'+str(d)+'kms_'+str(yr)+'-'+band_names[m-1][:2]+sfx),
                                  fileNamePrefix=str('albedo_bls_daily_lc-'+str(lc)+'_'+str(d)+'kms_'+str(yr)+'-'+band_names[m-1][:2]+sfx),
                                  # fileFormat='SHP',
                                  selectors=['sum_mean','sum_variance', 'lat', 'lon']
                                  ).start()
    
#%% Define other parameters needed for calculation
# Define sampling distance
d = 25

# Define tile scaling
tile_scale = 16

# Define geometries
if output_scale == 0.5:
    geom = albedo_geom_05
    sfx = ''
elif output_scale == 0.05:
    geom = albedo_geom_005
    sfx = '_005'
else:
    print('Invalid scale selected')
    
# Split geometry into ten chunks
geom_a = geom.filter(ee.Filter.lt('lat', -17.5))#
geom_b = geom.filter(ee.Filter.lt('lat', 0)).filter(ee.Filter.gte('lat', -17.5))
geom_c = geom.filter(ee.Filter.lt('lat', 17.5)).filter(ee.Filter.gte('lat', 0))
geom_d = geom.filter(ee.Filter.lt('lat', 29.5)).filter(ee.Filter.gte('lat', 17.5))
geom_e = geom.filter(ee.Filter.lt('lat', 40)).filter(ee.Filter.gte('lat', 29.5))
geom_f = geom.filter(ee.Filter.lt('lat', 49)).filter(ee.Filter.gte('lat', 40))
geom_g = geom.filter(ee.Filter.lt('lat', 57)).filter(ee.Filter.gte('lat', 49))#
geom_h = geom.filter(ee.Filter.lt('lat', 64)).filter(ee.Filter.gte('lat', 57))#
geom_i = geom.filter(ee.Filter.lt('lat', 70)).filter(ee.Filter.gte('lat', 64))#
geom_j = geom.filter(ee.Filter.gte('lat', 70))#

geom_list = [geom_a, geom_b, geom_c, geom_d,geom_e, geom_f, geom_g, geom_h, geom_i, geom_j]
sfx_list = [sfx+'_a-10', sfx+'_b-10', sfx+'_c-10', sfx+'_d-10', sfx+'_e-10',sfx+'_f-10', sfx+'_g-10', sfx+'_h-10', sfx+'_i-10', sfx+'_j-10']

#%% Calculate results for all months in year specified by input parameters
lc_list = list(range(1,15)) + [16,11.27,11.50,12.7,12.8,14.25,14.35,'forest']
if lc in lc_list:
    for m in range(1, 13):
        for geom, sfx in zip(geom_list, sfx_list):
            albedo_bls_export(yr, d, lc, geom, m, tile_scale, sfx)
            print(m)
    else:
        print('Invalid land cover type selected')
