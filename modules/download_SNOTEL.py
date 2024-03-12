"""
Filename:    download_SNOTEL.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Functions to download SNOTEL data from hydroportal based on https://snowex-2021.hackweek.io/tutorials/geospatial/SNOTEL_query.html
"""

import os
from datetime import datetime
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import ulmo

def pull_snotel_sitecode(snotel_station_name):
    #This is the latest CUAHSI API endpoint
    wsdlurl = 'https://hydroportal.cuahsi.org/Snotel/cuahsi_1_1.asmx?WSDL'
    sites = ulmo.cuahsi.wof.get_sites(wsdlurl)
    ## put sites into pandas dataframe
    sites_df = pd.DataFrame.from_dict(sites, orient='index').dropna()
    ## create point geometry from lat and lon values
    sites_df['geometry'] = [Point(float(loc['longitude']), float(loc['latitude'])) for loc in sites_df['location']]
    
    ## clean the dataframe 
    sites_df = sites_df.drop(columns='location')
    sites_df = sites_df.astype({"elevation_m":float})
    
    ## create a geopandas dataframe from the pandas dataframe
    sites_gdf_all = gpd.GeoDataFrame(sites_df, crs='EPSG:4326')
    
    ## get the sitecode of the station name we are interested in
    idx = (sites_gdf_all['name'] == snotel_station_name)
    sitecode = sites_gdf_all.loc[idx].index.values[0]
    
    return sitecode

def snotel_fetch(sitecode, variablecode='SNOTEL:SNWD_D', start_date='1950-10-01', end_date='2023-12-31'):
    #print(sitecode, variablecode, start_date, end_date)
    wsdlurl = 'https://hydroportal.cuahsi.org/Snotel/cuahsi_1_1.asmx?WSDL'
    values_df = None
    try:
        #Request data from the server
        site_values = ulmo.cuahsi.wof.get_values(wsdlurl, sitecode, variablecode, start=start_date, end=end_date)
        #Convert to a Pandas DataFrame   
        values_df = pd.DataFrame.from_dict(site_values['values'])
        #Parse the datetime values to Pandas Timestamp objects
        values_df['datetime'] = pd.to_datetime(values_df['datetime'], utc=True)
        #Set the DataFrame index to the Timestamps
        values_df = values_df.set_index('datetime')
        #Convert values to float and replace -9999 nodata values with NaN
        values_df['value'] = pd.to_numeric(values_df['value']).replace(-9999, np.nan)
        #Remove any records flagged with lower quality
        values_df = values_df[values_df['quality_control_level_code'] == '1']
    except:
        print("Unable to fetch %s" % variablecode)

    return values_df