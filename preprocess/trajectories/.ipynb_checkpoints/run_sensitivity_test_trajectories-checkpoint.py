######################################################################
# Filename:    calculate_trajectories.py
# Author:      Deanna Nash dnash@ucsd.edu
# Description: Script to run backwards trajectories for CO PRISM watershed >90th percentile precipitation events
#
######################################################################

import os, sys
import yaml
import xarray as xr
import pandas as pd

path_to_repo = '/home/dnash/repos/UCRB_WY2023_ARs/'
sys.path.append(path_to_repo+'modules')
from trajectory import calculate_backward_trajectory

path_to_data = '/expanse/lustre/scratch/dnash/temp_project/'

## get HUC8 from config file
config_file = str(sys.argv[1]) # this is the config file name
job_info = str(sys.argv[2]) # this is the job name

config = yaml.load(open(config_file), Loader=yaml.SafeLoader) # read the file
ddict = config[job_info] # pull the job info from the dict

HUC8_ID = str(ddict['HUC8_ID'])
event_date = ddict['date']
if ddict['grid'] == 'center':
    lat_offset = 0
    lon_offset = 0
elif ddict['grid'] == 'north':
    lat_offset = 0.25
    lon_offset = 0
elif ddict['grid'] == 'south':
    lat_offset = -0.25
    lon_offset = 0
elif ddict['grid'] == 'east':
    lat_offset = 0
    lon_offset = 0.25
elif ddict['grid'] == 'west':
    lat_offset = 0
    lon_offset = -0.25

h = ddict['hour']
lev_t = ddict['lev']

## set starting lat/lon
## choose this based on extreme precip days
fname = path_to_data + 'preprocessed/PRISM/PRISM_HUC8_CO.nc'
ds = xr.open_dataset(fname)
# start with single event from single watershed
ds = ds.sel(HUC8=HUC8_ID)

## calculate trajectory for current 
s = calculate_backward_trajectory(ds=ds, event_date=event_date, start_lev=lev_t, start_time=h, lat_offset=lat_offset, lon_offset=lon_offset)
df = s.compute_trajectory()
new_ds = df.to_xarray()

## add in the coordinate information so we can combine the trajectories easily
## pull start_date from first time step
start_date = new_ds.time.isel(index=0).values

new_coords = {"start_lev": [lev_t],
              "grid": [ddict['grid']],
              "start_date": [start_date]}
new_ds = new_ds.expand_dims(dim=new_coords)

## save trajectory data to netCDF file
print('Writing {0} {1} to netCDF ....'.format(HUC8_ID, event_date))
out_fname = path_to_data + 'preprocessed/UCRB_trajectories/PRISM_HUC8_{0}_{1}_{2}_{3}Z_{4}hPa.nc'.format(HUC8_ID, event_date, ddict['grid'], h, lev_t) 
new_ds.load().to_netcdf(path=out_fname, mode = 'w', format='NETCDF4')
