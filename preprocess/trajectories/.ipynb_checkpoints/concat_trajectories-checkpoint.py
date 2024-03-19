######################################################################
# Filename:    concat_trajectories.py
# Author:      Deanna Nash dnash@ucsd.edu
# Description: Script to concatenate trajectories from each watershed
#
######################################################################

import sys
import xarray as xr
import numpy as np
import pandas as pd
import glob


path_to_data = '/expanse/lustre/scratch/dnash/temp_project/'
HUC8_lst = ['14010001', '14080101']

for j, HUCID in enumerate(HUC8_lst):
    fname_pattern = '/expanse/lustre/scratch/dnash/temp_project/preprocessed/ERA5_sensitivity_test_trajectories/PRISM_HUC8_{0}*.nc'.format(HUCID)
    fname_lst = glob.glob(fname_pattern)
    # print(fname_lst)
    ds_lst = []
    for i, fname in enumerate(fname_lst):
        ds = xr.open_dataset(fname)
        ds_lst.append(ds)
    
    ds = xr.merge(ds_lst)
    
    out_fname = '/expanse/nfs/cw3e/cwp140/preprocessed/UCRB_trajectories/PRISM_HUC8_{0}.nc'.format(HUCID)
    ds.to_netcdf(path=out_fname, mode = 'w', format='NETCDF4')