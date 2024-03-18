######################################################################
# Filename:    create_job_configs.py
# Author:      Deanna Nash dnash@ucsd.edu
# Description: Script to create .yaml configuration file to run job_array with slurm for running ERA5 CO PRISM HUC8 trajectories
#
######################################################################

import yaml
from itertools import chain
import itertools
import xarray as xr
import pandas as pd

## get list of HUC8 trajectories
path_to_data = '/expanse/lustre/scratch/dnash/temp_project/'
fname = path_to_data + 'preprocessed/PRISM/PRISM_HUC8_CO.nc'
ds = xr.open_dataset(fname)

HUC8_lst = ['14010001', '14080101']
name_lst = ['colorado_headwaters', 'upper_san_juan']
HUC8_final = []
time_lst = []
for i, HUC8 in enumerate(HUC8_lst):
    ## load list of dates
    df = pd.read_csv('../../data/SNOTEL_80th_perc_WY2023_{0}'.format(name_lst[i]))
    time = df.time.values
    for j, t in enumerate(time):
        HUC8_final.append(HUC8)
        time_lst.append(t)


jobcounter = 0
filecounter = 0
## loop through to create dictionary for each job
d_lst = []
dest_lst = []
njob_lst = []
for i, (HUC8, time) in enumerate(zip(HUC8_final, time_lst)):
    jobcounter += 1
    d = {'job_{0}'.format(jobcounter):
         {'HUC8_ID': HUC8,
          'date': time,
          'lev': 650.,
          'hour': 12,
          'grid': 'center'}}
    d_lst.append(d)
    
    if (jobcounter == 999):
        filecounter += 1
        ## merge all the dictionaries to one
        dest = dict(chain.from_iterable(map(dict.items, d_lst)))
        njob_lst.append(len(d_lst))
        ## write to .yaml file and close
        file=open("config_{0}.yaml".format(str(filecounter)),"w")
        yaml.dump(dest,file, allow_unicode=True, default_flow_style=None)
        file.close()
        
        ## reset jobcounter and d_lst
        jobcounter = 0
        d_lst = []
        
## now save the final config
filecounter += 1
## merge all the dictionaries to one
dest = dict(chain.from_iterable(map(dict.items, d_lst)))
njob_lst.append(len(d_lst))
## write to .yaml file and close
file=open("config_{0}.yaml".format(str(filecounter)),"w")
yaml.dump(dest,file, allow_unicode=True, default_flow_style=None)
file.close()

## create calls.txt for config_1(-8)

for i, njobs in enumerate(njob_lst):
    call_str_lst = []
    for j, job in enumerate(range(1, njobs+1, 1)):
        call_string = "python run_sensitivity_test_trajectories.py config_{0}.yaml 'job_{1}'".format(i+1, j+1)
        call_str_lst.append(call_string)
        
    ## now write those lines to a text file
    with open('calls_{0}.txt'.format(i+1), 'w',encoding='utf-8') as f:
        for line in call_str_lst:
            f.write(line)
            f.write('\n')
        f.close()