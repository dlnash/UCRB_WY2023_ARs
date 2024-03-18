#!/bin/bash
######################################################################
# Filename:    download_ERA5.sh
# Author:      Deanna Nash dnash@ucsd.edu
# Description: Script to download all the necessary ERA5 files
# To run, activate conda env with "conda activate cds", then run the script with "bash download_ERA5.sh"
######################################################################

### Activate python

# source /home/dnash/miniconda3/etc/profile.d/conda.sh 

### Activate conda env

conda activate cds

# names of configuration dictionaries to loop through
array=(
sp # 12 hourly surface pressure over CONUS for WY 2023
)

# now loop through each configuration dictionary to download the ERA5 data

##TODO figure out how to input argument into python run file

for i in ${!array[*]}
do 
    inconfig="${array[$i]}"
    python getERA5_batch.py ${inconfig}
done