"""
Filename:    trajectory.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: class for calculating backwards trajectories with u, v, and w data
"""

import math
import xarray as xr
import pandas as pd
import numpy as np
import dask
from utils import  find_closest_MERRA2_lon_df, find_closest_MERRA2_lon, MERRA2_range, roundPartial


dask.config.set(**{'array.slicing.split_large_chunks': True})

class calculate_backward_trajectory:
    '''
    Returns an array 
    Parameters
    ----------
    ds : xarray
        xarray with list of dates, lats, and lons for extreme precip days
    idx : int
        integer of the index that is currently being processed
    start_lev: float
        starting level (in hPa) for trajectory
        
    Returns
    -------
    traj : pandas dataframe
        pandas dataframe of the backward trajectory for 72 hours
    
    '''
    
    def __init__(self, ds, event_date, start_lev=700., start_time=12, lat_offset=0, lon_offset=0):
    
        ## get center_date, start_lat, and start_lon
        ## center in the middle of the day
        self.center_date = ds.sel(date=event_date).date.values + np.timedelta64(start_time,'h')
        self.start_lat = ds.sel(date=event_date).lat.values + lat_offset
        self.start_lon = ds.sel(date=event_date).lon.values + lon_offset
        self.start_lev = start_lev
        print(self.center_date, self.start_lat, self.start_lon)
        
        self.varlst = ['time', 'latitude', 'longitude', 'level', 'q', 'u', 'v', 'w']     
        self.calc_vars = ['drying_ratio', 'dq']
        
        self.date_lst = pd.date_range(end=self.center_date, periods=72, freq='H')

        ## create list of dates based on start date
        start_date = ds.sel(date=event_date).date.values - np.timedelta64(3,'D')
        end_date = ds.sel(date=event_date).date.values
        self.date_lst_era = pd.date_range(start_date, end_date, freq='1D')

    def preprocess(self, ds):
        return ds.sel(latitude=slice(50., 15.), longitude=slice(-180., -80.))
    
    def read_data(self):    
        # read ERA5 data
        # create list of daily ERA5 files for each AR
        # ds_lst = []
        fname_lst = []
        for j, date in enumerate(self.date_lst_era):
            year = date.year
            month = date.strftime("%m")
            day = date.strftime("%d")
            
            path_to_data = '/expanse/nfs/cw3e/cwp140/downloads/ERA5/ERA5/{0}/'.format(year)
            fname = "era5_nhemi_025dg_1hr_uvwq_{0}{1}{2}.nc".format(year, month, day)
            fname_lst.append(path_to_data+fname)
            # ds2 = xr.open_dataset(path_to_data+fname)
            # ds_lst.append(ds2)
        
        # concatenate ds_lst
        # self.ds1 = xr.concat(ds_lst, dim='time')

        
        self.ds1 = xr.open_mfdataset(fname_lst, engine='netcdf4', combine='by_coords')

        ## read MERRA2 data - not sure if we need this anymore
        # fname = '/data/downloaded/Reanalysis/MERRA2/M2I3NPASM.5.12.4_raw/1980/MERRA2_100.inst3_3d_asm_Np.19801231.nc4'
        # calculate vertical velocity (w) if MERRA2
        # w = mpcalc.vertical_velocity(ds.OMEGA, ds.lev, ds.T, mixing_ratio=0)

    def create_empty_array(self):   
   
        # initial conditions
        t0 = self.ds1.interp(latitude=self.start_lat, longitude=self.start_lon, level=self.start_lev, time=self.center_date)
        
        ## append initial conditions to empty DataFrame
        t0_vals = []
        for i in self.varlst:
            if i == 'time':
                t0_vals.append(t0[i].values)
            else:
                t0_vals.append(float(t0[i].values))
        df = pd.DataFrame(columns = self.varlst, index=np.arange(0, 72, 1))
        df.iloc[0] = t0_vals
        self.df = df
        
    def find_distance_travelled(self, t0):
    
        ## not sure what this is
        nRes = 1.
        nHour = 24.
        nXDay = 4
        timeRatio = ((nHour/nXDay)/nRes)

        # t1 is where the parcel was 1 hour ago if it had been travelling at the speed of t0
        del_x = ((0-t0.u)*3600.)/1000. # convert to km
        del_y = ((0-t0.v)*3600.)/1000. # convert to km 
        del_z = ((0-t0.w)*3600.)/100. # convert to hPa
    
        return del_x, del_y, del_z
    
    def km_to_decimal_degrees(self, distx, disty, latitude):
        '''
        # Example usage:
        disty = 100 ## distance along latitude in km
        distx = 100 ## distance along longitude in km
        latitude = 52.5200 ## latitude to calculate along

        xdeg, ydeg = km_to_decimal_degrees(distx, disty, latitude)
        '''
        ## convert distance north-south to decimal degrees
        ## 111.2 km per degree
        ydeg = disty/111.2 # convert km to degree
        
        ## convert distance east-west to decimal degrees
        # Radius of the Earth in kilometers
        earth_radius = 6371

        # Convert latitude from degrees to radians
        lat_radians = math.radians(latitude)

        # Calculate decimal degrees
        xdeg = distx / (earth_radius * (math.pi / 180))

        return xdeg, ydeg

    
    def get_values_at_current_timestep(self, idx):
        '''
        Given u, v, and w at current lat/lon/lev, 
        what is u/v/w at previous time step
        
        idx: int
            index of current time step
        '''
        # get values of previous time step
        t0 = self.df.iloc[idx-1]
        
        ## find distance travelled between this hour and previous hour
        del_x, del_y, del_z = self.find_distance_travelled(t0)

        ## convert km to degrees
        del_x, del_y = self.km_to_decimal_degrees(del_x, del_y, t0.latitude)

        ## now the new location is the old location plus the distances
        new_lat = t0.latitude + del_y
        new_lon = t0.longitude + del_x
        new_lev = t0.level + del_z
        new_date = t0.time - np.timedelta64(1,'h')

        ## interpolate to new point
        t1 = self.ds1.interp(latitude=new_lat, longitude=new_lon, level=new_lev, time=new_date)
        ## put new values in dataframe
        t1_vals = []
        for i in self.varlst:
            if i == 'time':
                t1_vals.append(t1[i].values)
            else:
                t1_vals.append(float(t1[i].values))
        self.df.iloc[idx] = t1_vals

        return self.df
        
    def compute_trajectory(self):
        
        ## read data TODO flexible for date input
        print('Reading ERA5 data...')
        self.read_data()
        
        ## build dataframe with initial conditions
        print('Calculating initial conditions...')
        self.create_empty_array()
        
        ## loop through the remaining 72 hours
        print('Calculating trajectory ...')
        for i, idx in enumerate(np.arange(1, 72, 1)):
            print('... time step {0}'.format(idx))
            self.df = self.get_values_at_current_timestep(idx)
            
            ## write break statements
            if self.df.iloc[idx]['latitude'] <= 0.:
                break
            elif self.df.iloc[idx]['latitude'] >= 90.:
                break
            if self.df.iloc[idx]['longitude'] < -180.:
                self.df.iloc[idx]['longitude'] = self.df.iloc[idx]['longitude'] + 360.
            if self.df.iloc[idx]['level'] < 1.:
                break
            elif self.df.iloc[idx]['level'] > 1000.:
                self.df.iloc[idx]['level'] = 1000.
        
        ## convert specific humidity to g kg-1
        self.df['q'] = self.df['q']* 1000
        ## calculate change in q between this step and previous step
        self.df['dq'] = self.df['q'].diff()
        ## using del_q, calculate "drying ratio"
        ## how much water vapor has been lost/gained in between time steps
        self.df['drying_ratio'] = (self.df['dq']/self.df['q'])*100.# convert to %    
        
        return self.df
    
def combine_IVT_and_trajectory(ERA5):
    ## load ERA5 IVT data
    start_date = ERA5.start_date.values - np.timedelta64(3,'D')
    end_date = ERA5.start_date.values
    print(start_date, end_date)

    dates = pd.date_range(start=start_date, end=end_date, freq='1D')
    # put into pandas df
    d ={"date": dates}
    df = pd.DataFrame(data=d)
    df['day']= df['date'].dt.day.map("{:02}".format)
    df['month']= df['date'].dt.month.map("{:02}".format)
    df['year']= df['date'].dt.year

    # create list of daily ERA5 files
    filenames = []
    for j, row in df.iterrows():
        filenames.append('/data/downloaded/Reanalysis/ERA5/IVT/{0}/ERA5_IVT_{0}{1}{2}.nc'.format(row['year'], row['month'], row['day']))
        # open all files within the AR period

    ivt = xr.open_mfdataset(filenames, combine='by_coords', parallel=False)


    ## interpolate IVT to trajectory points
    ivt = ivt.interp(lat=ERA5.lat, lon=ERA5.lon, time=ERA5.time)
    ivt = ivt.compute()

    ## merge IVT, uIVT, vIVT, and IWV to trajectory ds
    ERA5 = xr.merge([ivt, ERA5])
    
    return ERA5

def combine_arscale_and_trajectory(ERA5, arscale, ar):
    ## create a list of lat/lons that match MERRA2 spacing
    ## lat and lon points from trajectory
    new_lst = []
    for lon in ERA5.lon.values:
        new_lst.append(find_closest_MERRA2_lon(lon))

    t = xr.DataArray(ERA5.time.values, dims=['location'], name='time') 
    x = xr.DataArray(new_lst, dims=['location'])
    y = xr.DataArray(roundPartial(ERA5.lat.values, 0.5), dims=['location'])

    x = xr.DataArray(ERA5.lon.values, dims=("location"), coords={"lon": x}, name='traj_lons')
    y = xr.DataArray(ERA5.lat.values, dims=("location"), coords={"lat": y}, name='traj_lats')

    # create a new dataset that has the trajectory lat and lons and the closest MERRA2 lat/lons as coords
    z = xr.merge([x, y, t])
    
    ## Open text file with coordinates of coastal region along N. America West Coast
    textpts_fname = '../data/latlon_coast-modified.txt'
    df = pd.read_csv(textpts_fname, header=None, sep=' ', names=['latitude', 'longitude'], engine='python')
    df['longitude'] = df['longitude']*-1

    ## create column with closest MERRA2 lons
    df['MERRA2_lon'] = df.apply(lambda row: find_closest_MERRA2_lon_df(row), axis=1)

    d = {'lat' : df['latitude'],
        'lon' : df['MERRA2_lon']}

    txtpts = pd.DataFrame(d)
    txtpts = txtpts.drop_duplicates(subset=['lat', 'lon'])

    ## Now loop through the lat/lon pairs and see where they match
    idx_lst = []
    for i, (x, y) in enumerate(zip(z.lon.values, z.lat.values)):
        for j, (lon, lat) in enumerate(zip(txtpts.lon.values, txtpts.lat.values)):
            ## test if lat/lon pair matches
            result_variable = (x == lon) & (y == lat)

            if (result_variable == True):
                idx = (i, j)
                idx_lst.append(idx)

    if len(idx_lst) > 0:
        ## take first time the trajectory crosses the coast
        idx = idx_lst[0]
        ## this is the time of the trajectory when it crosses west coast
        time_match = z.sel(location=idx[0]).time.values
        
        # get the location index value where the lat/lon matches the coastal intersection value
        idx_ds = int(arscale.location.where((arscale.lat==txtpts.iloc[idx[1]].lat) & (arscale.lon==txtpts.iloc[idx[1]].lon), drop=True).values)
        
        #####################
        ### STRICT METHOD ###
        #####################
        
        ## Gather arscale of closest grid and time value
        tmp1 = arscale.sel(location=idx_ds)
        arscale_val = tmp1.sel(time=time_match, method='nearest').ar_scale.values
        ERA5 = ERA5.assign(ar_scale_strict=arscale_val)
        
        ## Gather Rutz AR value of closest grid and time value
        tmp1 = ar.sel(location=idx_ds)
        ar_val = tmp1.sel(time=time_match, method='nearest').AR.values
        ERA5 = ERA5.assign(ar_strict=ar_val)
        
        coastal_IVT_val = tmp1.sel(time=time_match, method='nearest').IVT.values
        ERA5 = ERA5.assign(coastal_IVT_strict=coastal_IVT_val)
        
        
        #######################
        ### FLEXIBLE METHOD ###
        #######################
        
        ## select the surrounding grid points
        tmp = arscale.sel(location=slice(idx_ds-2, idx_ds+3))

        ## select the 12 hours on each side of the time step
        sta = time_match - np.timedelta64(12,'h')
        sto = time_match + np.timedelta64(12,'h')
        arscale_val = tmp.sel(time=slice(sta, sto)).ar_scale.values.max()
        ## now put those values into the trajectory dataset
        ERA5 = ERA5.assign(ar_scale=arscale_val)
        
        ## lets also grab whether rutz et al AR was detected
        tmp1 = ar.sel(location=slice(idx_ds-2, idx_ds+3))
        ar_val = tmp1.sel(time=slice(sta, sto)).AR.values.max()
        
        coastal_IVT_val = tmp1.sel(time=slice(sta, sto)).IVT.values.max()

        ## assign value to trajectory dataset
        ERA5 = ERA5.assign(ar=ar_val)
        ERA5 = ERA5.assign(coastal_IVT=coastal_IVT_val)
        
    else:
        ## since the trajectory didn't cross the west coast, set ar_scale to nan
        ERA5 = ERA5.assign(ar_scale=np.nan)
        ERA5 = ERA5.assign(ar=np.nan)
        ERA5 = ERA5.assign(coastal_IVT=np.nan)
        
        ERA5 = ERA5.assign(ar_scale_strict=np.nan)
        ERA5 = ERA5.assign(ar_strict=np.nan)
        ERA5 = ERA5.assign(coastal_IVT_strict=np.nan)

    return ERA5