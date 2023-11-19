"""
This is a script that uses xarray to calculate the linear trend of the EWS timeseries. There's an extra "set_time" option to set the years as integers, because some netcdf files have the time in a format that confuses xarray and it outputs a trend that isn't year^-1. 

If desired there is also a second function to fit a gaussian distribution to the trends across some chosen dimension, and output a file with the mu and sigma of that distribution.  
"""

import xarray as xr
from scipy.stats import linregress
import scipy.stats as st
import numpy as np

def get_mu(x):
    if not np.isnan(np.sum(x.values)):
        (mu, sigma) = st.norm.fit(x)
        return xr.DataArray(mu)
    else:
        return xr.DataArray(np.nan)

def get_sigma(x):
    if not np.isnan(np.sum(x.values)):
        (mu, sigma) = st.norm.fit(x)
        return xr.DataArray(sigma)
    else:
        return xr.DataArray(np.nan)

def make_trends(path, im=None, set_time=False,times=np.arange(1870,2020)):
    data = xr.open_dataset(path +'.nc').samples
    if set_time:
        data['time']=times
    mid_times = data.time.values[30:-30]
    data = data.where(data.time.isin(mid_times))
    result = data.polyfit(dim = 'time', deg = 1)
    if im:
        ## sometimes this is a good step to add a member dimension to the file for easy concatenation later
        result = result.expand_dims({'member':[im]})

    result.to_netcdf(path+'_trends.nc')
    print('done '+path+'_trends.nc')

def make_stats(path,mdims=['latitude','longitude']):
    ds=xr.open_dataset(path+'_trends.nc')
    trends = ds.polyfit_coefficients.sel(degree=1)
    stacked = trends.where(trends,drop=True).stack(allpoints=mdims)
    mu = stacked.groupby('allpoints').apply(get_mu)
    mu_unstacked = mu.unstack('allpoints')
    stacked = trends.where(trends,drop=True).stack(allpoints=mdims)
    sigma = stacked.groupby('allpoints').apply(get_sigma)
    sigma_unstacked = sigma.unstack('allpoints')

    sigma = sigma_unstacked.rename('sigma')
    mu = mu_unstacked.rename('mu')
    array = xr.merge((sigma, mu))

    array.to_netcdf(path+'_stats.nc')
    print('done '+path+'_stats.nc')

## EXAMPLES
path = 'HadCRUT.5.0.1.0.analysis.anomalies_full_ar1s2_w60_runmean'
make_trends(path)
make_stats(path,mdims=['latitude','longitude'])

path='sst2d.ano.1854.2017.ensemble_full_at2_stds_w60_rmean'
make_trends(path, set_time=True,times=np.arange(1854,2018))
make_stats(path,mdims=['lat','lon'])