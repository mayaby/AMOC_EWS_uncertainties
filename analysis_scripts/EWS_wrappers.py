## Because of the large number of surrogates/members and for a large number of gridcells, the EWS were calculated using parallel codes. An example of such a code is attached, but the exact details of these scripts are not important - all the actual analysis is in these wrappers and the functions they use from EWS_functions.py. Just apply the appropriate wrapper to every timeseries. ##

import xarray as xr
from EWS_functions import *
import numpy as np

## For true SST and salinity data the running mean is removed before calculating the EWS ("rmean" functions)

def lambda_wrapper_rmean(amoc,ws=60,rws=50):
    if np.count_nonzero(np.isnan(amoc))==0:
        if len(np.where(amoc==0)[0])<=100:
            try:
                amoc = np.nan_to_num(amoc)
                amoc_low = runmean(amoc, rws)
                if amoc.sum() != 0:
                    lamb = run_fit_a_ar1((amoc-amoc_low),ws)
                else:
                    lamb = np.full(len(amoc),np.nan)
            except:
                print('failed to get lambda')
                print(amoc-amoc_low)
                lamb = np.full(len(amoc),np.nan)
        else:
            lamb = np.full(len(amoc),np.nan)
    else:
        lamb = np.full(len(amoc),np.nan)
    return lamb

def ar1_wrapper_rmean(amoc,ws=60,rws=50):
    if np.count_nonzero(np.isnan(amoc))==0:
        try:
            amoc = np.nan_to_num(amoc)
            amoc_low = runmean(amoc, rws)
            if amoc.sum() != 0:
                ar1 = runac2((amoc - amoc_low), ws)
        except:
            print('failed to get ar1')
            print(amoc-amoc_low)
            ar1 = np.full(len(amoc),np.nan)
    else:
        ar1 = np.full(len(amoc),np.nan)
    return ar1

def std_wrapper_rmean(amoc,ws=60,rws=50):
    if np.count_nonzero(np.isnan(amoc))==0:
        try:
            amoc = np.nan_to_num(amoc)
            amoc_low = runmean(amoc, rws)
            if amoc.sum() != 0:
                std = runstd((amoc - amoc_low), ws)
        except:
            print('failed to get std')
            print(amoc-amoc_low)
            std = np.full(len(amoc),np.nan)
    else:
        std = np.full(len(amoc),np.nan)
    return std

## The surrogates are generated from the detrended timeseries, so there is no need to detrend them. The EWS are thus calculated straight from the timeseries.

def lambda_wrapper_ndt(amoc,ws=60):
    if np.count_nonzero(np.isnan(amoc))==0:
        if len(np.where(amoc==0)[0])<=100:
            try:
                amoc = np.nan_to_num(amoc)
                if amoc.sum() != 0:
                    lamb = run_fit_a_ar1((amoc),ws)
                else:
                    lamb = np.full(len(amoc),np.nan)
            except:
                print('failed to get lambda')
                print(amoc)
                lamb = np.full(len(amoc),np.nan)
        else:
            lamb = np.full(len(amoc),np.nan)
    else:
        lamb = np.full(len(amoc),np.nan)
    return lamb

def ar1_wrapper_ndt(amoc,ws=60):
    if np.count_nonzero(np.isnan(amoc))==0:
        try:
            amoc = np.nan_to_num(amoc)
            if amoc.sum() != 0:
                ar1 = runac2((amoc), ws)
            else:
                ar1 = np.full(len(amoc),np.nan)
        except:
            print('failed to get ar1')
            print(amoc)
            ar1 = np.full(len(amoc),np.nan)
    else:
        ar1 = np.full(len(amoc),np.nan)
    return ar1

def std_wrapper_ndt(amoc,ws=60):
    if np.count_nonzero(np.isnan(amoc))==0:
        try:
            amoc = np.nan_to_num(amoc)
            if amoc.sum() != 0:
                std = runstd((amoc), ws)
            else:
                std = np.full(len(amoc),np.nan)
        except:
            print('failed to get std')
            print(amoc)
            std = np.full(len(amoc),np.nan)
    else:
        std = np.full(len(amoc),np.nan)
    return std