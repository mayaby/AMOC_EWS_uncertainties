"""
This script contains the wrapper functions to:
- generate fourier surrogates
- generate AR2 surrogates
- modify a surrogate array using the salinity weights method described in the paper

In general surrogates can be generated in a parallel script similar to that in the example, using the wrappers and functions from this script. 

The process for the AR2 surrogates is described in detail in the paper - we basically first do a fit to get the correct a1 and a2 coefficients and then generate many AR2 surrogates using these (the ar2_wrapper needs to be run in a loop, as opposed to the fourier wrapper).

"""

import xarray as xr
from EWS_functions import *
from statsmodels.tsa.arima_process import ArmaProcess
from statsmodels.tsa.seasonal import seasonal_decompose
import scipy.optimize as optimize

def ymean(x):
    return np.mean(x.reshape(-1, 12), axis=1)

## FOURIER SURROGATES

def fourrier_surrogates(ts, ns):
    ## generate ns fourier surrogates from timeseries ts
    ts_fourier  = np.fft.rfft(ts)
    random_phases = np.exp(np.random.uniform(0, 2 * np.pi, (ns, ts.shape[0] // 2 + 1)) * 1.0j)
    ts_fourier_new = ts_fourier * random_phases
    new_ts = np.real(np.fft.irfft(ts_fourier_new))
    return new_ts

def fourier_wrapper(ts,nsmp):
    ts = ts - runmean(ts,w=50)
    return np.transpose(fourrier_surrogates(ts,nsmp))

## AR2 SURROGATES

def make_ar2(a1,a2,nt=2066):
    ar1 = np.array([1,-a1,-a2])
    ma1 = np.array([1])
    AR_object1 = ArmaProcess(ar1, ma1)
    return AR_object1.generate_sample(nsample=nt)

def g(A,a2):
    ## given the true ar2 and the estimated monthly AR(1) we can calculate the true ar1 analytically
    return A*(1-a2)

def F(a2,Am, Ay,ny=500*12):
    ## the function to minimize to get the right measured yearly AR(1) Ay (the measured monthly AR(1) Am is constrained by function g)
    a1 = g(Am,a2)
    x = make_ar2(a1[0],a2[0],nt=12*ny)
    Ay2 = sm.tsa.acf(ymean(x), nlags = 1)[1]
    return (Ay-Ay2)**2


def fit_ar2(ts):
    ## Function to fit an ar2 timeseries to ts. The criteria is that both the monthly and the yearly average measured AR(1) need to match that of ts. See paper for details.  

    ts = ts - ts.mean() ## this is NOT nanmean, so if there is any nan values the whole time series becomes nan
    ts = ts - runmean(ts,w=240) ## remove running mean with window of 20 years
    if np.count_nonzero(np.isnan(ts))==0:
        result=seasonal_decompose(ts,model='additive', period=12)
        ts2 = ts-result.seasonal ## remove the seasonal part
        Am2 = sm.tsa.acf(ts2, nlags = 1)[1] ## get the monthly AR(1)
        Ay2 = sm.tsa.acf(ymean(ts2), nlags = 1)[1] ## get the yearly mean AR(1)
        initial_guess = sm.tsa.pacf(ts, nlags = 2)[2]
        if initial_guess<0:
            initial_guess = 0
        results = optimize.minimize(F, initial_guess,args=(Am2,Ay2),method='Powell',bounds=[[0,1]]) 
        a2 = results.x
        a1 = g(Am2,a2)
        return a1, a2
    else:
        return np.nan, np.nan

def ar2_wrapper(ts):
    ## run this wrapper in a loop to generate the desired amount of AR2 surrogates
    a1, a2 = fit_ar2(ts)
    s = make_ar2(a1,a2)
    s = (s/np.nanstd(s))*np.nanstd(ts)
    return s

## MODIFY SALINITY SURROGATES

def modify_wrapper(smp):
    ## modification of salinity surrogates
    ## done on the whole surrogate array as a whole: smp should have time, lat, lon dimensions. 

    dmean_wei = xr.open_dataset('/p/tmp/mayayami/sst-uncertainties-data/ENS4_analyses/weights_dmean_at.nc').salinity_observation_weights.squeeze()[:-5]
    seas = xr.open_dataset('/p/tmp/mayayami/sst-uncertainties-data/ENS4_analyses/salinity_seasonal_cycle.nc').seas

    nt = smp.shape[0]

    wei_lim = 0.5
    fe = np.copy(smp)
    fe = fe + seas.values

    w = dmean_wei.values

    clim = np.transpose(np.tile(np.transpose(smp.values[0:12]),122))
    clim_an = clim - clim.mean(axis=0)

    ## iterate though times and modify the values where the weight is below the limit
    for it in np.arange(nt):
        idc = np.where(w[it]<=wei_lim) 
        if it==0:
            replc = clim_an[it]
        else: 
            replc=0.9*(fe[it-1]-clim_an[it-1])+clim_an[it]
        fe[it][idc]=replc[idc]

    ## make into xarray and take annual mean
    da = xr.DataArray(
        data=fe,
        dims=["time", "lat", "lon"],
        coords=dict(
            lon=dmean_wei.lon,
            lat=dmean_wei.lat,
            time=dmean_wei.time
        ))
    da = da.groupby('time.year').mean('time')
    return da

