"""
To plot the cross-hatched areas in the figures I use matplotlib's contourf. The data is simple: every coordinate that is e.g. statistically significant has a value, otherwise that data point is NaN. The problem is that contourf sees the lat/lon coordinates as the points on the intersection of grid lines, and so it will only fill in a square with a hatch if all four of its corners are not NaN. 

(Note that that's in the corner_mask=False case, if True it can fill in half of the gridcell. See https://matplotlib.org/2.0.2/examples/pylab_examples/contour_corner_mask.html for the difference). 

This is not what we want for the figure: because we are interested in all specific locations that are statistically significant we would want to plot a value even if it's the only one in the area. So we want to have hatches in a square around every lat/lon coordinate that isn't NaN. 

The way I get around the contourf method is by creating a new grid that is shifted half a lon/lat step from the true grid. Then I can run a loop through the data and make sure that all the non-NaN values have non-NaN values in all their four corners.

Below you can see the functional code for the regions in which the HadISST lambda trend is statistically significant. Some of the partial codes for the other datasets are commented out below that.
"""

import xarray as xr
import numpy as np
import scipy.stats as st

## GETTING THE SIGNIFICANCE AREAS

def f(x):
    if not np.isnan(x[-1]):
        x=x[~np.isnan(x)]
        return st.percentileofscore(x, x[-1])
    else:
        return np.nan

Hlams = xr.open_dataset('/p/tmp/mayayami/sst-uncertainties-data/obs_data/HadISST_lambdas_w60_rmean.nc').lams
Hlams['time']=np.arange(1870,2021)
Hmean_trend = Hlams[:,30:-30].polyfit(dim = "time", deg = 1).polyfit_coefficients.sel(degree=1)
ds = xr.open_dataset('/p/tmp/mayayami/sst-uncertainties-data/obs_data/fsurrogates_HadISST_dt50_lams_w60_ndt_100_trends.nc')
trends = ds.lams_polyfit_coefficients.sel(degree=1)
d = Hmean_trend.expand_dims(dim='smps', axis=0)
d['smps']=[100]
trend_arr = xr.concat((trends,d),dim='smps')

trend_pvs = 1- np.apply_along_axis(f, 2, trend_arr)/100
pv_array = xr.Dataset(
        data_vars = dict(pv=(['longitude','latitude'],trend_pvs)),
        coords = dict(
                lat = d.latitude,
        lon = d.longitude))
Hlam_significant = pv_array.pv.where(pv_array.pv<=0.05)

# SHIFTING THE GRID
Hlats = Hlam_significant.latitude
Hlons = Hlam_significant.longitude
shifted_sig = np.full((len(Hlats),len(Hlons)),np.nan)
for i, lat in enumerate(Hlats):
    for j, lon in enumerate(Hlons):
        x = Hlam_significant.sel(latitude=lat,longitude=lon).values
        if not np.isnan(x):
            if (i!=179) & (j!=359):
                shifted_sig[i,j]=x
                shifted_sig[i,j+1]=x
                shifted_sig[i+1,j]=x
                shifted_sig[i+1,j+1]=x
            else:
                if (i==179) & (j!=359):
                    shifted_sig[0,j]=x
                    shifted_sig[0,j+1]=x
                elif (j==359)&(i!=179):
                    shifted_sig[i,0]=x
                    shifted_sig[i+1,0]=x
                else:
                    print(i,j)
                    shifted_sig[0,0]=x
Hsig = xr.Dataset(
                data_vars = dict(pv=(['latitude','longitude'],shifted_sig)),
                coords = dict(
                        latitude = Hlats-0.5,
                longitude = Hlons-0.5))


## FOR OTHER DATASETS 
## (note the below code won't work without some of the code from Main Plots.ipynb)

## HADCRUT ENSEMBLE SPREAD
# Clats = esignificant.latitude
# Clons = esignificant.longitude
# shifted_sig = np.full((len(Clats),len(Clons)),np.nan)
# for i, lat in enumerate(Clats):
#     for j, lon in enumerate(Clons):
#         x = esignificant.sel(latitude=lat,longitude=lon).values
#         if not np.isnan(x):
#             if (i!=35) & (j!=71):
#                 shifted_sig[i,j]=x
#                 shifted_sig[i,j+1]=x
#                 shifted_sig[i+1,j]=x
#                 shifted_sig[i+1,j+1]=x
#             else:
#                 if (i==35) & (j!=71):
#                     shifted_sig[0,j]=x
#                     shifted_sig[0,j+1]=x
#                 elif (j==71)&(i!=35):
#                     shifted_sig[i,0]=x
#                     shifted_sig[i+1,0]=x
#                 else:
#                     print(i,j)
#                     shifted_sig[0,0]=x
# shifted_esig = xr.Dataset(
#                 data_vars = dict(pv=(['latitude','longitude'],shifted_sig)),
#                 coords = dict(
#                         latitude = Clats-2.5,
#                 longitude = Clons-2.5))

## HADCRUT SURROGATES
# Clats = Csignificant.latitude
# Clons = Csignificant.longitude
# shifted_sig = np.full((len(Clats),len(Clons)),np.nan)
# for i, lat in enumerate(Clats):
#     for j, lon in enumerate(Clons):
#         x = Csignificant.sel(latitude=lat,longitude=lon).values
#         if not np.isnan(x):
#             shifted_sig[i,j]=x
#             shifted_sig[i,j+1]=x
#             shifted_sig[i+1,j]=x
#             shifted_sig[i+1,j+1]=x
# Csig = xr.Dataset(
#                 data_vars = dict(pv=(['latitude','longitude'],shifted_sig)),
#                 coords = dict(
#                         latitude = Clats-2.5,
#                 longitude = Clons-2.5))

## SALINITY SURROGATES
# Slats = ndt05_sig.lat
# Slons = ndt05_sig.lon
# shifted_sig = np.full((len(Slats),len(Slons)),np.nan)
# for i, lat in enumerate(Slats):
#     for j, lon in enumerate(Slons):
#         x = ndt05_sig.sel(lat=lat,lon=lon).values
#         if not np.isnan(x):
#             if (i!=172) & (j!=120):
#                 shifted_sig[i,j]=x
#                 shifted_sig[i,j+1]=x
#                 shifted_sig[i+1,j]=x
#                 shifted_sig[i+1,j+1]=x
#             else:
#                 if (i==172) & (j!=120):
#                     shifted_sig[0,j]=x
#                     shifted_sig[0,j+1]=x
#                 elif (j==120)&(i!=172):
#                     shifted_sig[i,0]=x
#                     shifted_sig[i+1,0]=x
#                 else:
#                     print(i,j)
#                     shifted_sig[0,0]=x
# Ssig = xr.Dataset(
#                 data_vars = dict(pv=(['lat','lon'],shifted_sig)),
#                 coords = dict(
#                         lat = Slats-0.5,
#                 lon = Slons-0.5))