#### file that goes through the monthly salinity profile files in EN4 and counts the number of "good" profiles in a given lat-lon gridcell in the analysis grid. Only used for visualisation in figure S5. #####

import numpy as np
import xarray as xr

dmeanm = xr.open_dataset('/p/tmp/mayayami/sst-uncertainties-data/ENS4_analyses/salinity_dmean.nc').salinity.squeeze()[:-1] # open analysis netcdf file to get lat-lon grid

lats = dmeanm.lat.values
lons = dmeanm.lon.values
nlon = len(lons)
nlat = len(lats)
times = dmeanm.time[:-4]

yrs = np.arange(1900,2022) # don't include 2022 because incomplete
months = ['01','02','03','04','05','06','07','08','09','10','11','12']

nt = len(yrs)*12

density_grid = np.zeros((nt,nlat,nlon))


for iyr, yr in enumerate(yrs):
    for im, month in enumerate(months):
        time = '{}'.format(yr)+month
        print(time)

        it = iyr*12 + im
        data = xr.open_dataset('/p/tmp/mayayami/sst-uncertainties-data/ENS4_profiles/EN.4.2.2.f.profiles.g10.{}.nc'.format(time))

        ## count all profiles
        # plats = data['LATITUDE'][:].values
        # plons = data['LONGITUDE'][:].values
        # nprof = int(len(data.N_PROF.values))

        ## count only "good" profiles
        psal_qc = data['PROFILE_PSAL_QC'].astype('unicode')
        good_profs = np.where(psal_qc=='1')
        plats = data['LATITUDE'][:].values[good_profs]
        plons = data['LONGITUDE'][:].values[good_profs]
        nprof = int(len(data.N_PROF.values[good_profs]))

        nlev = int(len(data.N_LEVELS))

        idc = np.where(plons<0)
        plons[idc] = plons[idc]+360

        new_plats = np.array([round(lat) for lat in plats])
        new_plons = np.array([round(lon) for lon in plons])

        for i, plat in enumerate(new_plats):
            plon = new_plons[i]
            ilon = np.where(lons==plon)
            ilat = np.where(lats==plat)
            density_grid[it,ilat,ilon] = density_grid[it,ilat,ilon]+1

density_dataset = xr.Dataset(
            data_vars = dict(
                nprofiles = (['time','lat','lon'],density_grid)),
            coords = dict(
                    time = dmeanm.time[:-4],
                    lat = dmeanm.lat,
                    lon = dmeanm.lon)
        )

density_dataset.to_netcdf('profile_density_good.nc')
