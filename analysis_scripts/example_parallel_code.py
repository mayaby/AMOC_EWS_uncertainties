"""
This is an **example** of the mpi4py parallel codes that were used to run most of the analysis in the papers. It takes an SST dataset, divides it into slices along the longitude axis and assigns those slices to the different processors. In each slice you still need to use a loop to run the wrapper function on the timeseries in each gridcell.

The axis along which to slice can be adapted to the datasets - in cases where I had a surrogate/member axis I usually sliced along that. 

This is a very simple script. But because it uses comm.scatter() and comm.gather() the dimension along which you slice (here longitude) needs to be divisible by the number of processors. In this paper this was generally not a problem as my dimensions were nice round numbers. For a more general but more complicated code using comm.Scatterv() and comm.Gatherv() you can email maya.ben-yami@tum.de.  

Also note that comm.scatter() and comm.gather() have a 2GB limit, so if the file is larger you might need to divide it into smaller slices inside a loop. 
"""

import xarray as xr
from EWS_functions import *
from EWS_wrappers import *
import time
import sys
from mpi4py import MPI

ws = 60 ## window size
start_time = time.time()
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
print('Processor {}/{}'.format(rank,size))

if rank==0:
    ## in processor 0 we open the file, get the size parameters and divide the data into a list of longitude slices for the processors
    ## (to divide a difference coordinate, e.g. the ensemble member coordinate, just add the coordinate and substitute ncoord for nlon and the dimension name in lon_list)

    ds = xr.open_dataset('ERSSTv5/sst.mnmean_ymean.nc')
    y = ds.sst[:-1]

    ##########
    nt, nlat, nlon = y.shape[0],y.shape[1],y.shape[2]
    ##########

    ## this is a very basic code, in which the number of longitudes needs to be divisible by the number of processors (for a code that adapts to any number please contact me)
    if not nlon % size == 0: 
        sys.exit('Wrong number of processors')
    num_per_rank = nlon /size
    print('{} tasks per processor'.format(num_per_rank))

    ## create list of slices of the data by longitude and divide the list between the processors
    ## (for different data structures it's sometimes easier to divide the list of longitudes instead of data slices, but note that that doesn't help with size issues because comm.gather() will still be gathering data slices)

    lon_list = [y.isel(lon=i).values for i in np.arange(0,nlon)]
    print(len(lon_list))
    in_list = [lon_list[i:i+int(num_per_rank)] for i in np.arange(0,nlon,int(num_per_rank))]
    print('List for scatter: size {}'.format(len(in_list)))
    if len(in_list)!=size:
        print('Mismatch between number of processors and list for scatter')
        sys.exit('Wrong number of processors')
else:
    in_list = None

## scatter the list of lists from processor 0 to all processors
v = comm.scatter(in_list,root=0)

num_per_rank=len(v)
dat = v[0]
nt, nlat = dat.shape[0],dat.shape[1]
print(nt, nlat)
print('Lists in scatter: size {}'.format(num_per_rank))

res = []

## in each processor we go through the number of slices per processor
for l in np.arange(0,num_per_rank):
    print("Sample number {}/{} being done by processor {}/{}".format(l+1,num_per_rank,rank,size))
    smp = v[int(l)] ## get a longitude slice (dimenstion nt, nlat)
    lambdas = np.full((nt, nlat),np.nan)

    for i in np.arange(0,nlat):
        ts = smp[:,i]
        lambdas[:,i]= lambda_wrapper_rmean(ts,ws=ws) ## or whichever other wrapper is appropriate

    res.append(lambdas)

print('made lambdas')
res = np.array(res) ## make all the slices into an array on each processor

g = comm.gather(res,root=0) ## gather all the arrays from the different processors


if comm.rank==0:
    lams = np.concatenate(g)
    lams = np.transpose(lams, (1,2,0)) ## we added up the slices with lon as a first dim, so swap axes around so longitude is last dim

    ds = xr.open_dataset('/p/tmp/mayayami/sst-uncertainties-data/ERSSTv5/sst.mnmean_ymean.nc')
    data = ds.sst[:-1]

    array = data.copy(deep=False, data=lams)

    ## (if you don't swap axes around then you can't copy the shape of the original array and need to create a new one, see below)
    # array = xr.Dataset(
    #     data_vars = dict(lams=(['longitude','time','latitude'],lams)),
    #     coords = dict(
    #         time = data.time,
    #         latitude = data.latitude,
    #         longitude      = data.longitude)
    # )

    print('made dataset')
    print("--- %s seconds ---" % (time.time() - start_time))

    opath = 'ERSSTv5/ERSSTv5_lambdas_w60_rmean.nc'
    array.to_netcdf(opath)

    print(opath)
    print('saved')
    print("--- %s seconds ---" % (time.time() - start_time))