This is the github repository for the Nature Communications paper "Uncertainties in critical slowing down indicators of observation-based fingerprints of the Atlantic Overturning Circulation". https://www.nature.com/articles/s41467-023-44046-9

If you have any questions please email maya.ben-yami@tum.de

--------------------------

## Plotting notebooks

Two jupyter notebooks contain the python code for plotting the figures from the paper. Main_Plots.ipynb has the code for the main text figures, and SI_Plots.ipynb for the SI figures. 

## Analysis scripts

The folder analysis_scripts contains python and bash scripts which introduce all the code used to conduct the analysis in the paper. Note - most of these scripts are not operational, but rather introduce functions/ show examples of the code used. To get the results in the paper these scripts need to be modified to fit individual file structures, which differ across datasets. All the information needed to reproduce the analysis is contained in these scripts. 

* EWS_functions.py - functions for calculation CSD indicators (EWS)
* EWS_wrappers.py - wrappers for the above functions, with or without detrending
* count_profiles.py - goes through the monthly salinity profile files in EN4 and counts the number of "good" profiles in a given lat-lon gridcell in the analysis grid. Only used for visualisation in figure S5
* example_parallel_code.py - an example of the mpi4py parallel codes that were used to run most of the CSD analysis in the paper
* example_shift_grid.py - an example of the grid pre-processing necesseary for correct cross hatching in plots like figure 6
* generate_surrogates.py -  contains the wrapper functions to: generate fourier surrogates, generate AR2 surrogates and to modify a surrogate array using the salinity weights method described in the paper
* get_trends.py - calculates the linear trend of input timeseries using xarray
* process_files.sh - all the various cdo scripts used to pre-process files (select areas, get annual average, etc.)

## Examples

Here are two examples of the work processes:

1. To get the file with the number of salinity profiles, you need to:

* Download EN.4.2.2.p.analysis.g10.202404.nc.gz and EN.4.2.2.p.profiles.g10.202404.nc.gz from https://www.metoffice.gov.uk/hadobs/en4/download-en4-2-2.html
* Extract the files from the zip
* Run `cdo mergetime -selvar,salinity EN.4.2.2.f.analysis.g10*.nc salinity.nc` to get the merged salinity analysis file
* Run `cdo -L vertmean -sellevidx,1/19 salinity.nc salinity_dmean.nc` to get the file of the mean salinities in the top 300m
* Change the file locations in count_profiles.py to the locations on your computer
* Adjust the `[:-1]` in line 6 to get only complete years, and the maximum year in line 14 to match your data
* Run count_profiles.py

2. To get the lambda time series of the HadSST4 AMOC SST index from Figure 1a) and 1d), you need to:

* Download HadSST.4.0.1.0_median.nc  from https://www.metoffice.gov.uk/hadobs/hadsst4/data/download.html
* Run `cdo -L sub -fldmean -sellonlatbox,-55,-20,46,61 HadSST.4.0.1.0_median.nc -fldmean HadSST.4.0.1.0_median.nc HadSST.4.0.1.0_median_amoc.nc` to get the AMOC index
* Run `cdo -L divdpy -yearsum -muldpm HadSST.4.0.1.0_median_amoc.nc HadSST.4.0.1.0_median_amoc_ymean.nc` to get the annual average
* And then you can write your own python code to apply lambda_wrapper_rmean() from EWS_wrapper.py to the time series from HadSST.4.0.1.0_median_amoc_ymean.nc

